#!/usr/bin/env python3
"""Prepare a fail-fast replay of the deposited GSE302339 author workflow.

The Zenodo record does not contain an annotated AnnData object or a cell-level
metadata table.  It does, however, contain enough information to reproduce the
annotation run: exact core package versions, the 65-sample processing order,
recorded matrix shapes, and the object-specific cluster-to-label dictionaries.

This script extracts clean copies of notebooks 1 and 2, reconstructs the exact
sample order from the notebook output and the H5 barcode counts, writes the
185-gene ribosomal helper directly from the authors' deposited code, and adds
guards that stop before any author label is assigned if the replay diverges.
"""

from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import importlib.metadata
import importlib.util
import json
import re
import shutil
import sys
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


ZIP_CANDIDATES = (
    "data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip",
    "data/raw/downloads/zenodo/16341197/scanpy_workflow.zip",
)
PREPROCESSING_MEMBER = "scanpy_workflow/1_preprocessing_doublet_detection.ipynb"
ANNOTATION_MEMBER = "scanpy_workflow/2_celltype_annotation.ipynb"
DEFAULT_INPUT_DIR = "input/data_cellranger8"
DEFAULT_CODE_DIR = "intermediate/gse302339_author_exact_code"
DEFAULT_ORDER_FILE = "input/gse302339_author_exact/sample_order.txt"
DEFAULT_RIBOSOMAL_FILE = (
    "input/gse302339_author_exact/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv"
)
DEFAULT_OUTPUT_DIR = "output/gse302339_author_exact"
DEFAULT_META_FILE = "input/meta_cr8.csv"
DEFAULT_MIN_FREE_GIB = 25.0
EXPECTED_AUTHOR_SAMPLES = 65
EXPECTED_FULL_CLUSTER_COUNT = 62
EXPECTED_NOTEBOOK_ZIP_MD5 = "386d4d6cb2e7f813a13bb7aee2f866a7"

INFERRED_RUNTIME = {
    "doubletdetection": "4.2",
    "harmonypy": "0.0.10",
    "leidenalg": "0.10.2",
}
PACKAGE_DISTRIBUTIONS = {
    "scanpy": "scanpy",
    "anndata": "anndata",
    "umap": "umap-learn",
    "numpy": "numpy",
    "scipy": "scipy",
    "pandas": "pandas",
    "scikit-learn": "scikit-learn",
    "statsmodels": "statsmodels",
    "igraph": "igraph",
    "louvain": "louvain",
    "pynndescent": "pynndescent",
    "doubletdetection": "doubletdetection",
    "harmonypy": "harmonypy",
    "leidenalg": "leidenalg",
}


def resolve_zip(explicit: str | None) -> Path:
    if explicit:
        return project_path(explicit)
    for candidate in ZIP_CANDIDATES:
        path = project_path(candidate)
        if path.exists():
            return path
    return project_path(ZIP_CANDIDATES[0])


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_notebook(zip_path: Path, member: str) -> dict:
    with zipfile.ZipFile(zip_path) as archive:
        return json.loads(archive.read(member).decode("utf-8", errors="replace"))


def cell_source(cell: dict) -> str:
    source = cell.get("source", "")
    return "".join(source) if isinstance(source, list) else str(source)


def output_text(cell: dict) -> str:
    chunks: list[str] = []
    for output in cell.get("outputs", []):
        value = output.get("text", "")
        if isinstance(value, list):
            chunks.append("".join(str(item) for item in value))
        elif value:
            chunks.append(str(value))
        plain = output.get("data", {}).get("text/plain", "")
        if isinstance(plain, list):
            chunks.append("".join(str(item) for item in plain))
        elif plain:
            chunks.append(str(plain))
    return "".join(chunks)


def notebook_contract(preprocessing: dict, annotation: dict) -> dict:
    python_version = str(
        preprocessing.get("metadata", {}).get("language_info", {}).get("version", "")
    )
    header = "\n".join(output_text(cell) for cell in preprocessing.get("cells", [])[:3])
    runtime = dict(re.findall(r"([A-Za-z0-9_.-]+)==([^\s]+)", header))
    runtime.update(INFERRED_RUNTIME)

    sample_cell = next(
        (
            cell
            for cell in preprocessing.get("cells", [])
            if "Total number of cells:" in output_text(cell)
        ),
        None,
    )
    if sample_cell is None:
        raise ValueError("Could not recover per-sample checkpoints from notebook output.")
    sample_output = output_text(sample_cell)
    raw_totals = [
        int(value)
        for value in re.findall(r"Total number of cells:\s*(\d+)", sample_output)
    ]
    post_low_quality_counts = [
        int(value)
        for value in re.findall(
            r"Number of cells after filtering of low quality cells:\s*(\d+)",
            sample_output,
        )
    ]
    post_doublet_shapes = [
        (int(rows), int(columns))
        for rows, columns in re.findall(
            r"(?<![\w.])\((\d+),\s*(\d+)\)", sample_output
        )
    ]
    checkpoint_lengths = {
        len(raw_totals), len(post_low_quality_counts), len(post_doublet_shapes)
    }
    if checkpoint_lengths != {EXPECTED_AUTHOR_SAMPLES}:
        raise ValueError(
            "Expected 65 raw/low-quality/post-doublet sample checkpoints; observed "
            f"raw={len(raw_totals)}, low_quality={len(post_low_quality_counts)}, "
            f"post_doublet={len(post_doublet_shapes)}."
        )

    preprocessing_text = "\n".join(
        output_text(cell) for cell in preprocessing.get("cells", [])
    )
    shapes = [
        (int(rows), int(columns))
        for rows, columns in re.findall(r"\((\d+),\s*(\d+)\)", preprocessing_text)
    ]
    concat_shape = (160620, 18941)
    hvg_shape = (160620, 2323)
    if concat_shape not in shapes or hvg_shape not in shapes:
        raise ValueError(
            "Deposited notebook outputs did not contain the expected matrix checkpoints."
        )
    if sum(rows for rows, _ in post_doublet_shapes) != concat_shape[0]:
        raise ValueError(
            "Per-sample post-doublet cell counts do not sum to the notebook concat checkpoint."
        )

    annotation_text = "\n".join(output_text(cell) for cell in annotation.get("cells", []))
    subcluster_counts = [
        int(value)
        for value in re.findall(r"finished: found\s+(\d+)\s+clusters", annotation_text)
    ]
    if len(subcluster_counts) < 2:
        raise ValueError("Could not recover parenchyma/immune cluster counts from notebook.")

    return {
        "python": python_version,
        "runtime": runtime,
        "raw_cell_totals_in_author_order": raw_totals,
        "post_low_quality_cell_counts_in_author_order": post_low_quality_counts,
        "post_doublet_shapes_in_author_order": [
            list(shape) for shape in post_doublet_shapes
        ],
        "n_author_samples": len(raw_totals),
        "raw_cell_totals_are_unique": len(raw_totals) == len(set(raw_totals)),
        "concat_shape": list(concat_shape),
        "post_hvg_shape": list(hvg_shape),
        "full_cluster_count": EXPECTED_FULL_CLUSTER_COUNT,
        "full_cluster_ids": [str(value) for value in range(EXPECTED_FULL_CLUSTER_COUNT)],
        "parenchyma_cluster_count": subcluster_counts[0],
        "immune_cluster_count": subcluster_counts[1],
        "harmony_iterations_in_notebook_log": 26,
    }


def installed_runtime(contract: dict) -> list[dict]:
    rows: list[dict] = []
    expected_python = str(contract["python"])
    observed_python = ".".join(map(str, sys.version_info[:3]))
    rows.append(
        {
            "component": "python",
            "expected": expected_python,
            "observed": observed_python,
            "match": observed_python == expected_python,
            "evidence": "notebook metadata",
        }
    )
    printed = set(contract["runtime"]) - set(INFERRED_RUNTIME)
    for component, expected in contract["runtime"].items():
        distribution = PACKAGE_DISTRIBUTIONS.get(component, component)
        try:
            observed = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            observed = "missing"
        rows.append(
            {
                "component": component,
                "distribution": distribution,
                "expected": str(expected),
                "observed": observed,
                "match": observed == str(expected),
                "evidence": (
                    "notebook print_header"
                    if component in printed
                    else "paper/traceback or author-era release"
                ),
            }
        )
    return rows


def h5_barcode_count(path: Path) -> int:
    try:
        import h5py
    except ImportError as exc:  # pragma: no cover - environment diagnostic
        raise RuntimeError("h5py is required to inspect Cell Ranger H5 files") from exc
    with h5py.File(path, "r") as handle:
        if "matrix/barcodes" in handle:
            return int(handle["matrix/barcodes"].shape[0])
        if "matrix/shape" in handle:
            shape = handle["matrix/shape"][()]
            return int(shape[1])
    raise ValueError(f"No Cell Ranger barcode dataset found in {path}")


def reconstruct_sample_order(
    input_dir: Path,
    expected_totals: list[int],
    expected_low_quality_counts: list[int],
    expected_post_doublet_shapes: list[list[int]],
) -> tuple[list[dict], list[str]]:
    h5_paths = sorted(input_dir.glob("*_filtered_feature_bc_matrix.h5"))
    observed: dict[int, list[Path]] = {}
    for path in h5_paths:
        observed.setdefault(h5_barcode_count(path), []).append(path)

    expected_counter = Counter(expected_totals)
    rows: list[dict] = []
    failures: list[str] = []
    for position, total in enumerate(expected_totals, start=1):
        candidates = observed.get(total, [])
        filename = candidates[0].name if len(candidates) == 1 else ""
        if len(candidates) != 1:
            failures.append(
                f"author position {position} raw count {total}: found {len(candidates)} H5 candidates"
            )
        rows.append(
            {
                "author_position": position,
                "raw_cell_count": total,
                "expected_post_low_quality_cells": expected_low_quality_counts[
                    position - 1
                ],
                "expected_post_doublet_cells": expected_post_doublet_shapes[
                    position - 1
                ][0],
                "expected_post_doublet_genes": expected_post_doublet_shapes[
                    position - 1
                ][1],
                "h5_filename": filename,
                "n_h5_candidates": len(candidates),
            }
        )

    extra_counts = sorted(set(observed) - set(expected_counter))
    if extra_counts:
        failures.append(f"H5 raw-cell counts absent from notebook outputs: {extra_counts}")
    if len(h5_paths) != len(expected_totals):
        failures.append(
            f"expected {len(expected_totals)} H5 files from notebook; found {len(h5_paths)}"
        )
    if any(value != 1 for value in expected_counter.values()):
        failures.append("Notebook raw-cell totals are not unique; exact order is ambiguous")
    return rows, failures


def audit_metadata(meta_path: Path, order_rows: list[dict]) -> tuple[dict, list[str]]:
    required_columns = {
        "batch",
        "lobe_emphysema_simple",
        "patient",
        "total_emphysema",
        "lobe_emphysema",
    }
    summary = {
        "path": str(meta_path),
        "exists": meta_path.exists(),
        "required_columns": sorted(required_columns),
        "n_rows": 0,
        "n_expected_batches": len(order_rows),
        "n_expected_batches_with_one_match": 0,
    }
    failures: list[str] = []
    if not meta_path.exists():
        failures.append(f"author metadata file missing: {meta_path}")
        return summary, failures

    with meta_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    summary["n_rows"] = len(rows)
    columns = set(rows[0]) if rows else set()
    missing_columns = sorted(required_columns - columns)
    summary["missing_required_columns"] = missing_columns
    if missing_columns:
        failures.append(f"metadata missing required columns: {missing_columns}")
        return summary, failures

    batch_counts = Counter(str(row.get("batch", "")) for row in rows)
    expected_batches = [
        f"{DEFAULT_INPUT_DIR}/{str(row['h5_filename']).removesuffix('_filtered_feature_bc_matrix.h5')}"
        for row in order_rows
        if row.get("h5_filename")
    ]
    matched = sum(batch_counts[batch] == 1 for batch in expected_batches)
    summary["n_expected_batches_with_one_match"] = matched
    bad_batches = [batch for batch in expected_batches if batch_counts[batch] != 1]
    summary["first_bad_batches"] = bad_batches[:10]
    if bad_batches:
        failures.append(
            f"metadata must contain exactly one row for every exact author batch; "
            f"{len(bad_batches)} of {len(expected_batches)} failed"
        )
    rows_by_batch = {
        str(row.get("batch", "")): row
        for row in rows
        if batch_counts[str(row.get("batch", ""))] == 1
    }
    phenotype_columns = sorted(required_columns - {"batch"})
    incomplete_batches = [
        batch
        for batch in expected_batches
        if batch in rows_by_batch
        and any(
            str(rows_by_batch[batch].get(column, "")).strip().casefold()
            in {"", "nan", "none", "na"}
            for column in phenotype_columns
        )
    ]
    summary["first_batches_with_incomplete_phenotype"] = incomplete_batches[:10]
    if incomplete_batches:
        failures.append(
            f"{len(incomplete_batches)} exact author metadata rows contain empty phenotype/donor values"
        )
    return summary, failures


def extract_ribosomal_genes(preprocessing_code: str) -> list[str]:
    tree = ast.parse(preprocessing_code)
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "ribosome" for target in node.targets):
            continue
        value = node.value
        if not isinstance(value, ast.Call) or not value.args:
            continue
        tuple_arg = value.args[0]
        if not isinstance(tuple_arg, (ast.Tuple, ast.List)):
            continue
        genes = [
            str(element.value)
            for element in tuple_arg.elts
            if isinstance(element, ast.Constant) and isinstance(element.value, str)
        ]
        if genes:
            return genes
    raise ValueError("Could not recover the deposited 185-gene ribosomal tuple")


def load_export_module():
    path = Path(__file__).with_name("03_extract_gse302339_annotation_code.py")
    spec = importlib.util.spec_from_file_location("gse302339_code_export", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def replace_once(text: str, old: str, new: str, label: str) -> str:
    if old not in text:
        raise ValueError(f"Exact-replay patch point missing: {label}")
    return text.replace(old, new, 1)


def patch_preprocessing_code(
    path: Path,
    *,
    order_file: Path,
    ribosomal_file: Path,
    output_dir: Path,
    contract: dict,
) -> None:
    text = path.read_text(encoding="utf-8")
    if "CODEX-AUTHOR-EXACT-REPLAY" in text:
        return

    text = replace_once(
        text,
        "# This is for audit/reconstruction; it is not a cleaned executable script.",
        "# This is for audit/reconstruction; it is not a cleaned executable script.\n"
        "# CODEX-AUTHOR-EXACT-REPLAY: deterministic order, seed, and fail-fast gates.",
        "replay marker",
    )
    text = text.replace("'output/", f"'{output_dir.as_posix()}/")
    text = text.replace('"output/', f'"{output_dir.as_posix()}/')
    concat_pickle_block = (
        f"with open('{output_dir.as_posix()}/adata_concatenated_cr8', \"wb\") as file:\n"
        "    pickle.dump(adata, file)"
    )
    text = replace_once(
        text,
        concat_pickle_block,
        "pass  # CODEX-AUTHOR-EXACT-REPLAY: unused 2.3-GB concat checkpoint not serialized",
        "skip unused concat pickle",
    )
    text = replace_once(
        text,
        "pass  # NOTEBOOK-OPTIONAL-SCVI skipped for current environment: scvi.settings.seed = 0",
        "# Equivalent to the notebook's scvi.settings.seed = 0 for NumPy-based steps.\n"
        "import random\nrandom.seed(0)\nnp.random.seed(0)",
        "author seed",
    )
    text = replace_once(
        text,
        "    pseudocount=0.1,\n"
        "    n_jobs=-1)",
        "    pseudocount=0.1,\n"
        "    random_state=(\n"
        "        None\n"
        "        if os.environ.get('GSE302339_DD_RANDOM_STATE', '0').lower() == 'none'\n"
        "        else int(os.environ.get('GSE302339_DD_RANDOM_STATE', '0'))\n"
        "    ),\n"
        "    # Keep n_jobs as in the deposited notebook. The explicit random_state\n"
        "    # records the library default and makes the replay contract visible.\n"
        "    n_jobs=int(os.environ.get('GSE302339_DD_N_JOBS', '-1')))",
        "explicit doubletdetection seed",
    )
    text = text.replace(
        '"input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv"',
        f'"{ribosomal_file.as_posix()}"',
    )
    text = replace_once(
        text,
        "def pp(h5_path):",
        "def pp(h5_path, expected_raw_count, expected_low_quality_count, expected_post_doublet_shape):",
        "per-sample checkpoint arguments",
    )
    text = replace_once(
        text,
        'print(f"Total number of cells: {adata.n_obs}")',
        'print(f"Total number of cells: {adata.n_obs}")\n'
        "    if adata.n_obs != expected_raw_count:\n"
        "        raise RuntimeError(f'Exact replay raw-cell mismatch for {h5_path}: expected {expected_raw_count}, got {adata.n_obs}')",
        "per-sample raw-cell checkpoint",
    )
    text = replace_once(
        text,
        'print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")',
        'print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")\n'
        "    if adata.n_obs != expected_low_quality_count:\n"
        "        raise RuntimeError(f'Exact replay low-quality-filter mismatch for {h5_path}: expected {expected_low_quality_count}, got {adata.n_obs}')",
        "per-sample low-quality checkpoint",
    )
    text = replace_once(
        text,
        "    print(adata.shape)\n    sc.pl.violin",
        "    print(f\"Doublets removed: {int(adata.uns['doublets_removed'])}\")\n"
        "    print(adata.shape)\n"
        "    if tuple(adata.shape) != tuple(expected_post_doublet_shape):\n"
        "        raise RuntimeError(f'Exact replay post-doublet mismatch for {h5_path}: expected {tuple(expected_post_doublet_shape)}, got {adata.shape}')\n"
        "    sc.pl.violin",
        "per-sample post-doublet checkpoint",
    )

    old_loop = """os.listdir('input/data_cellranger8/')
pass  # NOTEBOOK-ONLY: !rm -r \"input/data_cellranger8/.DS_Store\"

# %% [cell 7]
out = []
for file in os.listdir('input/data_cellranger8/'):
    out.append(pp(\"input/data_cellranger8/\" + file))"""
    new_loop = f"""author_order_path = {order_file.as_posix()!r}
with open(author_order_path, encoding='utf-8') as handle:
    author_files = [line.strip() for line in handle if line.strip()]
expected_raw_counts = {contract['raw_cell_totals_in_author_order']!r}
expected_low_quality_counts = {contract['post_low_quality_cell_counts_in_author_order']!r}
expected_post_doublet_shapes = {contract['post_doublet_shapes_in_author_order']!r}
discovered_files = sorted(
    file for file in os.listdir('input/data_cellranger8/')
    if file.endswith('_filtered_feature_bc_matrix.h5')
)
if sorted(author_files) != discovered_files:
    raise RuntimeError('Exact-replay sample manifest does not match input H5 files')
diagnostic_max_samples = int(os.environ.get('GSE302339_REPLAY_MAX_SAMPLES', '0') or '0')
print(
    'DOUBLETD_REPLAY_SETTINGS '
    f"random_state={{os.environ.get('GSE302339_DD_RANDOM_STATE', '0')}} "
    f"n_jobs={{os.environ.get('GSE302339_DD_N_JOBS', '-1')}} "
    f"max_samples={{diagnostic_max_samples}}"
)

# %% [cell 7]
out = []
for position, file in enumerate(author_files):
    print(f'AUTHOR_SAMPLE_ORDER {{position + 1}}/{{len(author_files)}}: {{file}}')
    out.append(
        pp(
            \"input/data_cellranger8/\" + file,
            expected_raw_counts[position],
            expected_low_quality_counts[position],
            expected_post_doublet_shapes[position],
        )
    )
    if diagnostic_max_samples and position + 1 >= diagnostic_max_samples:
        print(f'CODEX-AUTHOR-EXACT-REPLAY diagnostic stop after {{position + 1}} sample(s)')
        raise SystemExit(0)"""
    text = replace_once(text, old_loop, new_loop, "author sample order")

    concat_shape = tuple(contract["concat_shape"])
    text = replace_once(
        text,
        "print(adata.shape)\n\n# %% [cell 9]",
        f"print(adata.shape)\nif adata.shape != {concat_shape!r}:\n"
        f"    raise RuntimeError('Exact replay failed concat checkpoint: expected {concat_shape!r}, got ' + str(adata.shape))\n\n"
        "# %% [cell 9]",
        "concat checkpoint",
    )
    hvg_shape = tuple(contract["post_hvg_shape"])
    text = replace_once(
        text,
        "print(adata.shape)\nsc.pp.scale(adata, max_value=10)",
        f"print(adata.shape)\nif adata.shape != {hvg_shape!r}:\n"
        f"    raise RuntimeError('Exact replay failed HVG checkpoint: expected {hvg_shape!r}, got ' + str(adata.shape))\n"
        "sc.pp.scale(adata, max_value=10)",
        "HVG checkpoint",
    )
    full_cluster_count = int(contract["full_cluster_count"])
    full_cluster_ids = set(contract["full_cluster_ids"])
    text = replace_once(
        text,
        "sc.tl.leiden(adata,resolution=3,flavor='igraph',n_iterations=2)",
        "sc.tl.leiden(adata,resolution=3,flavor='igraph',n_iterations=2)\n"
        f"required_full_clusters = {full_cluster_ids!r}\n"
        "observed_full_clusters = set(adata.obs['leiden'].astype(str))\n"
        f"if adata.obs['leiden'].nunique() != {full_cluster_count} or observed_full_clusters != required_full_clusters:\n"
        f"    raise RuntimeError('Exact replay full Leiden mismatch: expected {full_cluster_count} clusters 0-{full_cluster_count - 1}; observed ' + str(sorted(observed_full_clusters)))",
        "full-cluster checkpoint",
    )
    text = replace_once(
        text,
        "import os\n",
        f"import os\nos.makedirs({output_dir.as_posix()!r}, exist_ok=True)\n",
        "exact output directory",
    )
    path.write_text(text, encoding="utf-8")


def patch_annotation_code(path: Path, *, output_dir: Path, contract: dict) -> None:
    text = path.read_text(encoding="utf-8")
    if "CODEX-AUTHOR-EXACT-REPLAY" in text:
        return
    text = replace_once(
        text,
        "# This is for audit/reconstruction; it is not a cleaned executable script.",
        "# This is for audit/reconstruction; it is not a cleaned executable script.\n"
        "# CODEX-AUTHOR-EXACT-REPLAY: object-local author labels with fail-fast gates.",
        "replay marker",
    )
    text = text.replace("'output/", f"'{output_dir.as_posix()}/")
    text = text.replace('"output/', f'"{output_dir.as_posix()}/')
    text = replace_once(
        text,
        "pass  # NOTEBOOK-OPTIONAL-SCVI skipped for current environment: scvi.settings.seed = 0",
        "# Equivalent to the notebook's scvi.settings.seed = 0 for the used code paths.\n"
        "import random\nrandom.seed(0)\nnp.random.seed(0)",
        "annotation seed",
    )

    expected_shape = tuple(contract["post_hvg_shape"])
    expected_full_cluster_ids = set(contract["full_cluster_ids"])
    expected_full_cluster_count = int(contract["full_cluster_count"])
    broad_guard = f"""if tuple(adata.shape) != {expected_shape!r}:
    raise RuntimeError(
        'Exact replay annotation input shape mismatch: expected {expected_shape!r}, got '
        + str(adata.shape)
    )
expected_full_cluster_ids = {expected_full_cluster_ids!r}
observed_full_cluster_ids = set(adata.obs['leiden'].astype(str))
if adata.obs['leiden'].nunique() != {expected_full_cluster_count} or observed_full_cluster_ids != expected_full_cluster_ids:
    raise RuntimeError(
        'Exact replay annotation input must contain author full clusters 0-{expected_full_cluster_count - 1}; observed '
        + str(sorted(observed_full_cluster_ids))
    )
required_broad_clusters = set(
    cluster for clusters in celldict_level1.values() for cluster in clusters
)
observed_broad_clusters = set(adata.obs['leiden'].astype(str))
missing_broad_clusters = sorted(required_broad_clusters - observed_broad_clusters)
if missing_broad_clusters:
    raise RuntimeError(
        'Exact replay broad annotation mismatch; missing author clusters: '
        + str(missing_broad_clusters)
    )

"""
    text = replace_once(
        text,
        "adata.obs['celltype_level1']=np.NaN",
        broad_guard + "adata.obs['celltype_level1']=np.NaN",
        "broad-label guard",
    )

    parenchyma_count = int(contract["parenchyma_cluster_count"])
    parenchyma_ids = {str(value) for value in range(parenchyma_count)}
    text = replace_once(
        text,
        "sc.tl.leiden(parenchyma,resolution=2,flavor='igraph',n_iterations=2)",
        "sc.tl.leiden(parenchyma,resolution=2,flavor='igraph',n_iterations=2)\n"
        "observed_parenchyma_cluster_ids = set(parenchyma.obs['leiden'].astype(str))\n"
        f"expected_parenchyma_cluster_ids = {parenchyma_ids!r}\n"
        f"if parenchyma.obs['leiden'].nunique() != {parenchyma_count} or observed_parenchyma_cluster_ids != expected_parenchyma_cluster_ids:\n"
        f"    raise RuntimeError('Exact replay parenchyma Leiden mismatch: expected IDs 0-{parenchyma_count - 1}; observed ' + str(sorted(observed_parenchyma_cluster_ids)))",
        "parenchyma checkpoint",
    )
    parenchyma_guard = """required_parenchyma_clusters = set(
    cluster for clusters in celldict_level1.values() for cluster in clusters
)
observed_parenchyma_cluster_ids = set(parenchyma.obs['leiden'].astype(str))
missing_parenchyma_clusters = sorted(
    required_parenchyma_clusters - observed_parenchyma_cluster_ids
)
if missing_parenchyma_clusters:
    raise RuntimeError(
        'Exact replay parenchyma annotation mismatch; missing author clusters: '
        + str(missing_parenchyma_clusters)
    )

"""
    text = replace_once(
        text,
        "parenchyma.obs['parenchyma_celltype_level1']=np.NaN",
        parenchyma_guard + "parenchyma.obs['parenchyma_celltype_level1']=np.NaN",
        "parenchyma-label guard",
    )

    immune_count = int(contract["immune_cluster_count"])
    immune_ids = {str(value) for value in range(immune_count)}
    text = replace_once(
        text,
        "sc.tl.leiden(immune,resolution=2,flavor='igraph',n_iterations=2)",
        "sc.tl.leiden(immune,resolution=2,flavor='igraph',n_iterations=2)\n"
        "observed_immune_cluster_ids = set(immune.obs['leiden'].astype(str))\n"
        f"expected_immune_cluster_ids = {immune_ids!r}\n"
        f"if immune.obs['leiden'].nunique() != {immune_count} or observed_immune_cluster_ids != expected_immune_cluster_ids:\n"
        f"    raise RuntimeError('Exact replay immune Leiden mismatch: expected IDs 0-{immune_count - 1}; observed ' + str(sorted(observed_immune_cluster_ids)))",
        "immune checkpoint",
    )

    optional_read = (
        f"with open('{output_dir.as_posix()}/adata_mergedmeta_abt_cr8', \"rb\") as file:"
    )
    clean_stop = (
        "print('Exact author major-lineage annotations completed and checkpointed.')\n"
        "raise SystemExit(0)\n\n"
    )
    text = replace_once(text, optional_read, clean_stop + optional_read, "clean stop")
    path.write_text(text, encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted({field for row in rows for field in row})
    preferred = [
        "author_position",
        "raw_cell_count",
        "expected_post_low_quality_cells",
        "expected_post_doublet_cells",
        "expected_post_doublet_genes",
        "h5_filename",
        "n_h5_candidates",
        "component",
        "distribution",
        "expected",
        "observed",
        "match",
        "evidence",
    ]
    order = [field for field in preferred if field in fields]
    order.extend(field for field in fields if field not in order)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=order, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare an exact, checkpointed replay of GSE302339 author notebooks 1 and 2."
    )
    parser.add_argument("--zip-path")
    parser.add_argument("--code-dir", default=DEFAULT_CODE_DIR)
    parser.add_argument("--order-file", default=DEFAULT_ORDER_FILE)
    parser.add_argument("--ribosomal-file", default=DEFAULT_RIBOSOMAL_FILE)
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--min-free-gib", type=float, default=DEFAULT_MIN_FREE_GIB)
    parser.add_argument("--prepare-code", action="store_true")
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results = ensure_results_dirs(load_config())
    zip_path = resolve_zip(args.zip_path)
    input_dir = project_path(DEFAULT_INPUT_DIR)
    meta_file = project_path(DEFAULT_META_FILE)
    code_dir = project_path(args.code_dir)
    order_file = project_path(args.order_file)
    ribosomal_file = project_path(args.ribosomal_file)
    output_dir = project_path(args.output_dir)

    disk_usage = shutil.disk_usage(project_path("."))
    free_gib = disk_usage.free / (1024**3)
    disk_summary = {
        "free_gib": round(free_gib, 2),
        "required_free_gib": float(args.min_free_gib),
        "ready": free_gib >= float(args.min_free_gib),
        "note": "The exact replay skips the unused concatenated pickle but retains the integrated and annotated objects.",
    }

    failures: list[str] = []
    if not disk_summary["ready"]:
        failures.append(
            f"insufficient free disk: {free_gib:.2f} GiB available; "
            f"{args.min_free_gib:.2f} GiB required"
        )
    order_failures: list[str] = []
    metadata_summary: dict = {
        "path": str(meta_file),
        "exists": meta_file.exists(),
        "ready": False,
    }
    preprocessing_path: Path | None = None
    annotation_path: Path | None = None
    if not zip_path.exists():
        failures.append(f"notebook archive missing: {zip_path}")
        contract: dict = {}
        runtime_rows: list[dict] = []
        order_rows: list[dict] = []
    else:
        zip_md5 = file_digest(zip_path, "md5")
        if zip_md5 != EXPECTED_NOTEBOOK_ZIP_MD5:
            failures.append(
                f"notebook ZIP MD5 expected {EXPECTED_NOTEBOOK_ZIP_MD5}; observed {zip_md5}"
            )
        preprocessing = read_notebook(zip_path, PREPROCESSING_MEMBER)
        annotation = read_notebook(zip_path, ANNOTATION_MEMBER)
        contract = notebook_contract(preprocessing, annotation)
        runtime_rows = installed_runtime(contract)
        mismatches = [
            f"{row['component']} expected {row['expected']} observed {row['observed']}"
            for row in runtime_rows
            if not row["match"]
        ]
        failures.extend(f"runtime mismatch: {message}" for message in mismatches)

        if not input_dir.exists():
            order_rows = []
            failures.append(f"Cell Ranger H5 input directory missing: {input_dir}")
        else:
            order_rows, order_failures = reconstruct_sample_order(
                input_dir,
                contract["raw_cell_totals_in_author_order"],
                contract["post_low_quality_cell_counts_in_author_order"],
                contract["post_doublet_shapes_in_author_order"],
            )
            failures.extend(order_failures)

        if order_rows and not order_failures:
            order_file.parent.mkdir(parents=True, exist_ok=True)
            order_file.write_text(
                "\n".join(str(row["h5_filename"]) for row in order_rows) + "\n",
                encoding="utf-8",
            )

        metadata_summary, metadata_failures = audit_metadata(meta_file, order_rows)
        metadata_summary["ready"] = not metadata_failures
        failures.extend(metadata_failures)

        if args.prepare_code:
            export_module = load_export_module()
            code_dir.mkdir(parents=True, exist_ok=True)
            preprocessing_path, _ = export_module.export_notebook_code(
                zip_path, PREPROCESSING_MEMBER, code_dir
            )
            annotation_path, _ = export_module.export_notebook_code(
                zip_path, ANNOTATION_MEMBER, code_dir
            )
            genes = extract_ribosomal_genes(
                preprocessing_path.read_text(encoding="utf-8")
            )
            if len(genes) != 185:
                failures.append(
                    f"deposited ribosomal tuple expected 185 genes; recovered {len(genes)}"
                )
            ribosomal_file.parent.mkdir(parents=True, exist_ok=True)
            ribosomal_file.write_text("\n".join(genes) + "\n", encoding="utf-8")
            output_dir.mkdir(parents=True, exist_ok=True)
            patch_preprocessing_code(
                preprocessing_path,
                order_file=order_file,
                ribosomal_file=ribosomal_file,
                output_dir=output_dir,
                contract=contract,
            )
            patch_annotation_code(
                annotation_path, output_dir=output_dir, contract=contract
            )

    runtime_path = results["tables"] / "gse302339_author_exact_runtime_manifest.csv"
    order_path = results["tables"] / "gse302339_author_exact_sample_order.csv"
    summary_path = results["meta"] / "gse302339_author_exact_replay_preflight.json"
    write_csv(runtime_path, runtime_rows)
    write_csv(order_path, order_rows)

    summary = {
        "dataset": "GSE302339",
        "strategy": "checkpoint-matched deposited-notebook author-workflow replay",
        "zip_path": str(zip_path),
        "zip_md5": file_digest(zip_path, "md5") if zip_path.exists() else "",
        "expected_zip_md5": EXPECTED_NOTEBOOK_ZIP_MD5,
        "input_dir": str(input_dir),
        "meta_file": str(meta_file),
        "code_dir": str(code_dir),
        "order_file": str(order_file),
        "ribosomal_file": str(ribosomal_file),
        "output_dir": str(output_dir),
        "contract": contract,
        "disk": disk_summary,
        "metadata": metadata_summary,
        "runtime_matches": bool(runtime_rows) and all(row["match"] for row in runtime_rows),
        "sample_order_ready": bool(order_rows) and all(
            int(row["n_h5_candidates"]) == 1 for row in order_rows
        ),
        "code_prepared": bool(args.prepare_code and zip_path.exists()),
        "prepared_code_paths": {
            "preprocessing": str(preprocessing_path) if preprocessing_path else "",
            "annotation": str(annotation_path) if annotation_path else "",
        },
        "prepared_code_sha256": {
            "preprocessing": (
                file_digest(preprocessing_path) if preprocessing_path else ""
            ),
            "annotation": file_digest(annotation_path) if annotation_path else "",
        },
        "ready_to_run_exact_replay": not failures and bool(args.prepare_code),
        "failures": failures,
        "outputs": {
            "runtime_manifest": str(runtime_path),
            "sample_order_manifest": str(order_path),
            "summary": str(summary_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    for path in (runtime_path, order_path, summary_path):
        print(path)

    if args.strict and failures:
        raise SystemExit("Strict exact-replay preflight failed: " + "; ".join(failures))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
