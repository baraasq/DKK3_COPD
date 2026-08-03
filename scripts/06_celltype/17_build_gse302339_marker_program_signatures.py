#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import pickle
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_OBJECT = "output/adata_harmony_annotated_cr8"
DEFAULT_OUTPUT_SIGNATURE = (
    "data/processed/gse292993/"
    "gse302339_marker_program_signatures_raw_logcpm.csv"
)
DEFAULT_SELECTED_CELLS = (
    "results/tables/gse302339_marker_program_selected_cells.csv.gz"
)
DEFAULT_CELL_TYPE_COLUMN = "marker_program_celltype"
DEFAULT_REQUIRED_CELL_TYPES = ["AT1", "AT2", "Fibroblast"]

MARKER_SETS = {
    "AT1": [
        "AGER",
        "PDPN",
        "CAV1",
        "CAV2",
        "CLDN18",
        "HOPX",
        "EMP2",
        "RTKN2",
        "CYP4B1",
        "AQP5",
    ],
    "AT2": [
        "SFTPC",
        "SFTPB",
        "SFTPA1",
        "SFTPA2",
        "ABCA3",
        "NAPSA",
        "SLC34A2",
        "LPCAT1",
        "MUC1",
    ],
    "Fibroblast": [
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "DCN",
        "LUM",
        "PDGFRA",
        "COL6A1",
        "COL6A2",
        "FBLN1",
        "MFAP4",
        "C7",
    ],
    "Endothelial": [
        "PECAM1",
        "VWF",
        "KDR",
        "CLDN5",
        "RAMP2",
        "ACKR1",
        "EMCN",
    ],
    "Smooth muscle": [
        "ACTA2",
        "TAGLN",
        "MYH11",
        "CNN1",
        "MYL9",
        "TPM2",
    ],
    "Airway epithelial": [
        "FOXJ1",
        "PIFO",
        "CAPS",
        "TPPP3",
        "KRT5",
        "KRT15",
        "MUC5B",
        "MUC5AC",
        "SCGB1A1",
        "SCGB3A1",
    ],
    "Immune": [
        "PTPRC",
        "LST1",
        "TYROBP",
        "CD3D",
        "CD3E",
        "MS4A1",
        "NKG7",
        "LYZ",
    ],
}


class RawExpressionView:
    def __init__(self, adata: Any):
        raw = getattr(adata, "raw", None)
        if raw is None:
            raise ValueError("AnnData object has no .raw slot")
        self.X = raw.X
        self.obs = adata.obs
        self.var = raw.var
        self.var_names = raw.var_names
        self.layers = {}
        self.n_obs = adata.n_obs
        self.n_vars = raw.shape[1]

    def __getitem__(self, item):
        obs_index, var_index = item
        return RawExpressionSlice(self, obs_index, var_index)


class RawExpressionSlice:
    def __init__(self, parent: RawExpressionView, obs_index: Any, var_index: Any):
        self.X = parent.X[obs_index, :][:, var_index]
        self.obs = parent.obs.iloc[obs_index].copy()
        self.var = parent.var.iloc[var_index].copy()
        self.var_names = self.var.index
        self.layers = {}
        self.n_obs = int(self.X.shape[0])
        self.n_vars = int(self.X.shape[1])

    def copy(self):
        return self


def load_scrna_audit_module():
    path = Path(__file__).resolve().parent / "00_audit_scrna_reference.py"
    spec = importlib.util.spec_from_file_location("scrna_reference_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load scRNA reference audit helpers.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def load_object(path: Path) -> tuple[Any | None, dict]:
    if not path.exists():
        return None, {"path": str(path), "exists": False, "status": "missing"}
    try:
        if "".join(path.suffixes).casefold().endswith(".h5ad"):
            import anndata as ad

            adata = ad.read_h5ad(path)
            loader = "anndata.read_h5ad"
        else:
            with path.open("rb") as handle:
                adata = pickle.load(handle)
            loader = "pickle.load"
        if not all(hasattr(adata, attr) for attr in ("obs", "X", "var_names")):
            raise TypeError("Loaded object is not AnnData-like")
        return adata, {
            "path": str(path),
            "exists": True,
            "status": "ok",
            "loader": loader,
            "shape": [int(adata.shape[0]), int(adata.shape[1])],
        }
    except Exception as exc:
        return None, {
            "path": str(path),
            "exists": True,
            "status": "load_failed",
            "error": str(exc),
        }


def expression_reference(adata: Any, source: str) -> tuple[Any | None, dict]:
    if source == "raw":
        try:
            ref = RawExpressionView(adata)
            return ref, {
                "requested": source,
                "selected": "raw.X",
                "shape": [int(ref.n_obs), int(ref.n_vars)],
            }
        except Exception as exc:
            return None, {
                "requested": source,
                "selected": "",
                "failure": str(exc),
            }
    return adata, {
        "requested": source,
        "selected": "X",
        "shape": [int(adata.shape[0]), int(adata.shape[1])],
    }


def gene_index(var_names: Any) -> dict[str, int]:
    output = {}
    for index, value in enumerate([str(item) for item in var_names]):
        text = value.strip()
        if text:
            output.setdefault(text.upper(), index)
    return output


def marker_indices(var_names: Any, marker_sets: dict[str, list[str]]) -> dict[str, list[int]]:
    index_by_gene = gene_index(var_names)
    return {
        label: [
            index_by_gene[gene.upper()]
            for gene in genes
            if gene.upper() in index_by_gene
        ]
        for label, genes in marker_sets.items()
    }


def matrix_logcpm_columns(matrix: Any, gene_indices: list[int]):
    import numpy as np

    subset = matrix[:, gene_indices]
    library = np.asarray(matrix.sum(axis=1)).ravel()
    keep = library > 0
    if hasattr(subset, "tocsr"):
        subset = subset.tocsr()
        normalized = subset.copy()
        normalized = normalized.multiply((1_000_000 / np.maximum(library, 1))[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        values = np.asarray(normalized.mean(axis=1)).ravel()
        detected = np.asarray((subset > 0).sum(axis=1)).ravel()
    else:
        subset = np.asarray(subset)
        normalized = subset / np.maximum(library, 1)[:, None] * 1_000_000
        values = np.log1p(normalized).mean(axis=1)
        detected = (subset > 0).sum(axis=1)
    values[~keep] = 0.0
    return values, detected


def zscore(values):
    import numpy as np

    array = np.asarray(values, dtype=float)
    center = float(np.nanmedian(array))
    mad = float(np.nanmedian(np.abs(array - center)))
    scale = 1.4826 * mad if mad > 0 else float(np.nanstd(array))
    if not math.isfinite(scale) or scale == 0:
        scale = 1.0
    return (array - center) / scale


def score_marker_programs(ref: Any, marker_sets: dict[str, list[str]]) -> tuple[dict, dict]:
    import numpy as np

    indices = marker_indices(ref.var_names, marker_sets)
    scores = {}
    detected = {}
    marker_summary = {}
    for label, gene_indices in indices.items():
        present = [
            gene
            for gene in marker_sets[label]
            if gene.upper() in gene_index(ref.var_names)
        ]
        marker_summary[label] = {
            "cell_type": label,
            "n_configured_markers": len(marker_sets[label]),
            "n_present_markers": len(gene_indices),
            "present_markers": ";".join(present),
            "missing_markers": ";".join(
                gene for gene in marker_sets[label] if gene not in present
            ),
        }
        if not gene_indices:
            scores[label] = np.zeros(ref.n_obs)
            detected[label] = np.zeros(ref.n_obs, dtype=int)
            continue
        score, n_detected = matrix_logcpm_columns(ref.X, gene_indices)
        scores[label] = score
        detected[label] = n_detected
    return {
        "scores": scores,
        "z_scores": {label: zscore(value) for label, value in scores.items()},
        "detected": detected,
    }, marker_summary


def assign_marker_labels(
    *,
    score_bundle: dict,
    min_z: float,
    min_margin_z: float,
    min_detected_markers: int,
    max_cells_per_type: int,
    include_cell_types: list[str],
) -> tuple[list[str], list[dict]]:
    import numpy as np

    labels = include_cell_types
    z_matrix = np.vstack([score_bundle["z_scores"][label] for label in labels]).T
    best_indices = np.argmax(z_matrix, axis=1)
    sorted_z = np.sort(z_matrix, axis=1)
    best_z = z_matrix[np.arange(z_matrix.shape[0]), best_indices]
    second_z = sorted_z[:, -2] if len(labels) > 1 else np.zeros_like(best_z)

    assigned = ["Unassigned"] * z_matrix.shape[0]
    rows = []
    for label_index, label in enumerate(labels):
        detected = score_bundle["detected"][label]
        candidates = [
            index
            for index in range(z_matrix.shape[0])
            if best_indices[index] == label_index
            and best_z[index] >= min_z
            and (best_z[index] - second_z[index]) >= min_margin_z
            and int(detected[index]) >= min_detected_markers
        ]
        candidates.sort(key=lambda index: best_z[index], reverse=True)
        if max_cells_per_type > 0:
            candidates = candidates[:max_cells_per_type]
        for index in candidates:
            assigned[index] = label
        rows.append(
            {
                "cell_type": label,
                "n_selected_cells": len(candidates),
                "selection_rule": (
                    f"winner_z>={min_z}; margin_z>={min_margin_z}; "
                    f"detected_markers>={min_detected_markers}; "
                    f"max_cells={max_cells_per_type or 'none'}"
                ),
            }
        )
    return assigned, rows


def write_csv(path: Path, rows: list[dict], preferred: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = sorted({key for row in rows for key in row})
    ordered = [field for field in (preferred or []) if field in fields]
    ordered.extend(field for field in fields if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_selected_cells(
    path: Path,
    *,
    obs_names: list[str],
    obs: Any,
    assigned: list[str],
    score_bundle: dict,
    include_cell_types: list[str],
) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        column
        for column in ("batch", "sample", "patient", "donor", "lobe_emphysema_simple")
        if column in obs.columns
    ]
    n_written = 0
    with gzip.open(path, "wt", newline="", encoding="utf-8") as handle:
        fieldnames = [
            "cell_barcode",
            "marker_program_celltype",
            *columns,
            *[f"score_{label}" for label in include_cell_types],
            *[f"z_{label}" for label in include_cell_types],
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for index, label in enumerate(assigned):
            if label == "Unassigned":
                continue
            row = {
                "cell_barcode": obs_names[index],
                "marker_program_celltype": label,
            }
            for column in columns:
                row[column] = str(obs.iloc[index][column])
            for cell_type in include_cell_types:
                row[f"score_{cell_type}"] = float(score_bundle["scores"][cell_type][index])
                row[f"z_{cell_type}"] = float(score_bundle["z_scores"][cell_type][index])
            writer.writerow(row)
            n_written += 1
    return n_written


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build GeoMx-compatible signatures from GSE302339 using transparent "
            "marker-program cell selection instead of author cluster-ID replay."
        )
    )
    parser.add_argument("--input-object", default=DEFAULT_INPUT_OBJECT)
    parser.add_argument("--expression-source", choices=["raw", "X"], default="raw")
    parser.add_argument("--cell-type-column", default=DEFAULT_CELL_TYPE_COLUMN)
    parser.add_argument("--min-z", type=float, default=1.5)
    parser.add_argument("--min-margin-z", type=float, default=0.25)
    parser.add_argument("--min-detected-markers", type=int, default=2)
    parser.add_argument("--max-cells-per-cell-type", type=int, default=5000)
    parser.add_argument(
        "--include-cell-type",
        action="append",
        dest="include_cell_types",
        help="Marker program to include. Defaults to all built-in marker sets.",
    )
    parser.add_argument(
        "--required-cell-type",
        action="append",
        dest="required_cell_types",
        help="Cell type that must be present in the emitted signature.",
    )
    parser.add_argument("--signature-output", default=DEFAULT_OUTPUT_SIGNATURE)
    parser.add_argument("--selected-cells-output", default=DEFAULT_SELECTED_CELLS)
    parser.add_argument("--min-cells-per-cell-type", type=int, default=100)
    parser.add_argument("--min-overlap-genes", type=int, default=500)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    results = ensure_results_dirs(config)
    ref_helpers = load_scrna_audit_module()
    table_dir = results["tables"]
    summary_path = results["meta"] / "gse302339_marker_program_signature_summary.json"
    marker_path = table_dir / "gse302339_marker_program_marker_manifest.csv"
    selection_path = table_dir / "gse302339_marker_program_selection_summary.csv"
    cell_counts_path = table_dir / "gse302339_marker_program_celltype_counts.csv"
    signature_manifest_path = table_dir / "gse302339_marker_program_signature_manifest.csv"

    include_cell_types = args.include_cell_types or list(MARKER_SETS)
    unknown_marker_sets = [label for label in include_cell_types if label not in MARKER_SETS]
    if unknown_marker_sets:
        raise SystemExit(
            "Unknown marker-program cell type(s): " + ", ".join(unknown_marker_sets)
        )
    marker_sets = {label: MARKER_SETS[label] for label in include_cell_types}
    required_cell_types = args.required_cell_types or DEFAULT_REQUIRED_CELL_TYPES

    adata, load_summary = load_object(project_path(args.input_object))
    failures = []
    if adata is None:
        failures.append(f"Input object not loadable: {load_summary.get('status')}")
        summary = {"input_object": load_summary, "ready_for_deconvolution": False, "failures": failures}
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 2 if args.strict else 1

    ref, expression_summary = expression_reference(adata, args.expression_source)
    if ref is None:
        failures.append(f"Expression source not available: {expression_summary.get('failure')}")
        summary = {
            "input_object": load_summary,
            "expression_source": expression_summary,
            "ready_for_deconvolution": False,
            "failures": failures,
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 2 if args.strict else 1

    score_bundle, marker_summary = score_marker_programs(ref, marker_sets)
    assigned, selection_rows = assign_marker_labels(
        score_bundle=score_bundle,
        min_z=args.min_z,
        min_margin_z=args.min_margin_z,
        min_detected_markers=args.min_detected_markers,
        max_cells_per_type=args.max_cells_per_cell_type,
        include_cell_types=include_cell_types,
    )
    ref.obs[args.cell_type_column] = assigned
    selected_indices = [
        index for index, label in enumerate(assigned) if label != "Unassigned"
    ]
    selected_counts = Counter(label for label in assigned if label != "Unassigned")
    missing_required = [
        label
        for label in required_cell_types
        if selected_counts.get(label, 0) < args.min_cells_per_cell_type
    ]
    if missing_required:
        failures.append(
            "Required marker programs below min cell threshold: "
            + ", ".join(missing_required)
        )

    selected_cells_path = project_path(args.selected_cells_output)
    n_selected_rows = write_selected_cells(
        selected_cells_path,
        obs_names=[str(value) for value in ref.obs.index],
        obs=ref.obs,
        assigned=assigned,
        score_bundle=score_bundle,
        include_cell_types=include_cell_types,
    )

    geomx_feature_manifest = results["tables"] / "gse292993_geomx_feature_manifest.csv"
    geomx_genes = ref_helpers.read_geomx_genes(geomx_feature_manifest)
    overlap = ref_helpers.select_gene_overlap(ref, geomx_genes)
    if overlap["n_overlap"] < args.min_overlap_genes:
        failures.append(
            f"Too few GeoMx/scRNA overlap genes: {overlap['n_overlap']} "
            f"< {args.min_overlap_genes}"
        )

    signature_output = project_path(args.signature_output)
    signature_rows = []
    cell_count_rows = []
    if selected_indices and (not failures or overlap["n_overlap"] >= args.min_overlap_genes):
        selected_ref = ref[selected_indices, :].copy()
        signature_rows, cell_count_rows = ref_helpers.write_signature_matrix(
            adata=selected_ref,
            cell_type_column=args.cell_type_column,
            layer=None,
            overlap=overlap,
            output_path=signature_output,
            min_cells=args.min_cells_per_cell_type,
        )

    signature_cell_types = [row["cell_type"] for row in signature_rows]
    missing_emitted = [
        label for label in required_cell_types if label not in set(signature_cell_types)
    ]
    if missing_emitted:
        failures.append(
            "Required marker programs missing from emitted signature: "
            + ", ".join(missing_emitted)
        )

    marker_rows = list(marker_summary.values())
    for row in marker_rows:
        if int(row["n_present_markers"]) == 0:
            failures.append(f"No marker genes present for {row['cell_type']}")

    write_csv(
        marker_path,
        marker_rows,
        preferred=[
            "cell_type",
            "n_configured_markers",
            "n_present_markers",
            "present_markers",
            "missing_markers",
        ],
    )
    write_csv(
        selection_path,
        selection_rows,
        preferred=["cell_type", "n_selected_cells", "selection_rule"],
    )
    write_csv(
        cell_counts_path,
        cell_count_rows,
        preferred=["cell_type", "n_cells", "included_in_signature"],
    )
    write_csv(
        signature_manifest_path,
        signature_rows,
        preferred=["cell_type", "n_cells", "n_cells_used", "n_genes"],
    )

    summary = {
        "strategy": "marker-program high-confidence cell selection",
        "input_object": load_summary,
        "expression_source": expression_summary,
        "cell_type_column": args.cell_type_column,
        "selection_parameters": {
            "min_z": args.min_z,
            "min_margin_z": args.min_margin_z,
            "min_detected_markers": args.min_detected_markers,
            "max_cells_per_cell_type": args.max_cells_per_cell_type,
            "min_cells_per_cell_type": args.min_cells_per_cell_type,
        },
        "include_cell_types": include_cell_types,
        "required_cell_types": required_cell_types,
        "selected_cell_counts": dict(sorted(selected_counts.items())),
        "n_selected_cells_written": n_selected_rows,
        "geomx_feature_manifest": str(geomx_feature_manifest),
        "gene_overlap": {
            "source": overlap["source"],
            "n_overlap": overlap["n_overlap"],
            "first_overlap_genes": overlap["common_genes"][:20],
            "primary_gene_DKK3_in_overlap": "DKK3" in set(overlap["common_genes"]),
        },
        "signature_output": str(signature_output),
        "n_signature_cell_types": len(signature_rows),
        "signature_cell_types": signature_cell_types,
        "outputs": {
            "summary": str(summary_path),
            "marker_manifest": str(marker_path),
            "selection_summary": str(selection_path),
            "celltype_counts": str(cell_counts_path),
            "signature_manifest": str(signature_manifest_path),
            "selected_cells": str(selected_cells_path),
            "signature": str(signature_output),
        },
        "ready_for_deconvolution": not failures,
        "failures": failures,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    for path in (
        summary_path,
        marker_path,
        selection_path,
        cell_counts_path,
        signature_manifest_path,
        selected_cells_path,
        signature_output,
    ):
        print(path)
    if args.strict and failures:
        print("Strict marker-program signature build failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
