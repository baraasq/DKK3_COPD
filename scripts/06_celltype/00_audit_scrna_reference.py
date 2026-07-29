#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_ACCESSION = "GSE302339"
DEFAULT_CELL_TYPE_CANDIDATES = [
    "cell_type",
    "cell type",
    "celltype",
    "cell_type_major",
    "cell_type_fine",
    "cell_subtype",
    "cell subtype",
    "annotation",
    "annotations",
    "subclass",
    "class",
    "label",
]
DEFAULT_CONDITION_CANDIDATES = [
    "diagnosis",
    "condition",
    "disease",
    "Disease",
    "group",
    "phenotype",
]
GENE_SYMBOL_CANDIDATES = [
    "gene_symbols",
    "gene_symbol",
    "gene_ids",
    "feature_name",
    "features",
    "name",
    "symbol",
]
DEFAULT_MAX_CODES_PER_GEOMX_TARGET = 20


def raise_csv_field_limit() -> None:
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit = int(limit / 10)


def write_csv(path: Path, rows: list[dict], preferred: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    ordered = [field for field in (preferred or []) if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def resolve_column(columns: list[str], candidates: list[str]) -> str | None:
    exact = {column: column for column in columns}
    folded = {column.casefold(): column for column in columns}
    for candidate in candidates:
        if candidate in exact:
            return exact[candidate]
        if candidate.casefold() in folded:
            return folded[candidate.casefold()]
    return None


def file_size(path: Path) -> int | None:
    try:
        return path.stat().st_size
    except OSError:
        return None


def discover_reference_files(config: dict, reference_h5ad: Path | None) -> list[dict]:
    roots = [
        configured_path(config, "geomx_processed_dir"),
        project_path("data/external/scrna_reference"),
        project_path("data/external"),
        project_path(f"data/raw/downloads/geo/{DEFAULT_ACCESSION}"),
        project_path(f"data/raw/{DEFAULT_ACCESSION}"),
        project_path(f"data/processed/{DEFAULT_ACCESSION}"),
    ]
    rows = []
    seen = set()
    if reference_h5ad:
        roots.insert(0, reference_h5ad.parent)
    for root in roots:
        if not root.exists():
            continue
        candidates = [
            *root.rglob("*.h5ad"),
            *root.rglob("*.h5"),
            *root.rglob("*.hdf5"),
            *root.rglob("filtered_feature_bc_matrix.h5"),
            *root.rglob("raw_feature_bc_matrix.h5"),
            *root.rglob("matrix.mtx"),
            *root.rglob("matrix.mtx.gz"),
            *root.rglob("*.rds"),
            *root.rglob("*.RDS"),
            *root.rglob("*.rda"),
            *root.rglob("*.RData"),
            *root.rglob("*.csv"),
            *root.rglob("*.csv.gz"),
            *root.rglob("*.tsv"),
            *root.rglob("*.tsv.gz"),
            *root.rglob("*.tar"),
            *root.rglob("*.tar.gz"),
            *root.rglob("*.tgz"),
            *root.rglob("*.zip"),
        ]
        for path in sorted(candidates):
            if path in seen:
                continue
            seen.add(path)
            rows.append(
                {
                    "path": str(path),
                    "name": path.name,
                    "suffix": "".join(path.suffixes),
                    "size_bytes": file_size(path),
                    "kind": guess_reference_kind(path),
                }
            )
    return rows


def guess_reference_kind(path: Path) -> str:
    suffixes = "".join(path.suffixes).casefold()
    if suffixes.endswith(".h5ad"):
        return "annotated_h5ad_candidate"
    if path.name in {"filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5"}:
        return "cellranger_h5_counts_only"
    if suffixes.endswith(".h5") or suffixes.endswith(".hdf5"):
        return "hdf5_candidate"
    if path.name in {"matrix.mtx", "matrix.mtx.gz"}:
        return "cellranger_mtx_counts_only"
    if suffixes.endswith(".rds") or suffixes.endswith(".rda") or suffixes.endswith(".rdata"):
        return "r_object_candidate"
    if suffixes.endswith(".csv") or suffixes.endswith(".csv.gz") or suffixes.endswith(".tsv") or suffixes.endswith(".tsv.gz"):
        return "metadata_or_matrix_table_candidate"
    if suffixes.endswith(".tar") or suffixes.endswith(".tar.gz") or suffixes.endswith(".tgz") or suffixes.endswith(".zip"):
        return "archive_candidate"
    return "unknown"


def read_geomx_genes(feature_manifest: Path) -> list[str]:
    if not feature_manifest.exists():
        return []
    raise_csv_field_limit()
    with feature_manifest.open(newline="", encoding="utf-8-sig") as handle:
        return [
            row["target"]
            for row in csv.DictReader(handle)
            if row.get("target") and row.get("is_control_feature") != "True"
            and int(row.get("n_codes") or 0) <= DEFAULT_MAX_CODES_PER_GEOMX_TARGET
        ]


def candidate_gene_vectors(adata: Any) -> list[tuple[str, list[str]]]:
    output = [("var_names", [str(value) for value in adata.var_names])]
    for column in GENE_SYMBOL_CANDIDATES:
        if column in adata.var.columns:
            values = [str(value) for value in adata.var[column].tolist()]
            output.append((f"var.{column}", values))
    return output


def select_gene_overlap(adata: Any, geomx_genes: list[str]) -> dict:
    geomx_by_upper = {gene.upper(): gene for gene in geomx_genes}
    best = {
        "source": "var_names",
        "n_overlap": 0,
        "common_genes": [],
        "adata_gene_indices": [],
        "adata_gene_names": [],
    }
    for source, values in candidate_gene_vectors(adata):
        indices_by_upper: dict[str, int] = {}
        for index, value in enumerate(values):
            text = str(value).strip()
            if text:
                indices_by_upper.setdefault(text.upper(), index)
        common_upper = sorted(set(indices_by_upper) & set(geomx_by_upper))
        if len(common_upper) > best["n_overlap"]:
            best = {
                "source": source,
                "n_overlap": len(common_upper),
                "common_genes": [geomx_by_upper[item] for item in common_upper],
                "adata_gene_indices": [indices_by_upper[item] for item in common_upper],
                "adata_gene_names": [values[indices_by_upper[item]] for item in common_upper],
            }
    return best


def dense_or_sparse_logcpm_mean(matrix, cell_indices: list[int], gene_indices: list[int]):
    import numpy as np

    subset = matrix[cell_indices, :][:, gene_indices]
    if hasattr(subset, "tocsr"):
        subset = subset.tocsr()
        library = np.asarray(subset.sum(axis=1)).ravel()
        keep = library > 0
        if not keep.any():
            return None, 0
        subset = subset[keep, :]
        scale = 1_000_000 / library[keep]
        normalized = subset.multiply(scale[:, None]).tocsr()
        normalized.data = np.log1p(normalized.data)
        return np.asarray(normalized.mean(axis=0)).ravel(), int(keep.sum())

    subset = np.asarray(subset)
    library = subset.sum(axis=1)
    keep = library > 0
    if not keep.any():
        return None, 0
    normalized = subset[keep, :] / library[keep, None] * 1_000_000
    return np.log1p(normalized).mean(axis=0), int(keep.sum())


def write_signature_matrix(
    *,
    adata: Any,
    cell_type_column: str,
    layer: str | None,
    overlap: dict,
    output_path: Path,
    min_cells: int,
) -> tuple[list[dict], list[dict]]:
    import numpy as np

    matrix = adata.layers[layer] if layer and layer in adata.layers else adata.X
    labels = [str(value) for value in adata.obs[cell_type_column].tolist()]
    label_counts = Counter(labels)
    selected_labels = sorted(label for label, count in label_counts.items() if count >= min_cells)
    genes = overlap["common_genes"]
    gene_indices = overlap["adata_gene_indices"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    signature_rows = []
    cell_count_rows = []
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["cell_type", *genes])
        for label in selected_labels:
            cell_indices = [index for index, value in enumerate(labels) if value == label]
            means, n_cells_used = dense_or_sparse_logcpm_mean(
                matrix, cell_indices, gene_indices
            )
            if means is None:
                continue
            writer.writerow([label, *[float(value) for value in np.asarray(means).ravel()]])
            signature_rows.append(
                {
                    "cell_type": label,
                    "n_cells": len(cell_indices),
                    "n_cells_used": n_cells_used,
                    "n_genes": len(genes),
                }
            )
        for label, count in sorted(label_counts.items()):
            cell_count_rows.append(
                {
                    "cell_type": label,
                    "n_cells": count,
                    "included_in_signature": str(label in selected_labels),
                }
            )
    return signature_rows, cell_count_rows


def inspect_h5ad(
    *,
    path: Path,
    geomx_genes: list[str],
    cell_type_column_arg: str | None,
    condition_column_arg: str | None,
    layer_arg: str | None,
    signature_output: Path,
    min_cells: int,
    write_signatures: bool,
) -> dict:
    try:
        import anndata as ad
    except ImportError as exc:
        return {
            "path": str(path),
            "status": "anndata_import_failed",
            "error": str(exc),
            "ready_for_deconvolution": False,
        }

    if not path.exists():
        return {
            "path": str(path),
            "status": "missing",
            "ready_for_deconvolution": False,
        }

    adata = ad.read_h5ad(path)
    obs_columns = [str(column) for column in adata.obs.columns]
    cell_type_column = (
        cell_type_column_arg
        if cell_type_column_arg and cell_type_column_arg in adata.obs.columns
        else resolve_column(obs_columns, DEFAULT_CELL_TYPE_CANDIDATES)
    )
    condition_column = (
        condition_column_arg
        if condition_column_arg and condition_column_arg in adata.obs.columns
        else resolve_column(obs_columns, DEFAULT_CONDITION_CANDIDATES)
    )
    layer = layer_arg if layer_arg in adata.layers else None
    overlap = select_gene_overlap(adata, geomx_genes)

    cell_counts = []
    signature_rows = []
    if write_signatures and cell_type_column and overlap["n_overlap"] > 0:
        signature_rows, cell_counts = write_signature_matrix(
            adata=adata,
            cell_type_column=cell_type_column,
            layer=layer,
            overlap=overlap,
            output_path=signature_output,
            min_cells=min_cells,
        )
    elif cell_type_column:
        counts = Counter(str(value) for value in adata.obs[cell_type_column].tolist())
        cell_counts = [
            {
                "cell_type": label,
                "n_cells": count,
                "included_in_signature": "False",
            }
            for label, count in sorted(counts.items())
        ]

    return {
        "path": str(path),
        "status": "ok",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "obs_columns": obs_columns,
        "var_columns": [str(column) for column in adata.var.columns],
        "layers": [str(layer_name) for layer_name in adata.layers.keys()],
        "selected_cell_type_column": cell_type_column,
        "selected_condition_column": condition_column,
        "selected_expression_layer": layer or "X",
        "geomx_gene_overlap": {
            "source": overlap["source"],
            "n_geomx_genes": len(geomx_genes),
            "n_overlap": overlap["n_overlap"],
            "first_overlap_genes": overlap["common_genes"][:20],
        },
        "signature_output": str(signature_output) if signature_rows else None,
        "n_signature_cell_types": len(signature_rows),
        "ready_for_deconvolution": bool(
            cell_type_column and overlap["n_overlap"] >= 100 and len(signature_rows) >= 2
        ),
        "cell_counts": cell_counts,
        "signature_rows": signature_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit the scRNA-seq reference and optionally build cell-type signatures for GeoMx deconvolution."
    )
    parser.add_argument("--reference-h5ad")
    parser.add_argument("--cell-type-column")
    parser.add_argument("--condition-column")
    parser.add_argument("--layer")
    parser.add_argument("--min-cells-per-cell-type", type=int, default=25)
    parser.add_argument("--write-signatures", action="store_true", default=True)
    parser.add_argument("--no-write-signatures", action="store_false", dest="write_signatures")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    processed_dir = configured_path(config, "geomx_processed_dir")
    processed_dir.mkdir(parents=True, exist_ok=True)

    reference_h5ad = (
        Path(args.reference_h5ad).expanduser()
        if args.reference_h5ad
        else configured_path(config, "scrna_h5ad")
    )
    geomx_genes = read_geomx_genes(table_dir / "gse292993_geomx_feature_manifest.csv")
    discovered = discover_reference_files(config, reference_h5ad)
    signature_output = processed_dir / "gse302339_scrna_reference_signatures_logcpm.csv"
    h5ad_summary = inspect_h5ad(
        path=reference_h5ad,
        geomx_genes=geomx_genes,
        cell_type_column_arg=args.cell_type_column,
        condition_column_arg=args.condition_column,
        layer_arg=args.layer or config.get("scrna", {}).get("counts_layer"),
        signature_output=signature_output,
        min_cells=args.min_cells_per_cell_type,
        write_signatures=args.write_signatures,
    )

    summary = {
        "accession": DEFAULT_ACCESSION,
        "reference_h5ad": str(reference_h5ad),
        "geomx_feature_manifest_present": bool(geomx_genes),
        "n_geomx_genes": len(geomx_genes),
        "discovered_reference_files": discovered,
        "h5ad": {
            key: value
            for key, value in h5ad_summary.items()
            if key not in {"cell_counts", "signature_rows"}
        },
        "ready_for_deconvolution": bool(h5ad_summary.get("ready_for_deconvolution")),
    }
    summary_path = meta_dir / "gse302339_scrna_reference_audit_summary.json"
    cell_counts_path = table_dir / "gse302339_scrna_celltype_counts.csv"
    signature_manifest_path = table_dir / "gse302339_scrna_signature_manifest.csv"
    discovery_path = table_dir / "gse302339_scrna_reference_file_manifest.csv"

    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(discovery_path, discovered, preferred=["kind", "path", "size_bytes"])
    write_csv(
        cell_counts_path,
        h5ad_summary.get("cell_counts", []),
        preferred=["cell_type", "n_cells", "included_in_signature"],
    )
    write_csv(
        signature_manifest_path,
        h5ad_summary.get("signature_rows", []),
        preferred=["cell_type", "n_cells", "n_cells_used", "n_genes"],
    )

    print(json.dumps(summary, indent=2))
    print()
    for path in (summary_path, discovery_path, cell_counts_path, signature_manifest_path):
        print(path)
    if h5ad_summary.get("signature_output"):
        print(h5ad_summary["signature_output"])

    failures = []
    if not geomx_genes:
        failures.append("GeoMx feature manifest missing; run the WTA matrix builder first")
    if h5ad_summary.get("status") != "ok":
        failures.append(f"H5AD reference not ready: {h5ad_summary.get('status')}")
    if not h5ad_summary.get("selected_cell_type_column"):
        failures.append("No cell type annotation column resolved in scRNA reference")
    if h5ad_summary.get("geomx_gene_overlap", {}).get("n_overlap", 0) < 100:
        failures.append("Fewer than 100 overlapping genes between scRNA reference and GeoMx WTA matrix")
    if args.strict and failures:
        print("Strict scRNA reference audit failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
