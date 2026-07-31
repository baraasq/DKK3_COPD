#!/usr/bin/env python3
from __future__ import annotations

import argparse
import importlib.util
import json
import pickle
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_OBJECT = project_path("output/parenchyma_harmony_annotated_cr8")
DEFAULT_CELL_TYPE_COLUMN = "parenchyma_celltype_level1"
DEFAULT_REQUIRED_CELL_TYPES = ["AT1", "AT2", "Fibroblast"]
DEFAULT_REFERENCE_H5AD = project_path(
    "data/external/scrna_reference/gse302339_author_parenchyma_celltype_level1.h5ad"
)
DEFAULT_CLUSTER_MAP = project_path(
    "results/tables/gse302339_author_celltype_cluster_maps.csv"
)


def load_scrna_audit_module():
    path = Path(__file__).resolve().parent / "00_audit_scrna_reference.py"
    spec = importlib.util.spec_from_file_location("scrna_reference_audit", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("Could not load scRNA reference audit helpers.")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def valid_label(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return bool(text) and text.casefold() not in {"nan", "none", "na", "null", "<na>"}


def has_anndata_interface(value: object) -> bool:
    return all(hasattr(value, attr) for attr in ("obs", "var", "X")) and hasattr(value, "n_obs")


def coerce_anndata_like(value: object) -> Any:
    if has_anndata_interface(value):
        return value
    if isinstance(value, dict):
        for key in ("adata", "ann_data", "anndata", "rna"):
            candidate = value.get(key)
            if has_anndata_interface(candidate):
                return candidate
    raise TypeError(
        "Loaded object is not AnnData-like. Expected an object with obs, var, X, and n_obs."
    )


def load_author_object(path: Path) -> tuple[Any | None, dict]:
    if not path.exists():
        return None, {
            "path": str(path),
            "exists": False,
            "status": "missing",
        }
    try:
        if "".join(path.suffixes).casefold().endswith(".h5ad"):
            import anndata as ad

            adata = ad.read_h5ad(path)
            loader = "anndata.read_h5ad"
        else:
            with path.open("rb") as handle:
                adata = coerce_anndata_like(pickle.load(handle))
            loader = "pickle.load"
        return adata, {
            "path": str(path),
            "exists": True,
            "status": "ok",
            "loader": loader,
        }
    except Exception as exc:
        return None, {
            "path": str(path),
            "exists": True,
            "status": "load_failed",
            "error": str(exc),
        }


def read_csv(path: Path) -> list[dict]:
    import csv

    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def obs_columns(adata: Any) -> list[str]:
    return [str(column) for column in adata.obs.columns]


def available_layers(adata: Any) -> list[str]:
    layers = getattr(adata, "layers", {})
    try:
        return sorted(str(key) for key in layers.keys())
    except Exception:
        return []


def selected_layer(adata: Any, requested_layer: str | None) -> str | None:
    if not requested_layer:
        return None
    return requested_layer if requested_layer in available_layers(adata) else None


def filter_to_valid_labels(
    adata: Any,
    *,
    cell_type_column: str,
    include_cell_types: list[str] | None,
) -> tuple[Any, dict]:
    labels = list(adata.obs[cell_type_column].tolist())
    include_set = set(include_cell_types or [])
    mask = [
        valid_label(label) and (not include_set or str(label) in include_set)
        for label in labels
    ]
    n_keep = sum(mask)
    summary = {
        "n_cells_before_filter": int(adata.n_obs),
        "n_cells_after_filter": int(n_keep),
        "n_cells_dropped_invalid_or_unselected_label": int(adata.n_obs) - int(n_keep),
        "include_cell_types": sorted(include_set),
    }
    if n_keep == len(mask):
        return adata, summary
    return adata[mask, :].copy(), summary


def write_h5ad(adata: Any, path: Path) -> dict:
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        if hasattr(adata, "write_h5ad"):
            adata.write_h5ad(path)
            method = "write_h5ad"
        elif hasattr(adata, "write"):
            adata.write(path)
            method = "write"
        else:
            return {
                "path": str(path),
                "written": False,
                "status": "no_write_method",
            }
        return {
            "path": str(path),
            "written": True,
            "status": "ok",
            "method": method,
            "size_bytes": path.stat().st_size if path.exists() else None,
        }
    except Exception as exc:
        return {
            "path": str(path),
            "written": False,
            "status": "write_failed",
            "error": str(exc),
        }


def label_counts(adata: Any, cell_type_column: str) -> list[dict]:
    counts = Counter(str(value) for value in adata.obs[cell_type_column].tolist())
    return [
        {
            "cell_type": label,
            "n_cells": int(count),
        }
        for label, count in sorted(counts.items())
    ]


def cluster_label_map(
    cluster_map_path: Path,
    *,
    assigned_obs_column: str,
) -> dict[str, str]:
    rows = read_csv(cluster_map_path)
    mapping: dict[str, str] = {}
    for row in rows:
        if row.get("assigned_obs_column") != assigned_obs_column:
            continue
        label = str(row.get("celltype_label", "")).strip()
        for cluster in str(row.get("cluster_ids", "")).split(";"):
            cluster = cluster.strip()
            if cluster and label:
                mapping[cluster] = label
    return mapping


def apply_cluster_label_map(
    adata: Any,
    *,
    cluster_column: str,
    output_column: str,
    mapping: dict[str, str],
) -> dict:
    if cluster_column not in adata.obs.columns:
        return {
            "applied": False,
            "reason": f"cluster_column_missing:{cluster_column}",
            "cluster_column": cluster_column,
            "output_column": output_column,
            "n_cluster_labels": len(mapping),
        }
    labels = adata.obs[cluster_column].astype(str).map(mapping)
    n_assigned = int(labels.notna().sum())
    adata.obs[output_column] = labels
    return {
        "applied": True,
        "reason": "mapped_clusters",
        "cluster_column": cluster_column,
        "output_column": output_column,
        "n_cluster_labels": len(mapping),
        "n_cells_assigned": n_assigned,
        "n_cells_unassigned": int(len(labels) - n_assigned),
        "assigned_labels": sorted(str(value) for value in labels.dropna().unique()),
    }


def build_author_signatures(
    *,
    adata: Any,
    ref_helpers: Any,
    geomx_feature_manifest: Path,
    cell_type_column: str,
    layer: str | None,
    signature_output: Path,
    min_cells: int,
) -> dict:
    geomx_genes = ref_helpers.read_geomx_genes(geomx_feature_manifest)
    overlap = ref_helpers.select_gene_overlap(adata, geomx_genes)
    signature_rows, cell_count_rows = ref_helpers.write_signature_matrix(
        adata=adata,
        cell_type_column=cell_type_column,
        layer=layer,
        overlap=overlap,
        output_path=signature_output,
        min_cells=min_cells,
    )
    return {
        "geomx_feature_manifest": str(geomx_feature_manifest),
        "n_geomx_genes": len(geomx_genes),
        "gene_overlap": {
            "source": overlap["source"],
            "n_overlap": overlap["n_overlap"],
            "first_overlap_genes": overlap["common_genes"][:20],
            "primary_gene_DKK3_in_overlap": "DKK3" in set(overlap["common_genes"]),
        },
        "signature_output": str(signature_output),
        "n_signature_cell_types": len(signature_rows),
        "signature_cell_types": [row["cell_type"] for row in signature_rows],
        "signature_rows": signature_rows,
        "cell_count_rows": cell_count_rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Build GeoMx-compatible cell-type signature matrix from reconstructed "
            "GSE302339 author AnnData pickle/H5AD objects."
        )
    )
    parser.add_argument("--input-object", default=str(DEFAULT_INPUT_OBJECT))
    parser.add_argument("--cell-type-column", default=DEFAULT_CELL_TYPE_COLUMN)
    parser.add_argument("--layer", default="counts")
    parser.add_argument("--min-cells-per-cell-type", type=int, default=25)
    parser.add_argument("--min-overlap-genes", type=int, default=500)
    parser.add_argument(
        "--include-cell-type",
        action="append",
        dest="include_cell_types",
        help="Optional cell type to keep. Can be passed multiple times.",
    )
    parser.add_argument(
        "--required-cell-type",
        action="append",
        dest="required_cell_types",
        help=(
            "Cell type that must be present in the emitted signature. "
            "Can be passed multiple times. Defaults to AT1, AT2, Fibroblast."
        ),
    )
    parser.add_argument("--signature-output")
    parser.add_argument("--h5ad-output", default=str(DEFAULT_REFERENCE_H5AD))
    parser.add_argument("--cluster-map", default=str(DEFAULT_CLUSTER_MAP))
    parser.add_argument("--cluster-column", default="leiden")
    parser.add_argument(
        "--rebuild-labels-from-cluster-map",
        action="store_true",
        help=(
            "Recreate --cell-type-column from the extracted author cluster map "
            "before building signatures. Use only after confirming the cluster "
            "IDs in the object match the extracted map."
        ),
    )
    parser.add_argument("--no-write-h5ad", action="store_true")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    results = ensure_results_dirs(config)
    processed_dir = configured_path(config, "geomx_processed_dir")
    processed_dir.mkdir(parents=True, exist_ok=True)
    input_object = project_path(args.input_object)
    cluster_map_path = project_path(args.cluster_map)
    signature_output = (
        project_path(args.signature_output)
        if args.signature_output
        else processed_dir / "gse302339_author_parenchyma_signatures_logcpm.csv"
    )
    h5ad_output = project_path(args.h5ad_output)
    summary_path = results["meta"] / "gse302339_author_signature_reference_summary.json"
    cell_counts_path = results["tables"] / "gse302339_author_signature_celltype_counts.csv"
    signature_manifest_path = results["tables"] / "gse302339_author_signature_manifest.csv"

    adata, load_summary = load_author_object(input_object)
    failures = []
    if adata is None:
        failures.append(f"Author object not loadable: {load_summary.get('status')}")
        summary = {
            "input_object": load_summary,
            "ready_for_deconvolution": False,
            "failures": failures,
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        print()
        print(summary_path)
        if args.strict:
            print(
                "Strict author signature build failed: " + "; ".join(failures),
                file=sys.stderr,
            )
            return 2
        return 1

    required_cell_types = args.required_cell_types or DEFAULT_REQUIRED_CELL_TYPES
    label_rebuild_summary = {
        "requested": args.rebuild_labels_from_cluster_map,
        "applied": False,
        "cluster_map": str(cluster_map_path),
    }
    if args.rebuild_labels_from_cluster_map:
        if not cluster_map_path.exists():
            label_rebuild_summary["reason"] = "cluster_map_missing"
        else:
            mapping = cluster_label_map(
                cluster_map_path,
                assigned_obs_column=args.cell_type_column,
            )
            label_rebuild_summary = {
                "requested": True,
                "cluster_map": str(cluster_map_path),
                **apply_cluster_label_map(
                    adata,
                    cluster_column=args.cluster_column,
                    output_column=args.cell_type_column,
                    mapping=mapping,
                ),
            }

    if args.cell_type_column not in adata.obs.columns:
        failures.append(f"Cell type column not found: {args.cell_type_column}")

    chosen_layer = selected_layer(adata, args.layer)
    if args.layer and chosen_layer is None:
        layer_warning = f"Requested layer '{args.layer}' not found; using X"
    else:
        layer_warning = ""

    filtered = adata
    filter_summary = {}
    if not failures:
        filtered, filter_summary = filter_to_valid_labels(
            adata,
            cell_type_column=args.cell_type_column,
            include_cell_types=args.include_cell_types,
        )

    h5ad_summary = {
        "path": str(h5ad_output),
        "written": False,
        "status": "skipped",
    }
    if not failures and not args.no_write_h5ad:
        h5ad_summary = write_h5ad(filtered, h5ad_output)

    signature_summary = {}
    if not failures:
        ref_helpers = load_scrna_audit_module()
        signature_summary = build_author_signatures(
            adata=filtered,
            ref_helpers=ref_helpers,
            geomx_feature_manifest=results["tables"] / "gse292993_geomx_feature_manifest.csv",
            cell_type_column=args.cell_type_column,
            layer=chosen_layer,
            signature_output=signature_output,
            min_cells=args.min_cells_per_cell_type,
        )

    signature_cell_types = set(signature_summary.get("signature_cell_types", []))
    missing_required = [
        label for label in required_cell_types if label not in signature_cell_types
    ]
    if not failures:
        if signature_summary.get("gene_overlap", {}).get("n_overlap", 0) < args.min_overlap_genes:
            failures.append(
                "Fewer than "
                f"{args.min_overlap_genes} overlapping genes between author reference and GeoMx WTA"
            )
        if len(signature_cell_types) < 2:
            failures.append("Fewer than two cell types emitted in signature matrix")
        if missing_required:
            failures.append(
                "Required cell types missing from emitted signature: "
                + ", ".join(missing_required)
            )

    summary = {
        "dataset": "GSE302339",
        "input_object": load_summary,
        "n_cells": int(getattr(filtered, "n_obs", 0)),
        "n_genes": int(getattr(filtered, "n_vars", 0)),
        "obs_columns": obs_columns(filtered),
        "cell_type_column": args.cell_type_column,
        "available_layers": available_layers(filtered),
        "requested_layer": args.layer,
        "selected_expression_source": f"layer:{chosen_layer}" if chosen_layer else "X",
        "layer_warning": layer_warning,
        "filter_summary": filter_summary,
        "label_rebuild_summary": label_rebuild_summary,
        "required_cell_types": required_cell_types,
        "missing_required_cell_types": missing_required,
        "h5ad_output": h5ad_summary,
        **signature_summary,
        "ready_for_deconvolution": not failures,
        "failures": failures,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    ref_helpers = load_scrna_audit_module()
    ref_helpers.write_csv(
        cell_counts_path,
        signature_summary.get("cell_count_rows", label_counts(filtered, args.cell_type_column))
        if not failures or args.cell_type_column in getattr(filtered, "obs", {})
        else [],
        preferred=["cell_type", "n_cells", "included_in_signature"],
    )
    ref_helpers.write_csv(
        signature_manifest_path,
        signature_summary.get("signature_rows", []),
        preferred=["cell_type", "n_cells", "n_cells_used", "n_genes"],
    )

    print(json.dumps({key: value for key, value in summary.items() if key not in {"obs_columns"}}, indent=2))
    print()
    for path in (summary_path, cell_counts_path, signature_manifest_path):
        print(path)
    if signature_summary.get("signature_output"):
        print(signature_summary["signature_output"])
    if h5ad_summary.get("written"):
        print(h5ad_summary["path"])

    if args.strict and failures:
        print(
            "Strict author signature build failed: " + "; ".join(failures),
            file=sys.stderr,
        )
        return 2
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
