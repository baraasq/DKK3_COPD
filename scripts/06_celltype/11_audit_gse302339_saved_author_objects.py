#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gc
import json
import pickle
import sys
from collections import Counter
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_OBJECTS = [
    ("full", project_path("output/adata_harmony_annotated_cr8")),
    ("parenchyma", project_path("output/parenchyma_harmony_annotated_cr8")),
    ("immune", project_path("output/immune_harmony_annotated_cr8")),
]
DEFAULT_CLUSTER_MAP = project_path(
    "results/tables/gse302339_author_celltype_cluster_maps.csv"
)
DEFAULT_LABEL_COLUMNS = [
    "celltype_level1",
    "parenchyma_celltype_level1",
    "immune_celltype_level1",
]
DEFAULT_REQUIRED_CELL_TYPES = ["AT1", "AT2", "Fibroblast"]


def valid_label(value: object) -> bool:
    if value is None:
        return False
    text = str(value).strip()
    return bool(text) and text.casefold() not in {"nan", "none", "na", "null", "<na>"}


def has_anndata_interface(value: object) -> bool:
    return all(hasattr(value, attr) for attr in ("obs", "var", "X")) and hasattr(
        value, "n_obs"
    )


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
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_cluster_ids(value: object) -> set[str]:
    clusters: set[str] = set()
    for item in str(value or "").replace(",", ";").split(";"):
        item = item.strip()
        if item:
            clusters.add(item)
    return clusters


def cluster_maps(rows: list[dict]) -> dict[str, dict[str, set[str]]]:
    maps: dict[str, dict[str, set[str]]] = {}
    for row in rows:
        assigned = str(row.get("assigned_obs_column", "")).strip()
        label = str(row.get("celltype_label", "")).strip()
        clusters = parse_cluster_ids(row.get("cluster_ids"))
        if not assigned or not label or not clusters:
            continue
        maps.setdefault(assigned, {}).setdefault(label, set()).update(clusters)
    return maps


def sorted_cluster_values(values: set[str]) -> list[str]:
    def key(value: str) -> tuple[int, int | str]:
        return (0, int(value)) if value.isdigit() else (1, value)

    return sorted(values, key=key)


def label_counts(adata: Any, column: str) -> Counter:
    if column not in adata.obs.columns:
        return Counter()
    return Counter(
        str(value).strip()
        for value in adata.obs[column].astype(str)
        if valid_label(value)
    )


def cluster_counts(adata: Any, cluster_column: str) -> Counter:
    if cluster_column not in adata.obs.columns:
        return Counter()
    return Counter(
        str(value).strip()
        for value in adata.obs[cluster_column].astype(str)
        if valid_label(value)
    )


def object_summary(
    *,
    object_name: str,
    path: Path,
    load_summary: dict,
    adata: Any | None,
    label_columns: list[str],
    cluster_column: str,
) -> dict:
    summary = {
        "object_name": object_name,
        **load_summary,
        "n_cells": None,
        "n_genes": None,
        "obs_columns": [],
        "label_columns_present": [],
        "n_leiden_clusters": None,
        "leiden_clusters": [],
        "observed_label_counts": {},
    }
    if adata is None:
        return summary

    obs_columns = [str(column) for column in adata.obs.columns]
    labels_present = [column for column in label_columns if column in obs_columns]
    clusters = cluster_counts(adata, cluster_column)
    observed_labels = {
        column: dict(label_counts(adata, column).most_common())
        for column in labels_present
    }

    summary.update(
        {
            "n_cells": int(getattr(adata, "n_obs")),
            "n_genes": int(getattr(adata, "n_vars", len(getattr(adata, "var", [])))),
            "obs_columns": obs_columns,
            "label_columns_present": labels_present,
            "n_leiden_clusters": len(clusters),
            "leiden_clusters": sorted_cluster_values(set(clusters)),
            "observed_label_counts": observed_labels,
        }
    )
    return summary


def compatibility_rows(
    *,
    object_name: str,
    adata: Any,
    maps: dict[str, dict[str, set[str]]],
    label_columns: list[str],
    cluster_column: str,
    required_cell_types: list[str],
) -> tuple[list[dict], list[dict]]:
    clusters = cluster_counts(adata, cluster_column)
    observed_clusters = set(clusters)
    map_rows: list[dict] = []
    required_rows: list[dict] = []

    for assigned_column, label_map in sorted(maps.items()):
        mapped_clusters = set().union(*label_map.values()) if label_map else set()
        overlap_clusters = mapped_clusters & observed_clusters
        natural_column_present = assigned_column in adata.obs.columns
        direct_counts = label_counts(adata, assigned_column)
        labels_recoverable_by_cluster = {
            label
            for label, label_clusters in label_map.items()
            if label_clusters & observed_clusters
        }

        map_rows.append(
            {
                "object_name": object_name,
                "assigned_obs_column": assigned_column,
                "natural_label_column_present": natural_column_present,
                "n_observed_clusters": len(observed_clusters),
                "n_mapped_clusters": len(mapped_clusters),
                "n_overlap_clusters": len(overlap_clusters),
                "overlap_fraction_of_map": (
                    len(overlap_clusters) / len(mapped_clusters)
                    if mapped_clusters
                    else ""
                ),
                "labels_in_map": ";".join(sorted(label_map)),
                "labels_directly_observed": ";".join(sorted(direct_counts)),
                "labels_recoverable_by_cluster_id": ";".join(
                    sorted(labels_recoverable_by_cluster)
                ),
                "required_directly_observed": ";".join(
                    label for label in required_cell_types if label in direct_counts
                ),
                "required_recoverable_by_cluster_id": ";".join(
                    label
                    for label in required_cell_types
                    if label in labels_recoverable_by_cluster
                ),
                "safe_default_for_signature": natural_column_present
                and all(label in direct_counts for label in required_cell_types),
                "warning": (
                    ""
                    if natural_column_present
                    else "mapped column is not present in this saved object; cluster-ID relabeling would need manual verification"
                ),
            }
        )

        for label in required_cell_types:
            label_clusters = label_map.get(label, set())
            overlap = label_clusters & observed_clusters
            required_rows.append(
                {
                    "object_name": object_name,
                    "assigned_obs_column": assigned_column,
                    "required_cell_type": label,
                    "natural_label_column_present": natural_column_present,
                    "direct_cell_count": int(direct_counts.get(label, 0)),
                    "map_cluster_ids": ";".join(sorted_cluster_values(label_clusters)),
                    "overlap_cluster_ids": ";".join(sorted_cluster_values(overlap)),
                    "n_overlap_clusters": len(overlap),
                    "n_cells_if_cluster_map_forced": int(
                        sum(clusters[cluster] for cluster in overlap)
                    ),
                    "safe_default_for_signature": natural_column_present
                    and int(direct_counts.get(label, 0)) > 0,
                }
            )

    return map_rows, required_rows


def parse_object_args(values: list[str] | None) -> list[tuple[str, Path]]:
    if not values:
        return DEFAULT_OBJECTS
    objects: list[tuple[str, Path]] = []
    for value in values:
        if "=" in value:
            name, path = value.split("=", 1)
        else:
            path = value
            name = Path(value).name
        objects.append((name.strip(), project_path(path.strip())))
    return objects


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit saved GSE302339 author AnnData objects against extracted "
            "cluster-to-cell-type maps before building GeoMx deconvolution signatures."
        )
    )
    parser.add_argument(
        "--object",
        action="append",
        dest="objects",
        help=(
            "Object to audit as name=path or path. May be repeated. "
            "Defaults to full, parenchyma, and immune author outputs."
        ),
    )
    parser.add_argument("--cluster-map", default=str(DEFAULT_CLUSTER_MAP))
    parser.add_argument("--cluster-column", default="leiden")
    parser.add_argument(
        "--label-column",
        action="append",
        dest="label_columns",
        help="Label column to summarize. May be repeated.",
    )
    parser.add_argument(
        "--required-cell-type",
        action="append",
        dest="required_cell_types",
        help=(
            "Cell type required for the planned signature. May be repeated. "
            "Defaults to AT1, AT2, Fibroblast."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    config = load_config()
    result_dirs = ensure_results_dirs(config)
    objects = parse_object_args(args.objects)
    cluster_map_path = project_path(args.cluster_map)
    label_columns = args.label_columns or DEFAULT_LABEL_COLUMNS
    required_cell_types = args.required_cell_types or DEFAULT_REQUIRED_CELL_TYPES

    cluster_map_rows = read_csv(cluster_map_path)
    maps = cluster_maps(cluster_map_rows)

    object_summaries: list[dict] = []
    map_compatibility: list[dict] = []
    required_compatibility: list[dict] = []
    failures: list[str] = []

    if not cluster_map_path.exists():
        failures.append(f"Cluster map file missing: {cluster_map_path}")
    if not maps:
        failures.append("No cluster maps parsed from extracted author map file")

    for object_name, object_path in objects:
        adata, load_summary = load_author_object(object_path)
        summary = object_summary(
            object_name=object_name,
            path=object_path,
            load_summary=load_summary,
            adata=adata,
            label_columns=label_columns,
            cluster_column=args.cluster_column,
        )
        object_summaries.append(summary)

        if adata is None:
            failures.append(f"{object_name} object not loadable: {load_summary['status']}")
            continue

        obj_map_rows, obj_required_rows = compatibility_rows(
            object_name=object_name,
            adata=adata,
            maps=maps,
            label_columns=label_columns,
            cluster_column=args.cluster_column,
            required_cell_types=required_cell_types,
        )
        map_compatibility.extend(obj_map_rows)
        required_compatibility.extend(obj_required_rows)
        del adata
        gc.collect()

    safe_signature_candidates = [
        row
        for row in map_compatibility
        if str(row.get("safe_default_for_signature")) == "True"
        or row.get("safe_default_for_signature") is True
    ]
    if not safe_signature_candidates:
        failures.append(
            "No saved object naturally contains all required cell types in its own label column"
        )

    map_path = result_dirs["tables"] / "gse302339_saved_object_cluster_map_compatibility.csv"
    required_path = (
        result_dirs["tables"]
        / "gse302339_saved_object_required_celltype_compatibility.csv"
    )
    object_path = result_dirs["tables"] / "gse302339_saved_author_object_summary.csv"
    summary_path = result_dirs["meta"] / "gse302339_saved_author_object_audit_summary.json"

    write_csv(map_path, map_compatibility)
    write_csv(required_path, required_compatibility)
    write_csv(
        object_path,
        [
            {
                "object_name": row["object_name"],
                "path": row["path"],
                "exists": row["exists"],
                "status": row["status"],
                "loader": row.get("loader", ""),
                "n_cells": row["n_cells"],
                "n_genes": row["n_genes"],
                "label_columns_present": ";".join(row["label_columns_present"]),
                "n_leiden_clusters": row["n_leiden_clusters"],
                "leiden_clusters": ";".join(row["leiden_clusters"]),
            }
            for row in object_summaries
        ],
    )

    summary = {
        "objects": [
            {
                "object_name": row["object_name"],
                "path": row["path"],
                "exists": row["exists"],
                "status": row["status"],
                "n_cells": row["n_cells"],
                "n_genes": row["n_genes"],
                "label_columns_present": row["label_columns_present"],
                "n_leiden_clusters": row["n_leiden_clusters"],
                "observed_label_counts": row["observed_label_counts"],
            }
            for row in object_summaries
        ],
        "cluster_map": {
            "path": str(cluster_map_path),
            "exists": cluster_map_path.exists(),
            "n_rows": len(cluster_map_rows),
            "assigned_obs_columns": sorted(maps),
        },
        "required_cell_types": required_cell_types,
        "safe_signature_candidates": safe_signature_candidates,
        "ready_for_required_signature": not failures,
        "failures": failures,
        "outputs": {
            "object_summary": str(object_path),
            "map_compatibility": str(map_path),
            "required_celltype_compatibility": str(required_path),
            "summary": str(summary_path),
        },
    }

    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    for path in [summary_path, object_path, map_path, required_path]:
        print(path)

    if args.strict and failures:
        raise SystemExit("Strict saved-object audit failed: " + "; ".join(failures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
