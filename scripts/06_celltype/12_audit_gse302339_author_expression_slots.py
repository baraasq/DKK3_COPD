#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import pickle
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_OBJECTS = [
    ("full", project_path("output/adata_harmony_annotated_cr8")),
    ("parenchyma", project_path("output/parenchyma_harmony_annotated_cr8")),
]
DEFAULT_PRIMARY_GENE = "DKK3"


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


def matrix_shape(matrix: Any) -> tuple[int | None, int | None]:
    shape = getattr(matrix, "shape", None)
    if shape is None or len(shape) != 2:
        return None, None
    return int(shape[0]), int(shape[1])


def summarize_matrix(matrix: Any, *, max_dense_rows: int, max_dense_cols: int) -> dict:
    import numpy as np

    n_rows, n_cols = matrix_shape(matrix)
    summary = {
        "n_rows": n_rows,
        "n_cols": n_cols,
        "matrix_type": type(matrix).__name__,
        "sparse": hasattr(matrix, "tocsr"),
        "n_values_checked": 0,
        "min_checked": None,
        "max_checked": None,
        "n_negative_checked": 0,
        "n_nonfinite_checked": 0,
        "looks_nonnegative": None,
    }
    if n_rows is None or n_cols is None or n_rows == 0 or n_cols == 0:
        return summary

    if hasattr(matrix, "tocsr"):
        data = matrix.tocsr().data
        checked = np.asarray(data)
        if checked.size:
            summary.update(
                {
                    "n_values_checked": int(checked.size),
                    "min_checked": float(np.nanmin(checked)),
                    "max_checked": float(np.nanmax(checked)),
                    "n_negative_checked": int(np.sum(checked < 0)),
                    "n_nonfinite_checked": int(np.sum(~np.isfinite(checked))),
                }
            )
        else:
            summary.update(
                {
                    "n_values_checked": 0,
                    "min_checked": 0.0,
                    "max_checked": 0.0,
                    "n_negative_checked": 0,
                    "n_nonfinite_checked": 0,
                }
            )
    else:
        row_stop = min(n_rows, max_dense_rows)
        col_stop = min(n_cols, max_dense_cols)
        checked = np.asarray(matrix[:row_stop, :col_stop])
        summary.update(
            {
                "n_values_checked": int(checked.size),
                "min_checked": float(np.nanmin(checked)) if checked.size else None,
                "max_checked": float(np.nanmax(checked)) if checked.size else None,
                "n_negative_checked": int(np.sum(checked < 0)),
                "n_nonfinite_checked": int(np.sum(~np.isfinite(checked))),
            }
        )

    summary["looks_nonnegative"] = (
        summary["n_negative_checked"] == 0 and summary["n_nonfinite_checked"] == 0
    )
    return summary


def gene_names_from_var(var: Any) -> list[str]:
    try:
        return [str(value) for value in var.index.tolist()]
    except Exception:
        return []


def raw_summary(adata: Any, *, primary_gene: str, max_dense_rows: int, max_dense_cols: int) -> dict:
    raw = getattr(adata, "raw", None)
    if raw is None:
        return {
            "exists": False,
        }
    var_names = gene_names_from_var(raw.var)
    return {
        "exists": True,
        "shape": [int(raw.shape[0]), int(raw.shape[1])],
        "n_var_names": len(var_names),
        "primary_gene_present": primary_gene.upper()
        in {gene.upper() for gene in var_names},
        "first_var_names": var_names[:20],
        "matrix": summarize_matrix(
            raw.X,
            max_dense_rows=max_dense_rows,
            max_dense_cols=max_dense_cols,
        ),
    }


def layer_summaries(
    adata: Any,
    *,
    max_dense_rows: int,
    max_dense_cols: int,
) -> dict[str, dict]:
    output = {}
    layers = getattr(adata, "layers", {})
    try:
        layer_names = sorted(str(name) for name in layers.keys())
    except Exception:
        layer_names = []
    for layer_name in layer_names:
        output[layer_name] = summarize_matrix(
            layers[layer_name],
            max_dense_rows=max_dense_rows,
            max_dense_cols=max_dense_cols,
        )
    return output


def object_summary(
    *,
    object_name: str,
    path: Path,
    load_summary: dict,
    adata: Any | None,
    primary_gene: str,
    max_dense_rows: int,
    max_dense_cols: int,
) -> dict:
    summary = {
        "object_name": object_name,
        **load_summary,
        "n_cells": None,
        "n_genes": None,
        "primary_gene_present_in_var_names": None,
        "first_var_names": [],
        "layers": {},
        "raw": {},
        "X": {},
        "interpretation": "object not loaded",
    }
    if adata is None:
        return summary

    var_names = gene_names_from_var(adata.var)
    x_summary = summarize_matrix(
        adata.X,
        max_dense_rows=max_dense_rows,
        max_dense_cols=max_dense_cols,
    )
    layers = layer_summaries(
        adata,
        max_dense_rows=max_dense_rows,
        max_dense_cols=max_dense_cols,
    )
    raw = raw_summary(
        adata,
        primary_gene=primary_gene,
        max_dense_rows=max_dense_rows,
        max_dense_cols=max_dense_cols,
    )

    usable_sources = []
    if x_summary.get("looks_nonnegative"):
        usable_sources.append("X")
    usable_sources.extend(
        f"layer:{name}"
        for name, layer in layers.items()
        if layer.get("looks_nonnegative")
    )
    if raw.get("exists") and raw.get("matrix", {}).get("looks_nonnegative"):
        usable_sources.append("raw.X")

    interpretation = (
        "nonnegative expression source detected: " + ", ".join(usable_sources)
        if usable_sources
        else "no nonnegative expression source detected in checked values; avoid building logCPM signatures from this object until counts/raw data are recovered"
    )

    summary.update(
        {
            "n_cells": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "primary_gene_present_in_var_names": primary_gene.upper()
            in {gene.upper() for gene in var_names},
            "first_var_names": var_names[:20],
            "layers": layers,
            "raw": raw,
            "X": x_summary,
            "usable_nonnegative_sources": usable_sources,
            "interpretation": interpretation,
        }
    )
    return summary


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


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit expression matrices/layers/raw slots in saved GSE302339 author objects "
            "before building logCPM signatures."
        )
    )
    parser.add_argument(
        "--object",
        action="append",
        dest="objects",
        help=(
            "Object to audit as name=path or path. May be repeated. "
            "Defaults to full and parenchyma author outputs."
        ),
    )
    parser.add_argument("--primary-gene", default=DEFAULT_PRIMARY_GENE)
    parser.add_argument("--max-dense-rows", type=int, default=2000)
    parser.add_argument("--max-dense-cols", type=int, default=2000)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    results = ensure_results_dirs(load_config())
    summaries = []
    failures = []
    for object_name, object_path in parse_object_args(args.objects):
        adata, load_summary = load_author_object(object_path)
        summary = object_summary(
            object_name=object_name,
            path=object_path,
            load_summary=load_summary,
            adata=adata,
            primary_gene=args.primary_gene,
            max_dense_rows=args.max_dense_rows,
            max_dense_cols=args.max_dense_cols,
        )
        summaries.append(summary)
        if summary["status"] != "ok":
            failures.append(f"{object_name} object not loadable: {summary['status']}")
        elif not summary.get("usable_nonnegative_sources"):
            failures.append(
                f"{object_name} has no checked nonnegative expression source"
            )

    table_rows = []
    for summary in summaries:
        table_rows.append(
            {
                "object_name": summary["object_name"],
                "path": summary["path"],
                "status": summary["status"],
                "n_cells": summary["n_cells"],
                "n_genes": summary["n_genes"],
                "primary_gene_present_in_var_names": summary[
                    "primary_gene_present_in_var_names"
                ],
                "x_min_checked": summary.get("X", {}).get("min_checked"),
                "x_max_checked": summary.get("X", {}).get("max_checked"),
                "x_n_negative_checked": summary.get("X", {}).get(
                    "n_negative_checked"
                ),
                "x_n_nonfinite_checked": summary.get("X", {}).get(
                    "n_nonfinite_checked"
                ),
                "layer_names": ";".join(sorted(summary.get("layers", {}))),
                "raw_exists": summary.get("raw", {}).get("exists"),
                "raw_shape": ";".join(
                    str(value) for value in summary.get("raw", {}).get("shape", [])
                ),
                "raw_primary_gene_present": summary.get("raw", {}).get(
                    "primary_gene_present"
                ),
                "raw_min_checked": summary.get("raw", {})
                .get("matrix", {})
                .get("min_checked"),
                "raw_max_checked": summary.get("raw", {})
                .get("matrix", {})
                .get("max_checked"),
                "raw_n_negative_checked": summary.get("raw", {})
                .get("matrix", {})
                .get("n_negative_checked"),
                "usable_nonnegative_sources": ";".join(
                    summary.get("usable_nonnegative_sources", [])
                ),
                "interpretation": summary["interpretation"],
            }
        )

    table_path = results["tables"] / "gse302339_author_expression_slot_audit.csv"
    summary_path = results["meta"] / "gse302339_author_expression_slot_audit_summary.json"
    write_csv(table_path, table_rows)
    payload = {
        "objects": summaries,
        "ready_for_signature_expression_source": not failures,
        "failures": failures,
        "outputs": {
            "summary": str(summary_path),
            "table": str(table_path),
        },
    }
    summary_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps(payload, indent=2))
    print()
    print(summary_path)
    print(table_path)

    if args.strict and failures:
        raise SystemExit(
            "Strict expression-slot audit failed: " + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
