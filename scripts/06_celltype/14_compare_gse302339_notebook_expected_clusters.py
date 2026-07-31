#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import pickle
import re
import sys
import zipfile
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_ZIP = project_path("data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip")
DEFAULT_OBJECTS = {
    "adata": project_path("output/adata_harmony_annotated_cr8"),
    "parenchyma": project_path("output/parenchyma_harmony_annotated_cr8"),
    "immune": project_path("output/immune_harmony_annotated_cr8"),
}
DEFAULT_MEMBERS = [
    "scanpy_workflow/1_preprocessing_doublet_detection.ipynb",
    "scanpy_workflow/2_celltype_annotation.ipynb",
]

SHAPE_RE = re.compile(r"^\((?P<n_obs>\d+),\s*(?P<n_vars>\d+)\)\s*$", re.MULTILINE)
FOUND_CLUSTERS_RE = re.compile(r"finished:\s+found\s+(?P<n>\d+)\s+clusters", re.IGNORECASE)
LEIDEN_CALL_RE = re.compile(r"sc\.tl\.leiden\(\s*(?P<object>[A-Za-z_]\w*)")


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
    raise TypeError("Loaded object is not AnnData-like.")


def output_text(cell: dict) -> str:
    parts: list[str] = []
    for out in cell.get("outputs", []):
        if "text" in out:
            value = out["text"]
            parts.append("".join(value) if isinstance(value, list) else str(value))
        data = out.get("data", {})
        if "text/plain" in data:
            value = data["text/plain"]
            parts.append("".join(value) if isinstance(value, list) else str(value))
    return "\n".join(parts)


def read_notebook(zip_path: Path, member: str) -> dict | None:
    try:
        with zipfile.ZipFile(zip_path) as archive:
            return json.loads(archive.read(member).decode("utf-8", errors="replace"))
    except Exception:
        return None


def parse_notebook_expectations(zip_path: Path, members: list[str]) -> tuple[list[dict], dict]:
    events: list[dict] = []
    summary = {
        "zip_path": str(zip_path),
        "zip_exists": zip_path.exists(),
        "members_requested": members,
        "members_read": [],
        "read_failures": [],
    }
    if not zip_path.exists():
        return events, summary

    for member in members:
        notebook = read_notebook(zip_path, member)
        if notebook is None:
            summary["read_failures"].append(member)
            continue
        summary["members_read"].append(member)

        for cell_index, cell in enumerate(notebook.get("cells", [])):
            source = "".join(cell.get("source", ""))
            outputs = output_text(cell)

            if "sc.pp.highly_variable_genes" in source:
                shapes = [
                    {
                        "n_obs": int(match.group("n_obs")),
                        "n_vars": int(match.group("n_vars")),
                    }
                    for match in SHAPE_RE.finditer(outputs)
                ]
                if shapes:
                    events.append(
                        {
                            "kind": "post_hvg_shape",
                            "object": "adata",
                            "notebook": member,
                            "cell_index": cell_index,
                            **shapes[-1],
                        }
                    )

            leiden_objects = [
                match.group("object") for match in LEIDEN_CALL_RE.finditer(source)
            ]
            found_clusters = [
                int(match.group("n")) for match in FOUND_CLUSTERS_RE.finditer(outputs)
            ]
            for object_name, n_clusters in zip(leiden_objects, found_clusters):
                events.append(
                    {
                        "kind": "leiden_cluster_count",
                        "object": object_name,
                        "notebook": member,
                        "cell_index": cell_index,
                        "n_clusters": n_clusters,
                    }
                )

    return events, summary


def load_object_summary(object_name: str, path: Path, cluster_column: str) -> dict:
    summary = {
        "object": object_name,
        "path": str(path),
        "exists": path.exists(),
        "status": "missing" if not path.exists() else "unknown",
        "n_obs": None,
        "n_vars": None,
        "n_leiden_clusters": None,
        "neighbors_params": None,
        "leiden_params": None,
        "label_counts": {},
    }
    if not path.exists():
        return summary

    try:
        with path.open("rb") as handle:
            adata = coerce_anndata_like(pickle.load(handle))
        summary["status"] = "ok"
    except Exception as exc:
        summary["status"] = "load_failed"
        summary["error"] = str(exc)
        return summary

    summary["n_obs"] = int(getattr(adata, "n_obs"))
    summary["n_vars"] = int(getattr(adata, "n_vars", len(getattr(adata, "var", []))))
    if cluster_column in adata.obs.columns:
        summary["n_leiden_clusters"] = int(adata.obs[cluster_column].astype(str).nunique())
    summary["neighbors_params"] = adata.uns.get("neighbors", {}).get("params")
    leiden_uns = adata.uns.get("leiden", {})
    summary["leiden_params"] = (
        leiden_uns.get("params", leiden_uns) if isinstance(leiden_uns, dict) else leiden_uns
    )

    for column in (
        "celltype_level1",
        "parenchyma_celltype_level1",
        "immune_celltype_level1",
        "lobe_emphysema_simple",
    ):
        if column in adata.obs.columns:
            counts = adata.obs[column].astype(str).value_counts().to_dict()
            summary["label_counts"][column] = {str(key): int(value) for key, value in counts.items()}
    return summary


def build_comparison_rows(expectations: list[dict], object_summaries: dict[str, dict]) -> list[dict]:
    rows: list[dict] = []
    for event in expectations:
        object_name = event["object"]
        observed = object_summaries.get(object_name, {})
        if event["kind"] == "post_hvg_shape":
            rows.append(
                {
                    "object": object_name,
                    "metric": "post_hvg_n_obs",
                    "expected": event["n_obs"],
                    "observed": observed.get("n_obs"),
                    "matches": observed.get("n_obs") == event["n_obs"],
                    "notebook": event["notebook"],
                    "cell_index": event["cell_index"],
                }
            )
            rows.append(
                {
                    "object": object_name,
                    "metric": "post_hvg_n_vars",
                    "expected": event["n_vars"],
                    "observed": observed.get("n_vars"),
                    "matches": observed.get("n_vars") == event["n_vars"],
                    "notebook": event["notebook"],
                    "cell_index": event["cell_index"],
                }
            )
        elif event["kind"] == "leiden_cluster_count":
            rows.append(
                {
                    "object": object_name,
                    "metric": "leiden_cluster_count",
                    "expected": event["n_clusters"],
                    "observed": observed.get("n_leiden_clusters"),
                    "matches": observed.get("n_leiden_clusters") == event["n_clusters"],
                    "notebook": event["notebook"],
                    "cell_index": event["cell_index"],
                }
            )
    return rows


def write_csv(path: Path, rows: list[dict]) -> None:
    import csv

    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Compare GSE302339 author notebook-recorded clustering outputs against "
            "the locally reconstructed saved AnnData objects."
        )
    )
    parser.add_argument("--zip-path", type=Path, default=DEFAULT_ZIP)
    parser.add_argument("--member", action="append", dest="members")
    parser.add_argument("--cluster-column", default="leiden")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    result_dirs = ensure_results_dirs(config)
    members = args.members or DEFAULT_MEMBERS

    expectations, notebook_summary = parse_notebook_expectations(
        project_path(args.zip_path), members
    )
    object_summaries = {
        name: load_object_summary(name, path, args.cluster_column)
        for name, path in DEFAULT_OBJECTS.items()
    }
    comparison_rows = build_comparison_rows(expectations, object_summaries)

    failures = [
        f"{row['object']} {row['metric']}: expected {row['expected']} observed {row['observed']}"
        for row in comparison_rows
        if not row["matches"]
    ]
    if not expectations:
        failures.append("No notebook-recorded shape or Leiden cluster-count expectations parsed.")

    summary = {
        "dataset": "GSE302339",
        "notebook_summary": notebook_summary,
        "n_expectation_events": len(expectations),
        "expectations": expectations,
        "object_summaries": object_summaries,
        "comparison_rows": comparison_rows,
        "all_expected_values_match": not failures,
        "failures": failures,
        "interpretation": (
            "If parenchyma expected cluster count differs from the saved object, "
            "author cluster-ID dictionaries should not be used directly for signatures."
        ),
    }

    summary_path = result_dirs["meta"] / "gse302339_notebook_expected_cluster_comparison_summary.json"
    table_path = result_dirs["tables"] / "gse302339_notebook_expected_cluster_comparison.csv"

    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(table_path, comparison_rows)

    print(
        json.dumps(
            {
                "dataset": summary["dataset"],
                "n_expectation_events": summary["n_expectation_events"],
                "all_expected_values_match": summary["all_expected_values_match"],
                "failures": summary["failures"],
                "outputs": {
                    "summary": str(summary_path),
                    "comparison_table": str(table_path),
                },
            },
            indent=2,
        )
    )
    print()
    print(summary_path)
    print(table_path)

    if args.strict and failures:
        raise SystemExit(
            "Strict notebook/object cluster comparison failed: " + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
