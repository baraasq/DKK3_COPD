#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_NOTEBOOK_ZIP = project_path(
    "data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip"
)
DEFAULT_TERMS = [
    "read_h5ad",
    "write_h5ad",
    ".h5ad",
    "read_10x_h5",
    "read_10x_mtx",
    "sc.read",
    ".write",
    "to_csv",
    "obs[",
    ".obs",
    "cell_type",
    "celltype",
    "annotation",
    "annotat",
    "marker_gene",
    "rank_genes_groups",
    "leiden",
    "dendrogram",
]
PATH_PATTERN = re.compile(
    r"""["']([^"']+\.(?:h5ad|h5|csv|tsv|txt|xlsx|rds|RDS|pkl|pickle))["']"""
)


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


def notebook_members(zip_path: Path) -> list[str]:
    with zipfile.ZipFile(zip_path) as archive:
        return sorted(
            name
            for name in archive.namelist()
            if name.casefold().endswith(".ipynb") and not name.endswith("/")
        )


def read_notebook(zip_path: Path, member: str) -> dict:
    with zipfile.ZipFile(zip_path) as archive:
        return json.loads(archive.read(member).decode("utf-8", errors="replace"))


def cell_source(cell: dict) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return str(source)


def matching_terms(text: str, terms: list[str]) -> list[str]:
    folded = text.casefold()
    return [term for term in terms if term.casefold() in folded]


def relevant_lines(text: str, terms: list[str], *, max_lines: int = 12) -> list[str]:
    output = []
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if any(term.casefold() in stripped.casefold() for term in terms):
            output.append(stripped[:500])
        if len(output) >= max_lines:
            break
    return output


def inspect_notebook(zip_path: Path, member: str, terms: list[str]) -> tuple[dict, list[dict]]:
    notebook = read_notebook(zip_path, member)
    cells = notebook.get("cells", [])
    code_cells = [cell for cell in cells if cell.get("cell_type") == "code"]
    rows = []
    term_counter: Counter[str] = Counter()
    paths = set()
    for index, cell in enumerate(cells):
        if cell.get("cell_type") != "code":
            continue
        text = cell_source(cell)
        terms_found = matching_terms(text, terms)
        if not terms_found:
            continue
        term_counter.update(terms_found)
        paths.update(PATH_PATTERN.findall(text))
        lines = relevant_lines(text, terms)
        rows.append(
            {
                "notebook": member,
                "cell_index": index,
                "matched_terms": ";".join(terms_found),
                "mentioned_paths": ";".join(PATH_PATTERN.findall(text)),
                "relevant_lines": " || ".join(lines),
            }
        )
    summary = {
        "notebook": member,
        "n_cells": len(cells),
        "n_code_cells": len(code_cells),
        "n_matching_code_cells": len(rows),
        "matched_term_counts": dict(sorted(term_counter.items())),
        "mentioned_paths": sorted(paths),
        "has_read_h5ad": bool(term_counter.get("read_h5ad") or term_counter.get("sc.read")),
        "has_write_h5ad": bool(term_counter.get("write_h5ad") or term_counter.get(".write")),
        "has_celltype_or_annotation_terms": any(
            term_counter.get(term) for term in ("cell_type", "celltype", "annotation", "annotat")
        ),
    }
    return summary, rows


def summarize_archive(zip_path: Path, terms: list[str]) -> tuple[dict, list[dict], list[dict]]:
    members = notebook_members(zip_path)
    notebook_summaries = []
    match_rows = []
    for member in members:
        summary, rows = inspect_notebook(zip_path, member, terms)
        notebook_summaries.append(summary)
        match_rows.extend(rows)
    overall = {
        "zip_path": str(zip_path),
        "zip_exists": zip_path.exists(),
        "n_notebooks": len(members),
        "notebooks": members,
        "notebooks_with_write_h5ad_or_write": [
            row["notebook"] for row in notebook_summaries if row["has_write_h5ad"]
        ],
        "notebooks_with_celltype_or_annotation_terms": [
            row["notebook"]
            for row in notebook_summaries
            if row["has_celltype_or_annotation_terms"]
        ],
        "interpretation": (
            "A write_h5ad/.write mention suggests the notebooks may create an annotated object. "
            "If only read_10x_h5/leiden/marker terms are present, the deposited code likely "
            "contains an annotation workflow but not a ready reference object."
        ),
    }
    return overall, notebook_summaries, match_rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Inspect GSE302339 authors' Scanpy notebooks for annotated-object and cell-type-signature clues."
    )
    parser.add_argument("--zip", default=str(DEFAULT_NOTEBOOK_ZIP))
    parser.add_argument(
        "--terms",
        default=",".join(DEFAULT_TERMS),
        help="Comma-separated case-insensitive terms to search in notebook code cells.",
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    zip_path = Path(args.zip).expanduser()
    terms = [term.strip() for term in args.terms.split(",") if term.strip()]
    output = ensure_results_dirs(load_config())
    meta_dir = output["meta"]
    table_dir = output["tables"]

    if not zip_path.exists():
        summary = {
            "zip_path": str(zip_path),
            "zip_exists": False,
            "ready": False,
            "failure": "Notebook zip not found",
        }
        summary_path = meta_dir / "gse302339_scanpy_notebook_inspection_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        if args.strict:
            return 2
        return 1

    overall, notebook_summaries, match_rows = summarize_archive(zip_path, terms)
    summary = {
        **overall,
        "notebook_summaries": notebook_summaries,
    }
    summary_path = meta_dir / "gse302339_scanpy_notebook_inspection_summary.json"
    notebook_table = table_dir / "gse302339_scanpy_notebook_summary.csv"
    match_table = table_dir / "gse302339_scanpy_notebook_code_matches.csv"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        notebook_table,
        notebook_summaries,
        preferred=[
            "notebook",
            "n_cells",
            "n_code_cells",
            "n_matching_code_cells",
            "has_read_h5ad",
            "has_write_h5ad",
            "has_celltype_or_annotation_terms",
            "mentioned_paths",
        ],
    )
    write_csv(
        match_table,
        match_rows,
        preferred=["notebook", "cell_index", "matched_terms", "mentioned_paths", "relevant_lines"],
    )

    print(json.dumps(overall, indent=2))
    print()
    for path in (summary_path, notebook_table, match_table):
        print(path)

    failures = []
    if not notebook_summaries:
        failures.append("No notebooks found")
    if args.strict and failures:
        print("Strict notebook inspection failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
