#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_NOTEBOOK_ZIP = project_path(
    "data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip"
)
DEFAULT_OUTPUT_DIR = project_path("intermediate/gse302339_scanpy_workflow_code")
DEFAULT_NOTEBOOK_PATTERNS = [
    "2_celltype_annotation",
    "8_meta_merge",
]


def cell_source(cell: dict) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return str(source)


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


def safe_stem(member: str) -> str:
    stem = Path(member).stem
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", stem).strip("_")


def selected_members(members: list[str], patterns: list[str]) -> list[str]:
    folded_patterns = [pattern.casefold() for pattern in patterns]
    return [
        member
        for member in members
        if any(pattern in member.casefold() for pattern in folded_patterns)
    ]


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "notebook",
        "cell_index",
        "line_count",
        "has_celldict",
        "has_marker_dict",
        "has_open_output",
        "has_celltype_assignment",
        "code_file",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def export_notebook_code(zip_path: Path, member: str, output_dir: Path) -> tuple[Path, list[dict]]:
    notebook = read_notebook(zip_path, member)
    output_dir.mkdir(parents=True, exist_ok=True)
    code_file = output_dir / f"{safe_stem(member)}.py"
    rows = []
    chunks = [
        "# Exported from deposited GSE302339 Scanpy workflow notebook.",
        f"# Notebook: {member}",
        "# This is for audit/reconstruction; it is not a cleaned executable script.",
        "",
    ]
    for index, cell in enumerate(notebook.get("cells", [])):
        if cell.get("cell_type") != "code":
            continue
        source = cell_source(cell).rstrip()
        if not source:
            continue
        lowered = source.casefold()
        rows.append(
            {
                "notebook": member,
                "cell_index": index,
                "line_count": len(source.splitlines()),
                "has_celldict": "celldict" in lowered,
                "has_marker_dict": "marker_gene_dict" in lowered,
                "has_open_output": "open(" in lowered and "output/" in lowered,
                "has_celltype_assignment": "celltype" in lowered and ".obs" in lowered,
                "code_file": str(code_file),
            }
        )
        chunks.extend(
            [
                f"# %% [cell {index}]",
                source,
                "",
            ]
        )
    code_file.write_text("\n".join(chunks), encoding="utf-8")
    return code_file, rows


def export_annotation_code(
    zip_path: Path,
    output_dir: Path,
    patterns: list[str],
) -> tuple[dict, list[dict]]:
    members = notebook_members(zip_path)
    selected = selected_members(members, patterns)
    rows = []
    code_files = []
    for member in selected:
        code_file, member_rows = export_notebook_code(zip_path, member, output_dir)
        code_files.append(str(code_file))
        rows.extend(member_rows)
    summary = {
        "zip_path": str(zip_path),
        "zip_exists": zip_path.exists(),
        "patterns": patterns,
        "n_notebooks_in_zip": len(members),
        "selected_notebooks": selected,
        "n_selected_notebooks": len(selected),
        "code_files": code_files,
        "n_code_cells_exported": len(rows),
        "cells_with_celldict": sum(bool(row["has_celldict"]) for row in rows),
        "cells_with_marker_gene_dict": sum(bool(row["has_marker_dict"]) for row in rows),
        "cells_with_open_output": sum(bool(row["has_open_output"]) for row in rows),
        "cells_with_celltype_assignment": sum(
            bool(row["has_celltype_assignment"]) for row in rows
        ),
    }
    return summary, rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Export full code cells from GSE302339 annotation notebooks for reconstruction/audit."
    )
    parser.add_argument("--zip", default=str(DEFAULT_NOTEBOOK_ZIP))
    parser.add_argument("--output-dir", default=str(DEFAULT_OUTPUT_DIR))
    parser.add_argument(
        "--notebook-pattern",
        action="append",
        dest="patterns",
        help=(
            "Case-insensitive substring of notebook paths to export. "
            "Can be passed multiple times. Defaults to celltype annotation and meta merge."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    zip_path = Path(args.zip).expanduser()
    output_dir = Path(args.output_dir).expanduser()
    patterns = args.patterns or DEFAULT_NOTEBOOK_PATTERNS
    results = ensure_results_dirs(load_config())
    summary_path = results["meta"] / "gse302339_annotation_code_export_summary.json"
    table_path = results["tables"] / "gse302339_annotation_code_cells.csv"

    if not zip_path.exists():
        summary = {
            "zip_path": str(zip_path),
            "zip_exists": False,
            "failure": "Notebook zip not found",
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        return 2 if args.strict else 1

    summary, rows = export_annotation_code(zip_path, output_dir, patterns)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(table_path, rows)

    print(json.dumps(summary, indent=2))
    print()
    for path in [summary_path, table_path, *map(Path, summary["code_files"])]:
        print(path)

    if args.strict and not summary["selected_notebooks"]:
        print("Strict annotation code export failed: no matching notebooks", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
