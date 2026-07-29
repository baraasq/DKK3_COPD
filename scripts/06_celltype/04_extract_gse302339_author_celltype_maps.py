#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_CODE_DIR = project_path("intermediate/gse302339_scanpy_workflow_code")
DEFAULT_CODE_FILES = [
    DEFAULT_CODE_DIR / "2_celltype_annotation.py",
]
DICT_ASSIGNMENT_PATTERN = re.compile(r"(?P<name>\w*celldict\w*)\s*=\s*\{", re.IGNORECASE)
OBS_COLUMN_PATTERN = re.compile(
    r"(?P<object>[A-Za-z_]\w*)\.obs(?:\[['\"](?P<bracket>[^'\"]*celltype[^'\"]*)['\"]\]|"
    r"\.(?P<attribute>[A-Za-z_]\w*celltype[A-Za-z0-9_]*))\s*=",
    re.IGNORECASE,
)


def line_number_at(text: str, index: int) -> int:
    return text.count("\n", 0, index) + 1


def extract_balanced_brace_text(text: str, start_index: int) -> tuple[str, int]:
    depth = 0
    in_string: str | None = None
    escaped = False
    for index in range(start_index, len(text)):
        char = text[index]
        if in_string:
            if escaped:
                escaped = False
            elif char == "\\":
                escaped = True
            elif char == in_string:
                in_string = None
            continue
        if char in {"'", '"'}:
            in_string = char
            continue
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0:
                return text[start_index : index + 1], index + 1
    raise ValueError("Could not find matching closing brace")


def infer_obs_assignment(text: str, end_index: int, *, window_lines: int = 35) -> dict[str, str]:
    following = "\n".join(text[end_index:].splitlines()[:window_lines])
    match = OBS_COLUMN_PATTERN.search(following)
    if not match:
        return {
            "assigned_object": "",
            "assigned_obs_column": "",
        }
    return {
        "assigned_object": match.group("object") or "",
        "assigned_obs_column": match.group("bracket") or match.group("attribute") or "",
    }


def parse_cluster_values(value: object) -> list[str]:
    if isinstance(value, (list, tuple, set)):
        return [str(item) for item in value]
    return [str(value)]


def extract_celltype_maps_from_text(text: str, *, source: str) -> tuple[list[dict], list[dict]]:
    rows = []
    blocks = []
    for block_index, match in enumerate(DICT_ASSIGNMENT_PATTERN.finditer(text), start=1):
        dict_start = match.end() - 1
        name = match.group("name")
        start_line = line_number_at(text, match.start())
        try:
            block_text, end_index = extract_balanced_brace_text(text, dict_start)
            parsed = ast.literal_eval(block_text)
            parse_error = ""
        except Exception as exc:  # pragma: no cover - defensive audit path
            block_text = ""
            end_index = match.end()
            parsed = {}
            parse_error = str(exc)

        assignment = infer_obs_assignment(text, end_index)
        block = {
            "source_file": source,
            "dictionary_index": block_index,
            "dictionary_name": name,
            "start_line": start_line,
            "assigned_object": assignment["assigned_object"],
            "assigned_obs_column": assignment["assigned_obs_column"],
            "n_celltype_labels": len(parsed) if isinstance(parsed, dict) else 0,
            "parse_error": parse_error,
        }
        blocks.append(block)

        if not isinstance(parsed, dict):
            continue
        for label, clusters in parsed.items():
            cluster_values = parse_cluster_values(clusters)
            rows.append(
                {
                    **block,
                    "celltype_label": str(label),
                    "cluster_ids": ";".join(cluster_values),
                    "n_clusters": len(cluster_values),
                }
            )
    return rows, blocks


def write_csv(path: Path, rows: list[dict], preferred: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def extract_maps(code_files: list[Path]) -> tuple[dict, list[dict], list[dict]]:
    rows = []
    blocks = []
    file_summaries = []
    for code_file in code_files:
        exists = code_file.exists()
        summary = {
            "code_file": str(code_file),
            "exists": exists,
            "n_dictionary_blocks": 0,
            "n_mapping_rows": 0,
        }
        if exists:
            text = code_file.read_text(encoding="utf-8", errors="replace")
            file_rows, file_blocks = extract_celltype_maps_from_text(
                text,
                source=str(code_file),
            )
            rows.extend(file_rows)
            blocks.extend(file_blocks)
            summary["n_dictionary_blocks"] = len(file_blocks)
            summary["n_mapping_rows"] = len(file_rows)
        file_summaries.append(summary)

    summary = {
        "code_files": [str(path) for path in code_files],
        "file_summaries": file_summaries,
        "n_dictionary_blocks": len(blocks),
        "n_mapping_rows": len(rows),
        "obs_columns": sorted(
            {
                row["assigned_obs_column"]
                for row in rows
                if row.get("assigned_obs_column")
            }
        ),
        "celltype_labels": sorted({row["celltype_label"] for row in rows}),
        "ready_to_reconstruct_author_labels": bool(rows)
        and all(not block.get("parse_error") for block in blocks),
        "interpretation": (
            "These are the authors' notebook cluster-to-cell-type dictionaries. "
            "For pneumocyte-fibroblast deconvolution, prioritize mappings assigned "
            "to parenchyma_celltype_level1."
        ),
    }
    return summary, rows, blocks


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Parse author celldict cluster-to-celltype mappings from exported GSE302339 notebook code."
    )
    parser.add_argument(
        "--code-file",
        action="append",
        dest="code_files",
        help=(
            "Exported notebook .py file to parse. Can be passed multiple times. "
            "Defaults to intermediate/gse302339_scanpy_workflow_code/2_celltype_annotation.py."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    code_files = [Path(path).expanduser() for path in args.code_files] if args.code_files else DEFAULT_CODE_FILES
    results = ensure_results_dirs(load_config())
    summary_path = results["meta"] / "gse302339_author_celltype_map_summary.json"
    map_table = results["tables"] / "gse302339_author_celltype_cluster_maps.csv"
    block_table = results["tables"] / "gse302339_author_celltype_dictionary_blocks.csv"

    summary, rows, blocks = extract_maps(code_files)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_csv(
        map_table,
        rows,
        preferred=[
            "source_file",
            "dictionary_index",
            "dictionary_name",
            "start_line",
            "assigned_object",
            "assigned_obs_column",
            "celltype_label",
            "cluster_ids",
            "n_clusters",
            "parse_error",
        ],
    )
    write_csv(
        block_table,
        blocks,
        preferred=[
            "source_file",
            "dictionary_index",
            "dictionary_name",
            "start_line",
            "assigned_object",
            "assigned_obs_column",
            "n_celltype_labels",
            "parse_error",
        ],
    )

    print(json.dumps(summary, indent=2))
    print()
    for path in (summary_path, map_table, block_table):
        print(path)

    failures = []
    if any(not row["exists"] for row in summary["file_summaries"]):
        failures.append("One or more exported code files are missing")
    if not rows:
        failures.append("No celldict mapping rows parsed")
    if not summary["ready_to_reconstruct_author_labels"]:
        failures.append("One or more celldict blocks had parse errors")
    if args.strict and failures:
        print(
            "Strict author celltype map extraction failed: " + "; ".join(failures),
            file=sys.stderr,
        )
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
