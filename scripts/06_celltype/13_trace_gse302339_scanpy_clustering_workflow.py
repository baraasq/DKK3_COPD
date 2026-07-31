#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_CODE_DIR = project_path("intermediate/gse302339_scanpy_workflow_code")

CELL_MARKER_PATTERN = re.compile(r"^# %% \[cell (?P<cell_index>\d+)\]")
NOTEBOOK_HEADER_PATTERN = re.compile(r"^# Notebook:\s*(?P<notebook>.+)")
OPEN_PATTERN = re.compile(
    r"""open\s*\(\s*["'](?P<path>[^"']+)["']\s*,\s*["'](?P<mode>[^"']+)["']""",
    re.IGNORECASE,
)
SCANPY_CALL_PATTERN = re.compile(
    r"""(?P<prefix>sc|sce)\.(?P<section>pp|tl|pl)\.(?P<function>[A-Za-z_]\w*)\s*\(\s*(?P<object>[A-Za-z_]\w*)\s*(?P<args>.*)"""
)
SUBSET_ASSIGNMENT_PATTERN = re.compile(
    r"""^(?P<target>[A-Za-z_]\w*)\s*=\s*(?P<source>[A-Za-z_]\w*)\s*\[(?P<selector>.+)"""
)
OBS_ASSIGNMENT_PATTERN = re.compile(
    r"""(?P<object>[A-Za-z_]\w*)\.obs(?:\[['"](?P<bracket>[^'"]+)['"]\]|\.(?P<attr>[A-Za-z_]\w*))\s*="""
)
DICT_ASSIGNMENT_PATTERN = re.compile(
    r"""(?P<dict_name>\w*celldict\w*)\s*=\s*\{""",
    re.IGNORECASE,
)
CLUSTERING_FUNCTIONS = {
    "highly_variable_genes",
    "regress_out",
    "scale",
    "pca",
    "neighbors",
    "harmony_integrate",
    "umap",
    "leiden",
    "dendrogram",
}
PLOT_FUNCTIONS = {
    "umap",
    "dotplot",
}


def code_files(code_dir: Path, selected: list[str] | None) -> list[Path]:
    if selected:
        return [project_path(path) for path in selected]
    return sorted(code_dir.glob("*.py"))


def parse_call_argument(args: str, name: str) -> str:
    match = re.search(rf"{re.escape(name)}\s*=\s*([^,\)]+)", args)
    return match.group(1).strip().strip("\"'") if match else ""


def compact_code(line: str) -> str:
    return re.sub(r"\s+", " ", line.strip())[:500]


def open_direction(mode: str) -> str:
    folded = mode.casefold()
    if any(flag in folded for flag in ("w", "a", "x")):
        return "write"
    if "r" in folded:
        return "read"
    return "unknown"


def detect_steps_in_line(
    *,
    line: str,
    code_file: Path,
    notebook: str,
    cell_index: str,
    source_line: int,
    step_order: int,
) -> list[dict]:
    stripped = line.strip()
    if not stripped or stripped.startswith("#"):
        return []

    rows: list[dict] = []
    common = {
        "code_file": str(code_file),
        "notebook": notebook,
        "cell_index": cell_index,
        "source_line": source_line,
        "step_order": step_order,
        "code": compact_code(line),
    }

    open_match = OPEN_PATTERN.search(line)
    if open_match:
        mode = open_match.group("mode")
        rows.append(
            {
                **common,
                "kind": "artifact_open",
                "object": "",
                "function": "open",
                "parameters": f"mode={mode}",
                "resolution": "",
                "flavor": "",
                "n_iterations": "",
                "use_rep": "",
                "n_pcs": "",
                "artifact_path": open_match.group("path"),
                "direction": open_direction(mode),
                "target_object": "",
                "source_object": "",
                "assigned_obs_column": "",
                "dictionary_name": "",
            }
        )

    scanpy_match = SCANPY_CALL_PATTERN.search(line)
    if scanpy_match:
        function = scanpy_match.group("function")
        section = scanpy_match.group("section")
        kind = "scanpy_call"
        if function in CLUSTERING_FUNCTIONS:
            kind = "clustering_step"
        elif section == "pl" or function in PLOT_FUNCTIONS:
            kind = "plot_step"
        args = scanpy_match.group("args")
        rows.append(
            {
                **common,
                "kind": kind,
                "object": scanpy_match.group("object"),
                "function": f"{scanpy_match.group('prefix')}.{section}.{function}",
                "parameters": args.strip()[:300],
                "resolution": parse_call_argument(args, "resolution"),
                "flavor": parse_call_argument(args, "flavor"),
                "n_iterations": parse_call_argument(args, "n_iterations"),
                "use_rep": parse_call_argument(args, "use_rep"),
                "n_pcs": parse_call_argument(args, "n_pcs"),
                "artifact_path": "",
                "direction": "",
                "target_object": "",
                "source_object": "",
                "assigned_obs_column": "",
                "dictionary_name": "",
            }
        )

    subset_match = SUBSET_ASSIGNMENT_PATTERN.search(stripped)
    if subset_match:
        rows.append(
            {
                **common,
                "kind": "object_subset_assignment",
                "object": subset_match.group("target"),
                "function": "subset",
                "parameters": subset_match.group("selector").strip()[:300],
                "resolution": "",
                "flavor": "",
                "n_iterations": "",
                "use_rep": "",
                "n_pcs": "",
                "artifact_path": "",
                "direction": "",
                "target_object": subset_match.group("target"),
                "source_object": subset_match.group("source"),
                "assigned_obs_column": "",
                "dictionary_name": "",
            }
        )

    obs_match = OBS_ASSIGNMENT_PATTERN.search(line)
    if obs_match:
        rows.append(
            {
                **common,
                "kind": "obs_assignment",
                "object": obs_match.group("object"),
                "function": "obs_assignment",
                "parameters": "",
                "resolution": "",
                "flavor": "",
                "n_iterations": "",
                "use_rep": "",
                "n_pcs": "",
                "artifact_path": "",
                "direction": "",
                "target_object": obs_match.group("object"),
                "source_object": "",
                "assigned_obs_column": obs_match.group("bracket")
                or obs_match.group("attr")
                or "",
                "dictionary_name": "",
            }
        )

    dict_match = DICT_ASSIGNMENT_PATTERN.search(line)
    if dict_match:
        rows.append(
            {
                **common,
                "kind": "celltype_dictionary_start",
                "object": "",
                "function": "dictionary_assignment",
                "parameters": "",
                "resolution": "",
                "flavor": "",
                "n_iterations": "",
                "use_rep": "",
                "n_pcs": "",
                "artifact_path": "",
                "direction": "",
                "target_object": "",
                "source_object": "",
                "assigned_obs_column": "",
                "dictionary_name": dict_match.group("dict_name"),
            }
        )

    return rows


def trace_code_file(code_file: Path, *, start_order: int) -> tuple[list[dict], int]:
    rows = []
    notebook = ""
    cell_index = ""
    step_order = start_order
    text = code_file.read_text(encoding="utf-8", errors="replace")
    for line_number, line in enumerate(text.splitlines(), start=1):
        notebook_match = NOTEBOOK_HEADER_PATTERN.match(line)
        if notebook_match:
            notebook = notebook_match.group("notebook").strip()
        cell_match = CELL_MARKER_PATTERN.match(line)
        if cell_match:
            cell_index = cell_match.group("cell_index")
        detected = detect_steps_in_line(
            line=line,
            code_file=code_file,
            notebook=notebook,
            cell_index=cell_index,
            source_line=line_number,
            step_order=step_order,
        )
        if detected:
            rows.extend(detected)
            step_order += len(detected)
    return rows, step_order


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


def summarize(rows: list[dict], input_files: list[Path]) -> dict:
    by_object: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("object"):
            by_object[row["object"]].append(row)

    clustering_by_object = []
    for object_name, object_rows in sorted(by_object.items()):
        clustering_rows = [
            row
            for row in object_rows
            if row.get("kind") == "clustering_step"
        ]
        if not clustering_rows:
            continue
        clustering_by_object.append(
            {
                "object": object_name,
                "n_clustering_steps": len(clustering_rows),
                "steps": [
                    {
                        "function": row["function"],
                        "source": f"{Path(row['code_file']).name}:{row['source_line']}",
                        "cell_index": row["cell_index"],
                        "resolution": row["resolution"],
                        "flavor": row["flavor"],
                        "n_iterations": row["n_iterations"],
                        "use_rep": row["use_rep"],
                        "n_pcs": row["n_pcs"],
                        "code": row["code"],
                    }
                    for row in clustering_rows
                ],
            }
        )

    artifact_edges = [
        {
            "direction": row["direction"],
            "artifact_path": row["artifact_path"],
            "source": f"{Path(row['code_file']).name}:{row['source_line']}",
            "cell_index": row["cell_index"],
            "code": row["code"],
        }
        for row in rows
        if row.get("kind") == "artifact_open"
    ]

    return {
        "code_dir": str(DEFAULT_CODE_DIR),
        "input_files": [str(path) for path in input_files],
        "n_input_files": len(input_files),
        "n_steps": len(rows),
        "n_clustering_steps": sum(
            1 for row in rows if row.get("kind") == "clustering_step"
        ),
        "objects_with_clustering": [row["object"] for row in clustering_by_object],
        "clustering_by_object": clustering_by_object,
        "artifact_edges": artifact_edges,
        "interpretation": (
            "Leiden cluster IDs are object-local. A cell-type dictionary should be "
            "applied only to the object variable that was clustered immediately "
            "before that dictionary assignment, or to a saved object known to be "
            "that exact variable state."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Trace clustering, subsetting, annotation, and pickle read/write steps "
            "in exported GSE302339 Scanpy workflow code."
        )
    )
    parser.add_argument("--code-dir", default=str(DEFAULT_CODE_DIR))
    parser.add_argument(
        "--code-file",
        action="append",
        dest="code_files",
        help="Specific exported notebook .py file to trace. May be repeated.",
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args(argv)

    code_dir = project_path(args.code_dir)
    selected_files = code_files(code_dir, args.code_files)
    results = ensure_results_dirs(load_config())
    steps_path = results["tables"] / "gse302339_scanpy_clustering_workflow_steps.csv"
    summary_path = results["meta"] / "gse302339_scanpy_clustering_workflow_summary.json"

    rows = []
    next_order = 1
    for code_file in selected_files:
        if not code_file.exists():
            continue
        file_rows, next_order = trace_code_file(code_file, start_order=next_order)
        rows.extend(file_rows)

    summary = summarize(rows, selected_files)
    summary.update(
        {
            "code_dir_exists": code_dir.exists(),
            "output_paths": {
                "steps": str(steps_path),
                "summary": str(summary_path),
            },
        }
    )

    write_csv(
        steps_path,
        rows,
        preferred=[
            "step_order",
            "notebook",
            "cell_index",
            "code_file",
            "source_line",
            "kind",
            "object",
            "function",
            "resolution",
            "flavor",
            "n_iterations",
            "use_rep",
            "n_pcs",
            "target_object",
            "source_object",
            "assigned_obs_column",
            "dictionary_name",
            "artifact_path",
            "direction",
            "parameters",
            "code",
        ],
    )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(steps_path)

    failures = []
    if not code_dir.exists():
        failures.append(f"Code directory missing: {code_dir}")
    if not selected_files:
        failures.append("No exported notebook .py files selected")
    if not rows:
        failures.append("No workflow steps detected")
    if args.strict and failures:
        raise SystemExit("Strict clustering workflow trace failed: " + "; ".join(failures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
