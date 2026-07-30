#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import csv
import importlib.util
import json
import re
import sys
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_CODE_DIR = "intermediate/gse302339_scanpy_workflow_code"
DEFAULT_NOTEBOOK_ZIP = "data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip"
DEFAULT_OPTIONAL_MODULES = {"scvi", "mudata"}
IMPORT_RE = re.compile(r"^\s*import\s+([A-Za-z_][\w.]*)(?:\s+as\s+\w+)?", re.MULTILINE)
FROM_RE = re.compile(r"^\s*from\s+([A-Za-z_][\w.]*)\s+import\s+", re.MULTILINE)

PACKAGE_BY_MODULE = {
    "PIL": "pillow",
    "anndata": "anndata",
    "cv2": "opencv-python",
    "doubletdetection": "doubletdetection",
    "gseapy": "gseapy",
    "harmonypy": "harmonypy",
    "igraph": "igraph",
    "leidenalg": "leidenalg",
    "louvain": "louvain",
    "matplotlib": "matplotlib",
    "muon": "muon",
    "mudata": "mudata",
    "numpy": "numpy",
    "pandas": "pandas",
    "pertpy": "pertpy",
    "scanpy": "scanpy",
    "scipy": "scipy",
    "scvi": "scvi-tools",
    "seaborn": "seaborn",
    "sklearn": "scikit-learn",
    "statsmodels": "statsmodels",
}

INDIRECT_DEPENDENCY_RULES = [
    {
        "module": "harmonypy",
        "package": "harmonypy",
        "pattern": "harmony_integrate",
        "reason": "scanpy.external.pp.harmony_integrate imports harmonypy at runtime",
    },
    {
        "module": "leidenalg",
        "package": "leidenalg",
        "pattern": ".tl.leiden",
        "reason": "Scanpy Leiden clustering commonly requires leidenalg",
    },
    {
        "module": "igraph",
        "package": "igraph",
        "pattern": ".tl.leiden",
        "reason": "Scanpy Leiden clustering commonly requires igraph",
    },
]


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


def source_from_notebook_cell(cell: dict) -> str:
    source = cell.get("source", "")
    if isinstance(source, list):
        return "".join(source)
    return str(source)


def notebook_sources(zip_path: Path) -> list[dict]:
    rows: list[dict] = []
    if not zip_path.exists():
        return rows
    with zipfile.ZipFile(zip_path) as archive:
        for member in sorted(archive.namelist()):
            if not member.casefold().endswith(".ipynb"):
                continue
            notebook = json.loads(archive.read(member).decode("utf-8", errors="replace"))
            code_chunks = [
                source_from_notebook_cell(cell)
                for cell in notebook.get("cells", [])
                if cell.get("cell_type") == "code"
            ]
            rows.append(
                {
                    "source_id": member,
                    "source_kind": "notebook",
                    "text": "\n".join(code_chunks),
                }
            )
    return rows


def exported_code_sources(code_dir: Path) -> list[dict]:
    if not code_dir.exists():
        return []
    rows: list[dict] = []
    for path in sorted(code_dir.glob("*.py")):
        rows.append(
            {
                "source_id": str(path),
                "source_kind": "exported_py",
                "text": path.read_text(encoding="utf-8", errors="replace"),
            }
        )
    return rows


def imported_modules_from_source(text: str) -> set[str]:
    modules: set[str] = set()
    try:
        tree = ast.parse(text)
    except SyntaxError:
        tree = None

    if tree is not None:
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                modules.update(alias.name.split(".", 1)[0] for alias in node.names)
            elif isinstance(node, ast.ImportFrom) and node.module:
                modules.add(node.module.split(".", 1)[0])
        return modules

    modules.update(match.group(1).split(".", 1)[0] for match in IMPORT_RE.finditer(text))
    modules.update(match.group(1).split(".", 1)[0] for match in FROM_RE.finditer(text))
    return modules


def indirect_modules_from_source(text: str) -> list[dict]:
    rows = []
    for rule in INDIRECT_DEPENDENCY_RULES:
        if rule["pattern"] in text:
            rows.append(
                {
                    "module": rule["module"],
                    "package": rule["package"],
                    "source": "indirect_runtime",
                    "reason": rule["reason"],
                }
            )
    return rows


def is_stdlib_module(module: str) -> bool:
    if module in {"__future__"}:
        return True
    stdlib = getattr(sys, "stdlib_module_names", set())
    if module in stdlib:
        return True
    return module in sys.builtin_module_names


def module_available(module: str) -> bool:
    try:
        return importlib.util.find_spec(module) is not None
    except (ImportError, AttributeError, ValueError):
        return False


def dependency_rows(sources: list[dict], optional_modules: set[str]) -> tuple[list[dict], list[dict]]:
    module_sources: dict[str, set[str]] = {}
    source_rows: list[dict] = []
    indirect_rows: list[dict] = []
    for source in sources:
        modules = imported_modules_from_source(source["text"])
        for module in modules:
            module_sources.setdefault(module, set()).add(source["source_id"])
        for indirect in indirect_modules_from_source(source["text"]):
            indirect_rows.append({**indirect, "source_id": source["source_id"]})
            module_sources.setdefault(indirect["module"], set()).add(source["source_id"])
        source_rows.append(
            {
                "source_id": source["source_id"],
                "source_kind": source["source_kind"],
                "n_imported_modules": len(modules),
                "imported_modules": ";".join(sorted(modules)),
                "indirect_modules": ";".join(
                    sorted({row["module"] for row in indirect_rows if row["source_id"] == source["source_id"]})
                ),
            }
        )

    rows: list[dict] = []
    for module in sorted(module_sources):
        package = PACKAGE_BY_MODULE.get(module, module)
        optional = module in optional_modules
        stdlib = is_stdlib_module(module)
        available = True if stdlib else module_available(module)
        status = "ok"
        if stdlib:
            status = "stdlib"
        elif optional and not available:
            status = "optional_missing"
        elif not available:
            status = "missing"
        rows.append(
            {
                "module": module,
                "package": package,
                "status": status,
                "available": available,
                "optional": optional,
                "source_ids": ";".join(sorted(module_sources[module])),
                "source_count": len(module_sources[module]),
            }
        )
    return rows, source_rows


def install_command(missing_rows: list[dict]) -> str:
    packages = sorted({row["package"] for row in missing_rows if row.get("package")})
    if not packages:
        return ""
    return "python -m pip install " + " ".join(packages)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit Python dependencies needed by the deposited GSE302339 Scanpy "
            "workflow before rerunning long preprocessing/annotation jobs."
        )
    )
    parser.add_argument("--code-dir", default=DEFAULT_CODE_DIR)
    parser.add_argument("--zip", default=DEFAULT_NOTEBOOK_ZIP)
    parser.add_argument(
        "--scan-zip",
        action="store_true",
        help="Scan notebooks inside the Zenodo zip in addition to exported .py files.",
    )
    parser.add_argument(
        "--optional-module",
        action="append",
        dest="optional_modules",
        help=(
            "Module to treat as optional if missing. Can be repeated. "
            "Defaults to scvi and mudata, which are unused/skipped in our exported preprocessing."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dirs = ensure_results_dirs(load_config())
    code_dir = project_path(args.code_dir)
    zip_path = project_path(args.zip)
    optional_modules = set(args.optional_modules or DEFAULT_OPTIONAL_MODULES)

    sources = exported_code_sources(code_dir)
    if args.scan_zip or not sources:
        sources.extend(notebook_sources(zip_path))

    rows, source_rows = dependency_rows(sources, optional_modules)
    missing_required = [row for row in rows if row["status"] == "missing"]
    optional_missing = [row for row in rows if row["status"] == "optional_missing"]
    command = install_command(missing_required)

    dependency_path = results_dirs["tables"] / "gse302339_author_dependency_manifest.csv"
    source_path = results_dirs["tables"] / "gse302339_author_dependency_source_manifest.csv"
    summary_path = results_dirs["meta"] / "gse302339_author_dependency_audit_summary.json"
    write_csv(
        dependency_path,
        rows,
        preferred=[
            "module",
            "package",
            "status",
            "available",
            "optional",
            "source_count",
            "source_ids",
        ],
    )
    write_csv(
        source_path,
        source_rows,
        preferred=["source_id", "source_kind", "n_imported_modules", "imported_modules", "indirect_modules"],
    )

    summary = {
        "code_dir": str(code_dir),
        "code_dir_exists": code_dir.exists(),
        "zip_path": str(zip_path),
        "zip_exists": zip_path.exists(),
        "scan_zip": args.scan_zip,
        "n_sources_scanned": len(sources),
        "n_dependency_modules": len(rows),
        "missing_required_modules": [
            {"module": row["module"], "package": row["package"]} for row in missing_required
        ],
        "optional_missing_modules": [
            {"module": row["module"], "package": row["package"]} for row in optional_missing
        ],
        "install_missing_required_command": command,
        "ready_for_author_workflow": not missing_required and bool(sources),
        "dependency_manifest": str(dependency_path),
        "source_manifest": str(source_path),
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(dependency_path)
    print(source_path)
    if command:
        print()
        print(command)

    if args.strict and missing_required:
        raise SystemExit(
            "Strict dependency audit failed: missing required modules "
            + ", ".join(f"{row['module']} ({row['package']})" for row in missing_required)
        )
    if args.strict and not sources:
        raise SystemExit("Strict dependency audit failed: no notebook/code sources found")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
