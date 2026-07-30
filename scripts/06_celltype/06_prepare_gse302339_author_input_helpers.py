#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import zipfile
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_DIR = "input/data_cellranger8"
DEFAULT_RIBOSOMAL_OUTPUT = "input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv"
DEFAULT_META_OUTPUT = "input/meta_cr8.csv"
DEFAULT_SCANPY_ZIP = "data/raw/downloads/zenodo/16341197/files/scanpy_workflow.zip"
DEFAULT_PREFIXES = ("RPL", "RPS", "MRPL", "MRPS")
SIDECAR_EXTENSIONS = {".csv", ".tsv", ".txt", ".xlsx", ".xls", ".pkl", ".pickle"}
READ_CALL_PATTERNS = [
    (
        "pandas_read",
        re.compile(
            r"""(?:pd\.)?(read_csv|read_table|read_excel)\s*\(\s*["']([^"']+)["']""",
            re.IGNORECASE,
        ),
    ),
    (
        "scanpy_read",
        re.compile(
            r"""(?:sc|ad|anndata|mu)\.(read(?:_h5ad|_10x_h5|_10x_mtx)?)\s*\(\s*["']([^"']+)["']""",
            re.IGNORECASE,
        ),
    ),
]
OPEN_PATTERN = re.compile(
    r"""open\s*\(\s*["']([^"']+)["']\s*,\s*["']([^"']+)["']""",
    re.IGNORECASE,
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


def decode_feature_names(values: Iterable[object]) -> list[str]:
    names: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            names.append(value.decode())
        else:
            names.append(str(value))
    return names


def read_10x_h5_feature_names(path: Path) -> list[str]:
    try:
        import h5py
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "h5py is required to inspect Cell Ranger h5 feature names. "
            "Install h5py in the active environment or run inside spatial_omics."
        ) from exc

    with h5py.File(path, "r") as handle:
        candidates = [
            "matrix/features/name",
            "matrix/features/gene_names",
            "matrix/gene_names",
        ]
        for key in candidates:
            if key in handle:
                return decode_feature_names(handle[key][:])
        available = []
        handle.visit(available.append)
        raise KeyError(
            f"Could not find a Cell Ranger feature-name dataset in {path}. "
            f"Tried {candidates}. First available h5 keys: {available[:20]}"
        )


def ribosomal_gene_pattern(prefixes: tuple[str, ...]) -> re.Pattern:
    escaped = "|".join(re.escape(prefix.upper()) for prefix in prefixes)
    return re.compile(rf"^({escaped})[A-Z0-9-]*$")


def collect_ribosomal_genes(
    h5_paths: list[Path],
    prefixes: tuple[str, ...],
) -> tuple[set[str], list[dict]]:
    pattern = ribosomal_gene_pattern(prefixes)
    all_genes: set[str] = set()
    rows: list[dict] = []
    for h5_path in h5_paths:
        names = read_10x_h5_feature_names(h5_path)
        ribosomal = sorted({name for name in names if pattern.match(name.upper())})
        all_genes.update(ribosomal)
        rows.append(
            {
                "path": str(h5_path),
                "filename": h5_path.name,
                "n_features": len(names),
                "n_unique_features": len(set(names)),
                "n_duplicate_feature_names": len(names) - len(set(names)),
                "n_ribosomal_prefix_genes": len(ribosomal),
                "first_ribosomal_genes": ";".join(ribosomal[:20]),
            }
        )
    return all_genes, rows


def normalize_notebook_path(path: str) -> str:
    return path.replace("\\", "/").lstrip("./")


def is_literal_relative_path(path: str) -> bool:
    return not any(token in path for token in ("{", "}", "*", "?")) and "://" not in path


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


def line_for_offset(text: str, offset: int) -> str:
    running = 0
    for line in text.splitlines():
        next_running = running + len(line) + 1
        if running <= offset < next_running:
            return line.strip()[:500]
        running = next_running
    return ""


def discover_notebook_file_references(zip_path: Path) -> list[dict]:
    if not zip_path.exists():
        return []

    rows: list[dict] = []
    for notebook in notebook_members(zip_path):
        notebook_json = read_notebook(zip_path, notebook)
        for cell_index, cell in enumerate(notebook_json.get("cells", [])):
            if cell.get("cell_type") != "code":
                continue
            text = cell_source(cell)
            for family, pattern in READ_CALL_PATTERNS:
                for match in pattern.finditer(text):
                    operation = match.group(1)
                    path = match.group(2)
                    rows.append(
                        {
                            "notebook": notebook,
                            "cell_index": cell_index,
                            "operation_family": family,
                            "operation": operation,
                            "path": path,
                            "normalized_path": normalize_notebook_path(path),
                            "direction": "read",
                            "line": line_for_offset(text, match.start()),
                        }
                    )
            for match in OPEN_PATTERN.finditer(text):
                path = match.group(1)
                mode = match.group(2)
                folded_mode = mode.casefold()
                direction = "unknown"
                if any(flag in folded_mode for flag in ("w", "a", "x")):
                    direction = "write"
                elif "r" in folded_mode:
                    direction = "read"
                rows.append(
                    {
                        "notebook": notebook,
                        "cell_index": cell_index,
                        "operation_family": "open",
                        "operation": f"open:{mode}",
                        "path": path,
                        "normalized_path": normalize_notebook_path(path),
                        "direction": direction,
                        "line": line_for_offset(text, match.start()),
                    }
                )
    return rows


def is_input_sidecar(path: str) -> bool:
    normalized = normalize_notebook_path(path)
    if not normalized.startswith("input/"):
        return False
    if normalized.startswith("input/data_cellranger8/"):
        return False
    return Path(normalized).suffix.casefold() in SIDECAR_EXTENSIONS


def unique_read_sidecars(reference_rows: list[dict]) -> list[str]:
    return sorted(
        {
            str(row["normalized_path"])
            for row in reference_rows
            if row.get("direction") == "read"
            and is_input_sidecar(str(row.get("normalized_path", "")))
            and is_literal_relative_path(str(row.get("normalized_path", "")))
        }
    )


def find_zip_member_for_path(zip_path: Path, path: str) -> str | None:
    if not zip_path.exists():
        return None
    normalized_path = normalize_notebook_path(path)
    basename = Path(normalized_path).name
    with zipfile.ZipFile(zip_path) as archive:
        names = archive.namelist()
        for member in names:
            if normalize_notebook_path(member).endswith(normalized_path):
                return member
        for member in archive.namelist():
            if Path(member).name == basename:
                return member
    return None


def find_zip_member(zip_path: Path, filename: str) -> str | None:
    return find_zip_member_for_path(zip_path, filename)


def extract_zip_member(zip_path: Path, member: str, destination: Path) -> int:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path) as archive:
        data = archive.read(member)
    destination.write_bytes(data)
    return len(data)


def prepare_zip_sidecars(zip_path: Path, sidecar_paths: list[str]) -> tuple[list[dict], list[dict]]:
    prepared: list[dict] = []
    unresolved: list[dict] = []
    for sidecar in sidecar_paths:
        destination = project_path(sidecar)
        if destination.exists():
            prepared.append(
                {
                    "path": sidecar,
                    "destination": str(destination),
                    "status": "already_exists",
                    "source": "local",
                }
            )
            continue
        member = find_zip_member_for_path(zip_path, sidecar)
        if member:
            n_bytes = extract_zip_member(zip_path, member, destination)
            prepared.append(
                {
                    "path": sidecar,
                    "destination": str(destination),
                    "status": "extracted",
                    "source": str(zip_path),
                    "zip_member": member,
                    "bytes": n_bytes,
                }
            )
            continue
        unresolved.append(
            {
                "path": sidecar,
                "destination": str(destination),
                "status": "missing",
                "source": "not_found_in_scanpy_zip",
            }
        )
    return prepared, unresolved


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare small input helper files expected by the deposited GSE302339 "
            "Scanpy notebooks."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=DEFAULT_INPUT_DIR,
        help="Directory containing extracted Cell Ranger *_filtered_feature_bc_matrix.h5 files.",
    )
    parser.add_argument(
        "--ribosomal-output",
        default=DEFAULT_RIBOSOMAL_OUTPUT,
        help=(
            "Path to write the one-gene-per-line ribosomal gene list expected by "
            "the author preprocessing notebook."
        ),
    )
    parser.add_argument(
        "--meta-output",
        default=DEFAULT_META_OUTPUT,
        help="Path to write the meta_cr8.csv sidecar expected by the author notebook.",
    )
    parser.add_argument(
        "--scanpy-zip",
        default=DEFAULT_SCANPY_ZIP,
        help=(
            "Zenodo scanpy_workflow.zip path. If it contains meta_cr8.csv, the "
            "file is extracted into --meta-output."
        ),
    )
    parser.add_argument(
        "--prefix",
        action="append",
        dest="prefixes",
        help=(
            "Ribosomal gene-symbol prefix to include. Can be repeated. "
            "Defaults to RPL, RPS, MRPL, and MRPS."
        ),
    )
    parser.add_argument("--strict", action="store_true", help="Exit non-zero if helpers are not ready.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config()
    results_dirs = ensure_results_dirs(config)

    input_dir = project_path(args.input_dir)
    output_path = project_path(args.ribosomal_output)
    meta_output_path = project_path(args.meta_output)
    scanpy_zip_path = project_path(args.scanpy_zip)
    prefixes = tuple(args.prefixes or DEFAULT_PREFIXES)

    h5_paths = sorted(input_dir.glob("*filtered_feature_bc_matrix.h5"))
    failures: list[str] = []
    if not input_dir.exists():
        failures.append(f"Input directory missing: {input_dir}")
    if not h5_paths:
        failures.append(f"No *_filtered_feature_bc_matrix.h5 files found in {input_dir}")

    ribosomal_genes: set[str] = set()
    manifest_rows: list[dict] = []
    if h5_paths:
        ribosomal_genes, manifest_rows = collect_ribosomal_genes(h5_paths, prefixes)
        if not ribosomal_genes:
            failures.append(
                f"No ribosomal genes found with prefixes {list(prefixes)} in {input_dir}"
            )

    if ribosomal_genes:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(sorted(ribosomal_genes)) + "\n", encoding="utf-8")

    meta_member = find_zip_member(scanpy_zip_path, Path(args.meta_output).name)
    meta_extracted_bytes: int | None = None
    if meta_member:
        meta_extracted_bytes = extract_zip_member(
            scanpy_zip_path, meta_member, meta_output_path
        )
    if not meta_output_path.exists():
        failures.append(
            f"Author metadata helper missing: {meta_output_path}. "
            f"No {Path(args.meta_output).name} member found in {scanpy_zip_path}."
        )

    manifest_path = results_dirs["tables"] / "gse302339_ribosomal_gene_helper_manifest.csv"
    sidecar_manifest_path = (
        results_dirs["tables"] / "gse302339_author_notebook_sidecar_manifest.csv"
    )
    summary_path = results_dirs["meta"] / "gse302339_author_input_helper_summary.json"

    notebook_reference_rows = discover_notebook_file_references(scanpy_zip_path)
    read_sidecars = unique_read_sidecars(notebook_reference_rows)
    prepared_sidecars, unresolved_sidecars = prepare_zip_sidecars(
        scanpy_zip_path,
        [
            sidecar
            for sidecar in read_sidecars
            if normalize_notebook_path(sidecar)
            != normalize_notebook_path(str(Path(args.ribosomal_output)))
        ],
    )
    if unresolved_sidecars:
        failures.append(
            "Unresolved author notebook sidecars: "
            + ", ".join(row["path"] for row in unresolved_sidecars)
        )

    write_csv(
        manifest_path,
        manifest_rows,
        preferred=[
            "filename",
            "path",
            "n_features",
            "n_unique_features",
            "n_duplicate_feature_names",
            "n_ribosomal_prefix_genes",
            "first_ribosomal_genes",
        ],
    )
    write_csv(
        sidecar_manifest_path,
        [
            {
                **row,
                "local_exists": project_path(str(row["normalized_path"])).exists()
                if is_literal_relative_path(str(row["normalized_path"]))
                else "",
                "is_input_sidecar": is_input_sidecar(str(row["normalized_path"])),
            }
            for row in notebook_reference_rows
        ],
        preferred=[
            "notebook",
            "cell_index",
            "direction",
            "operation_family",
            "operation",
            "path",
            "normalized_path",
            "is_input_sidecar",
            "local_exists",
            "line",
        ],
    )

    summary = {
        "input_dir": str(input_dir),
        "input_dir_exists": input_dir.exists(),
        "n_cellranger_h5": len(h5_paths),
        "ribosomal_output": str(output_path),
        "ribosomal_output_exists": output_path.exists(),
        "meta_output": str(meta_output_path),
        "meta_output_exists": meta_output_path.exists(),
        "scanpy_zip": str(scanpy_zip_path),
        "scanpy_zip_exists": scanpy_zip_path.exists(),
        "meta_cr8_zip_member": meta_member,
        "meta_cr8_extracted_bytes": meta_extracted_bytes,
        "n_notebook_file_references": len(notebook_reference_rows),
        "n_required_input_sidecars": len(read_sidecars),
        "required_input_sidecars": read_sidecars,
        "prepared_sidecars": prepared_sidecars,
        "unresolved_sidecars": unresolved_sidecars,
        "prefixes": list(prefixes),
        "n_ribosomal_genes": len(ribosomal_genes),
        "first_ribosomal_genes": sorted(ribosomal_genes)[:30],
        "manifest_path": str(manifest_path),
        "sidecar_manifest_path": str(sidecar_manifest_path),
        "ready_for_author_preprocessing": (
            output_path.exists()
            and bool(ribosomal_genes)
            and meta_output_path.exists()
            and not unresolved_sidecars
        ),
        "failures": failures,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(manifest_path)
    print(sidecar_manifest_path)

    if args.strict and failures:
        raise SystemExit(
            "Strict helper preparation failed: " + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
