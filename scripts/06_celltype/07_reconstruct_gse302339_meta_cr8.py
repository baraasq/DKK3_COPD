#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
import urllib.request
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_DIR = "input/data_cellranger8"
DEFAULT_OUTPUT = "input/meta_cr8.csv"
DEFAULT_CODE_GLOB = "intermediate/gse302339_scanpy_workflow_code/*.py"
DEFAULT_SOFT_OUTPUT = (
    "data/raw/downloads/geo/GSE302339/metadata/GSE302339_family.soft.gz"
)
DEFAULT_SOFT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE302nnn/"
    "GSE302339/soft/GSE302339_family.soft.gz"
)
DEFAULT_SOFT_CANDIDATES = [
    DEFAULT_SOFT_OUTPUT,
    "data/raw/gse302339/geo_metadata/GSE302339_family.soft.gz",
    "data/raw/downloads/geo/GSE302339/metadata/GSE302339_family.soft.gz",
    "data/raw/downloads/geo/GSE302339/GSE302339_family.soft.gz",
]

H5_FILENAME_PATTERN = re.compile(
    r"^(?P<geo_accession>GSM\d+)_(?P<sample>.+?)_filtered_feature_bc_matrix\.h5$"
)
META_SINGLE_COLUMN_PATTERN = re.compile(
    r"""meta\s*\[\s*["']([^"']+)["']\s*\]"""
)
META_LIST_PATTERN = re.compile(
    r"""meta\s*\[\s*\[([^\]]+)\]\s*\]""",
    re.DOTALL,
)
QUOTED_PATTERN = re.compile(r"""["']([^"']+)["']""")


def sanitize_column(value: str) -> str:
    text = value.strip().casefold()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return text.strip("_")


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


def read_csv_rows(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def open_text(path: Path):
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def first_existing(paths: Iterable[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def parse_h5_samples(input_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(input_dir.glob("*filtered_feature_bc_matrix.h5")):
        match = H5_FILENAME_PATTERN.match(path.name)
        sample = path.stem.replace("_filtered_feature_bc_matrix", "")
        geo_accession = ""
        if match:
            geo_accession = match.group("geo_accession")
            sample = match.group("sample")
        row = {
            "sample": sample,
            "Sample": sample,
            "sample_id": sample,
            "sample_name": sample,
            "orig.ident": sample,
            "orig_ident": sample,
            "batch": sample,
            "geo_accession": geo_accession,
            "gsm": geo_accession,
            "cellranger_h5": str(path),
            "cellranger_h5_filename": path.name,
        }
        number_match = re.search(r"(\d+)$", sample)
        if number_match:
            row["sample_number"] = number_match.group(1)
        rows.append(row)
    return rows


def parse_soft_samples(path: Path) -> dict[str, dict]:
    samples: dict[str, dict] = {}
    current: dict | None = None
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current and current.get("geo_accession"):
                    samples[str(current["geo_accession"])] = current
                accession = line.split(" = ", 1)[1].strip()
                current = {"geo_accession": accession, "gsm": accession}
                continue
            if current is None or not line.startswith("!"):
                continue
            if " = " not in line:
                continue
            key, value = line[1:].split(" = ", 1)
            key = key.strip()
            value = value.strip()
            column = sanitize_column(key)
            if column.startswith("sample_"):
                column = column.removeprefix("sample_")
            if key == "Sample_characteristics_ch1" and ":" in value:
                characteristic_key, characteristic_value = value.split(":", 1)
                characteristic_column = "characteristics_" + sanitize_column(
                    characteristic_key
                )
                current[characteristic_column] = characteristic_value.strip()
            else:
                if column in current and current[column] != value:
                    suffix = 2
                    while f"{column}_{suffix}" in current:
                        suffix += 1
                    column = f"{column}_{suffix}"
                current[column] = value
    if current and current.get("geo_accession"):
        samples[str(current["geo_accession"])] = current
    return samples


def add_alias_columns(row: dict) -> dict:
    output = dict(row)
    title = output.get("title") or output.get("sample_title") or output.get("sample")
    if title:
        output.setdefault("title", title)
    condition = (
        output.get("condition")
        or output.get("characteristics_condition")
        or output.get("characteristics_disease")
        or output.get("characteristics_disease_state")
        or output.get("characteristics_status")
        or output.get("characteristics_group")
    )
    if condition:
        output.setdefault("condition", condition)
        output.setdefault("disease", condition)
        output.setdefault("diagnosis", condition)
        output.setdefault("group", condition)
    patient = (
        output.get("patient")
        or output.get("characteristics_patient")
        or output.get("characteristics_donor")
        or output.get("characteristics_individual")
        or output.get("characteristics_subject_identity")
        or output.get("characteristics_lab_id")
    )
    if patient:
        output.setdefault("patient", patient)
        output.setdefault("donor", patient)
        output.setdefault("donor_id", patient)
    lobe_emphysema = (
        output.get("lobe_emphysema_simple")
        or output.get("characteristics_lobe_emphysema")
        or output.get("characteristics_lobe_emphysema_severity_laa")
    )
    if lobe_emphysema:
        output.setdefault("lobe_emphysema_simple", lobe_emphysema)
        output.setdefault("lobe_emphysema", lobe_emphysema)
    total_emphysema = (
        output.get("total_emphysema")
        or output.get("characteristics_subject_group")
        or output.get("characteristics_copd")
        or output.get("condition")
        or lobe_emphysema
    )
    if total_emphysema:
        output.setdefault("total_emphysema", total_emphysema)
    return output


def sample_batch_aliases(row: dict, *, row_index: int | None = None) -> list[str]:
    aliases: list[str] = []
    for key in [
        "sample",
        "sample_id",
        "sample_name",
        "orig.ident",
        "orig_ident",
        "batch",
        "geo_accession",
        "gsm",
        "cellranger_h5_filename",
    ]:
        value = str(row.get(key, "")).strip()
        if value and value not in aliases:
            aliases.append(value)

    filename = str(row.get("cellranger_h5_filename", "")).strip()
    if filename:
        stem = filename.replace("_filtered_feature_bc_matrix.h5", "")
        for value in [
            stem,
            filename.removesuffix(".h5"),
            f"input/data_cellranger8/{stem}",
            f"input/data_cellranger8/{filename}",
        ]:
            if value and value not in aliases:
                aliases.append(value)

    sample_number = str(row.get("sample_number", "")).strip()
    if sample_number:
        for value in [sample_number, f"s{sample_number}"]:
            if value and value not in aliases:
                aliases.append(value)
    if row_index is not None:
        for value in [str(row_index), f"{row_index}"]:
            if value and value not in aliases:
                aliases.append(value)

    return aliases


def expand_sample_name_alias_rows(rows: list[dict]) -> list[dict]:
    alias_counts: dict[str, int] = {}
    row_aliases: list[tuple[dict, list[str]]] = []
    for row_index, row in enumerate(rows):
        aliases = sample_batch_aliases(row, row_index=row_index)
        row_aliases.append((row, aliases))
        for alias in set(aliases):
            alias_counts[alias] = alias_counts.get(alias, 0) + 1

    expanded: list[dict] = []
    seen: set[str] = set()
    for row, aliases in row_aliases:
        canonical = str(row.get("sample", row.get("sample_name", ""))).strip()
        for alias in aliases:
            if alias_counts.get(alias, 0) != 1:
                continue
            if alias in seen:
                continue
            seen.add(alias)
            alias_row = dict(row)
            alias_row["batch"] = alias
            alias_row["batch_alias"] = alias
            alias_row["batch_canonical"] = canonical
            alias_row["sample_name"] = alias
            alias_row["sample_name_alias"] = alias
            alias_row["sample_name_canonical"] = canonical
            expanded.append(alias_row)
    return expanded


def infer_expected_meta_columns(code_paths: list[Path]) -> list[str]:
    columns: set[str] = set()
    for path in code_paths:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8", errors="replace")
        if "meta_cr8.csv" not in text and "meta" not in text:
            continue
        columns.update(META_SINGLE_COLUMN_PATTERN.findall(text))
        for match in META_LIST_PATTERN.findall(text):
            columns.update(QUOTED_PATTERN.findall(match))
    return sorted(columns)


def merge_metadata(h5_rows: list[dict], soft_rows: dict[str, dict]) -> list[dict]:
    merged: list[dict] = []
    for h5_row in h5_rows:
        geo_accession = str(h5_row.get("geo_accession", ""))
        row = dict(h5_row)
        if geo_accession and geo_accession in soft_rows:
            for key, value in soft_rows[geo_accession].items():
                row.setdefault(key, value)
        merged.append(add_alias_columns(row))
    return merged


def duplicated_values(rows: list[dict], column: str) -> list[str]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row.get(column, "")).strip()
        if value:
            counts[value] = counts.get(value, 0) + 1
    return sorted(value for value, count in counts.items() if count > 1)


def download_soft(url: str, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as response:
        destination.write_bytes(response.read())
    return destination


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct the missing input/meta_cr8.csv sidecar expected by the "
            "deposited GSE302339 Scanpy preprocessing notebook."
        )
    )
    parser.add_argument("--input-dir", default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output", default=DEFAULT_OUTPUT)
    parser.add_argument("--code-glob", default=DEFAULT_CODE_GLOB)
    parser.add_argument("--soft", action="append", dest="soft_candidates")
    parser.add_argument("--soft-output", default=DEFAULT_SOFT_OUTPUT)
    parser.add_argument("--soft-url", default=DEFAULT_SOFT_URL)
    parser.add_argument(
        "--download-soft",
        action="store_true",
        help="Download GSE302339 family SOFT from NCBI if no local SOFT file is found.",
    )
    parser.add_argument(
        "--allow-minimal-without-soft",
        action="store_true",
        help=(
            "Allow a filename-derived meta_cr8.csv if GEO SOFT metadata is absent. "
            "This is useful only when the notebook needs sample IDs, not phenotype labels."
        ),
    )
    parser.add_argument(
        "--no-sample-name-alias-rows",
        action="store_true",
        help=(
            "Do not expand meta_cr8.csv with additional rows for plausible sample_name "
            "forms such as GSM accession, h5 filename, and sample stem."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dirs = ensure_results_dirs(load_config())

    input_dir = project_path(args.input_dir)
    output_path = project_path(args.output)
    soft_output = project_path(args.soft_output)
    soft_candidates = [
        project_path(path)
        for path in (args.soft_candidates or DEFAULT_SOFT_CANDIDATES)
    ]
    soft_path = first_existing(soft_candidates)
    failures: list[str] = []

    if soft_path is None and args.download_soft:
        try:
            soft_path = download_soft(args.soft_url, soft_output)
        except Exception as exc:  # pragma: no cover - network path
            failures.append(f"SOFT download failed from {args.soft_url}: {exc}")

    h5_rows = parse_h5_samples(input_dir)
    if not input_dir.exists():
        failures.append(f"Input directory missing: {input_dir}")
    if not h5_rows:
        failures.append(f"No Cell Ranger h5 files found in {input_dir}")

    soft_rows: dict[str, dict] = {}
    if soft_path:
        soft_rows = parse_soft_samples(soft_path)
    elif not args.allow_minimal_without_soft:
        failures.append(
            "GSE302339 family SOFT metadata is missing. Download it or rerun with "
            "--download-soft. Use --allow-minimal-without-soft only for a sample-ID-only scaffold."
        )

    base_output_rows = merge_metadata(h5_rows, soft_rows) if h5_rows else []
    output_rows = (
        base_output_rows
        if args.no_sample_name_alias_rows
        else expand_sample_name_alias_rows(base_output_rows)
    )
    output_columns = sorted({key for row in output_rows for key in row})
    duplicated_batch_values = duplicated_values(output_rows, "batch")
    if duplicated_batch_values:
        failures.append(
            "Reconstructed meta_cr8.csv has duplicated batch keys; this would "
            f"make pd.merge(old_meta, meta, on='batch') unsafe. First duplicates: "
            f"{duplicated_batch_values[:20]}"
        )

    code_paths = sorted(project_path(".").glob(args.code_glob))
    expected_columns = infer_expected_meta_columns(code_paths)
    missing_expected_columns = [
        column for column in expected_columns if column not in output_columns
    ]
    if missing_expected_columns:
        failures.append(
            "Reconstructed meta_cr8.csv is missing columns referenced by exported "
            f"author code: {missing_expected_columns}"
        )

    matched_soft = sum(
        1
        for row in output_rows
        if row.get("geo_accession") and row.get("geo_accession") in soft_rows
    )
    if soft_rows and matched_soft != len(output_rows):
        failures.append(
            f"Only {matched_soft}/{len(output_rows)} h5 files matched GEO SOFT samples"
        )

    manifest_path = results_dirs["tables"] / "gse302339_meta_cr8_reconstruction_manifest.csv"
    summary_path = results_dirs["meta"] / "gse302339_meta_cr8_reconstruction_summary.json"
    if output_rows and (not failures or not args.strict):
        write_csv(
            output_path,
            output_rows,
            preferred=[
                "sample",
                "Sample",
                "sample_id",
                "orig.ident",
                "orig_ident",
                "batch",
                "geo_accession",
                "gsm",
                "condition",
                "disease",
                "diagnosis",
                "group",
                "patient",
                "donor",
                "donor_id",
                "title",
                "cellranger_h5_filename",
                "cellranger_h5",
            ],
        )
    write_csv(
        manifest_path,
        output_rows,
        preferred=[
            "sample",
            "geo_accession",
            "condition",
            "patient",
            "title",
            "cellranger_h5_filename",
        ],
    )

    summary = {
        "input_dir": str(input_dir),
        "input_dir_exists": input_dir.exists(),
        "n_cellranger_h5": len(h5_rows),
        "n_base_meta_rows": len(base_output_rows),
        "sample_name_alias_rows_enabled": not args.no_sample_name_alias_rows,
        "n_meta_rows_written": len(output_rows),
        "n_unique_batch_values": len(
            {str(row.get("batch", "")).strip() for row in output_rows if row.get("batch")}
        ),
        "duplicated_batch_values": duplicated_batch_values,
        "output": str(output_path),
        "output_exists": output_path.exists(),
        "soft_path": str(soft_path) if soft_path else None,
        "soft_candidates": [str(path) for path in soft_candidates],
        "soft_url": args.soft_url,
        "n_soft_samples": len(soft_rows),
        "n_h5_rows_matched_to_soft": matched_soft,
        "code_paths": [str(path) for path in code_paths],
        "expected_meta_columns_from_code": expected_columns,
        "missing_expected_meta_columns": missing_expected_columns,
        "output_columns": output_columns,
        "manifest_path": str(manifest_path),
        "ready_for_author_preprocessing": output_path.exists()
        and not missing_expected_columns
        and not duplicated_batch_values
        and bool(output_rows),
        "failures": failures,
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(manifest_path)
    if output_path.exists():
        print(output_path)

    if args.strict and failures:
        raise SystemExit(
            "Strict meta_cr8 reconstruction failed: " + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
