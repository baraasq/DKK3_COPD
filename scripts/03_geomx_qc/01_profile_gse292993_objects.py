#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_GSE = "GSE292993"


def open_text(path: Path):
    if path.name.casefold().endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def dcc_id_from_name(name: str) -> str:
    filename = Path(name).name
    for suffix in (".dcc.gz", ".DCC.gz", ".dcc", ".DCC"):
        if filename.endswith(suffix):
            return filename[: -len(suffix)]
    return Path(filename).stem


def gsm_from_name(name: str) -> str | None:
    match = re.search(r"(GSM\d+)", Path(name).name)
    return match.group(1) if match else None


def sanitize_key(value: str) -> str:
    cleaned = re.sub(r"[^0-9A-Za-z]+", "_", value.strip().casefold())
    return cleaned.strip("_") or "unknown"


def find_soft_file(config: dict) -> Path | None:
    metadata_dir = configured_path(config, "geomx_metadata_dir")
    candidates = [
        metadata_dir / f"{DEFAULT_GSE}_family.soft.gz",
        project_path(
            f"data/raw/downloads/geo/{DEFAULT_GSE}/metadata/{DEFAULT_GSE}_family.soft.gz"
        ),
    ]
    return next((path for path in candidates if path.exists()), None)


def dcc_paths(config: dict) -> list[Path]:
    dcc_dir = configured_path(config, "geomx_dcc_dir")
    return sorted([*dcc_dir.rglob("*.dcc"), *dcc_dir.rglob("*.dcc.gz")])


def split_tabular_line(line: str) -> tuple[str, list[str]]:
    stripped = line.rstrip("\n\r")
    if "\t" in stripped:
        return "tab", stripped.split("\t")
    if "," in stripped:
        return "comma", next(csv.reader([stripped]))
    return "none", [stripped]


def likely_key_value(line: str) -> tuple[str, str] | None:
    stripped = line.strip()
    if not stripped:
        return None
    if "=" in stripped:
        key, value = stripped.split("=", 1)
    elif "\t" in stripped:
        parts = stripped.split("\t", 1)
        key, value = parts[0], parts[1]
    elif "," in stripped:
        parts = next(csv.reader([stripped]))
        if len(parts) < 2:
            return None
        key, value = parts[0], ",".join(parts[1:])
    else:
        return None

    key = key.strip().strip('"')
    value = value.strip().strip('"')
    if not key or len(key) > 80:
        return None
    if not re.search(r"[A-Za-z]", key):
        return None
    return key, value


def profile_dcc(path: Path, *, preview_lines: int = 25, max_scan_lines: int = 50000) -> dict:
    first_nonempty: list[str] = []
    section_headers: Counter[str] = Counter()
    delimiter_counts: Counter[str] = Counter()
    table_header_candidates: list[dict] = []
    key_values: dict[str, str] = {}
    total_lines = 0

    with open_text(path) as handle:
        for raw_line in handle:
            total_lines += 1
            line = raw_line.rstrip("\n\r")
            stripped = line.strip()
            if stripped and len(first_nonempty) < preview_lines:
                first_nonempty.append(stripped[:500])

            bracket = re.match(r"^\s*\[([^\]]+)\]\s*$", stripped)
            xml = re.match(r"^\s*<([^/!?][^>\s]*)[^>]*>\s*$", stripped)
            if bracket:
                section_headers[bracket.group(1).strip()] += 1
            elif xml:
                section_headers[xml.group(1).strip()] += 1

            delimiter, fields = split_tabular_line(line)
            delimiter_counts[delimiter] += 1
            if delimiter != "none" and len(fields) >= 3 and len(table_header_candidates) < 30:
                if any(re.search(r"[A-Za-z]", field) for field in fields):
                    table_header_candidates.append(
                        {
                            "line_number": total_lines,
                            "delimiter": delimiter,
                            "n_fields": len(fields),
                            "fields": fields[:20],
                        }
                    )

            pair = likely_key_value(line)
            if pair:
                key, value = pair
                lowered = key.casefold()
                if any(
                    token in lowered
                    for token in (
                        "sample",
                        "segment",
                        "roi",
                        "aoi",
                        "slide",
                        "scan",
                        "nuclei",
                        "area",
                        "read",
                        "count",
                        "aligned",
                        "saturation",
                    )
                ):
                    key_values.setdefault(key, value)

            if total_lines >= max_scan_lines:
                break

    return {
        "dcc_file": str(path),
        "dcc_filename": path.name,
        "dcc_id": dcc_id_from_name(path.name),
        "geo_accession": gsm_from_name(path.name),
        "size_bytes": path.stat().st_size,
        "scanned_lines": total_lines,
        "first_nonempty_lines": first_nonempty,
        "section_headers": dict(section_headers),
        "delimiter_counts": dict(delimiter_counts),
        "table_header_candidates": table_header_candidates,
        "key_values": key_values,
    }


def parse_soft_samples(path: Path | None) -> list[dict]:
    if path is None:
        return []

    samples: list[dict] = []
    current: dict[str, str] | None = None
    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n\r")
            if line.startswith("^SAMPLE = "):
                if current:
                    samples.append(current)
                current = {"geo_accession": line.split("=", 1)[1].strip()}
                continue
            if current is None or not line.startswith("!Sample_"):
                continue

            key, _, value = line[1:].partition(" = ")
            key = key.removeprefix("Sample_")
            value = value.strip()
            if key == "characteristics_ch1" and ":" in value:
                char_key, char_value = value.split(":", 1)
                current[f"characteristics_{sanitize_key(char_key)}"] = char_value.strip()
            else:
                column = sanitize_key(key)
                if column in current:
                    suffix = 2
                    while f"{column}_{suffix}" in current:
                        suffix += 1
                    column = f"{column}_{suffix}"
                current[column] = value

    if current:
        samples.append(current)
    return samples


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "geo_accession",
        "dcc_filename",
        "dcc_id",
        "title",
        "source_name_ch1",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def merge_dcc_and_geo(dcc_rows: list[dict], geo_rows: list[dict]) -> list[dict]:
    by_gsm = {row.get("geo_accession"): row for row in geo_rows if row.get("geo_accession")}
    merged = []
    for row in dcc_rows:
        gsm = row.get("geo_accession")
        merged_row = dict(by_gsm.get(gsm, {}))
        merged_row.update(row)
        merged_row["metadata_matched"] = bool(gsm and gsm in by_gsm)
        merged.append(merged_row)
    return merged


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Profile copied GSE292993 DCC objects and GEO sample metadata."
    )
    parser.add_argument("--preview-dccs", type=int, default=3)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    soft_file = find_soft_file(config)
    paths = dcc_paths(config)
    dcc_manifest_rows = [
        {
            "dcc_filename": path.name,
            "dcc_id": dcc_id_from_name(path.name),
            "geo_accession": gsm_from_name(path.name),
            "size_bytes": path.stat().st_size,
        }
        for path in paths
    ]
    geo_rows = parse_soft_samples(soft_file)
    merged_rows = merge_dcc_and_geo(dcc_manifest_rows, geo_rows)
    profiles = [profile_dcc(path) for path in paths[: args.preview_dccs]]

    matched = sum(1 for row in merged_rows if row["metadata_matched"])
    summary = {
        "dataset": DEFAULT_GSE,
        "dcc_dir": str(configured_path(config, "geomx_dcc_dir")),
        "dcc_count": len(paths),
        "geo_soft_file": str(soft_file) if soft_file else None,
        "geo_sample_count": len(geo_rows),
        "dccs_matched_to_geo_samples": matched,
        "dccs_without_geo_sample": len(merged_rows) - matched,
        "previewed_dccs": len(profiles),
        "first_dcc_profile": profiles[0] if profiles else None,
    }

    meta_dir = output["meta"]
    table_dir = output["tables"]
    (meta_dir / "gse292993_object_profile.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(table_dir / "gse292993_geo_sample_metadata.csv", geo_rows)
    write_csv(table_dir / "gse292993_roi_metadata_initial.csv", merged_rows)

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {meta_dir / 'gse292993_object_profile.json'}")
    print(f"Wrote {table_dir / 'gse292993_geo_sample_metadata.csv'}")
    print(f"Wrote {table_dir / 'gse292993_roi_metadata_initial.csv'}")

    failures = []
    if not paths:
        failures.append("No copied DCC files found")
    if not geo_rows:
        failures.append("No GEO sample metadata parsed from SOFT file")
    if paths and matched != len(paths):
        failures.append(f"Only {matched}/{len(paths)} DCC files matched GEO samples")
    if args.strict and failures:
        print("Strict profile failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
