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
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_GSE = "GSE292993"
CODE_PATTERN = re.compile(r"\bRTS\d+\b")
DISPLAY_KEYS = (
    "DisplayName",
    "displayName",
    "TargetName",
    "targetName",
    "Name",
    "name",
    "Gene",
    "gene",
    "GeneSymbol",
    "geneSymbol",
    "Symbol",
    "symbol",
)
CLASS_KEYS = (
    "CodeClass",
    "codeClass",
    "CodeClassName",
    "codeClassName",
    "AnalyteType",
    "analyteType",
)


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


def find_pkc_file(config: dict) -> Path | None:
    pkc_dir = configured_path(config, "geomx_pkc_dir")
    candidates = [
        *sorted(pkc_dir.glob("*.pkc.gz")),
        *sorted(pkc_dir.glob("*.pkc")),
        project_path(
            f"data/raw/downloads/geo/{DEFAULT_GSE}/supplementary/"
            f"{DEFAULT_GSE}_Hs_R_NGS_WTA_v1.0.pkc.gz"
        ),
    ]
    return next((path for path in candidates if path.exists()), None)


def dcc_paths(config: dict) -> list[Path]:
    dcc_dir = configured_path(config, "geomx_dcc_dir")
    return sorted([*dcc_dir.rglob("*.dcc"), *dcc_dir.rglob("*.dcc.gz")])


def scalar_text(value: Any) -> str | None:
    if isinstance(value, (str, int, float)):
        return str(value)
    return None


def value_contains_gene(value: Any, gene: str) -> bool:
    if isinstance(value, dict):
        return any(value_contains_gene(item, gene) for item in value.values())
    if isinstance(value, list):
        return any(value_contains_gene(item, gene) for item in value)
    text = scalar_text(value)
    return bool(text and text.casefold() == gene.casefold())


def direct_record_contains_gene(record: dict, gene: str) -> bool:
    for value in record.values():
        if isinstance(value, dict):
            continue
        if isinstance(value, list):
            if any(
                scalar_text(item)
                and scalar_text(item).casefold() == gene.casefold()
                for item in value
                if not isinstance(item, (dict, list))
            ):
                return True
            continue
        text = scalar_text(value)
        if text and text.casefold() == gene.casefold():
            return True
    return False


def code_ids_from_value(value: Any) -> set[str]:
    codes: set[str] = set()
    if isinstance(value, dict):
        for item in value.values():
            codes.update(code_ids_from_value(item))
    elif isinstance(value, list):
        for item in value:
            codes.update(code_ids_from_value(item))
    else:
        text = scalar_text(value)
        if text:
            codes.update(CODE_PATTERN.findall(text))
    return codes


def first_present(record: dict, keys: tuple[str, ...]) -> str | None:
    for key in keys:
        value = record.get(key)
        text = scalar_text(value)
        if text:
            return text
    return None


def iter_dicts(value: Any, path: str = "$"):
    if isinstance(value, dict):
        yield path, value
        for key, item in value.items():
            yield from iter_dicts(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            yield from iter_dicts(item, f"{path}[{index}]")


def parse_pkc_code_map(path: Path | None, gene: str) -> tuple[dict, list[dict]]:
    summary = {
        "path": str(path) if path else None,
        "status": "missing",
        "code_count": 0,
        "target_count": 0,
        "gene": gene,
        "gene_code_count": 0,
    }
    if path is None or not path.exists():
        return summary, []

    with open_text(path) as handle:
        pkc = json.load(handle)

    rows: list[dict] = []
    seen: set[tuple[str, str, str | None]] = set()
    for path_string, record in iter_dicts(pkc):
        code_ids = sorted(code_ids_from_value(record))
        if not code_ids:
            continue
        target = first_present(record, DISPLAY_KEYS)
        code_class = first_present(record, CLASS_KEYS)
        is_gene_record = direct_record_contains_gene(record, gene)
        if not target and not code_class and not is_gene_record:
            continue
        if target is None and is_gene_record:
            target = gene
        for code_id in code_ids:
            key = (code_id, str(target or ""), code_class)
            if key in seen:
                continue
            seen.add(key)
            rows.append(
                {
                    "code_id": code_id,
                    "target": target,
                    "code_class": code_class,
                    "pkc_path": path_string,
                    "is_primary_gene": bool(
                        target and target.casefold() == gene.casefold()
                    )
                    or is_gene_record,
                }
            )

    targets = {row["target"] for row in rows if row.get("target")}
    gene_codes = {row["code_id"] for row in rows if row["is_primary_gene"]}
    summary.update(
        status="ok",
        code_count=len({row["code_id"] for row in rows}),
        target_count=len(targets),
        gene_code_count=len(gene_codes),
        gene_codes=sorted(gene_codes),
    )
    return summary, rows


def parse_numeric(value: str | None) -> int | float | str | None:
    if value is None:
        return None
    text = value.strip().strip('"')
    if not text:
        return None
    try:
        number = float(text)
    except ValueError:
        return text
    if number.is_integer():
        return int(number)
    return number


def numeric_float(value: int | float | str | None) -> float | None:
    if isinstance(value, (int, float)):
        return float(value)
    return None


def safe_fraction(numerator: int | float | str | None, denominator: int | float | str | None) -> float | None:
    top = numeric_float(numerator)
    bottom = numeric_float(denominator)
    if top is None or bottom is None or bottom == 0:
        return None
    return top / bottom


def summarize_code_subset(code_counts: dict[str, int], codes: set[str], prefix: str) -> dict:
    values = [count for code, count in code_counts.items() if code in codes]
    if not values:
        return {
            f"{prefix}_n_codes": 0,
            f"{prefix}_total_counts": 0,
            f"{prefix}_mean_counts": None,
            f"{prefix}_max_counts": None,
        }
    return {
        f"{prefix}_n_codes": len(values),
        f"{prefix}_total_counts": sum(values),
        f"{prefix}_mean_counts": sum(values) / len(values),
        f"{prefix}_max_counts": max(values),
    }


def parse_dcc(
    path: Path,
    gene_codes: set[str],
    negative_codes: set[str] | None = None,
) -> tuple[dict, dict[str, int]]:
    attrs: dict[str, str] = {}
    code_counts: dict[str, int] = {}
    section = None

    with open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            open_tag = re.match(r"^<([^/][^>]*)>$", line)
            close_tag = re.match(r"^</([^>]*)>$", line)
            if open_tag:
                section = open_tag.group(1)
                continue
            if close_tag:
                section = None
                continue
            fields = next(csv.reader([line]))
            if len(fields) < 2:
                continue
            key, value = fields[0].strip(), fields[1].strip()
            if section == "Code_Summary" and CODE_PATTERN.fullmatch(key):
                numeric = parse_numeric(value)
                if isinstance(numeric, int):
                    code_counts[key] = numeric
                continue
            if section in {"Scan_Attributes", "NGS_Processing_Attributes"}:
                attrs[key] = value

    total_counts = sum(code_counts.values())
    gene_counts = sum(count for code, count in code_counts.items() if code in gene_codes)
    raw_reads = parse_numeric(attrs.get("Raw"))
    trimmed_reads = parse_numeric(attrs.get("Trimmed"))
    stitched_reads = parse_numeric(attrs.get("Stitched"))
    aligned_reads = parse_numeric(attrs.get("Aligned"))
    negative_summary = summarize_code_subset(
        code_counts, negative_codes or set(), "negative_probe"
    )
    row = {
        "geo_accession": gsm_from_name(path.name),
        "dcc_filename": path.name,
        "dcc_id": dcc_id_from_name(path.name),
        "scan_id": attrs.get("ID"),
        "plate_id": attrs.get("Plate_ID"),
        "well": attrs.get("Well"),
        "seq_set_id": attrs.get("SeqSetId"),
        "raw_reads": raw_reads,
        "trimmed_reads": trimmed_reads,
        "stitched_reads": stitched_reads,
        "aligned_reads": aligned_reads,
        "trimmed_fraction": safe_fraction(trimmed_reads, raw_reads),
        "stitched_fraction": safe_fraction(stitched_reads, trimmed_reads),
        "aligned_fraction_raw": safe_fraction(aligned_reads, raw_reads),
        "aligned_fraction_trimmed": safe_fraction(aligned_reads, trimmed_reads),
        "aligned_fraction_stitched": safe_fraction(aligned_reads, stitched_reads),
        "umi_q30": parse_numeric(attrs.get("umiQ30")),
        "rts_q30": parse_numeric(attrs.get("rtsQ30")),
        "n_code_counts": len(code_counts),
        "total_code_counts": total_counts,
        "primary_gene_counts": gene_counts,
        **negative_summary,
    }
    return row, code_counts


def summarize_numeric(rows: list[dict], column: str) -> dict:
    values = sorted(
        float(row[column])
        for row in rows
        if isinstance(row.get(column), (int, float))
    )
    if not values:
        return {"n": 0}
    indexes = {
        "min": 0,
        "q25": int((len(values) - 1) * 0.25),
        "median": int((len(values) - 1) * 0.50),
        "q75": int((len(values) - 1) * 0.75),
        "max": len(values) - 1,
    }
    summary = {"n": len(values)}
    summary.update({key: values[index] for key, index in indexes.items()})
    return summary


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
        "scan_id",
        "plate_id",
        "well",
        "raw_reads",
        "trimmed_reads",
        "stitched_reads",
        "aligned_reads",
        "trimmed_fraction",
        "stitched_fraction",
        "aligned_fraction_raw",
        "aligned_fraction_trimmed",
        "aligned_fraction_stitched",
        "umi_q30",
        "rts_q30",
        "n_code_counts",
        "total_code_counts",
        "primary_gene_counts",
        "negative_probe_n_codes",
        "negative_probe_total_counts",
        "negative_probe_mean_counts",
        "negative_probe_max_counts",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Extract GSE292993 DCC QC metrics and PKC code mappings."
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    gene = config["project"]["gene"]
    output = ensure_results_dirs(config)
    paths = dcc_paths(config)
    pkc_summary, pkc_rows = parse_pkc_code_map(find_pkc_file(config), gene)
    gene_codes = {
        row["code_id"]
        for row in pkc_rows
        if row.get("is_primary_gene") and row.get("code_id")
    }
    negative_codes = {
        row["code_id"]
        for row in pkc_rows
        if row.get("code_id")
        and (
            "negative" in str(row.get("code_class") or "").casefold()
            or "negative" in str(row.get("target") or "").casefold()
        )
    }

    qc_rows: list[dict] = []
    n_codes: Counter[int] = Counter()
    for path in paths:
        row, code_counts = parse_dcc(path, gene_codes, negative_codes)
        qc_rows.append(row)
        n_codes[len(code_counts)] += 1

    summary = {
        "dataset": DEFAULT_GSE,
        "dcc_count": len(paths),
        "pkc": pkc_summary,
        "primary_gene": gene,
        "primary_gene_codes_used": sorted(gene_codes),
        "negative_probe_codes_used": sorted(negative_codes),
        "negative_probe_code_count": len(negative_codes),
        "dccs_with_primary_gene_counts_gt0": sum(
            1 for row in qc_rows if row["primary_gene_counts"] > 0
        ),
        "n_code_counts_distribution": dict(sorted(n_codes.items())),
        "qc_summaries": {
            column: summarize_numeric(qc_rows, column)
            for column in (
                "raw_reads",
                "trimmed_reads",
                "stitched_reads",
                "aligned_reads",
                "trimmed_fraction",
                "stitched_fraction",
                "aligned_fraction_raw",
                "aligned_fraction_trimmed",
                "aligned_fraction_stitched",
                "umi_q30",
                "rts_q30",
                "n_code_counts",
                "total_code_counts",
                "primary_gene_counts",
                "negative_probe_n_codes",
                "negative_probe_total_counts",
                "negative_probe_mean_counts",
                "negative_probe_max_counts",
            )
        },
    }

    meta_dir = output["meta"]
    table_dir = output["tables"]
    (meta_dir / "gse292993_dcc_qc_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    write_csv(table_dir / "gse292993_pkc_code_map.csv", pkc_rows)
    write_csv(table_dir / "gse292993_dcc_qc_metrics.csv", qc_rows)

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {meta_dir / 'gse292993_dcc_qc_summary.json'}")
    print(f"Wrote {table_dir / 'gse292993_pkc_code_map.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dcc_qc_metrics.csv'}")

    failures = []
    if not paths:
        failures.append("No DCC files found")
    if pkc_summary["status"] != "ok":
        failures.append("PKC file missing or unreadable")
    if pkc_summary.get("gene_code_count", 0) == 0:
        failures.append(f"No PKC code IDs resolved for {gene}")
    if args.strict and failures:
        print("Strict DCC QC failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
