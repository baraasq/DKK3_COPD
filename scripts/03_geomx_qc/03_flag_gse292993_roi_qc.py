#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


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


def as_float(row: dict, column: str) -> float | None:
    value = row.get(column)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def first_existing(row: dict, candidates: list[str]) -> str | None:
    for column in candidates:
        value = row.get(column)
        if value not in (None, ""):
            return str(value)
    return None


def normalize_label(value: str | None) -> str:
    if value is None:
        return ""
    return re.sub(r"\s+", " ", str(value).strip())


def disease_group(value: str | None) -> str:
    label = normalize_label(value)
    folded = label.casefold().replace("-", " ")
    if folded == "copd":
        return "COPD"
    if folded in {"non smoker", "nonsmoker", "never smoker", "smoker"}:
        return "Control"
    return label or "unknown"


def smoking_group(row: dict) -> str:
    explicit = normalize_label(first_existing(row, ["characteristics_former_or_current_smoker"]))
    if explicit:
        return explicit
    condition = normalize_label(row.get("diagnosis_guess"))
    folded = condition.casefold().replace("-", " ")
    if folded in {"non smoker", "nonsmoker", "never smoker"}:
        return "Never smoker"
    if folded == "smoker":
        return "Smoker control"
    return "unknown"


def merge_by_key(left: list[dict], right: list[dict], key: str) -> list[dict]:
    right_by_key = {row.get(key): row for row in right if row.get(key)}
    merged = []
    for row in left:
        merged_row = dict(row)
        match = right_by_key.get(row.get(key))
        if match:
            for column, value in match.items():
                if column == key:
                    continue
                if column in merged_row and merged_row[column] not in ("", None):
                    merged_row[f"qc_{column}"] = value
                else:
                    merged_row[column] = value
        merged_row["qc_metrics_matched"] = bool(match)
        merged.append(merged_row)
    return merged


def qc_flags(row: dict, args: argparse.Namespace) -> tuple[bool, list[str]]:
    reasons = []
    checks = [
        ("aligned_reads", args.min_aligned_reads, "aligned_reads_below_min"),
        ("n_code_counts", args.min_code_counts, "n_code_counts_below_min"),
        ("total_code_counts", args.min_total_code_counts, "total_code_counts_below_min"),
        ("trimmed_fraction", args.min_trimmed_fraction, "trimmed_fraction_below_min"),
        ("stitched_fraction", args.min_stitched_fraction, "stitched_fraction_below_min"),
        (
            "aligned_fraction_stitched",
            args.min_aligned_fraction_stitched,
            "aligned_fraction_stitched_below_min",
        ),
        ("umi_q30", args.min_umi_q30, "umi_q30_below_min"),
        ("rts_q30", args.min_rts_q30, "rts_q30_below_min"),
    ]
    for column, minimum, reason in checks:
        value = as_float(row, column)
        if value is None:
            reasons.append(f"{column}_missing")
        elif value < minimum:
            reasons.append(reason)
    if not row.get("metadata_matched") in (True, "True", "true", "1", 1):
        reasons.append("geo_metadata_unmatched")
    if not row.get("qc_metrics_matched") in (True, "True", "true", "1", 1):
        reasons.append("dcc_qc_unmatched")
    return not reasons, reasons


def summarize_counts(rows: list[dict], columns: list[str]) -> list[dict]:
    counter: Counter[tuple[str, ...]] = Counter()
    for row in rows:
        key = tuple(normalize_label(row.get(column)) or "unknown" for column in columns)
        counter[key] += 1
    output = []
    for key, count in sorted(counter.items()):
        item = dict(zip(columns, key))
        item["n_rois"] = count
        output.append(item)
    return output


def summarize_qc_by_group(rows: list[dict], group_columns: list[str]) -> list[dict]:
    groups: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    for row in rows:
        key = tuple(
            normalize_label(first_existing(row, [column])) or "unknown"
            for column in group_columns
        )
        groups[key].append(row)

    summaries = []
    for key, group_rows in sorted(groups.items()):
        passed = sum(row.get("include_qc") == "True" for row in group_rows)
        item = dict(zip(group_columns, key))
        item["n_rois"] = len(group_rows)
        item["n_pass_qc"] = passed
        item["n_fail_qc"] = len(group_rows) - passed
        summaries.append(item)
    return summaries


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Merge GSE292993 ROI metadata with DCC QC metrics and flag low-quality ROIs."
    )
    parser.add_argument("--min-aligned-reads", type=float, default=100000)
    parser.add_argument("--min-code-counts", type=float, default=10000)
    parser.add_argument("--min-total-code-counts", type=float, default=10000)
    parser.add_argument("--min-trimmed-fraction", type=float, default=0.90)
    parser.add_argument("--min-stitched-fraction", type=float, default=0.80)
    parser.add_argument("--min-aligned-fraction-stitched", type=float, default=0.80)
    parser.add_argument("--min-umi-q30", type=float, default=0.98)
    parser.add_argument("--min-rts-q30", type=float, default=0.98)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    roi_metadata_path = table_dir / "gse292993_roi_metadata_initial.csv"
    dcc_qc_path = table_dir / "gse292993_dcc_qc_metrics.csv"

    roi_rows = read_csv(roi_metadata_path)
    qc_rows = read_csv(dcc_qc_path)
    rows = merge_by_key(roi_rows, qc_rows, "geo_accession")

    for row in rows:
        include, reasons = qc_flags(row, args)
        row["include_qc"] = str(include)
        row["exclusion_reason"] = ";".join(reasons)
        row["diagnosis_guess"] = normalize_label(
            first_existing(
                row,
                [
                    "characteristics_diagnosis",
                    "characteristics_disease",
                    "characteristics_condition",
                    "diagnosis",
                    "disease",
                ],
            )
        )
        row["diagnosis_group"] = disease_group(row["diagnosis_guess"])
        row["smoking_status_guess"] = smoking_group(row)
        row["compartment_guess"] = normalize_label(
            first_existing(
                row,
                [
                    "characteristics_compartment",
                    "characteristics_region",
                    "characteristics_segment",
                    "characteristics_tissue_region",
                    "characteristics_location",
                    "compartment",
                ],
            )
        )
        row["donor_guess"] = normalize_label(
            first_existing(
                row,
                [
                    "characteristics_patient",
                    "characteristics_patient_id",
                    "characteristics_donor",
                    "characteristics_subject",
                    "characteristics_case",
                    "source_name_ch1",
                ],
            )
        )

    summary = {
        "dataset": config["project"]["primary_dataset"],
        "thresholds": {
            "min_aligned_reads": args.min_aligned_reads,
            "min_code_counts": args.min_code_counts,
            "min_total_code_counts": args.min_total_code_counts,
            "min_trimmed_fraction": args.min_trimmed_fraction,
            "min_stitched_fraction": args.min_stitched_fraction,
            "min_aligned_fraction_stitched": args.min_aligned_fraction_stitched,
            "min_umi_q30": args.min_umi_q30,
            "min_rts_q30": args.min_rts_q30,
        },
        "n_rois": len(rows),
        "n_include_qc": sum(row["include_qc"] == "True" for row in rows),
        "n_exclude_qc": sum(row["include_qc"] != "True" for row in rows),
        "metadata_columns": sorted({key for row in roi_rows for key in row}),
        "qc_columns": sorted({key for row in qc_rows for key in row}),
        "exclusion_reasons": dict(
            Counter(
                reason
                for row in rows
                for reason in row["exclusion_reason"].split(";")
                if reason
            )
        ),
        "roi_counts_by_diagnosis_guess": summarize_counts(rows, ["diagnosis_guess"]),
        "roi_counts_by_diagnosis_group": summarize_counts(rows, ["diagnosis_group"]),
        "roi_counts_by_smoking_status_guess": summarize_counts(
            rows, ["smoking_status_guess"]
        ),
        "roi_counts_by_compartment_guess": summarize_counts(rows, ["compartment_guess"]),
        "qc_by_diagnosis_compartment_guess": summarize_qc_by_group(
            rows, ["diagnosis_guess", "compartment_guess"]
        ),
        "qc_by_diagnosis_group_compartment_guess": summarize_qc_by_group(
            rows, ["diagnosis_group", "compartment_guess"]
        ),
    }

    write_csv(
        table_dir / "gse292993_roi_qc_flags.csv",
        rows,
        preferred=[
            "geo_accession",
            "dcc_filename",
            "dcc_id",
            "diagnosis_guess",
            "diagnosis_group",
            "smoking_status_guess",
            "compartment_guess",
            "donor_guess",
            "include_qc",
            "exclusion_reason",
            "aligned_reads",
            "trimmed_fraction",
            "stitched_fraction",
            "aligned_fraction_stitched",
            "n_code_counts",
            "total_code_counts",
            "primary_gene_counts",
            "negative_probe_mean_counts",
            "negative_probe_max_counts",
        ],
    )
    (meta_dir / "gse292993_roi_qc_flag_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {table_dir / 'gse292993_roi_qc_flags.csv'}")
    print(f"Wrote {meta_dir / 'gse292993_roi_qc_flag_summary.json'}")

    failures = []
    if len(rows) != 794:
        failures.append(f"Expected 794 ROIs, found {len(rows)}")
    if not rows:
        failures.append("No merged ROI/QC rows")
    if args.strict and failures:
        print("Strict ROI QC flagging failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
