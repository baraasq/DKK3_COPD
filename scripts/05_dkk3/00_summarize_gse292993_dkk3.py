#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


COMPARTMENTS_FOR_DIAGNOSTIC_SUMMARY = {"airway", "parenchyma", "vessel", "unknown"}


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = sorted({key for row in rows for key in row})
    preferred = [
        "donor_guess",
        "diagnosis_group",
        "diagnosis_guess",
        "compartment_guess",
        "n_rois",
        "n_dkk3_gt0",
        "n_dkk3_above_geometric_loq",
        "fraction_dkk3_above_geometric_loq",
        "median_dkk3_count",
        "median_dkk3_cpm",
        "median_log1p_dkk3_cpm",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: object) -> bool:
    return value in (True, "True", "true", "1", 1)


def as_float(row: dict, column: str) -> float | None:
    value = row.get(column)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def median(values: list[float]) -> float | None:
    if not values:
        return None
    return statistics.median(values)


def quantile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    sorted_values = sorted(values)
    index = int((len(sorted_values) - 1) * fraction)
    return sorted_values[index]


def enrich_roi_rows(rows: list[dict]) -> list[dict]:
    enriched = []
    for row in rows:
        item = dict(row)
        dkk3_count = as_float(item, "dkk3_count")
        total_counts = as_float(item, "total_code_counts")
        if dkk3_count is not None and total_counts and total_counts > 0:
            dkk3_cpm = dkk3_count / total_counts * 1_000_000
            item["dkk3_cpm"] = dkk3_cpm
            item["log1p_dkk3_cpm"] = math.log1p(dkk3_cpm)
        else:
            item["dkk3_cpm"] = None
            item["log1p_dkk3_cpm"] = None
        item["dkk3_gt0"] = bool(dkk3_count is not None and dkk3_count > 0)
        enriched.append(item)
    return enriched


def summarize_group(rows: list[dict], columns: list[str]) -> list[dict]:
    groups: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    for row in rows:
        key = tuple(str(row.get(column) or "unknown") for column in columns)
        groups[key].append(row)

    output = []
    for key, group_rows in sorted(groups.items()):
        dkk3_counts = [
            value for row in group_rows if (value := as_float(row, "dkk3_count")) is not None
        ]
        dkk3_cpms = [
            value for row in group_rows if (value := as_float(row, "dkk3_cpm")) is not None
        ]
        log_cpms = [
            value
            for row in group_rows
            if (value := as_float(row, "log1p_dkk3_cpm")) is not None
        ]
        above_geometric = [
            row for row in group_rows if as_bool(row.get("dkk3_above_geometric_loq"))
        ]
        above_arithmetic = [
            row for row in group_rows if as_bool(row.get("dkk3_above_arithmetic_loq"))
        ]
        item = dict(zip(columns, key))
        item.update(
            n_rois=len(group_rows),
            n_dkk3_gt0=sum(as_bool(row.get("dkk3_gt0")) for row in group_rows),
            n_dkk3_above_geometric_loq=len(above_geometric),
            fraction_dkk3_above_geometric_loq=(
                len(above_geometric) / len(group_rows) if group_rows else None
            ),
            n_dkk3_above_arithmetic_loq=len(above_arithmetic),
            fraction_dkk3_above_arithmetic_loq=(
                len(above_arithmetic) / len(group_rows) if group_rows else None
            ),
            median_dkk3_count=median(dkk3_counts),
            q25_dkk3_count=quantile(dkk3_counts, 0.25),
            q75_dkk3_count=quantile(dkk3_counts, 0.75),
            median_dkk3_cpm=median(dkk3_cpms),
            median_log1p_dkk3_cpm=median(log_cpms),
        )
        output.append(item)
    return output


def primary_rows(rows: list[dict]) -> list[dict]:
    return [
        row
        for row in rows
        if as_bool(row.get("include_qc"))
        and row.get("diagnosis_group") in {"COPD", "Control"}
        and row.get("compartment_guess") in COMPARTMENTS_FOR_DIAGNOSTIC_SUMMARY
        and row.get("donor_guess") not in (None, "", "unknown")
    ]


def donor_counts_by_column(rows: list[dict], column: str) -> list[dict]:
    pairs = {
        (row["donor_guess"], row[column])
        for row in rows
        if row.get("donor_guess") and row.get(column)
    }
    counts = Counter(label for _, label in pairs)
    return [
        {column: label, "n_donors": count}
        for label, count in sorted(counts.items())
    ]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize GSE292993 DKK3 ROI and donor-compartment signal."
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    input_path = table_dir / "gse292993_dkk3_loq_flags.csv"

    roi_rows = enrich_roi_rows(read_csv(input_path))
    primary = primary_rows(roi_rows)
    roi_by_group = summarize_group(primary, ["diagnosis_group", "compartment_guess"])
    donor_compartment = summarize_group(
        primary, ["donor_guess", "diagnosis_group", "compartment_guess"]
    )
    donor_diagnosis_compartment = summarize_group(
        primary,
        ["donor_guess", "diagnosis_guess", "diagnosis_group", "compartment_guess"],
    )
    donor_overall = summarize_group(primary, ["donor_guess", "diagnosis_group"])
    donor_diagnosis_overall = summarize_group(
        primary, ["donor_guess", "diagnosis_guess", "diagnosis_group"]
    )

    summary = {
        "dataset": config["project"]["primary_dataset"],
        "primary_gene": config["project"]["gene"],
        "n_rois_total": len(roi_rows),
        "n_rois_primary": len(primary),
        "n_donors_primary": len({row["donor_guess"] for row in primary}),
        "primary_roi_counts_by_diagnosis_compartment": roi_by_group,
        "primary_donor_counts_by_diagnosis_group": donor_counts_by_column(
            primary, "diagnosis_group"
        ),
        "primary_donor_counts_by_diagnosis_guess": donor_counts_by_column(
            primary, "diagnosis_guess"
        ),
    }

    write_csv(table_dir / "gse292993_dkk3_roi_signal.csv", roi_rows)
    write_csv(table_dir / "gse292993_dkk3_primary_roi_by_group.csv", roi_by_group)
    write_csv(
        table_dir / "gse292993_dkk3_donor_compartment_summary.csv",
        donor_compartment,
    )
    write_csv(
        table_dir / "gse292993_dkk3_donor_diagnosis_compartment_summary.csv",
        donor_diagnosis_compartment,
    )
    write_csv(table_dir / "gse292993_dkk3_donor_overall_summary.csv", donor_overall)
    write_csv(
        table_dir / "gse292993_dkk3_donor_diagnosis_overall_summary.csv",
        donor_diagnosis_overall,
    )
    (meta_dir / "gse292993_dkk3_signal_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {table_dir / 'gse292993_dkk3_roi_signal.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dkk3_primary_roi_by_group.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dkk3_donor_compartment_summary.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dkk3_donor_diagnosis_compartment_summary.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dkk3_donor_overall_summary.csv'}")
    print(f"Wrote {table_dir / 'gse292993_dkk3_donor_diagnosis_overall_summary.csv'}")
    print(f"Wrote {meta_dir / 'gse292993_dkk3_signal_summary.json'}")

    failures = []
    if not primary:
        failures.append("No primary QC-passing DKK3 ROI rows")
    if len({row["donor_guess"] for row in primary}) < 2:
        failures.append("Fewer than two donors in primary DKK3 summary")
    if args.strict and failures:
        print("Strict DKK3 summary failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
