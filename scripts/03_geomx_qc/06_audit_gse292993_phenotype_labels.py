#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_INPUT = "results/tables/gse292993_roi_qc_flags.csv"
PATTERN_GROUPS = {
    "emphysema_or_laa": ["emphy", "laa", "attenuation"],
    "gold": ["gold"],
    "spirometry": ["fev", "fvc"],
    "smoking": ["smok", "pack"],
    "condition": ["condition", "diagnosis", "disease", "group"],
    "donor": ["patient", "donor", "subject", "case", "source_name"],
    "compartment": ["compartment", "location", "region", "tissue"],
}
DEFAULT_GROUP_COLUMNS = [
    "diagnosis_guess",
    "diagnosis_group",
    "smoking_status_guess",
    "characteristics_condition",
    "characteristics_laa950",
    "characteristics_gold",
    "characteristics_fev1_pred",
    "compartment_guess",
]


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


def normalize(value: object) -> str:
    return str(value).strip() if value not in (None, "") else ""


def as_float(value: object) -> float | None:
    text = normalize(value)
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def matching_columns(columns: list[str]) -> dict[str, list[str]]:
    output: dict[str, list[str]] = {}
    for group, patterns in PATTERN_GROUPS.items():
        output[group] = [
            column
            for column in columns
            if any(pattern in column.casefold() for pattern in patterns)
        ]
    return output


def column_summary(rows: list[dict], column: str) -> dict:
    values = [normalize(row.get(column)) for row in rows]
    nonempty = [value for value in values if value]
    unique_values = sorted(set(nonempty))
    numeric_values = [as_float(value) for value in nonempty]
    numeric_values = [value for value in numeric_values if value is not None]
    output = {
        "column": column,
        "n_rows": len(rows),
        "n_nonempty": len(nonempty),
        "n_unique": len(unique_values),
        "unique_preview": ";".join(unique_values[:30]),
        "is_all_nonempty_numeric": bool(nonempty) and len(numeric_values) == len(nonempty),
    }
    if numeric_values:
        output.update(
            {
                "numeric_min": min(numeric_values),
                "numeric_median": statistics.median(numeric_values),
                "numeric_max": max(numeric_values),
            }
        )
    return output


def count_rows(rows: list[dict], columns: list[str], donor_column: str | None = None) -> list[dict]:
    counter: Counter[tuple[str, ...]] = Counter()
    donors: dict[tuple[str, ...], set[str]] = defaultdict(set)
    for row in rows:
        key = tuple(normalize(row.get(column)) or "unknown" for column in columns)
        counter[key] += 1
        if donor_column:
            donor = normalize(row.get(donor_column))
            if donor:
                donors[key].add(donor)
    output = []
    for key, n_rows in sorted(counter.items()):
        item = dict(zip(columns, key))
        item["n_rois"] = n_rows
        if donor_column:
            item["n_donors"] = len(donors[key])
        output.append(item)
    return output


def threshold_counts(rows: list[dict], column: str, thresholds: list[float]) -> list[dict]:
    values_by_donor: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        donor = normalize(row.get("donor_guess"))
        value = as_float(row.get(column))
        if donor and value is not None:
            values_by_donor[donor].append(value)

    output = []
    for threshold in thresholds:
        roi_ge = sum(
            1
            for row in rows
            if (value := as_float(row.get(column))) is not None and value >= threshold
        )
        donor_ge = sum(
            1
            for values in values_by_donor.values()
            if values and max(values) >= threshold
        )
        output.append(
            {
                "column": column,
                "threshold_ge": threshold,
                "n_rois_ge_threshold": roi_ge,
                "n_rois_with_numeric_value": sum(
                    as_float(row.get(column)) is not None for row in rows
                ),
                "n_donors_ge_threshold": donor_ge,
                "n_donors_with_numeric_value": len(values_by_donor),
            }
        )
    return output


def classify_label_support(matches: dict[str, list[str]]) -> dict:
    emphysema_columns = [
        column
        for column in matches.get("emphysema_or_laa", [])
        if "emphy" in column.casefold()
    ]
    laa_columns = [
        column
        for column in matches.get("emphysema_or_laa", [])
        if "laa" in column.casefold() or "attenuation" in column.casefold()
    ]
    if emphysema_columns:
        status = "direct_emphysema_label_candidate_present"
    elif laa_columns:
        status = "numeric_laa_proxy_candidate_present"
    else:
        status = "no_emphysema_or_laa_label_detected"
    return {
        "status": status,
        "direct_emphysema_label_candidate_columns": emphysema_columns,
        "laa_proxy_candidate_columns": laa_columns,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Audit GSE292993 GeoMx ROI metadata for emphysema, LAA950, GOLD, "
            "smoking, and condition labels."
        )
    )
    parser.add_argument("--input", default=DEFAULT_INPUT)
    parser.add_argument(
        "--include-qc-only",
        action="store_true",
        help="Restrict summaries to ROIs passing include_qc == True.",
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results_dirs = ensure_results_dirs(load_config())
    input_path = project_path(args.input)

    if not input_path.exists():
        summary = {
            "input": str(input_path),
            "input_exists": False,
            "failure": "Input table missing; run scripts/03_geomx_qc/03_flag_gse292993_roi_qc.py first.",
        }
        print(json.dumps(summary, indent=2))
        return 2 if args.strict else 1

    all_rows = read_csv(input_path)
    rows = [
        row for row in all_rows if row.get("include_qc") == "True"
    ] if args.include_qc_only else all_rows
    columns = sorted({key for row in all_rows for key in row})
    matches = matching_columns(columns)
    label_support = classify_label_support(matches)

    candidate_columns = sorted(
        {
            column
            for group_columns in matches.values()
            for column in group_columns
        }
    )
    candidate_summary = [column_summary(rows, column) for column in candidate_columns]

    group_columns = [column for column in DEFAULT_GROUP_COLUMNS if column in columns]
    group_summaries: list[dict] = []
    for column in group_columns:
        group_summaries.extend(
            {
                "grouping_column": column,
                **row,
            }
            for row in count_rows(rows, [column], donor_column="donor_guess" if "donor_guess" in columns else None)
        )

    threshold_rows: list[dict] = []
    for column in label_support["laa_proxy_candidate_columns"]:
        threshold_rows.extend(threshold_counts(rows, column, [1.0, 5.0, 10.0, 15.0]))

    column_path = results_dirs["tables"] / "gse292993_phenotype_label_column_summary.csv"
    group_path = results_dirs["tables"] / "gse292993_phenotype_label_group_counts.csv"
    threshold_path = results_dirs["tables"] / "gse292993_laa_threshold_sensitivity_counts.csv"
    summary_path = results_dirs["meta"] / "gse292993_phenotype_label_audit_summary.json"

    write_csv(
        column_path,
        candidate_summary,
        preferred=[
            "column",
            "n_rows",
            "n_nonempty",
            "n_unique",
            "unique_preview",
            "is_all_nonempty_numeric",
            "numeric_min",
            "numeric_median",
            "numeric_max",
        ],
    )
    write_csv(group_path, group_summaries)
    write_csv(
        threshold_path,
        threshold_rows,
        preferred=[
            "column",
            "threshold_ge",
            "n_rois_ge_threshold",
            "n_rois_with_numeric_value",
            "n_donors_ge_threshold",
            "n_donors_with_numeric_value",
        ],
    )

    summary = {
        "input": str(input_path),
        "input_exists": True,
        "include_qc_only": args.include_qc_only,
        "n_rows_total": len(all_rows),
        "n_rows_summarized": len(rows),
        "n_donors_summarized": len(
            {normalize(row.get("donor_guess")) for row in rows if normalize(row.get("donor_guess"))}
        ),
        "matching_columns": matches,
        "label_support": label_support,
        "interpretation": (
            "Direct emphysema columns can support Control/NE-COPD/E-COPD grouping. "
            "LAA/GOLD/FEV fields are proxies unless the paper's exact threshold or "
            "author-provided label is confirmed."
        ),
        "outputs": {
            "column_summary": str(column_path),
            "group_counts": str(group_path),
            "laa_threshold_sensitivity_counts": str(threshold_path),
            "summary": str(summary_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(column_path)
    print(group_path)
    print(threshold_path)

    if args.strict and not rows:
        raise SystemExit("Strict phenotype label audit failed: no rows to summarize")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
