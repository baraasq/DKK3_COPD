#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import random
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


METRICS = [
    "median_log1p_dkk3_cpm",
    "median_dkk3_cpm",
    "median_dkk3_count",
    "fraction_dkk3_above_geometric_loq",
    "fraction_dkk3_above_arithmetic_loq",
]
INFERENCE_COMPARTMENTS = {"airway", "parenchyma", "vessel"}


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
        "compartment_guess",
        "metric",
        "n_copd_donors",
        "n_control_donors",
        "copd_mean",
        "control_mean",
        "mean_difference_copd_minus_control",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "permutation_p_two_sided",
        "copd_median",
        "control_median",
        "median_difference_copd_minus_control",
    ]
    ordered = [field for field in preferred if field in fieldnames]
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


def mean(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def mean_difference(copd: list[float], control: list[float]) -> float | None:
    if not copd or not control:
        return None
    return statistics.mean(copd) - statistics.mean(control)


def median_difference(copd: list[float], control: list[float]) -> float | None:
    if not copd or not control:
        return None
    return statistics.median(copd) - statistics.median(control)


def percentile(values: list[float], fraction: float) -> float:
    sorted_values = sorted(values)
    index = int((len(sorted_values) - 1) * fraction)
    return sorted_values[index]


def bootstrap_ci(
    copd: list[float],
    control: list[float],
    *,
    iterations: int,
    rng: random.Random,
) -> tuple[float | None, float | None]:
    if not copd or not control or iterations <= 0:
        return None, None
    estimates = []
    for _ in range(iterations):
        sample_copd = [rng.choice(copd) for _ in copd]
        sample_control = [rng.choice(control) for _ in control]
        estimate = mean_difference(sample_copd, sample_control)
        if estimate is not None:
            estimates.append(estimate)
    if not estimates:
        return None, None
    return percentile(estimates, 0.025), percentile(estimates, 0.975)


def permutation_p_value(
    copd: list[float],
    control: list[float],
    *,
    iterations: int,
    rng: random.Random,
) -> float | None:
    observed = mean_difference(copd, control)
    if observed is None or iterations <= 0:
        return None
    pooled = copd + control
    n_copd = len(copd)
    exceed = 0
    for _ in range(iterations):
        shuffled = list(pooled)
        rng.shuffle(shuffled)
        permuted = mean_difference(shuffled[:n_copd], shuffled[n_copd:])
        if permuted is not None and abs(permuted) >= abs(observed):
            exceed += 1
    return (exceed + 1) / (iterations + 1)


def values_by_diagnosis(rows: list[dict], metric: str) -> tuple[list[float], list[float]]:
    copd = []
    control = []
    for row in rows:
        value = as_float(row, metric)
        if value is None:
            continue
        if row.get("diagnosis_group") == "COPD":
            copd.append(value)
        elif row.get("diagnosis_group") == "Control":
            control.append(value)
    return copd, control


def test_compartment_metric(
    rows: list[dict],
    compartment: str,
    metric: str,
    *,
    permutations: int,
    bootstraps: int,
    rng: random.Random,
) -> dict:
    compartment_rows = [row for row in rows if row.get("compartment_guess") == compartment]
    copd, control = values_by_diagnosis(compartment_rows, metric)
    ci_low, ci_high = bootstrap_ci(copd, control, iterations=bootstraps, rng=rng)
    return {
        "compartment_guess": compartment,
        "metric": metric,
        "n_copd_donors": len(copd),
        "n_control_donors": len(control),
        "copd_mean": mean(copd),
        "control_mean": mean(control),
        "mean_difference_copd_minus_control": mean_difference(copd, control),
        "bootstrap_ci_low": ci_low,
        "bootstrap_ci_high": ci_high,
        "permutation_p_two_sided": permutation_p_value(
            copd, control, iterations=permutations, rng=rng
        ),
        "copd_median": median(copd),
        "control_median": median(control),
        "median_difference_copd_minus_control": median_difference(copd, control),
    }


def donor_counts(rows: list[dict]) -> list[dict]:
    groups: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        diagnosis = row.get("diagnosis_group")
        donor = row.get("donor_guess")
        if diagnosis and donor:
            groups[diagnosis].add(donor)
    return [
        {"diagnosis_group": diagnosis, "n_donors": len(donors)}
        for diagnosis, donors in sorted(groups.items())
    ]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Donor-aware COPD vs control tests for GSE292993 DKK3 summaries."
    )
    parser.add_argument("--permutations", type=int, default=10000)
    parser.add_argument("--bootstraps", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    rng = random.Random(args.seed)
    rows = read_csv(table_dir / "gse292993_dkk3_donor_compartment_summary.csv")
    compartments = sorted({row["compartment_guess"] for row in rows})

    results = [
        test_compartment_metric(
            rows,
            compartment,
            metric,
            permutations=args.permutations,
            bootstraps=args.bootstraps,
            rng=rng,
        )
        for compartment in compartments
        for metric in METRICS
    ]
    primary_results = [
        row for row in results if row["metric"] == "median_log1p_dkk3_cpm"
    ]
    summary = {
        "dataset": config["project"]["primary_dataset"],
        "primary_gene": config["project"]["gene"],
        "unit_of_inference": "donor_compartment",
        "permutations": args.permutations,
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "n_donor_compartment_rows": len(rows),
        "donor_counts_by_diagnosis": donor_counts(rows),
        "primary_metric": "median_log1p_dkk3_cpm",
        "primary_metric_results": primary_results,
    }

    write_csv(table_dir / "gse292993_dkk3_donor_effect_tests.csv", results)
    (meta_dir / "gse292993_dkk3_donor_effect_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {table_dir / 'gse292993_dkk3_donor_effect_tests.csv'}")
    print(f"Wrote {meta_dir / 'gse292993_dkk3_donor_effect_summary.json'}")

    failures = []
    for result in primary_results:
        if result["compartment_guess"] not in INFERENCE_COMPARTMENTS:
            continue
        if result["n_copd_donors"] < 2 or result["n_control_donors"] < 2:
            failures.append(
                f"Too few donors for {result['compartment_guess']} "
                "primary comparison"
            )
    if args.strict and failures:
        print("Strict donor DKK3 tests failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
