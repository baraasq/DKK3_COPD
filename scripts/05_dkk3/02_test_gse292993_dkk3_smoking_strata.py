#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import random
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


PRIMARY_METRIC = "median_log1p_dkk3_cpm"
DEFAULT_PAIRS = [
    ("COPD", "Non Smoker"),
    ("COPD", "Smoker"),
    ("Smoker", "Non Smoker"),
]


def load_effect_module():
    path = Path(__file__).resolve().with_name("01_test_gse292993_dkk3_donor_effects.py")
    spec = importlib.util.spec_from_file_location("gse292993_dkk3_effects", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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
        "label_a",
        "label_b",
        "n_label_a_donors",
        "n_label_b_donors",
        "label_a_mean",
        "label_b_mean",
        "mean_difference_label_a_minus_label_b",
        "bootstrap_ci_low",
        "bootstrap_ci_high",
        "permutation_p_two_sided",
        "label_a_median",
        "label_b_median",
        "median_difference_label_a_minus_label_b",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_pairs(text: str | None, labels: list[str]) -> list[tuple[str, str]]:
    if text:
        pairs = []
        for item in text.split(","):
            left, separator, right = item.partition(":")
            if not separator:
                raise ValueError(f"Pair must be formatted as A:B, got {item}")
            pairs.append((left.strip(), right.strip()))
        return pairs
    default = [pair for pair in DEFAULT_PAIRS if pair[0] in labels and pair[1] in labels]
    return default or list(combinations(labels, 2))


def values_by_label(
    rows: list[dict],
    *,
    label_column: str,
    label: str,
    metric: str,
    effects,
) -> list[float]:
    values = []
    for row in rows:
        if row.get(label_column) != label:
            continue
        value = effects.as_float(row, metric)
        if value is not None:
            values.append(value)
    return values


def test_pair_metric(
    rows: list[dict],
    *,
    label_column: str,
    label_a: str,
    label_b: str,
    compartment: str,
    metric: str,
    permutations: int,
    bootstraps: int,
    rng: random.Random,
    effects,
) -> dict:
    compartment_rows = [row for row in rows if row.get("compartment_guess") == compartment]
    values_a = values_by_label(
        compartment_rows,
        label_column=label_column,
        label=label_a,
        metric=metric,
        effects=effects,
    )
    values_b = values_by_label(
        compartment_rows,
        label_column=label_column,
        label=label_b,
        metric=metric,
        effects=effects,
    )
    ci_low, ci_high = effects.bootstrap_ci(
        values_a, values_b, iterations=bootstraps, rng=rng
    )
    return {
        "compartment_guess": compartment,
        "metric": metric,
        "label_column": label_column,
        "label_a": label_a,
        "label_b": label_b,
        "n_label_a_donors": len(values_a),
        "n_label_b_donors": len(values_b),
        "label_a_mean": effects.mean(values_a),
        "label_b_mean": effects.mean(values_b),
        "mean_difference_label_a_minus_label_b": effects.mean_difference(
            values_a, values_b
        ),
        "bootstrap_ci_low": ci_low,
        "bootstrap_ci_high": ci_high,
        "permutation_p_two_sided": effects.permutation_p_value(
            values_a, values_b, iterations=permutations, rng=rng
        ),
        "label_a_median": effects.median(values_a),
        "label_b_median": effects.median(values_b),
        "median_difference_label_a_minus_label_b": effects.median_difference(
            values_a, values_b
        ),
    }


def donor_counts(rows: list[dict], label_column: str) -> list[dict]:
    groups: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        label = row.get(label_column)
        donor = row.get("donor_guess")
        if label and donor:
            groups[label].add(donor)
    return [
        {label_column: label, "n_donors": len(donors)}
        for label, donors in sorted(groups.items())
    ]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Donor-aware DKK3 tests stratified by COPD, smoker, and non-smoker labels."
    )
    parser.add_argument("--label-column", default="diagnosis_guess")
    parser.add_argument(
        "--pairs",
        help="Comma-separated label pairs formatted as 'A:B,C:D'. Defaults to COPD/non-smoker/smoker comparisons.",
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
    effects = load_effect_module()
    rng = random.Random(args.seed)
    rows = read_csv(table_dir / "gse292993_dkk3_donor_diagnosis_compartment_summary.csv")
    labels = sorted(
        {
            row[args.label_column]
            for row in rows
            if row.get(args.label_column) not in (None, "", "unknown")
        }
    )
    compartments = sorted({row["compartment_guess"] for row in rows})
    pairs = parse_pairs(args.pairs, labels)

    results = [
        test_pair_metric(
            rows,
            label_column=args.label_column,
            label_a=label_a,
            label_b=label_b,
            compartment=compartment,
            metric=metric,
            permutations=args.permutations,
            bootstraps=args.bootstraps,
            rng=rng,
            effects=effects,
        )
        for compartment in compartments
        for label_a, label_b in pairs
        for metric in effects.METRICS
    ]
    primary_results = [
        row for row in results if row["metric"] == PRIMARY_METRIC
    ]
    summary = {
        "dataset": config["project"]["primary_dataset"],
        "primary_gene": config["project"]["gene"],
        "unit_of_inference": "donor_diagnosis_compartment",
        "label_column": args.label_column,
        "labels": labels,
        "pairs": [{"label_a": left, "label_b": right} for left, right in pairs],
        "permutations": args.permutations,
        "bootstraps": args.bootstraps,
        "seed": args.seed,
        "donor_counts_by_label": donor_counts(rows, args.label_column),
        "primary_metric": PRIMARY_METRIC,
        "primary_metric_results": primary_results,
    }

    write_csv(table_dir / "gse292993_dkk3_smoking_strata_effect_tests.csv", results)
    (meta_dir / "gse292993_dkk3_smoking_strata_effect_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {table_dir / 'gse292993_dkk3_smoking_strata_effect_tests.csv'}")
    print(f"Wrote {meta_dir / 'gse292993_dkk3_smoking_strata_effect_summary.json'}")

    failures = []
    for result in primary_results:
        if result["n_label_a_donors"] < 2 or result["n_label_b_donors"] < 2:
            failures.append(
                f"Too few donors for {result['label_a']} vs {result['label_b']} "
                f"in {result['compartment_guess']}"
            )
    if args.strict and failures:
        print("Strict stratified DKK3 tests failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
