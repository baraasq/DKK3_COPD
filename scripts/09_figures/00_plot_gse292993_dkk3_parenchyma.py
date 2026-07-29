#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import math
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


LABEL_ORDER = ["Non Smoker", "Smoker", "COPD"]
METRICS = [
    {
        "column": "median_log1p_dkk3_cpm",
        "title": "Abundance",
        "ylabel": "median log1p(DKK3 CPM)",
    },
    {
        "column": "fraction_dkk3_above_geometric_loq",
        "title": "Above-background frequency",
        "ylabel": "fraction above geometric LOQ",
    },
    {
        "column": "median_dkk3_count",
        "title": "Raw signal",
        "ylabel": "median DKK3 count",
    },
]
PAIR_ORDER = [
    ("COPD", "Non Smoker"),
    ("COPD", "Smoker"),
    ("Smoker", "Non Smoker"),
]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict, column: str) -> float | None:
    value = row.get(column)
    if value in (None, ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def stable_jitter(label: str, *, width: float = 0.18) -> float:
    digest = hashlib.sha1(label.encode("utf-8")).hexdigest()
    value = int(digest[:8], 16) / 0xFFFFFFFF
    return (value - 0.5) * 2 * width


def parenchyma_values(rows: list[dict], metric: str, compartment: str) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        if row.get("compartment_guess") != compartment:
            continue
        label = row.get("diagnosis_guess")
        value = as_float(row, metric)
        if label in LABEL_ORDER and value is not None:
            item = dict(row)
            item["plot_value"] = value
            grouped[label].append(item)
    return {label: grouped.get(label, []) for label in LABEL_ORDER}


def effect_lookup(rows: list[dict], compartment: str, metric: str) -> dict[tuple[str, str], dict]:
    output = {}
    for row in rows:
        if row.get("compartment_guess") != compartment or row.get("metric") != metric:
            continue
        output[(row.get("label_a"), row.get("label_b"))] = row
    return output


def format_p(value: float | None) -> str:
    if value is None:
        return "p=NA"
    if value < 0.001:
        return "p<0.001"
    return f"p={value:.3f}"


def metric_effect_text(effects: dict[tuple[str, str], dict]) -> str:
    parts = []
    for label_a, label_b in PAIR_ORDER:
        row = effects.get((label_a, label_b))
        if not row:
            continue
        difference = as_float(row, "mean_difference_label_a_minus_label_b")
        p_value = as_float(row, "permutation_p_two_sided")
        if difference is None:
            continue
        parts.append(f"{label_a}-{label_b}: {difference:+.2f}, {format_p(p_value)}")
    return "\n".join(parts)


def finite_values(grouped: dict[str, list[dict]]) -> list[float]:
    return [
        float(row["plot_value"])
        for label in LABEL_ORDER
        for row in grouped.get(label, [])
        if math.isfinite(float(row["plot_value"]))
    ]


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot donor-level parenchymal DKK3 signal in GSE292993."
    )
    parser.add_argument("--compartment", default="parenchyma")
    parser.add_argument("--formats", default="png,svg,pdf")
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    figure_dir = output["figures"] / "gse292993_dkk3"
    figure_dir.mkdir(parents=True, exist_ok=True)

    donor_rows = read_csv(
        table_dir / "gse292993_dkk3_donor_diagnosis_compartment_summary.csv"
    )
    effect_rows = read_csv(table_dir / "gse292993_dkk3_smoking_strata_effect_tests.csv")

    fig, axes = plt.subplots(1, len(METRICS), figsize=(12.8, 4.4), constrained_layout=True)
    colors = {
        "Non Smoker": "#5B8FF9",
        "Smoker": "#F6BD16",
        "COPD": "#E8684A",
    }
    for axis, metric in zip(axes, METRICS):
        grouped = parenchyma_values(donor_rows, metric["column"], args.compartment)
        data = [
            [row["plot_value"] for row in grouped[label]]
            for label in LABEL_ORDER
        ]
        box = axis.boxplot(
            data,
            positions=range(1, len(LABEL_ORDER) + 1),
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.4},
        )
        for patch, label in zip(box["boxes"], LABEL_ORDER):
            patch.set_facecolor(colors[label])
            patch.set_alpha(0.28)
            patch.set_edgecolor(colors[label])
        for whisker in box["whiskers"]:
            whisker.set_color("#666666")
        for cap in box["caps"]:
            cap.set_color("#666666")

        for index, label in enumerate(LABEL_ORDER, start=1):
            for row in grouped[label]:
                donor = row.get("donor_guess", "")
                x_value = index + stable_jitter(f"{metric['column']}:{donor}")
                axis.scatter(
                    x_value,
                    row["plot_value"],
                    s=34,
                    color=colors[label],
                    edgecolor="black",
                    linewidth=0.35,
                    alpha=0.82,
                    zorder=3,
                )
            axis.text(
                index,
                0.02,
                f"n={len(grouped[label])}",
                transform=axis.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=8,
                color="#555555",
            )

        values = finite_values(grouped)
        if metric["column"].startswith("fraction"):
            axis.set_ylim(-0.04, 1.04)
        elif values:
            span = max(values) - min(values)
            padding = span * 0.12 if span else 0.5
            axis.set_ylim(min(values) - padding, max(values) + padding)
        effects = effect_lookup(effect_rows, args.compartment, metric["column"])
        text = metric_effect_text(effects)
        if text:
            axis.text(
                0.02,
                0.98,
                text,
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "alpha": 0.75, "edgecolor": "none"},
            )
        axis.set_title(metric["title"], fontsize=11, fontweight="bold")
        axis.set_ylabel(metric["ylabel"])
        axis.set_xticks(range(1, len(LABEL_ORDER) + 1), LABEL_ORDER, rotation=18)
        axis.grid(axis="y", color="#dddddd", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    fig.suptitle("GSE292993 donor-level parenchymal DKK3 signal", fontsize=13, fontweight="bold")
    output_paths = []
    for fmt in [item.strip().lower() for item in args.formats.split(",") if item.strip()]:
        path = figure_dir / f"gse292993_dkk3_parenchyma_donor_signal.{fmt}"
        fig.savefig(path, dpi=300 if fmt == "png" else None)
        output_paths.append(str(path))
    plt.close(fig)
    print("\n".join(output_paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
