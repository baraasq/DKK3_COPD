#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


DIAGNOSIS_ORDER = ["Non Smoker", "Smoker", "COPD", "unknown"]
COMPARTMENT_ORDER = ["airway", "parenchyma", "vessel", "unknown"]
QC_METRICS = [
    {
        "column": "total_code_counts",
        "title": "Total code counts",
        "ylabel": "log10(total code counts + 1)",
        "transform": "log10p1",
    },
    {
        "column": "n_code_counts",
        "title": "Detected codes",
        "ylabel": "detected nonzero codes",
        "transform": "identity",
    },
    {
        "column": "aligned_reads",
        "title": "Aligned reads",
        "ylabel": "log10(aligned reads + 1)",
        "transform": "log10p1",
    },
    {
        "column": "trimmed_fraction",
        "title": "Trimmed fraction",
        "ylabel": "trimmed / raw reads",
        "transform": "identity",
    },
    {
        "column": "stitched_fraction",
        "title": "Stitched fraction",
        "ylabel": "stitched / trimmed reads",
        "transform": "identity",
    },
    {
        "column": "aligned_fraction_stitched",
        "title": "Aligned fraction",
        "ylabel": "aligned / stitched reads",
        "transform": "identity",
    },
    {
        "column": "umi_q30",
        "title": "UMI Q30",
        "ylabel": "fraction UMI bases Q30",
        "transform": "identity",
    },
    {
        "column": "rts_q30",
        "title": "RTS Q30",
        "ylabel": "fraction RTS bases Q30",
        "transform": "identity",
    },
    {
        "column": "negative_probe_mean_counts",
        "title": "Negative-probe background",
        "ylabel": "mean negative-probe counts",
        "transform": "identity",
    },
]


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


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


def ordered_labels(rows: list[dict], column: str, preferred: list[str]) -> list[str]:
    labels = {str(row.get(column) or "unknown") for row in rows}
    ordered = [label for label in preferred if label in labels]
    ordered.extend(sorted(label for label in labels if label not in ordered))
    return ordered


def stable_jitter(label: str, *, width: float = 0.18) -> float:
    digest = hashlib.sha1(label.encode("utf-8")).hexdigest()
    value = int(digest[:8], 16) / 0xFFFFFFFF
    return (value - 0.5) * 2 * width


def transform_value(value: float, transform: str) -> float:
    if transform == "log10p1":
        return math.log10(value + 1)
    return value


def metric_values(rows: list[dict], metric: dict, group_column: str) -> dict[str, list[dict]]:
    grouped: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        value = as_float(row, metric["column"])
        if value is None:
            continue
        label = str(row.get(group_column) or "unknown")
        item = dict(row)
        item["plot_value"] = transform_value(value, metric.get("transform", "identity"))
        grouped[label].append(item)
    return grouped


def count_table(rows: list[dict], row_column: str, stack_column: str) -> dict[str, Counter]:
    table: dict[str, Counter] = defaultdict(Counter)
    for row in rows:
        row_label = str(row.get(row_column) or "unknown")
        stack_label = str(row.get(stack_column) or "unknown")
        table[row_label][stack_label] += 1
    return table


def donor_count_table(rows: list[dict], row_column: str, stack_column: str) -> dict[str, Counter]:
    seen = set()
    table: dict[str, Counter] = defaultdict(Counter)
    for row in rows:
        donor = row.get("donor_guess")
        if donor in (None, "", "unknown"):
            continue
        row_label = str(row.get(row_column) or "unknown")
        stack_label = str(row.get(stack_column) or "unknown")
        key = (donor, row_label, stack_label)
        if key in seen:
            continue
        seen.add(key)
        table[row_label][stack_label] += 1
    return table


def plot_metric_distributions(rows: list[dict], path_stem: Path, formats: list[str]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = ordered_labels(rows, "diagnosis_guess", DIAGNOSIS_ORDER)
    colors = {
        True: "#3BA272",
        False: "#D14A61",
    }
    fig, axes = plt.subplots(3, 3, figsize=(14.2, 11.2), constrained_layout=True)
    for axis, metric in zip(axes.reshape(-1), QC_METRICS):
        grouped = metric_values(rows, metric, "diagnosis_guess")
        data = [
            [row["plot_value"] for row in grouped.get(label, [])]
            for label in labels
        ]
        box = axis.boxplot(
            data,
            positions=range(1, len(labels) + 1),
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "black", "linewidth": 1.25},
        )
        for patch in box["boxes"]:
            patch.set_facecolor("#B8D7F0")
            patch.set_alpha(0.35)
            patch.set_edgecolor("#5B8FF9")
        for index, label in enumerate(labels, start=1):
            for row in grouped.get(label, []):
                roi_id = row.get("geo_accession") or row.get("dcc_id") or ""
                x_value = index + stable_jitter(f"{metric['column']}:{roi_id}")
                passed = as_bool(row.get("include_qc"))
                axis.scatter(
                    x_value,
                    row["plot_value"],
                    s=12 if passed else 28,
                    marker="o" if passed else "x",
                    color=colors[passed],
                    alpha=0.55 if passed else 0.85,
                    linewidth=0.7,
                    zorder=3,
                )
            axis.text(
                index,
                0.02,
                f"n={len(grouped.get(label, []))}",
                transform=axis.get_xaxis_transform(),
                ha="center",
                va="bottom",
                fontsize=8,
                color="#555555",
            )
        axis.set_title(metric["title"], fontsize=10.5, fontweight="bold")
        axis.set_ylabel(metric["ylabel"])
        axis.set_xticks(range(1, len(labels) + 1), labels, rotation=20)
        axis.grid(axis="y", color="#dddddd", linewidth=0.8)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    fig.suptitle("GSE292993 GeoMx ROI QC metrics", fontsize=13, fontweight="bold")
    outputs = []
    for fmt in formats:
        path = path_stem.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=300 if fmt == "png" else None)
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_balance(rows: list[dict], path_stem: Path, formats: list[str]) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    diagnosis_labels = ordered_labels(rows, "diagnosis_guess", DIAGNOSIS_ORDER)
    compartment_labels = ordered_labels(rows, "compartment_guess", COMPARTMENT_ORDER)
    colors = {
        "airway": "#66C2A5",
        "parenchyma": "#FC8D62",
        "vessel": "#8DA0CB",
        "unknown": "#B3B3B3",
    }
    roi_counts = count_table(rows, "diagnosis_guess", "compartment_guess")
    donor_counts = donor_count_table(rows, "diagnosis_guess", "compartment_guess")
    qc_counts = count_table(rows, "diagnosis_guess", "include_qc")

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), constrained_layout=True)
    for axis, table, title, ylabel in [
        (axes[0], roi_counts, "ROI counts by compartment", "n ROIs"),
        (axes[1], donor_counts, "Donor-compartment counts", "n donors with compartment"),
    ]:
        bottoms = [0] * len(diagnosis_labels)
        for compartment in compartment_labels:
            values = [table[label][compartment] for label in diagnosis_labels]
            axis.bar(
                diagnosis_labels,
                values,
                bottom=bottoms,
                label=compartment,
                color=colors.get(compartment, "#999999"),
                edgecolor="white",
            )
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        axis.set_title(title, fontsize=10.5, fontweight="bold")
        axis.set_ylabel(ylabel)
        axis.tick_params(axis="x", rotation=20)
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)

    pass_labels = ["True", "False"]
    pass_colors = {"True": "#3BA272", "False": "#D14A61"}
    bottoms = [0] * len(diagnosis_labels)
    for passed in pass_labels:
        values = [qc_counts[label][passed] for label in diagnosis_labels]
        axes[2].bar(
            diagnosis_labels,
            values,
            bottom=bottoms,
            label="pass QC" if passed == "True" else "fail QC",
            color=pass_colors[passed],
            edgecolor="white",
        )
        bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
    axes[2].set_title("ROI QC inclusion", fontsize=10.5, fontweight="bold")
    axes[2].set_ylabel("n ROIs")
    axes[2].tick_params(axis="x", rotation=20)
    axes[2].spines["top"].set_visible(False)
    axes[2].spines["right"].set_visible(False)

    axes[0].legend(title="Compartment", fontsize=8, title_fontsize=8)
    axes[2].legend(fontsize=8)
    fig.suptitle("GSE292993 GeoMx sample/QC balance", fontsize=13, fontweight="bold")
    outputs = []
    for fmt in formats:
        path = path_stem.with_suffix(f".{fmt}")
        fig.savefig(path, dpi=300 if fmt == "png" else None)
        outputs.append(path)
    plt.close(fig)
    return outputs


def qc_summary(rows: list[dict]) -> dict:
    return {
        "n_rois": len(rows),
        "n_pass_qc": sum(as_bool(row.get("include_qc")) for row in rows),
        "n_fail_qc": sum(not as_bool(row.get("include_qc")) for row in rows),
        "diagnosis_counts": dict(Counter(str(row.get("diagnosis_guess") or "unknown") for row in rows)),
        "compartment_counts": dict(Counter(str(row.get("compartment_guess") or "unknown") for row in rows)),
        "standard_ngs_like_metrics_available": {
            "total_counts": "total_code_counts",
            "detected_features": "n_code_counts",
            "umi_quality": "umi_q30",
            "mitochondrial_percent": None,
        },
        "geomx_specific_metrics_available": [
            "aligned_reads",
            "trimmed_fraction",
            "stitched_fraction",
            "aligned_fraction_stitched",
            "umi_q30",
            "rts_q30",
            "negative_probe_mean_counts",
            "negative_probe_max_counts",
        ],
        "notes": [
            "GeoMx DCC QC does not provide a direct percent-mitochondrial metric like scRNA-seq or whole-transcriptome spot matrices.",
            "Nuclei count, ROI area, counts per nucleus, and counts per area require image/ROI annotation metadata if available outside the DCC files.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot GSE292993 GeoMx ROI QC distributions and sample balance."
    )
    parser.add_argument("--formats", default="png,svg,pdf")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    figure_dir = output["figures"] / "gse292993_qc"
    meta_dir = output["meta"]
    figure_dir.mkdir(parents=True, exist_ok=True)
    formats = [item.strip().lower() for item in args.formats.split(",") if item.strip()]
    rows = read_csv(table_dir / "gse292993_roi_qc_flags.csv")

    outputs = []
    outputs.extend(
        plot_metric_distributions(
            rows,
            figure_dir / "gse292993_roi_qc_metric_distributions",
            formats,
        )
    )
    outputs.extend(
        plot_balance(
            rows,
            figure_dir / "gse292993_roi_qc_sample_balance",
            formats,
        )
    )
    summary = qc_summary(rows)
    (meta_dir / "gse292993_roi_qc_plot_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print()
    for path in outputs:
        print(path)
    print(meta_dir / "gse292993_roi_qc_plot_summary.json")

    if args.strict and not rows:
        print("Strict QC plotting failed: no ROI QC rows found", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
