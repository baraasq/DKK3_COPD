#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


BIOLOGICAL_COMPARTMENTS = ["airway", "parenchyma", "vessel"]
LABEL_ORDER = ["Non Smoker", "Smoker", "COPD"]
PREFERRED_CELLTYPE_ORDER = [
    "Club cell",
    "AT1",
    "AT2",
    "Fibroblast",
    "Smooth muscle",
    "Vascular Endothelia",
    "Lymphatic Endothelia",
    "Mesothelia",
]


def read_csv(path: Path) -> list[dict]:
    if not path.exists():
        return []
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
        parsed = float(value)
    except ValueError:
        return None
    return parsed if math.isfinite(parsed) else None


def mean(values: list[float]) -> float | None:
    return statistics.mean(values) if values else None


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def parse_compartments(value: str) -> list[str]:
    normalized = value.strip().lower()
    if normalized == "all":
        return BIOLOGICAL_COMPARTMENTS
    requested = [item.strip().lower() for item in value.split(",") if item.strip()]
    unknown = [item for item in requested if item not in BIOLOGICAL_COMPARTMENTS]
    if unknown:
        raise ValueError(
            "Unknown compartment(s): "
            + ", ".join(unknown)
            + ". Expected all or one/more of: "
            + ", ".join(BIOLOGICAL_COMPARTMENTS)
        )
    return requested


def output_stem(compartments: list[str]) -> str:
    if compartments == BIOLOGICAL_COMPARTMENTS:
        return "gse292993_nnls_cell_composition_all_compartments"
    return "gse292993_nnls_cell_composition_" + "_".join(compartments)


def read_celltype_manifest(table_dir: Path, compartment: str) -> list[dict]:
    return read_csv(table_dir / f"gse292993_{compartment}_nnls_deconvolution_celltype_manifest.csv")


def read_roi_deconvolution(table_dir: Path, compartment: str) -> list[dict]:
    return read_csv(table_dir / f"gse292993_{compartment}_nnls_deconvolution_roi.csv")


def fraction_columns(manifest_rows: list[dict]) -> list[dict]:
    output = []
    for row in manifest_rows:
        cell_type = row.get("cell_type")
        column = row.get("fraction_column")
        if cell_type and column:
            output.append({"cell_type": cell_type, "fraction_column": column})
    return output


def celltype_order(manifest_by_compartment: dict[str, list[dict]]) -> list[str]:
    seen = []
    for preferred in PREFERRED_CELLTYPE_ORDER:
        if any(preferred == item.get("cell_type") for rows in manifest_by_compartment.values() for item in rows):
            seen.append(preferred)
    all_labels = sorted(
        {
            item.get("cell_type")
            for rows in manifest_by_compartment.values()
            for item in rows
            if item.get("cell_type")
        }
    )
    seen.extend(label for label in all_labels if label not in seen)
    return seen


def donor_fraction_rows(
    *,
    roi_rows: list[dict],
    manifest_rows: list[dict],
    compartment: str,
) -> list[dict]:
    manifest = fraction_columns(manifest_rows)
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in roi_rows:
        donor = row.get("donor_guess")
        diagnosis = row.get("diagnosis_guess")
        if not donor or donor == "unknown":
            continue
        if diagnosis not in LABEL_ORDER:
            continue
        groups[(donor, diagnosis)].append(row)

    output = []
    for (donor, diagnosis), rows in sorted(groups.items()):
        item = {
            "compartment": compartment,
            "donor_guess": donor,
            "diagnosis_guess": diagnosis,
            "n_rois": len(rows),
        }
        for entry in manifest:
            values = [
                value
                for row in rows
                if (value := as_float(row, entry["fraction_column"])) is not None
            ]
            item[f"mean_{entry['fraction_column']}"] = mean(values)
        output.append(item)
    return output


def composition_summary_rows(
    donor_rows: list[dict],
    *,
    manifest_by_compartment: dict[str, list[dict]],
    compartments: list[str],
    labels: list[str],
) -> list[dict]:
    rows = []
    grouped: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in donor_rows:
        grouped[(row["compartment"], row["diagnosis_guess"])].append(row)

    for compartment in compartments:
        manifest = fraction_columns(manifest_by_compartment.get(compartment, []))
        for diagnosis in labels:
            group_rows = grouped.get((compartment, diagnosis), [])
            n_donors = len({row["donor_guess"] for row in group_rows})
            n_rois = sum(int(row.get("n_rois") or 0) for row in group_rows)
            for entry in manifest:
                donor_values = [
                    value
                    for row in group_rows
                    if (value := as_float(row, f"mean_{entry['fraction_column']}")) is not None
                ]
                rows.append(
                    {
                        "compartment": compartment,
                        "diagnosis_guess": diagnosis,
                        "cell_type": entry["cell_type"],
                        "fraction_column": entry["fraction_column"],
                        "mean_donor_fraction": mean(donor_values),
                        "median_donor_fraction": median(donor_values),
                        "n_donors": n_donors,
                        "n_donors_with_fraction": len(donor_values),
                        "n_rois": n_rois,
                    }
                )
    return rows


def value_lookup(summary_rows: list[dict]) -> dict[tuple[str, str, str], float]:
    output = {}
    for row in summary_rows:
        value = row.get("mean_donor_fraction")
        if value in (None, ""):
            continue
        output[(row["compartment"], row["diagnosis_guess"], row["cell_type"])] = float(value)
    return output


def donor_count_lookup(summary_rows: list[dict]) -> dict[tuple[str, str], int]:
    output = {}
    for row in summary_rows:
        key = (row["compartment"], row["diagnosis_guess"])
        output[key] = max(output.get(key, 0), int(row.get("n_donors") or 0))
    return output


def stack_values(
    lookup: dict[tuple[str, str, str], float],
    *,
    compartment: str,
    diagnosis: str,
    cell_types: list[str],
    normalize: bool,
) -> list[float]:
    values = [lookup.get((compartment, diagnosis, cell_type), 0.0) for cell_type in cell_types]
    total = sum(values)
    if normalize and total > 0:
        return [value / total for value in values]
    return values


def plot_stacked_bars(
    *,
    axes,
    summary_rows: list[dict],
    compartments: list[str],
    labels: list[str],
    cell_types: list[str],
    colors: dict[str, object],
    normalize: bool,
) -> None:
    lookup = value_lookup(summary_rows)
    n_lookup = donor_count_lookup(summary_rows)
    for axis, compartment in zip(axes, compartments):
        bottoms = [0.0] * len(labels)
        for cell_type in cell_types:
            values = [
                stack_values(
                    lookup,
                    compartment=compartment,
                    diagnosis=label,
                    cell_types=cell_types,
                    normalize=normalize,
                )[cell_types.index(cell_type)]
                for label in labels
            ]
            axis.bar(
                range(len(labels)),
                values,
                bottom=bottoms,
                label=cell_type,
                color=colors[cell_type],
                edgecolor="white",
                linewidth=0.35,
            )
            bottoms = [bottom + value for bottom, value in zip(bottoms, values)]
        for index, label in enumerate(labels):
            axis.text(
                index,
                1.01 if normalize else bottoms[index] + 0.01,
                f"n={n_lookup.get((compartment, label), 0)}",
                ha="center",
                va="bottom",
                fontsize=8,
                color="#555555",
            )
        axis.set_title(compartment.title(), fontsize=11, fontweight="bold")
        axis.set_xticks(range(len(labels)), labels, rotation=20, ha="right")
        axis.set_ylim(0, 1.08 if normalize else max([1.0, *bottoms]) * 1.12)
        axis.set_ylabel("mean donor fraction" + (" (normalized stack)" if normalize else ""))
        axis.spines["top"].set_visible(False)
        axis.spines["right"].set_visible(False)
        axis.grid(axis="y", color="#dddddd", linewidth=0.7, alpha=0.7)


def plot_heatmap(
    *,
    axis,
    summary_rows: list[dict],
    compartments: list[str],
    labels: list[str],
    cell_types: list[str],
    normalize: bool,
    cmap: str,
) -> None:
    import numpy as np

    lookup = value_lookup(summary_rows)
    columns = [(compartment, label) for compartment in compartments for label in labels]
    matrix = []
    for cell_type in cell_types:
        row = []
        for compartment, label in columns:
            values = stack_values(
                lookup,
                compartment=compartment,
                diagnosis=label,
                cell_types=cell_types,
                normalize=normalize,
            )
            row.append(values[cell_types.index(cell_type)])
        matrix.append(row)
    data = np.asarray(matrix, dtype=float)
    image = axis.imshow(data, aspect="auto", cmap=cmap, vmin=0)
    axis.set_yticks(range(len(cell_types)), cell_types)
    axis.set_xticks(
        range(len(columns)),
        [f"{compartment}\n{label}" for compartment, label in columns],
        rotation=35,
        ha="right",
    )
    axis.set_title("Mean donor cell-fraction heatmap", fontsize=11, fontweight="bold")
    for row_index in range(data.shape[0]):
        for col_index in range(data.shape[1]):
            value = data[row_index, col_index]
            axis.text(
                col_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=6.5,
                color="white" if value > max(0.35, data.max() * 0.55) else "#222222",
            )
    return image


def collect_inputs(table_dir: Path, compartments: list[str]) -> tuple[list[dict], dict[str, list[dict]], dict]:
    donor_rows = []
    manifest_by_compartment = {}
    input_summary = {
        "compartments": compartments,
        "missing_roi_tables": [],
        "missing_manifest_tables": [],
        "n_roi_rows_by_compartment": {},
        "n_donor_rows_by_compartment": {},
        "n_cell_types_by_compartment": {},
    }
    for compartment in compartments:
        roi_path = table_dir / f"gse292993_{compartment}_nnls_deconvolution_roi.csv"
        manifest_path = table_dir / f"gse292993_{compartment}_nnls_deconvolution_celltype_manifest.csv"
        roi_rows = read_roi_deconvolution(table_dir, compartment)
        manifest_rows = read_celltype_manifest(table_dir, compartment)
        if not roi_path.exists():
            input_summary["missing_roi_tables"].append(str(roi_path))
        if not manifest_path.exists():
            input_summary["missing_manifest_tables"].append(str(manifest_path))
        manifest_by_compartment[compartment] = manifest_rows
        compartment_donor_rows = donor_fraction_rows(
            roi_rows=roi_rows,
            manifest_rows=manifest_rows,
            compartment=compartment,
        )
        donor_rows.extend(compartment_donor_rows)
        input_summary["n_roi_rows_by_compartment"][compartment] = len(roi_rows)
        input_summary["n_donor_rows_by_compartment"][compartment] = len(compartment_donor_rows)
        input_summary["n_cell_types_by_compartment"][compartment] = len(fraction_columns(manifest_rows))
    return donor_rows, manifest_by_compartment, input_summary


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Plot donor-balanced NNLS cell composition across GSE292993 GeoMx "
            "airway, parenchyma, and vessel compartments."
        )
    )
    parser.add_argument(
        "--compartment",
        default="all",
        help="all or comma-separated subset of airway, parenchyma, vessel.",
    )
    parser.add_argument("--formats", default="png,svg,pdf")
    parser.add_argument(
        "--no-normalize-stacks",
        action="store_true",
        help=(
            "Do not normalize each compartment/diagnosis stacked bar to 1. "
            "Default normalizes after donor-balanced averaging."
        ),
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    figure_dir = output["figures"] / "gse292993_cell_composition"
    figure_dir.mkdir(parents=True, exist_ok=True)
    compartments = parse_compartments(args.compartment)
    normalize = not args.no_normalize_stacks

    donor_rows, manifest_by_compartment, input_summary = collect_inputs(table_dir, compartments)
    cell_types = celltype_order(manifest_by_compartment)
    summary_rows = composition_summary_rows(
        donor_rows,
        manifest_by_compartment=manifest_by_compartment,
        compartments=compartments,
        labels=LABEL_ORDER,
    )

    composition_table = table_dir / f"{output_stem(compartments)}_summary.csv"
    donor_table = table_dir / f"{output_stem(compartments)}_donor_means.csv"
    write_csv(
        composition_table,
        summary_rows,
        preferred=[
            "compartment",
            "diagnosis_guess",
            "cell_type",
            "fraction_column",
            "mean_donor_fraction",
            "median_donor_fraction",
            "n_donors",
            "n_donors_with_fraction",
            "n_rois",
        ],
    )
    write_csv(
        donor_table,
        donor_rows,
        preferred=["compartment", "donor_guess", "diagnosis_guess", "n_rois"],
    )

    failures = []
    if input_summary["missing_roi_tables"]:
        failures.append("Missing deconvolution ROI tables")
    if input_summary["missing_manifest_tables"]:
        failures.append("Missing deconvolution celltype manifest tables")
    if not donor_rows:
        failures.append("No donor composition rows available")
    if not cell_types:
        failures.append("No cell types available to plot")

    output_paths = []
    heatmap_paths = []
    if failures:
        summary = {
            **input_summary,
            "labels": LABEL_ORDER,
            "cell_types": cell_types,
            "normalization": "stack values normalized to sum to 1" if normalize else "raw mean donor fractions",
            "n_donor_rows_total": len(donor_rows),
            "n_summary_rows": len(summary_rows),
            "output_paths": {
                "stacked": output_paths,
                "heatmap": heatmap_paths,
                "composition_table": str(composition_table),
                "donor_table": str(donor_table),
            },
            "failures": failures,
            "ready": False,
        }
        summary_path = meta_dir / f"{output_stem(compartments)}_plot_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        print()
        for path in [composition_table, donor_table, summary_path]:
            print(path)
        if args.strict:
            print(
                "Strict cell-composition plot failed: " + "; ".join(failures),
                file=sys.stderr,
            )
            return 2
        return 1

    cmap = plt.get_cmap("tab20")
    colors = {cell_type: cmap(index % cmap.N) for index, cell_type in enumerate(cell_types)}
    fig, axes = plt.subplots(
        1,
        len(compartments),
        figsize=(4.2 * len(compartments), 5.0),
        constrained_layout=True,
        squeeze=False,
    )
    plot_stacked_bars(
        axes=axes[0],
        summary_rows=summary_rows,
        compartments=compartments,
        labels=LABEL_ORDER,
        cell_types=cell_types,
        colors=colors,
        normalize=normalize,
    )
    handles, legend_labels = axes[0][-1].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.06),
        ncol=min(4, max(1, len(cell_types))),
        frameon=False,
        fontsize=8,
    )
    fig.suptitle(
        "GSE292993 donor-balanced NNLS cell composition by GeoMx compartment",
        fontsize=13,
        fontweight="bold",
        y=1.13,
    )

    heatmap_height = max(4.8, 0.42 * len(cell_types) + 2.2)
    heatmap_fig, heatmap_axis = plt.subplots(
        1,
        1,
        figsize=(1.25 * len(compartments) * len(LABEL_ORDER) + 2.5, heatmap_height),
        constrained_layout=True,
    )
    image = plot_heatmap(
        axis=heatmap_axis,
        summary_rows=summary_rows,
        compartments=compartments,
        labels=LABEL_ORDER,
        cell_types=cell_types,
        normalize=normalize,
        cmap="viridis",
    )
    heatmap_fig.colorbar(
        image,
        ax=heatmap_axis,
        label="mean donor fraction" + (" (normalized)" if normalize else ""),
        shrink=0.82,
    )

    for fmt in [item.strip().lower() for item in args.formats.split(",") if item.strip()]:
        stacked_path = figure_dir / f"{output_stem(compartments)}_stacked.{fmt}"
        heatmap_path = figure_dir / f"{output_stem(compartments)}_heatmap.{fmt}"
        fig.savefig(stacked_path, dpi=300 if fmt == "png" else None, bbox_inches="tight")
        heatmap_fig.savefig(heatmap_path, dpi=300 if fmt == "png" else None)
        output_paths.append(str(stacked_path))
        heatmap_paths.append(str(heatmap_path))
    plt.close(fig)
    plt.close(heatmap_fig)

    summary = {
        **input_summary,
        "labels": LABEL_ORDER,
        "cell_types": cell_types,
        "normalization": "stack values normalized to sum to 1" if normalize else "raw mean donor fractions",
        "n_donor_rows_total": len(donor_rows),
        "n_summary_rows": len(summary_rows),
        "output_paths": {
            "stacked": output_paths,
            "heatmap": heatmap_paths,
            "composition_table": str(composition_table),
            "donor_table": str(donor_table),
        },
        "failures": failures,
        "ready": not failures,
    }
    summary_path = meta_dir / f"{output_stem(compartments)}_plot_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    for path in [*output_paths, *heatmap_paths, composition_table, donor_table, summary_path]:
        print(path)

    if args.strict and failures:
        print(
            "Strict cell-composition plot failed: " + "; ".join(failures),
            file=sys.stderr,
        )
        return 2
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
