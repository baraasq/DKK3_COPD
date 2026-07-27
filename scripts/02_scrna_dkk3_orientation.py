#!/usr/bin/env python3
from __future__ import annotations

import json

from common import (
    configured_path,
    ensure_results_dirs,
    expression_vector,
    load_config,
    normalize_condition,
    resolve_column,
)


def require_dependencies():
    try:
        import anndata as ad
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import seaborn as sns
    except ImportError as exc:
        raise SystemExit(
            "Missing scRNA dependency. Install anndata, pandas, NumPy, "
            "matplotlib, and seaborn."
        ) from exc
    return ad, plt, np, pd, sns


def study_eligibility(frame, sample_col: str, study_col: str | None, pd):
    if study_col is None:
        return pd.DataFrame(
            [
                {
                    "study": "not available",
                    "Control": frame.loc[
                        frame["condition"] == "Control", sample_col
                    ].nunique(),
                    "COPD": frame.loc[
                        frame["condition"] == "COPD", sample_col
                    ].nunique(),
                    "contains_both": True,
                }
            ]
        )

    sample_metadata = frame[
        [sample_col, study_col, "condition"]
    ].drop_duplicates()
    table = (
        sample_metadata.groupby([study_col, "condition"], observed=True)[sample_col]
        .nunique()
        .unstack(fill_value=0)
        .reset_index()
        .rename(columns={study_col: "study"})
    )
    for condition in ["Control", "COPD"]:
        if condition not in table:
            table[condition] = 0
    table["contains_both"] = (table["Control"] > 0) & (table["COPD"] > 0)
    table["balanced_information"] = table[["Control", "COPD"]].min(axis=1)
    return table.sort_values(
        ["contains_both", "balanced_information", "COPD", "Control"],
        ascending=False,
    )


def choose_primary_study(eligibility):
    eligible = eligibility.loc[eligibility["contains_both"]]
    if eligible.empty or eligible.iloc[0]["study"] == "not available":
        return None
    return str(eligible.iloc[0]["study"])


def make_plot(summary, gene: str, path, plt, sns) -> None:
    ranking = (
        summary.groupby("cell_type", observed=True)["mean_expression"]
        .median()
        .sort_values(ascending=False)
    )
    top = list(ranking.head(20).index)
    plot_frame = summary.loc[summary["cell_type"].isin(top)].copy()
    plot_frame["cell_type"] = plot_frame["cell_type"].astype(str)
    order = [cell_type for cell_type in top if cell_type in set(plot_frame.cell_type)]

    height = max(6, 0.38 * len(order))
    figure, axis = plt.subplots(figsize=(11, height))
    sns.boxplot(
        data=plot_frame,
        x="mean_expression",
        y="cell_type",
        hue="condition",
        order=order,
        hue_order=["Control", "COPD"],
        showfliers=False,
        ax=axis,
    )
    sns.stripplot(
        data=plot_frame,
        x="mean_expression",
        y="cell_type",
        hue="condition",
        order=order,
        hue_order=["Control", "COPD"],
        dodge=True,
        alpha=0.65,
        size=3,
        linewidth=0,
        ax=axis,
    )
    handles, labels = axis.get_legend_handles_labels()
    axis.legend(
        handles[:2],
        labels[:2],
        title="Condition",
        frameon=False,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
    )
    axis.set_xlabel(f"Donor-level mean {gene} expression")
    axis.set_ylabel("")
    axis.set_title(
        f"{gene} expression by donor and cell type\n"
        "Descriptive orientation; no cell-level disease test"
    )
    figure.tight_layout()
    figure.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> int:
    ad, plt, np, pd, sns = require_dependencies()
    config = load_config()
    output = ensure_results_dirs(config)
    path = configured_path(
        config, "scrna_h5ad", environment_variable="COPD_SCRNA_H5AD"
    )
    if not path.exists():
        raise SystemExit(
            f"scRNA h5ad not found: {path}\n"
            "Set COPD_SCRNA_H5AD or update config/project.toml."
        )

    gene = config["project"]["gene"]
    backed = ad.read_h5ad(path, backed="r")
    try:
        if gene not in backed.var_names:
            raise SystemExit(f"{gene} is absent from {path}.")
        gene_data = backed[:, [gene]].to_memory()
    finally:
        if getattr(backed, "file", None) is not None:
            backed.file.close()

    obs_columns = list(map(str, gene_data.obs.columns))
    sample_col = resolve_column(
        obs_columns,
        config["metadata"]["sample_candidates"],
        label="sample/donor column",
    )
    condition_col = resolve_column(
        obs_columns,
        config["metadata"]["condition_candidates"],
        label="condition column",
    )
    cell_type_col = resolve_column(
        obs_columns,
        config["metadata"]["cell_type_candidates"],
        label="cell-type column",
    )
    study_col = resolve_column(
        obs_columns,
        config["metadata"]["study_candidates"],
        label="study column",
        required=False,
    )

    counts_layer = config["scrna"]["counts_layer"]
    expression_layer = config["scrna"]["expression_layer"]
    if counts_layer not in gene_data.layers:
        raise SystemExit(
            f"Raw-count layer {counts_layer!r} is required but absent. "
            f"Available layers: {list(gene_data.layers.keys())}"
        )
    if expression_layer not in gene_data.layers:
        raise SystemExit(
            f"Expression layer {expression_layer!r} is required but absent. "
            f"Available layers: {list(gene_data.layers.keys())}"
        )

    frame = gene_data.obs.copy()
    frame["condition"] = frame[condition_col].map(
        lambda value: normalize_condition(value, config)
    )
    frame = frame.loc[frame["condition"].isin(["Control", "COPD"])].copy()
    frame["cell_type"] = (
        frame[cell_type_col].astype("string").fillna("not annotated").astype(str)
    )
    frame["gene_counts"] = expression_vector(gene_data.layers[counts_layer])
    frame["gene_expression"] = expression_vector(
        gene_data.layers[expression_layer]
    )
    frame["gene_positive"] = frame["gene_counts"] > 0

    eligibility = study_eligibility(frame, sample_col, study_col, pd)
    eligibility.to_csv(
        output["tables"] / "scrna_copd_control_study_eligibility.csv",
        index=False,
    )
    primary_study = (
        choose_primary_study(eligibility)
        if config["scrna"]["primary_same_study"]
        else None
    )
    if primary_study and study_col:
        frame["primary_cohort"] = frame[study_col].astype(str) == primary_study
    else:
        frame["primary_cohort"] = True

    group_columns = [sample_col, "condition", "cell_type", "primary_cohort"]
    if study_col:
        group_columns.insert(1, study_col)
    summary = (
        frame.groupby(group_columns, observed=True)
        .agg(
            n_cells=("gene_expression", "size"),
            mean_expression=("gene_expression", "mean"),
            median_expression=("gene_expression", "median"),
            fraction_positive=("gene_positive", "mean"),
            pseudobulk_gene_counts=("gene_counts", "sum"),
        )
        .reset_index()
        .rename(columns={sample_col: "sample"})
    )
    summary.to_csv(
        output["tables"] / "scrna_dkk3_sample_celltype_summary.csv",
        index=False,
    )

    primary_summary = summary.loc[summary["primary_cohort"]].copy()
    thresholds = list(config["scrna"]["minimum_cells_per_sample_cell_type"])
    coverage_rows = []
    for threshold in thresholds:
        eligible = primary_summary.loc[primary_summary["n_cells"] >= threshold]
        counts = (
            eligible.groupby(["cell_type", "condition"], observed=True)["sample"]
            .nunique()
            .reset_index(name="n_samples")
        )
        counts["minimum_cells"] = int(threshold)
        coverage_rows.append(counts)
    coverage = pd.concat(coverage_rows, ignore_index=True)
    coverage.to_csv(
        output["tables"] / "scrna_dkk3_pseudobulk_coverage.csv",
        index=False,
    )

    make_plot(
        primary_summary,
        gene,
        output["figures"] / "scrna_dkk3_donor_celltype_orientation.png",
        plt,
        sns,
    )

    status = {
        "gene": gene,
        "n_cells": int(len(frame)),
        "n_samples": int(summary["sample"].nunique()),
        "n_cell_types": int(summary["cell_type"].nunique()),
        "primary_same_study": bool(config["scrna"]["primary_same_study"]),
        "primary_study": primary_study,
        "primary_samples": int(primary_summary["sample"].nunique()),
        "formal_inference": "not performed by this descriptive script",
    }
    (output["meta"] / "scrna_dkk3_orientation_status.json").write_text(
        json.dumps(status, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

