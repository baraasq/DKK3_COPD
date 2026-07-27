#!/usr/bin/env python3
from __future__ import annotations

import json
import zipfile
from pathlib import Path

from common import (
    configured_path,
    ensure_results_dirs,
    load_config,
    normalize_condition,
    resolve_column,
    scrinshot_cell_map_members,
)


def require_dependencies():
    try:
        import anndata as ad
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from sklearn.neighbors import NearestNeighbors
    except ImportError as exc:
        raise SystemExit(
            "Missing spatial-analysis dependency. Install anndata, scanpy, "
            "pandas, matplotlib, scikit-learn, igraph, and leidenalg."
        ) from exc
    return ad, plt, np, pd, sc, NearestNeighbors


def first_available(frame, candidates: list[str], *, label: str) -> str:
    return resolve_column(frame.columns, candidates, label=label, required=True)


def load_cell_maps(path: Path, config: dict, pd) -> "object":
    if not path.exists():
        raise SystemExit(
            f"SCRINSHOT archive not found: {path}\n"
            "Set COPD_SCRINSHOT_ZIP or update config/project.toml."
        )
    if not zipfile.is_zipfile(path):
        raise SystemExit(f"Configured SCRINSHOT path is not a ZIP archive: {path}")

    frames = []
    with zipfile.ZipFile(path) as archive:
        members = scrinshot_cell_map_members(archive.namelist())
        if not members:
            raise SystemExit("No processed cell-map CSVs found in the archive.")

        for member in members:
            with archive.open(member) as handle:
                frame = pd.read_csv(handle, low_memory=False)

            x_col, y_col = [
                resolve_column(
                    frame.columns,
                    [candidate],
                    label=f"{candidate} coordinate in {member}",
                )
                for candidate in config["scrinshot"]["coordinate_columns"]
            ]
            donor_col = first_available(
                frame,
                config["scrinshot"]["donor_candidates"],
                label=f"donor in {member}",
            )
            condition_col = first_available(
                frame,
                config["scrinshot"]["condition_candidates"],
                label=f"condition in {member}",
            )
            sample_col = first_available(
                frame,
                config["scrinshot"]["section_candidates"],
                label=f"sample/section in {member}",
            )

            cell_type_columns = [
                resolve_column(
                    frame.columns,
                    [candidate],
                    label=f"{candidate} in {member}",
                    required=False,
                )
                for candidate in config["scrinshot"]["cell_type_candidates"]
            ]
            cell_type_columns = [
                column for column in cell_type_columns if column is not None
            ]
            if not cell_type_columns:
                raise SystemExit(f"No cell-type annotation found in {member}.")

            annotation = frame[cell_type_columns[0]].astype("string")
            for column in cell_type_columns[1:]:
                fallback = frame[column].astype("string")
                annotation = annotation.mask(
                    annotation.isna()
                    | annotation.str.strip().isin(["", "nan", "NA", "None"]),
                    fallback,
                )

            section_id = Path(member).stem
            standardized = pd.DataFrame(
                {
                    "section_id": section_id,
                    "source_member": member,
                    "sample": frame[sample_col].astype(str),
                    "donor": frame[donor_col].astype(str),
                    "condition_original": frame[condition_col].astype(str),
                    "x": pd.to_numeric(frame[x_col], errors="coerce"),
                    "y": pd.to_numeric(frame[y_col], errors="coerce"),
                    "cell_type": annotation.fillna("not annotated").astype(str),
                }
            )
            standardized["condition"] = standardized["condition_original"].map(
                lambda value: normalize_condition(value, config)
            )
            standardized = standardized.dropna(subset=["x", "y"]).reset_index(
                drop=True
            )
            standardized["cell_id"] = [
                f"{section_id}:{index}" for index in range(len(standardized))
            ]
            frames.append(standardized)

    cells = pd.concat(frames, ignore_index=True)
    if cells["cell_id"].duplicated().any():
        raise RuntimeError("Cell identifiers are not unique after concatenation.")
    return cells


def build_neighbor_profiles(cells, k: int, np, pd, NearestNeighbors):
    cell_types = sorted(cells["cell_type"].unique())
    type_index = {cell_type: index for index, cell_type in enumerate(cell_types)}
    profiles = np.zeros((len(cells), len(cell_types)), dtype=np.float32)

    for section_id, section in cells.groupby("section_id", sort=False):
        if len(section) < 2:
            continue
        section_k = min(k, len(section) - 1)
        coordinates = section[["x", "y"]].to_numpy(dtype=float)
        model = NearestNeighbors(n_neighbors=section_k + 1)
        model.fit(coordinates)
        indices = model.kneighbors(
            coordinates, return_distance=False
        )[:, 1:]
        section_types = section["cell_type"].to_numpy()
        section_rows = section.index.to_numpy()

        for local_row, neighbor_rows in enumerate(indices):
            counts = np.zeros(len(cell_types), dtype=np.float32)
            for neighbor_type in section_types[neighbor_rows]:
                counts[type_index[neighbor_type]] += 1
            profiles[section_rows[local_row], :] = counts

    profile_frame = pd.DataFrame(
        profiles,
        index=cells["cell_id"],
        columns=cell_types,
    )
    return profile_frame


def plot_sections(cells, output_dir: Path, plt) -> None:
    labels = sorted(cells["neighborhood"].astype(str).unique())
    palette = {
        label: plt.get_cmap("tab20")(index % 20)
        for index, label in enumerate(labels)
    }

    for section_id, frame in cells.groupby("section_id", sort=False):
        figure, axis = plt.subplots(figsize=(8, 8))
        for label, group in frame.groupby("neighborhood", sort=False):
            axis.scatter(
                group["x"],
                group["y"],
                s=2,
                linewidths=0,
                alpha=0.8,
                color=palette[str(label)],
                label=str(label),
            )
        axis.invert_yaxis()
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        condition = ", ".join(sorted(frame["condition"].unique()))
        donor = ", ".join(sorted(frame["donor"].unique()))
        axis.set_title(f"{section_id} | {condition} | donor {donor}")
        axis.legend(
            title="Neighborhood",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            markerscale=4,
            frameon=False,
            fontsize=7,
        )
        figure.tight_layout()
        figure.savefig(
            output_dir / f"scrinshot_neighborhoods_{section_id}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close(figure)


def main() -> int:
    ad, plt, np, pd, sc, NearestNeighbors = require_dependencies()
    config = load_config()
    output = ensure_results_dirs(config)
    archive_path = configured_path(
        config,
        "scrinshot_zip",
        environment_variable="COPD_SCRINSHOT_ZIP",
    )

    cells = load_cell_maps(archive_path, config, pd)
    profiles = build_neighbor_profiles(
        cells,
        int(config["scrinshot"]["neighbors_per_cell"]),
        np,
        pd,
        NearestNeighbors,
    )

    neighborhood = ad.AnnData(
        X=profiles.to_numpy(),
        obs=cells.set_index("cell_id"),
        var=pd.DataFrame(index=profiles.columns),
    )
    sc.pp.normalize_total(neighborhood, target_sum=1e4)
    sc.pp.log1p(neighborhood)
    n_components = min(30, neighborhood.n_vars - 1, neighborhood.n_obs - 1)
    if n_components < 2:
        raise SystemExit("Too few cells or cell types to cluster neighborhoods.")
    sc.pp.pca(neighborhood, n_comps=n_components)
    sc.pp.neighbors(neighborhood, n_neighbors=20, n_pcs=n_components)
    try:
        sc.tl.leiden(
            neighborhood,
            resolution=float(config["scrinshot"]["leiden_resolution"]),
            key_added="neighborhood",
        )
    except ImportError as exc:
        raise SystemExit(
            "Leiden clustering requires igraph and leidenalg."
        ) from exc
    sc.tl.umap(neighborhood)

    cluster_sizes = neighborhood.obs["neighborhood"].value_counts()
    minimum_size = int(config["scrinshot"]["minimum_neighborhood_size"])
    neighborhood.obs["neighborhood_retained"] = (
        neighborhood.obs["neighborhood"].map(cluster_sizes).astype(int)
        >= minimum_size
    )
    neighborhood.obs["neighborhood_plot"] = neighborhood.obs[
        "neighborhood"
    ].astype(str)
    rare = ~neighborhood.obs["neighborhood_retained"]
    neighborhood.obs.loc[rare, "neighborhood_plot"] = (
        "rare_" + neighborhood.obs.loc[rare, "neighborhood"].astype(str)
    )

    cells = neighborhood.obs.reset_index(names="cell_id")
    cells["neighborhood"] = cells["neighborhood_plot"].astype(str)
    cells.to_csv(
        output["tables"] / "scrinshot_cells_with_neighborhoods.csv.gz",
        index=False,
    )

    composition = (
        cells.groupby(["neighborhood", "cell_type"], observed=True)
        .size()
        .reset_index(name="n_cells")
    )
    composition["proportion_within_neighborhood"] = composition.groupby(
        "neighborhood"
    )["n_cells"].transform(lambda values: values / values.sum())
    composition.to_csv(
        output["tables"] / "scrinshot_neighborhood_celltype_composition.csv",
        index=False,
    )

    donor_counts = (
        cells.groupby(
            ["donor", "condition", "neighborhood"], observed=True
        )
        .size()
        .reset_index(name="n_cells")
    )
    donor_counts["proportion_within_donor"] = donor_counts.groupby(
        ["donor", "condition"]
    )["n_cells"].transform(lambda values: values / values.sum())
    donor_counts.to_csv(
        output["tables"] / "scrinshot_donor_neighborhood_proportions.csv",
        index=False,
    )

    neighborhood.write_h5ad(output["root"] / "scrinshot_neighborhoods.h5ad")
    plot_sections(cells, output["figures"], plt)

    status = {
        "n_cells": int(neighborhood.n_obs),
        "n_cell_types": int(neighborhood.n_vars),
        "n_sections": int(cells["section_id"].nunique()),
        "n_donors": int(cells["donor"].nunique()),
        "n_neighborhoods": int(cells["neighborhood"].nunique()),
        "neighbors_per_cell": int(config["scrinshot"]["neighbors_per_cell"]),
        "minimum_neighborhood_size": minimum_size,
        "dkk3_measured": False,
        "interpretation": "cellular-neighborhood context only",
    }
    (output["meta"] / "scrinshot_neighborhood_status.json").write_text(
        json.dumps(status, indent=2),
        encoding="utf-8",
    )
    print(json.dumps(status, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
