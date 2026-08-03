#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_OBJECT = "output/adata_harmony_annotated_cr8"
DEFAULT_FIGURE_DIR = "results/figures/gse302339_marker_program_features"
DEFAULT_FEATURE_GENES = [
    "AGER",
    "HOPX",
    "SFTPC",
    "ABCA3",
    "COL1A1",
    "PDGFRA",
    "PECAM1",
    "ACTA2",
    "FOXJ1",
    "SCGB1A1",
    "PTPRC",
    "NKG7",
]


def load_marker_module():
    path = Path(__file__).resolve().parent / "17_build_gse302339_marker_program_signatures.py"
    spec = importlib.util.spec_from_file_location("marker_program_signatures", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


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


def expression_vector(ref, gene: str):
    import numpy as np

    index_by_gene = {str(value).upper(): index for index, value in enumerate(ref.var_names)}
    index = index_by_gene.get(gene.upper())
    if index is None:
        return None
    values = ref.X[:, index]
    library = np.asarray(ref.X.sum(axis=1)).ravel()
    if hasattr(values, "toarray"):
        values = np.asarray(values.toarray()).ravel()
    else:
        values = np.asarray(values).ravel()
    cpm = values / np.maximum(library, 1) * 1_000_000
    return np.log1p(cpm)


def quantile_clip(values, upper_quantile: float):
    import numpy as np

    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return values, 0.0
    vmax = float(np.quantile(finite, upper_quantile))
    if vmax <= 0:
        vmax = float(finite.max())
    return np.clip(values, 0, vmax), vmax


def plot_label_umap(coords, labels: list[str], output_path: Path, *, point_size: float) -> dict:
    import matplotlib.pyplot as plt
    import numpy as np

    output_path.parent.mkdir(parents=True, exist_ok=True)
    labels_array = np.asarray(labels, dtype=object)
    preferred = [
        "Unassigned",
        "AT1",
        "AT2",
        "Fibroblast",
        "Endothelial",
        "Smooth muscle",
        "Airway epithelial",
        "Immune",
    ]
    present = [label for label in preferred if label in set(labels)]
    present.extend(sorted(label for label in set(labels) if label not in set(present)))
    colors = {
        "Unassigned": "#d3d3d3",
        "AT1": "#1f77b4",
        "AT2": "#aec7e8",
        "Fibroblast": "#ff7f0e",
        "Endothelial": "#98df8a",
        "Smooth muscle": "#ffbb78",
        "Airway epithelial": "#2ca02c",
        "Immune": "#d62728",
    }

    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    for label in present:
        mask = labels_array == label
        alpha = 0.12 if label == "Unassigned" else 0.75
        size = point_size * 0.7 if label == "Unassigned" else point_size
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=size,
            c=colors.get(label, "#333333"),
            alpha=alpha,
            linewidths=0,
            label=f"{label} (n={int(mask.sum())})",
        )
    ax.set_title("GSE302339 marker-program selected cells")
    ax.set_xlabel("UMAP1")
    ax.set_ylabel("UMAP2")
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), frameon=False, markerscale=4)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return {"path": str(output_path), "written": True}


def plot_feature_umaps(
    coords,
    ref,
    genes: list[str],
    output_path: Path,
    *,
    point_size: float,
    clip_quantile: float,
) -> tuple[dict, list[dict]]:
    import matplotlib.pyplot as plt
    import numpy as np

    rows = []
    present = []
    for gene in genes:
        vector = expression_vector(ref, gene)
        if vector is None:
            rows.append({"gene": gene, "present": "False", "max_log1p_cpm": ""})
            continue
        clipped, vmax = quantile_clip(vector, clip_quantile)
        present.append((gene, clipped, vmax))
        rows.append(
            {
                "gene": gene,
                "present": "True",
                "max_log1p_cpm": float(np.nanmax(vector)),
                f"q{int(clip_quantile * 100)}_log1p_cpm": vmax,
            }
        )
    if not present:
        return {"path": str(output_path), "written": False, "error": "no genes present"}, rows

    n_cols = 4
    n_rows = math.ceil(len(present) / n_cols)
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.2 * n_cols, 3.0 * n_rows),
        squeeze=False,
        constrained_layout=True,
    )
    for axis in axes.ravel():
        axis.axis("off")
    for axis, (gene, values, vmax) in zip(axes.ravel(), present):
        axis.axis("on")
        axis.scatter(coords[:, 0], coords[:, 1], s=point_size, c="#d7d7d7", alpha=0.08, linewidths=0)
        image = axis.scatter(
            coords[:, 0],
            coords[:, 1],
            s=point_size,
            c=values,
            cmap="magma",
            vmin=0,
            vmax=vmax if vmax > 0 else None,
            alpha=0.85,
            linewidths=0,
        )
        axis.set_title(gene)
        axis.set_xticks([])
        axis.set_yticks([])
        fig.colorbar(image, ax=axis, fraction=0.046, pad=0.02)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220)
    plt.close(fig)
    return {"path": str(output_path), "written": True, "n_genes_plotted": len(present)}, rows


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot UMAP marker-program labels and feature expression for GSE302339."
    )
    parser.add_argument("--input-object", default=DEFAULT_INPUT_OBJECT)
    parser.add_argument("--expression-source", choices=["raw", "X"], default="raw")
    parser.add_argument("--min-z", type=float, default=1.5)
    parser.add_argument("--min-margin-z", type=float, default=0.25)
    parser.add_argument("--min-detected-markers", type=int, default=2)
    parser.add_argument("--max-cells-per-cell-type", type=int, default=5000)
    parser.add_argument("--feature-gene", action="append", dest="feature_genes")
    parser.add_argument("--figure-dir", default=DEFAULT_FIGURE_DIR)
    parser.add_argument("--clip-quantile", type=float, default=0.99)
    parser.add_argument("--point-size", type=float, default=1.2)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    results = ensure_results_dirs(config)
    marker_module = load_marker_module()
    figure_dir = project_path(args.figure_dir)
    summary_path = results["meta"] / "gse302339_marker_program_feature_plot_summary.json"
    count_path = results["tables"] / "gse302339_marker_program_feature_plot_counts.csv"
    feature_path = results["tables"] / "gse302339_marker_program_feature_plot_genes.csv"

    adata, load_summary = marker_module.load_object(project_path(args.input_object))
    failures = []
    if adata is None:
        failures.append(f"Input object not loadable: {load_summary.get('status')}")
    elif "X_umap" not in getattr(adata, "obsm", {}):
        failures.append("Input object lacks adata.obsm['X_umap']")

    ref = None
    expression_summary = {}
    if not failures:
        ref, expression_summary = marker_module.expression_reference(adata, args.expression_source)
        if ref is None:
            failures.append(f"Expression source unavailable: {expression_summary.get('failure')}")

    if failures:
        summary = {
            "input_object": load_summary,
            "ready": False,
            "failures": failures,
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        if args.strict:
            print("Strict marker feature plotting failed: " + "; ".join(failures), file=sys.stderr)
            return 2
        return 1

    score_bundle, _ = marker_module.score_marker_programs(ref, marker_module.MARKER_SETS)
    assigned, selection_rows = marker_module.assign_marker_labels(
        score_bundle=score_bundle,
        min_z=args.min_z,
        min_margin_z=args.min_margin_z,
        min_detected_markers=args.min_detected_markers,
        max_cells_per_type=args.max_cells_per_cell_type,
        include_cell_types=list(marker_module.MARKER_SETS),
    )
    counts = Counter(assigned)
    n_total = int(len(assigned))
    n_selected = int(n_total - counts.get("Unassigned", 0))
    n_unassigned = int(counts.get("Unassigned", 0))
    count_rows = [
        {
            "cell_type": label,
            "n_cells": int(count),
            "fraction_of_total": int(count) / n_total if n_total else 0,
            "used_for_signature": str(label != "Unassigned"),
        }
        for label, count in sorted(counts.items())
    ]
    count_rows.append(
        {
            "cell_type": "__TOTAL_SELECTED__",
            "n_cells": n_selected,
            "fraction_of_total": n_selected / n_total if n_total else 0,
            "used_for_signature": "True",
        }
    )
    count_rows.append(
        {
            "cell_type": "__NEGLECTED_UNASSIGNED__",
            "n_cells": n_unassigned,
            "fraction_of_total": n_unassigned / n_total if n_total else 0,
            "used_for_signature": "False",
        }
    )
    write_csv(
        count_path,
        count_rows,
        preferred=["cell_type", "n_cells", "fraction_of_total", "used_for_signature"],
    )

    coords = adata.obsm["X_umap"]
    label_plot = plot_label_umap(
        coords,
        assigned,
        figure_dir / "gse302339_marker_program_selected_labels_umap.png",
        point_size=args.point_size,
    )
    feature_plot, feature_rows = plot_feature_umaps(
        coords,
        ref,
        args.feature_genes or DEFAULT_FEATURE_GENES,
        figure_dir / "gse302339_marker_program_feature_umaps.png",
        point_size=args.point_size,
        clip_quantile=args.clip_quantile,
    )
    write_csv(
        feature_path,
        feature_rows,
        preferred=["gene", "present", "max_log1p_cpm", f"q{int(args.clip_quantile * 100)}_log1p_cpm"],
    )

    summary = {
        "input_object": load_summary,
        "expression_source": expression_summary,
        "selection_parameters": {
            "min_z": args.min_z,
            "min_margin_z": args.min_margin_z,
            "min_detected_markers": args.min_detected_markers,
            "max_cells_per_cell_type": args.max_cells_per_cell_type,
        },
        "n_total_cells": n_total,
        "n_selected_cells": n_selected,
        "n_neglected_unassigned_cells": n_unassigned,
        "fraction_selected": n_selected / n_total if n_total else 0,
        "fraction_neglected_unassigned": n_unassigned / n_total if n_total else 0,
        "selection_rows": selection_rows,
        "outputs": {
            "summary": str(summary_path),
            "counts": str(count_path),
            "feature_genes": str(feature_path),
            "label_umap": label_plot["path"],
            "feature_umaps": feature_plot["path"],
        },
        "ready": True,
        "failures": [],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    for path in (
        summary_path,
        count_path,
        feature_path,
        Path(label_plot["path"]),
        Path(feature_plot["path"]),
    ):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
