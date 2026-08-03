#!/usr/bin/env python3
"""Export barcode-keyed author annotations only after exact replay checkpoints pass."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import pickle
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_OUTPUT_DIR = "output/gse302339_author_exact"
DEFAULT_PREPROCESSING_LOG = "logs/gse302339_author_exact_preprocessing.log"
DEFAULT_ANNOTATION_LOG = "logs/gse302339_author_exact_annotation.log"
EXPECTED_SHAPE = (160620, 2323)
EXPECTED_FULL_CLUSTERS = 62
EXPECTED_PARENCHYMA_CLUSTERS = 41
EXPECTED_IMMUNE_CLUSTERS = 38
REQUIRED_PARENCHYMA_LABELS = {
    "AT1",
    "AT2",
    "Airway Epithelia",
    "Club cell",
    "Fibroblast",
    "Lymphatic Endothelia",
    "Mesothelia",
    "Smooth muscle",
    "Vascular Endothelia",
}
BARCODE_PATTERN = re.compile(r"^([ACGTN]+-\d+)")


def load_pickle(path: Path):
    with path.open("rb") as handle:
        return pickle.load(handle)


def sample_barcode_key(cell_id: str, obs_row) -> str:
    batch = scalar(obs_row.get("batch", ""))
    return f"{batch}::{raw_barcode(cell_id)}"


def labels_by_sample_barcode(adata, column: str) -> tuple[dict[str, str], list[str]]:
    if column not in adata.obs.columns:
        return {}, [f"annotation column missing: {column}"]
    labels: dict[str, str] = {}
    failures: list[str] = []
    for index, obs_row in adata.obs.iterrows():
        value = scalar(obs_row.get(column, ""))
        if not value:
            continue
        key = sample_barcode_key(str(index), obs_row)
        if key in labels:
            failures.append(f"duplicate sample/barcode annotation key in {column}: {key}")
            continue
        labels[key] = value
    return labels, failures


def raw_barcode(cell_id: str) -> str:
    match = BARCODE_PATTERN.match(cell_id)
    return match.group(1) if match else cell_id


def scalar(value) -> str:
    if value is None:
        return ""
    text = str(value)
    return "" if text in {"nan", "None", "<NA>"} else text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate exact replay outputs and export author labels keyed by sample and barcode."
    )
    parser.add_argument("--output-dir", default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--preprocessing-log", default=DEFAULT_PREPROCESSING_LOG)
    parser.add_argument("--annotation-log", default=DEFAULT_ANNOTATION_LOG)
    parser.add_argument(
        "--annotation-output",
        default="data/external/scrna_reference/gse302339_author_exact_cell_annotations.tsv.gz",
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    output_dir = project_path(args.output_dir)
    full_path = output_dir / "adata_harmony_annotated_cr8"
    parenchyma_path = output_dir / "parenchyma_harmony_annotated_cr8"
    immune_path = output_dir / "immune_harmony_annotated_cr8"
    annotation_output = project_path(args.annotation_output)
    preprocessing_log = project_path(args.preprocessing_log)
    annotation_log = project_path(args.annotation_log)
    results = ensure_results_dirs(load_config())
    summary_path = results["meta"] / "gse302339_author_exact_annotation_export_summary.json"
    counts_path = results["tables"] / "gse302339_author_exact_celltype_counts.csv"

    failures: list[str] = []
    if not preprocessing_log.exists():
        failures.append(f"missing exact-replay preprocessing log: {preprocessing_log}")
    else:
        preprocessing_text = preprocessing_log.read_text(
            encoding="utf-8", errors="replace"
        )
        harmony_iterations = [
            int(value)
            for value in re.findall(
                r"Converged after\s+(\d+)\s+iterations", preprocessing_text
            )
        ]
        if harmony_iterations != [26]:
            failures.append(
                f"Harmony checkpoint expected one 26-iteration convergence; observed {harmony_iterations}"
            )
    if not annotation_log.exists():
        failures.append(f"missing exact-replay annotation log: {annotation_log}")
    else:
        annotation_text = annotation_log.read_text(encoding="utf-8", errors="replace")
        if "Exact author major-lineage annotations completed and checkpointed." not in annotation_text:
            failures.append("annotation log lacks the exact-replay completion checkpoint")
    for path in (full_path, parenchyma_path, immune_path):
        if not path.exists():
            failures.append(f"missing exact-replay object: {path}")

    if failures:
        summary = {
            "ready": False,
            "failures": failures,
            "annotation_output": str(annotation_output),
            "preprocessing_log": str(preprocessing_log),
            "annotation_log": str(annotation_log),
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        if args.strict:
            raise SystemExit("Strict author annotation export failed: " + "; ".join(failures))
        return 1

    full = load_pickle(full_path)
    parenchyma = load_pickle(parenchyma_path)
    immune = load_pickle(immune_path)

    if tuple(full.shape) != EXPECTED_SHAPE:
        failures.append(f"full shape expected {EXPECTED_SHAPE}; observed {tuple(full.shape)}")
    full_cluster_ids = set(full.obs["leiden"].astype(str))
    expected_full_cluster_ids = {str(value) for value in range(EXPECTED_FULL_CLUSTERS)}
    if (
        int(full.obs["leiden"].nunique()) != EXPECTED_FULL_CLUSTERS
        or full_cluster_ids != expected_full_cluster_ids
    ):
        failures.append(
            f"full clusters expected IDs 0-{EXPECTED_FULL_CLUSTERS - 1}; "
            f"observed {sorted(full_cluster_ids)}"
        )
    parenchyma_cluster_ids = set(parenchyma.obs["leiden"].astype(str))
    expected_parenchyma_cluster_ids = {
        str(value) for value in range(EXPECTED_PARENCHYMA_CLUSTERS)
    }
    if (
        int(parenchyma.obs["leiden"].nunique()) != EXPECTED_PARENCHYMA_CLUSTERS
        or parenchyma_cluster_ids != expected_parenchyma_cluster_ids
    ):
        failures.append(
            f"parenchyma clusters expected IDs 0-{EXPECTED_PARENCHYMA_CLUSTERS - 1}; "
            f"observed {sorted(parenchyma_cluster_ids)}"
        )
    immune_cluster_ids = set(immune.obs["leiden"].astype(str))
    expected_immune_cluster_ids = {
        str(value) for value in range(EXPECTED_IMMUNE_CLUSTERS)
    }
    if (
        int(immune.obs["leiden"].nunique()) != EXPECTED_IMMUNE_CLUSTERS
        or immune_cluster_ids != expected_immune_cluster_ids
    ):
        failures.append(
            f"immune clusters expected IDs 0-{EXPECTED_IMMUNE_CLUSTERS - 1}; "
            f"observed {sorted(immune_cluster_ids)}"
        )

    full_key_counts = Counter(
        sample_barcode_key(str(index), obs_row)
        for index, obs_row in full.obs.iterrows()
    )
    broad_by_key = {
        sample_barcode_key(str(index), obs_row): scalar(
            obs_row.get("celltype_level1", "")
        )
        for index, obs_row in full.obs.iterrows()
    }
    observed_broad_labels = set(broad_by_key.values())
    expected_broad_labels = {"Parenchyma", "Immune", "Megakarycyte"}
    if observed_broad_labels != expected_broad_labels:
        failures.append(
            f"broad author labels expected {sorted(expected_broad_labels)}; "
            f"observed {sorted(observed_broad_labels)}"
        )
    duplicate_full_keys = sorted(
        key for key, count in full_key_counts.items() if count != 1
    )
    if duplicate_full_keys:
        failures.append(
            "full object contains non-unique sample/barcode keys; first duplicates: "
            + str(duplicate_full_keys[:10])
        )
    if any(key.startswith("::") for key in full_key_counts):
        failures.append("full object contains cells without a batch value")

    parenchyma_labels, parenchyma_key_failures = labels_by_sample_barcode(
        parenchyma, "parenchyma_celltype_level1"
    )
    immune_labels, immune_key_failures = labels_by_sample_barcode(
        immune, "immune_celltype_level1"
    )
    failures.extend(parenchyma_key_failures)
    failures.extend(immune_key_failures)
    if len(parenchyma_labels) != int(parenchyma.n_obs):
        failures.append(
            f"parenchyma label keys expected {parenchyma.n_obs}; observed {len(parenchyma_labels)}"
        )
    if len(immune_labels) != int(immune.n_obs):
        failures.append(
            f"immune label keys expected {immune.n_obs}; observed {len(immune_labels)}"
        )
    expected_parenchyma_keys = {
        key for key, broad_label in broad_by_key.items() if broad_label == "Parenchyma"
    }
    expected_immune_keys = {
        key for key, broad_label in broad_by_key.items() if broad_label == "Immune"
    }
    if set(parenchyma_labels) != expected_parenchyma_keys:
        failures.append(
            "parenchyma detailed-label membership does not match full-object broad Parenchyma membership"
        )
    if set(immune_labels) != expected_immune_keys:
        failures.append(
            "immune detailed-label membership does not match full-object broad Immune membership"
        )
    overlapping_lineage_keys = set(parenchyma_labels) & set(immune_labels)
    if overlapping_lineage_keys:
        failures.append(
            "cells appear in both parenchyma and immune objects; first keys: "
            + str(sorted(overlapping_lineage_keys)[:10])
        )
    detailed_required_keys = {
        key
        for key, broad_label in broad_by_key.items()
        if broad_label in {"Parenchyma", "Immune"}
    }
    detailed_label_keys = set(parenchyma_labels) | set(immune_labels)
    missing_detailed_keys = detailed_required_keys - detailed_label_keys
    extra_detailed_keys = (
        detailed_label_keys - detailed_required_keys
    )
    if missing_detailed_keys:
        failures.append(
            f"{len(missing_detailed_keys)} full-object cells lack a detailed author label"
        )
    if extra_detailed_keys:
        failures.append(
            f"{len(extra_detailed_keys)} detailed labels do not map to the full object"
        )
    observed_parenchyma_labels = set(parenchyma_labels.values())
    missing_labels = sorted(REQUIRED_PARENCHYMA_LABELS - observed_parenchyma_labels)
    if missing_labels:
        failures.append(f"required author parenchyma labels missing: {missing_labels}")

    if failures:
        summary = {
            "ready": False,
            "full_shape": list(full.shape),
            "parenchyma_cluster_count": int(parenchyma.obs["leiden"].nunique()),
            "immune_cluster_count": int(immune.obs["leiden"].nunique()),
            "failures": failures,
            "annotation_output": str(annotation_output),
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        if args.strict:
            raise SystemExit("Strict author annotation export failed: " + "; ".join(failures))
        return 1

    preferred_metadata = [
        "batch",
        "sample",
        "sample_id",
        "geo_accession",
        "patient",
        "donor",
        "lobe_emphysema",
        "lobe_emphysema_simple",
    ]
    fieldnames = [
        "cell_id",
        "raw_barcode",
        "sample_barcode_key",
        *preferred_metadata,
        "celltype_level1",
        "parenchyma_celltype_level1",
        "immune_celltype_level1",
        "author_celltype",
        "annotation_provenance",
    ]
    annotation_output.parent.mkdir(parents=True, exist_ok=True)
    counts: Counter[str] = Counter()
    with gzip.open(annotation_output, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for index, obs_row in full.obs.iterrows():
            cell_id = str(index)
            key = sample_barcode_key(cell_id, obs_row)
            broad = scalar(obs_row.get("celltype_level1", ""))
            parenchyma_label = parenchyma_labels.get(key, "")
            immune_label = immune_labels.get(key, "")
            author_celltype = parenchyma_label or immune_label or broad
            counts[author_celltype] += 1
            row = {
                "cell_id": cell_id,
                "raw_barcode": raw_barcode(cell_id),
                "sample_barcode_key": key,
                "celltype_level1": broad,
                "parenchyma_celltype_level1": parenchyma_label,
                "immune_celltype_level1": immune_label,
                "author_celltype": author_celltype,
                "annotation_provenance": "GSE302339 checkpoint-matched deposited-notebook replay",
            }
            for column in preferred_metadata:
                row[column] = scalar(obs_row.get(column, ""))
            writer.writerow(row)

    count_rows = [
        {"author_celltype": label, "n_cells": count}
        for label, count in sorted(counts.items())
    ]
    with counts_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["author_celltype", "n_cells"])
        writer.writeheader()
        writer.writerows(count_rows)

    summary = {
        "ready": True,
        "strategy": "sample-plus-barcode checkpoint-matched deposited-notebook annotation table",
        "n_cells": int(full.n_obs),
        "n_author_celltypes": len(counts),
        "full_shape": list(full.shape),
        "full_cluster_count": int(full.obs["leiden"].nunique()),
        "parenchyma_cluster_count": int(parenchyma.obs["leiden"].nunique()),
        "immune_cluster_count": int(immune.obs["leiden"].nunique()),
        "parenchyma_labels": sorted(observed_parenchyma_labels),
        "annotation_output": str(annotation_output),
        "celltype_counts": str(counts_path),
        "preprocessing_log": str(preprocessing_log),
        "annotation_log": str(annotation_log),
        "failures": [],
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))
    print()
    for path in (annotation_output, counts_path, summary_path):
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
