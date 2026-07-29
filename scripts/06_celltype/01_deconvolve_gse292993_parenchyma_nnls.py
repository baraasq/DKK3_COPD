#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config


BIOLOGICAL_LABELS = {"COPD", "Control"}


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


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


def safe_name(value: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_").lower()
    return text or "cell_type"


def median(values: list[float]) -> float | None:
    return statistics.median(values) if values else None


def rank(values: list[float]) -> list[float]:
    indexed = sorted(enumerate(values), key=lambda item: item[1])
    output = [0.0] * len(values)
    index = 0
    while index < len(indexed):
        end = index + 1
        while end < len(indexed) and indexed[end][1] == indexed[index][1]:
            end += 1
        average_rank = (index + end - 1) / 2 + 1
        for original, _ in indexed[index:end]:
            output[original] = average_rank
        index = end
    return output


def pearson(x_values: list[float], y_values: list[float]) -> float | None:
    if len(x_values) < 3 or len(y_values) != len(x_values):
        return None
    x_mean = statistics.mean(x_values)
    y_mean = statistics.mean(y_values)
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    x_denominator = math.sqrt(sum((x - x_mean) ** 2 for x in x_values))
    y_denominator = math.sqrt(sum((y - y_mean) ** 2 for y in y_values))
    if x_denominator == 0 or y_denominator == 0:
        return None
    return numerator / (x_denominator * y_denominator)


def spearman(x_values: list[float], y_values: list[float]) -> float | None:
    return pearson(rank(x_values), rank(y_values))


def read_signature_matrix(path: Path) -> tuple[list[str], list[str], list[list[float]]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        genes = header[1:]
        cell_types = []
        matrix = []
        for row in reader:
            if not row:
                continue
            cell_types.append(row[0])
            matrix.append([float(value) if value not in ("", None) else 0.0 for value in row[1:]])
    return cell_types, genes, matrix


def roi_metadata(rows: list[dict], compartment: str) -> dict[str, dict]:
    output = {}
    for row in rows:
        if not as_bool(row.get("include_qc")):
            continue
        if row.get("diagnosis_group") not in BIOLOGICAL_LABELS:
            continue
        if row.get("compartment_guess") != compartment:
            continue
        if row.get("donor_guess") in (None, "", "unknown"):
            continue
        geo = row.get("geo_accession")
        if geo:
            output[geo] = row
    return output


def read_geomx_matrix(
    path: Path,
    *,
    roi_ids: set[str],
    signature_genes: list[str],
) -> tuple[list[str], dict[str, list[float]]]:
    signature_set = set(signature_genes)
    with open_maybe_gzip(path, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        feature_columns = header[2:]
        common_genes = [gene for gene in signature_genes if gene in feature_columns]
        column_by_gene = {gene: index + 2 for index, gene in enumerate(feature_columns)}
        indices = [column_by_gene[gene] for gene in common_genes]
        matrix = {}
        for row in reader:
            if not row:
                continue
            geo = row[0]
            if geo not in roi_ids:
                continue
            matrix[geo] = [float(row[index]) if row[index] else 0.0 for index in indices]
    return common_genes, matrix


def subset_signature(
    *,
    signature_genes: list[str],
    common_genes: list[str],
    signature_matrix: list[list[float]],
) -> list[list[float]]:
    index_by_gene = {gene: index for index, gene in enumerate(signature_genes)}
    indices = [index_by_gene[gene] for gene in common_genes]
    return [[row[index] for index in indices] for row in signature_matrix]


def nnls_coefficients(signature_by_celltype: list[list[float]], roi_vector: list[float]) -> tuple[list[float], float, str]:
    import numpy as np

    signature = np.asarray(signature_by_celltype, dtype=float)
    target = np.asarray(roi_vector, dtype=float)
    design = signature.T
    try:
        from scipy.optimize import nnls

        coefficients, residual = nnls(design, target)
        method = "scipy.optimize.nnls"
    except Exception:
        coefficients, *_ = np.linalg.lstsq(design, target, rcond=None)
        coefficients = np.clip(coefficients, 0, None)
        residual = float(np.linalg.norm(design @ coefficients - target))
        method = "numpy_lstsq_nonnegative_clip"
    total = float(coefficients.sum())
    if total > 0:
        fractions = coefficients / total
    else:
        fractions = coefficients
    return [float(value) for value in fractions], float(residual), method


def dkk3_by_roi(rows: list[dict]) -> dict[str, dict]:
    return {row["geo_accession"]: row for row in rows if row.get("geo_accession")}


def donor_summaries(rows: list[dict], fraction_columns: list[str]) -> list[dict]:
    groups: dict[tuple[str, str, str], list[dict]] = defaultdict(list)
    for row in rows:
        key = (
            row.get("donor_guess", "unknown"),
            row.get("diagnosis_guess", "unknown"),
            row.get("diagnosis_group", "unknown"),
        )
        groups[key].append(row)

    output = []
    for (donor, diagnosis, group), group_rows in sorted(groups.items()):
        item = {
            "donor_guess": donor,
            "diagnosis_guess": diagnosis,
            "diagnosis_group": group,
            "n_rois": len(group_rows),
        }
        for column in fraction_columns:
            values = [value for row in group_rows if (value := as_float(row, column)) is not None]
            item[f"median_{column}"] = median(values)
        for column in ("log1p_dkk3_cpm", "dkk3_cpm", "dkk3_count"):
            values = [value for row in group_rows if (value := as_float(row, column)) is not None]
            item[f"median_{column}"] = median(values)
        output.append(item)
    return output


def correlation_rows(rows: list[dict], fraction_columns: list[str], *, unit: str) -> list[dict]:
    output = []
    dkk3_column = "log1p_dkk3_cpm" if unit == "roi" else "median_log1p_dkk3_cpm"
    for column in fraction_columns:
        x_values = []
        y_values = []
        for row in rows:
            x = as_float(row, column)
            y = as_float(row, dkk3_column)
            if x is None or y is None:
                continue
            x_values.append(x)
            y_values.append(y)
        output.append(
            {
                "unit": unit,
                "cell_fraction_column": column,
                "dkk3_column": dkk3_column,
                "n_pairs": len(x_values),
                "spearman_rho": spearman(x_values, y_values),
            }
        )
    return output


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run baseline NNLS deconvolution for QC-passing GSE292993 parenchymal ROIs."
    )
    parser.add_argument("--compartment", default="parenchyma")
    parser.add_argument("--signature-csv")
    parser.add_argument("--min-overlap-genes", type=int, default=500)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    processed_dir = configured_path(config, "geomx_processed_dir")
    signature_path = Path(args.signature_csv).expanduser() if args.signature_csv else (
        processed_dir / "gse302339_scrna_reference_signatures_logcpm.csv"
    )
    geomx_matrix_path = processed_dir / "gse292993_geomx_logcpm_by_roi.tsv.gz"
    prefix = f"gse292993_{args.compartment}_nnls_deconvolution"
    summary_path = meta_dir / f"{prefix}_summary.json"

    initial_failures = []
    if not signature_path.exists():
        initial_failures.append(f"Signature matrix not found: {signature_path}")
    if not geomx_matrix_path.exists():
        initial_failures.append(f"GeoMx logCPM matrix not found: {geomx_matrix_path}")
    if initial_failures:
        summary = {
            "dataset": config["project"]["primary_dataset"],
            "compartment": args.compartment,
            "signature_csv": str(signature_path),
            "geomx_matrix": str(geomx_matrix_path),
            "ready_for_deconvolution": False,
            "failures": initial_failures,
        }
        summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
        print(json.dumps(summary, indent=2))
        print()
        print(summary_path)
        if args.strict:
            print("Strict NNLS deconvolution failed: " + "; ".join(initial_failures), file=sys.stderr)
            return 2
        return 1

    cell_types, signature_genes, signature_matrix = read_signature_matrix(signature_path)
    metadata_by_roi = roi_metadata(
        read_csv(table_dir / "gse292993_roi_qc_flags.csv"), args.compartment
    )
    common_genes, geomx_matrix = read_geomx_matrix(
        geomx_matrix_path,
        roi_ids=set(metadata_by_roi),
        signature_genes=signature_genes,
    )
    signature_subset = subset_signature(
        signature_genes=signature_genes,
        common_genes=common_genes,
        signature_matrix=signature_matrix,
    )
    dkk3_rows = dkk3_by_roi(read_csv(table_dir / "gse292993_dkk3_roi_signal.csv"))

    safe_columns = [f"fraction_{safe_name(cell_type)}" for cell_type in cell_types]
    cell_type_manifest = [
        {
            "cell_type": cell_type,
            "fraction_column": column,
        }
        for cell_type, column in zip(cell_types, safe_columns)
    ]
    roi_rows = []
    methods = set()
    for geo, metadata in sorted(metadata_by_roi.items()):
        vector = geomx_matrix.get(geo)
        if vector is None:
            continue
        fractions, residual, method = nnls_coefficients(signature_subset, vector)
        methods.add(method)
        row = {
            "geo_accession": geo,
            "donor_guess": metadata.get("donor_guess"),
            "diagnosis_guess": metadata.get("diagnosis_guess"),
            "diagnosis_group": metadata.get("diagnosis_group"),
            "compartment_guess": metadata.get("compartment_guess"),
            "n_genes_used": len(common_genes),
            "nnls_residual": residual,
            "nnls_method": method,
        }
        for column, value in zip(safe_columns, fractions):
            row[column] = value
        if geo in dkk3_rows:
            for column in ("dkk3_count", "dkk3_cpm", "log1p_dkk3_cpm"):
                row[column] = dkk3_rows[geo].get(column)
        roi_rows.append(row)

    donor_rows = donor_summaries(roi_rows, safe_columns)
    donor_fraction_columns = [f"median_{column}" for column in safe_columns]
    correlations = [
        *correlation_rows(roi_rows, safe_columns, unit="roi"),
        *correlation_rows(donor_rows, donor_fraction_columns, unit="donor"),
    ]

    roi_path = table_dir / f"{prefix}_roi.csv"
    donor_path = table_dir / f"{prefix}_donor.csv"
    correlations_path = table_dir / f"{prefix}_dkk3_correlations.csv"
    celltype_path = table_dir / f"{prefix}_celltype_manifest.csv"
    write_csv(
        roi_path,
        roi_rows,
        preferred=[
            "geo_accession",
            "donor_guess",
            "diagnosis_guess",
            "diagnosis_group",
            "compartment_guess",
            "n_genes_used",
            "nnls_residual",
            "dkk3_count",
            "dkk3_cpm",
            "log1p_dkk3_cpm",
            *safe_columns,
        ],
    )
    write_csv(
        donor_path,
        donor_rows,
        preferred=[
            "donor_guess",
            "diagnosis_guess",
            "diagnosis_group",
            "n_rois",
            "median_log1p_dkk3_cpm",
            *donor_fraction_columns,
        ],
    )
    write_csv(correlations_path, correlations)
    write_csv(celltype_path, cell_type_manifest, preferred=["cell_type", "fraction_column"])
    summary = {
        "dataset": config["project"]["primary_dataset"],
        "compartment": args.compartment,
        "signature_csv": str(signature_path),
        "geomx_matrix": str(geomx_matrix_path),
        "n_signature_cell_types": len(cell_types),
        "cell_types": cell_types,
        "n_signature_genes": len(signature_genes),
        "n_overlap_genes_used": len(common_genes),
        "first_overlap_genes": common_genes[:20],
        "n_candidate_rois": len(metadata_by_roi),
        "n_deconvolved_rois": len(roi_rows),
        "n_deconvolved_donors": len({row.get("donor_guess") for row in roi_rows}),
        "nnls_methods_used": sorted(methods),
        "output_paths": {
            "roi": str(roi_path),
            "donor": str(donor_path),
            "dkk3_correlations": str(correlations_path),
            "celltype_manifest": str(celltype_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    for path in (roi_path, donor_path, correlations_path, celltype_path, summary_path):
        print(path)

    failures = []
    if len(common_genes) < args.min_overlap_genes:
        failures.append(
            f"Only {len(common_genes)} overlapping genes available; expected at least {args.min_overlap_genes}"
        )
    if not roi_rows:
        failures.append("No ROIs deconvolved")
    if args.strict and failures:
        print("Strict NNLS deconvolution failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
