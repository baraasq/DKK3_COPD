#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config


DEFAULT_DATASET = "GSE292993"


def load_dcc_module():
    path = Path(__file__).resolve().parents[1] / "03_geomx_qc" / "02_extract_gse292993_dcc_qc.py"
    spec = importlib.util.spec_from_file_location("gse292993_dcc_qc", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load DCC parser from {path}")
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


def open_maybe_gzip(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, newline="")
    return path.open(mode, newline="")


def is_control_feature(row: dict) -> bool:
    text = " ".join(
        str(row.get(column) or "")
        for column in ("target", "code_class")
    ).casefold()
    return any(token in text for token in ("negative", "control", "no template"))


def clean_target(value: object) -> str | None:
    text = str(value or "").strip()
    if not text:
        return None
    return text


def feature_manifest(pkc_rows: list[dict], *, include_controls: bool) -> tuple[list[dict], dict[str, str]]:
    code_to_target: dict[str, str] = {}
    target_rows: dict[str, dict] = {}
    code_to_targets: defaultdict[str, set[str]] = defaultdict(set)

    for row in pkc_rows:
        code_id = str(row.get("code_id") or "").strip()
        target = clean_target(row.get("target"))
        if not code_id or not target:
            continue
        is_control = is_control_feature(row)
        if is_control and not include_controls:
            continue
        code_to_targets[code_id].add(target)
        code_to_target[code_id] = target
        item = target_rows.setdefault(
            target,
            {
                "target": target,
                "n_codes": 0,
                "code_ids": [],
                "code_classes": set(),
                "is_control_feature": False,
            },
        )
        item["n_codes"] += 1
        item["code_ids"].append(code_id)
        if row.get("code_class"):
            item["code_classes"].add(str(row["code_class"]))
        item["is_control_feature"] = item["is_control_feature"] or is_control

    rows = []
    for target in sorted(target_rows):
        item = target_rows[target]
        rows.append(
            {
                "target": target,
                "n_codes": item["n_codes"],
                "code_ids": ";".join(sorted(set(item["code_ids"]))),
                "code_classes": ";".join(sorted(item["code_classes"])),
                "is_control_feature": str(bool(item["is_control_feature"])),
            }
        )
    ambiguous = {
        code_id: sorted(targets)
        for code_id, targets in code_to_targets.items()
        if len(targets) > 1
    }
    for code_id in ambiguous:
        code_to_target.pop(code_id, None)
    return rows, code_to_target


def aggregate_counts(code_counts: dict[str, int], code_to_target: dict[str, str]) -> dict[str, int]:
    target_counts: defaultdict[str, int] = defaultdict(int)
    for code_id, count in code_counts.items():
        target = code_to_target.get(code_id)
        if target:
            target_counts[target] += int(count)
    return dict(target_counts)


def write_matrix(
    path: Path,
    roi_rows: list[dict],
    features: list[str],
    matrix: dict[str, dict[str, int | float]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_maybe_gzip(path, "wt") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["geo_accession", "dcc_id", *features])
        for row in roi_rows:
            key = str(row["geo_accession"])
            values = matrix.get(key, {})
            writer.writerow(
                [
                    key,
                    row.get("dcc_id", ""),
                    *[values.get(feature, 0) for feature in features],
                ]
            )


def logcpm_matrix(counts_matrix: dict[str, dict[str, int]], features: list[str]) -> dict[str, dict[str, float]]:
    output: dict[str, dict[str, float]] = {}
    for roi, counts in counts_matrix.items():
        total = sum(float(counts.get(feature, 0)) for feature in features)
        output[roi] = {}
        for feature in features:
            count = float(counts.get(feature, 0))
            cpm = count / total * 1_000_000 if total > 0 else 0.0
            output[roi][feature] = math.log1p(cpm)
    return output


def matrix_summary(
    *,
    dataset: str,
    dcc_count: int,
    features: list[str],
    roi_rows: list[dict],
    include_controls: bool,
    counts_matrix: dict[str, dict[str, int]],
) -> dict:
    nonzero_by_roi = [
        sum(1 for value in counts_matrix.get(str(row["geo_accession"]), {}).values() if value)
        for row in roi_rows
    ]
    selected_totals = [
        sum(counts_matrix.get(str(row["geo_accession"]), {}).values())
        for row in roi_rows
    ]
    def quantiles(values: list[int]) -> dict:
        if not values:
            return {"n": 0}
        sorted_values = sorted(values)
        return {
            "n": len(values),
            "min": sorted_values[0],
            "q25": sorted_values[int((len(values) - 1) * 0.25)],
            "median": sorted_values[int((len(values) - 1) * 0.50)],
            "q75": sorted_values[int((len(values) - 1) * 0.75)],
            "max": sorted_values[-1],
        }
    return {
        "dataset": dataset,
        "include_controls": include_controls,
        "dcc_count": dcc_count,
        "n_roi_rows": len(roi_rows),
        "n_features": len(features),
        "first_features": features[:10],
        "nonzero_features_per_roi": quantiles(nonzero_by_roi),
        "selected_feature_total_counts_per_roi": quantiles(selected_totals),
    }


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Build full GeoMx WTA target-by-ROI count and logCPM matrices from GSE292993 DCC files."
    )
    parser.add_argument("--include-controls", action="store_true")
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    processed_dir = configured_path(config, "geomx_processed_dir")
    processed_dir.mkdir(parents=True, exist_ok=True)

    dcc = load_dcc_module()
    pkc_summary, pkc_rows = dcc.parse_pkc_code_map(
        dcc.find_pkc_file(config), config["project"]["gene"]
    )
    manifest_rows, code_to_target = feature_manifest(
        pkc_rows, include_controls=args.include_controls
    )
    features = [row["target"] for row in manifest_rows]
    paths = dcc.dcc_paths(config)

    roi_rows = []
    counts_matrix: dict[str, dict[str, int]] = {}
    n_codes_observed: Counter[int] = Counter()
    for path in paths:
        qc_row, code_counts = dcc.parse_dcc(path, set(), set())
        geo_accession = qc_row.get("geo_accession") or qc_row.get("dcc_id")
        qc_row["geo_accession"] = geo_accession
        roi_rows.append(qc_row)
        target_counts = aggregate_counts(code_counts, code_to_target)
        counts_matrix[str(geo_accession)] = target_counts
        n_codes_observed[len(code_counts)] += 1

    log_matrix = logcpm_matrix(counts_matrix, features)
    counts_path = processed_dir / "gse292993_geomx_counts_by_roi.tsv.gz"
    logcpm_path = processed_dir / "gse292993_geomx_logcpm_by_roi.tsv.gz"
    manifest_path = table_dir / "gse292993_geomx_feature_manifest.csv"
    roi_manifest_path = table_dir / "gse292993_geomx_matrix_roi_manifest.csv"
    summary_path = meta_dir / "gse292993_geomx_matrix_summary.json"

    write_matrix(counts_path, roi_rows, features, counts_matrix)
    write_matrix(logcpm_path, roi_rows, features, log_matrix)
    write_csv(
        manifest_path,
        manifest_rows,
        preferred=["target", "n_codes", "code_classes", "is_control_feature", "code_ids"],
    )
    write_csv(
        roi_manifest_path,
        roi_rows,
        preferred=["geo_accession", "dcc_id", "dcc_filename", "total_code_counts", "n_code_counts"],
    )
    summary = matrix_summary(
        dataset=DEFAULT_DATASET,
        dcc_count=len(paths),
        features=features,
        roi_rows=roi_rows,
        include_controls=args.include_controls,
        counts_matrix=counts_matrix,
    )
    summary.update(
        pkc=pkc_summary,
        code_to_target_count=len(code_to_target),
        n_codes_observed_distribution=dict(sorted(n_codes_observed.items())),
        output_paths={
            "counts_by_roi": str(counts_path),
            "logcpm_by_roi": str(logcpm_path),
            "feature_manifest": str(manifest_path),
            "roi_manifest": str(roi_manifest_path),
        },
    )
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    for path in (counts_path, logcpm_path, manifest_path, roi_manifest_path, summary_path):
        print(path)

    failures = []
    if not paths:
        failures.append("No DCC files found")
    if pkc_summary.get("status") != "ok":
        failures.append("PKC file missing or unreadable")
    if len(features) < 1000:
        failures.append(f"Suspiciously few GeoMx features resolved: {len(features)}")
    if args.strict and failures:
        print("Strict GeoMx matrix build failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
