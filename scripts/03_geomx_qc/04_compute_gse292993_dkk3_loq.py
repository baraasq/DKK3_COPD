#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config


DEFAULT_GSE = "GSE292993"


def load_dcc_qc_module():
    path = Path(__file__).resolve().with_name("02_extract_gse292993_dcc_qc.py")
    spec = importlib.util.spec_from_file_location("gse292993_dcc_qc", path)
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
        "geo_accession",
        "dcc_filename",
        "diagnosis_group",
        "compartment_guess",
        "include_qc",
        "dkk3_count",
        "negative_probe_n",
        "negative_geomean",
        "negative_geosd",
        "dkk3_geometric_loq",
        "dkk3_above_geometric_loq",
        "negative_arithmetic_mean",
        "negative_arithmetic_sd",
        "dkk3_arithmetic_loq",
        "dkk3_above_arithmetic_loq",
    ]
    ordered = [field for field in preferred if field in fieldnames]
    ordered.extend(field for field in fieldnames if field not in ordered)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=ordered, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def as_bool(value: object) -> bool:
    return value in (True, "True", "true", "1", 1)


def mean(values: list[float]) -> float | None:
    if not values:
        return None
    return sum(values) / len(values)


def sample_sd(values: list[float]) -> float | None:
    if not values:
        return None
    if len(values) == 1:
        return 0.0
    center = mean(values)
    assert center is not None
    return math.sqrt(sum((value - center) ** 2 for value in values) / (len(values) - 1))


def loq_metrics(
    negative_counts: list[int],
    *,
    pseudocount: float,
    sd_multiplier: float,
) -> dict:
    values = [float(value) for value in negative_counts]
    arithmetic_mean = mean(values)
    arithmetic_sd = sample_sd(values)
    arithmetic_loq = None
    if arithmetic_mean is not None and arithmetic_sd is not None:
        arithmetic_loq = arithmetic_mean + sd_multiplier * arithmetic_sd

    log_values = [math.log(value + pseudocount) for value in values]
    log_mean = mean(log_values)
    log_sd = sample_sd(log_values)
    geomean = None
    geosd = None
    geometric_loq = None
    if log_mean is not None and log_sd is not None:
        geomean = math.exp(log_mean) - pseudocount
        geosd = math.exp(log_sd)
        geometric_loq = math.exp(log_mean + sd_multiplier * log_sd) - pseudocount

    return {
        "negative_probe_n": len(values),
        "negative_arithmetic_mean": arithmetic_mean,
        "negative_arithmetic_sd": arithmetic_sd,
        "dkk3_arithmetic_loq": arithmetic_loq,
        "negative_geomean": geomean,
        "negative_geosd": geosd,
        "dkk3_geometric_loq": geometric_loq,
    }


def summarize_by_group(rows: list[dict], columns: list[str]) -> list[dict]:
    grouped: dict[tuple[str, ...], list[dict]] = defaultdict(list)
    for row in rows:
        key = tuple(str(row.get(column) or "unknown") for column in columns)
        grouped[key].append(row)

    output = []
    for key, group_rows in sorted(grouped.items()):
        included = [row for row in group_rows if as_bool(row.get("include_qc"))]
        above_geo = [
            row for row in included if as_bool(row.get("dkk3_above_geometric_loq"))
        ]
        above_arith = [
            row for row in included if as_bool(row.get("dkk3_above_arithmetic_loq"))
        ]
        item = dict(zip(columns, key))
        item["n_rois"] = len(group_rows)
        item["n_include_qc"] = len(included)
        item["n_dkk3_above_geometric_loq"] = len(above_geo)
        item["fraction_dkk3_above_geometric_loq"] = (
            len(above_geo) / len(included) if included else None
        )
        item["n_dkk3_above_arithmetic_loq"] = len(above_arith)
        item["fraction_dkk3_above_arithmetic_loq"] = (
            len(above_arith) / len(included) if included else None
        )
        output.append(item)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Compute per-ROI DKK3 background/LOQ flags from GSE292993 negative probes."
    )
    parser.add_argument("--negative-pseudocount", type=float, default=1.0)
    parser.add_argument("--loq-sd-multiplier", type=float, default=2.0)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    gene = config["project"]["gene"]
    output = ensure_results_dirs(config)
    table_dir = output["tables"]
    meta_dir = output["meta"]
    dcc_qc = load_dcc_qc_module()

    pkc_summary, pkc_rows = dcc_qc.parse_pkc_code_map(dcc_qc.find_pkc_file(config), gene)
    gene_codes = {
        row["code_id"]
        for row in pkc_rows
        if row.get("is_primary_gene") and row.get("code_id")
    }
    negative_codes = {
        row["code_id"]
        for row in pkc_rows
        if row.get("code_id")
        and (
            "negative" in str(row.get("code_class") or "").casefold()
            or "negative" in str(row.get("target") or "").casefold()
        )
    }
    roi_qc_rows = read_csv(table_dir / "gse292993_roi_qc_flags.csv")
    roi_by_gsm = {row["geo_accession"]: row for row in roi_qc_rows}

    rows: list[dict] = []
    for path in dcc_qc.dcc_paths(config):
        dcc_row, code_counts = dcc_qc.parse_dcc(path, gene_codes, negative_codes)
        gsm = dcc_row["geo_accession"]
        roi_row = dict(roi_by_gsm.get(gsm, {}))
        dkk3_count = sum(code_counts.get(code, 0) for code in gene_codes)
        negative_counts = [code_counts[code] for code in negative_codes if code in code_counts]
        metrics = loq_metrics(
            negative_counts,
            pseudocount=args.negative_pseudocount,
            sd_multiplier=args.loq_sd_multiplier,
        )
        geometric_loq = metrics["dkk3_geometric_loq"]
        arithmetic_loq = metrics["dkk3_arithmetic_loq"]
        row = {
            **roi_row,
            "geo_accession": gsm,
            "dcc_filename": dcc_row["dcc_filename"],
            "dkk3_count": dkk3_count,
            **metrics,
            "dkk3_above_geometric_loq": bool(
                geometric_loq is not None and dkk3_count > geometric_loq
            ),
            "dkk3_above_arithmetic_loq": bool(
                arithmetic_loq is not None and dkk3_count > arithmetic_loq
            ),
        }
        rows.append(row)

    include_rows = [row for row in rows if as_bool(row.get("include_qc"))]
    summary = {
        "dataset": DEFAULT_GSE,
        "primary_gene": gene,
        "primary_gene_codes_used": sorted(gene_codes),
        "negative_probe_code_count": len(negative_codes),
        "negative_pseudocount": args.negative_pseudocount,
        "loq_sd_multiplier": args.loq_sd_multiplier,
        "n_rois": len(rows),
        "n_include_qc": len(include_rows),
        "n_include_qc_dkk3_count_gt0": sum(
            1 for row in include_rows if float(row["dkk3_count"]) > 0
        ),
        "n_include_qc_dkk3_above_geometric_loq": sum(
            1 for row in include_rows if as_bool(row["dkk3_above_geometric_loq"])
        ),
        "n_include_qc_dkk3_above_arithmetic_loq": sum(
            1 for row in include_rows if as_bool(row["dkk3_above_arithmetic_loq"])
        ),
        "negative_probe_n_distribution": dict(
            sorted(Counter(int(row["negative_probe_n"]) for row in rows).items())
        ),
        "by_diagnosis_group_compartment": summarize_by_group(
            rows, ["diagnosis_group", "compartment_guess"]
        ),
        "by_diagnosis_guess_compartment": summarize_by_group(
            rows, ["diagnosis_guess", "compartment_guess"]
        ),
    }

    write_csv(table_dir / "gse292993_dkk3_loq_flags.csv", rows)
    (meta_dir / "gse292993_dkk3_loq_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(json.dumps(summary, indent=2))
    print(f"\nWrote {table_dir / 'gse292993_dkk3_loq_flags.csv'}")
    print(f"Wrote {meta_dir / 'gse292993_dkk3_loq_summary.json'}")

    failures = []
    if not gene_codes:
        failures.append(f"No PKC code IDs resolved for {gene}")
    if not negative_codes:
        failures.append("No negative probe code IDs resolved")
    if len(rows) != 794:
        failures.append(f"Expected 794 ROIs, found {len(rows)}")
    if args.strict and failures:
        print("Strict DKK3 LOQ failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
