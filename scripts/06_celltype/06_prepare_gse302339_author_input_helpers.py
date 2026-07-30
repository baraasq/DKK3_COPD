#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import ensure_results_dirs, load_config, project_path


DEFAULT_INPUT_DIR = "input/data_cellranger8"
DEFAULT_RIBOSOMAL_OUTPUT = "input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv"
DEFAULT_PREFIXES = ("RPL", "RPS", "MRPL", "MRPS")


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


def decode_feature_names(values: Iterable[object]) -> list[str]:
    names: list[str] = []
    for value in values:
        if isinstance(value, bytes):
            names.append(value.decode())
        else:
            names.append(str(value))
    return names


def read_10x_h5_feature_names(path: Path) -> list[str]:
    try:
        import h5py
    except ModuleNotFoundError as exc:
        raise ModuleNotFoundError(
            "h5py is required to inspect Cell Ranger h5 feature names. "
            "Install h5py in the active environment or run inside spatial_omics."
        ) from exc

    with h5py.File(path, "r") as handle:
        candidates = [
            "matrix/features/name",
            "matrix/features/gene_names",
            "matrix/gene_names",
        ]
        for key in candidates:
            if key in handle:
                return decode_feature_names(handle[key][:])
        available = []
        handle.visit(available.append)
        raise KeyError(
            f"Could not find a Cell Ranger feature-name dataset in {path}. "
            f"Tried {candidates}. First available h5 keys: {available[:20]}"
        )


def ribosomal_gene_pattern(prefixes: tuple[str, ...]) -> re.Pattern:
    escaped = "|".join(re.escape(prefix.upper()) for prefix in prefixes)
    return re.compile(rf"^({escaped})[A-Z0-9-]*$")


def collect_ribosomal_genes(
    h5_paths: list[Path],
    prefixes: tuple[str, ...],
) -> tuple[set[str], list[dict]]:
    pattern = ribosomal_gene_pattern(prefixes)
    all_genes: set[str] = set()
    rows: list[dict] = []
    for h5_path in h5_paths:
        names = read_10x_h5_feature_names(h5_path)
        ribosomal = sorted({name for name in names if pattern.match(name.upper())})
        all_genes.update(ribosomal)
        rows.append(
            {
                "path": str(h5_path),
                "filename": h5_path.name,
                "n_features": len(names),
                "n_unique_features": len(set(names)),
                "n_duplicate_feature_names": len(names) - len(set(names)),
                "n_ribosomal_prefix_genes": len(ribosomal),
                "first_ribosomal_genes": ";".join(ribosomal[:20]),
            }
        )
    return all_genes, rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare small input helper files expected by the deposited GSE302339 "
            "Scanpy notebooks."
        )
    )
    parser.add_argument(
        "--input-dir",
        default=DEFAULT_INPUT_DIR,
        help="Directory containing extracted Cell Ranger *_filtered_feature_bc_matrix.h5 files.",
    )
    parser.add_argument(
        "--ribosomal-output",
        default=DEFAULT_RIBOSOMAL_OUTPUT,
        help=(
            "Path to write the one-gene-per-line ribosomal gene list expected by "
            "the author preprocessing notebook."
        ),
    )
    parser.add_argument(
        "--prefix",
        action="append",
        dest="prefixes",
        help=(
            "Ribosomal gene-symbol prefix to include. Can be repeated. "
            "Defaults to RPL, RPS, MRPL, and MRPS."
        ),
    )
    parser.add_argument("--strict", action="store_true", help="Exit non-zero if helpers are not ready.")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    config = load_config()
    results_dirs = ensure_results_dirs(config)

    input_dir = project_path(args.input_dir)
    output_path = project_path(args.ribosomal_output)
    prefixes = tuple(args.prefixes or DEFAULT_PREFIXES)

    h5_paths = sorted(input_dir.glob("*filtered_feature_bc_matrix.h5"))
    failures: list[str] = []
    if not input_dir.exists():
        failures.append(f"Input directory missing: {input_dir}")
    if not h5_paths:
        failures.append(f"No *_filtered_feature_bc_matrix.h5 files found in {input_dir}")

    ribosomal_genes: set[str] = set()
    manifest_rows: list[dict] = []
    if h5_paths:
        ribosomal_genes, manifest_rows = collect_ribosomal_genes(h5_paths, prefixes)
        if not ribosomal_genes:
            failures.append(
                f"No ribosomal genes found with prefixes {list(prefixes)} in {input_dir}"
            )

    if ribosomal_genes:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text("\n".join(sorted(ribosomal_genes)) + "\n", encoding="utf-8")

    manifest_path = results_dirs["tables"] / "gse302339_ribosomal_gene_helper_manifest.csv"
    summary_path = results_dirs["meta"] / "gse302339_author_input_helper_summary.json"
    write_csv(
        manifest_path,
        manifest_rows,
        preferred=[
            "filename",
            "path",
            "n_features",
            "n_unique_features",
            "n_duplicate_feature_names",
            "n_ribosomal_prefix_genes",
            "first_ribosomal_genes",
        ],
    )

    summary = {
        "input_dir": str(input_dir),
        "input_dir_exists": input_dir.exists(),
        "n_cellranger_h5": len(h5_paths),
        "ribosomal_output": str(output_path),
        "ribosomal_output_exists": output_path.exists(),
        "prefixes": list(prefixes),
        "n_ribosomal_genes": len(ribosomal_genes),
        "first_ribosomal_genes": sorted(ribosomal_genes)[:30],
        "manifest_path": str(manifest_path),
        "ready_for_author_preprocessing": output_path.exists() and bool(ribosomal_genes),
        "failures": failures,
    }
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(json.dumps(summary, indent=2))
    print()
    print(summary_path)
    print(manifest_path)

    if args.strict and failures:
        raise SystemExit(
            "Strict helper preparation failed: " + "; ".join(failures)
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
