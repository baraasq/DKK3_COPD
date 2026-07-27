#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import tarfile
import tempfile
from pathlib import Path, PurePosixPath
from typing import Iterator

from common import configured_path, ensure_results_dirs, load_config


def decode_feature_names(values) -> list[str]:
    return [
        value.decode("utf-8") if isinstance(value, bytes) else str(value)
        for value in values
    ]


def audit_h5(
    path: Path,
    sample: str,
    gene: str,
    *,
    matrix_label: str | None = None,
) -> dict:
    try:
        import h5py
        import numpy as np
    except ImportError as exc:
        raise RuntimeError(
            "This analysis requires h5py and numpy. Install environment.yml."
        ) from exc

    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        names = np.asarray(
            decode_feature_names(matrix["features"]["name"][:]),
            dtype=object,
        )
        shape = tuple(int(value) for value in matrix["shape"][:])
        hits = np.flatnonzero(names == gene)

        row = {
            "sample": sample,
            "matrix": matrix_label or str(path),
            "n_features": shape[0],
            "n_spots": shape[1],
            "gene": gene,
            "feature_hits": int(hits.size),
            "positive_spots": 0,
            "positive_spot_fraction": 0.0,
            "total_counts": 0,
            "max_spot_count": 0,
        }
        if hits.size != 1:
            return row

        gene_index = int(hits[0])
        indices = matrix["indices"][:]
        data = matrix["data"][:]
        values = data[indices == gene_index]
        row.update(
            positive_spots=int(values.size),
            positive_spot_fraction=(
                float(values.size / shape[1]) if shape[1] else 0.0
            ),
            total_counts=int(values.sum()) if values.size else 0,
            max_spot_count=int(values.max()) if values.size else 0,
        )
        return row


def directory_matrices(root: Path, matrix_name: str) -> Iterator[tuple[str, Path]]:
    for path in sorted(root.rglob(matrix_name)):
        yield path.parent.name, path


def archive_matrices(
    archive_path: Path,
    matrix_name: str,
    temporary_root: Path,
) -> Iterator[tuple[str, Path, str]]:
    with tarfile.open(archive_path, mode="r:*") as archive:
        members = sorted(
            (
                member
                for member in archive.getmembers()
                if member.isfile()
                and PurePosixPath(member.name).name == matrix_name
            ),
            key=lambda member: member.name,
        )
        for index, member in enumerate(members):
            sample = PurePosixPath(member.name).parent.name
            safe_name = f"{index:03d}_{sample}_{matrix_name}"
            destination = temporary_root / safe_name
            source = archive.extractfile(member)
            if source is None:
                continue
            with source, destination.open("wb") as target:
                shutil.copyfileobj(source, target)
            yield sample, destination, member.name


def audit_input(path: Path, matrix_name: str, gene: str) -> list[dict]:
    if path.is_dir():
        return [
            audit_h5(matrix_path, sample, gene)
            for sample, matrix_path in directory_matrices(path, matrix_name)
        ]

    if path.is_file() and (
        path.name.casefold().endswith(".tar.gz")
        or path.suffix.casefold() in {".tgz", ".tar"}
    ):
        with tempfile.TemporaryDirectory(prefix="dkk3_spatial_audit_") as directory:
            temporary_root = Path(directory)
            return [
                audit_h5(
                    matrix_path,
                    sample,
                    gene,
                    matrix_label=f"{path}!{member_name}",
                )
                for sample, matrix_path, member_name in archive_matrices(
                    path, matrix_name, temporary_root
                )
            ]

    if path.is_file() and path.name == matrix_name:
        return [audit_h5(path, path.parent.name, gene)]

    raise FileNotFoundError(
        f"Expected a Space Ranger directory, archive, or {matrix_name}: {path}"
    )


def write_csv(path: Path, rows: list[dict]) -> None:
    fieldnames = [
        "sample",
        "matrix",
        "n_features",
        "n_spots",
        "gene",
        "feature_hits",
        "positive_spots",
        "positive_spot_fraction",
        "total_counts",
        "max_spot_count",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    config = load_config()
    default_path = configured_path(
        config,
        "healthy_spatial_input",
        environment_variable="DISCOVAIR_LUNG_INPUT",
    )
    parser = argparse.ArgumentParser(
        description=(
            "Audit DKK3 in healthy DiscovAir/RRST Space Ranger matrices."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=default_path,
        help=(
            "A lung.tar.gz archive, extracted Space Ranger directory, or "
            "filtered_feature_bc_matrix.h5 file."
        ),
    )
    args = parser.parse_args()

    input_path = Path(
        os.path.expandvars(os.path.expanduser(str(args.input)))
    ).resolve()
    gene = config["project"]["gene"]
    matrix_name = config["healthy_spatial"]["matrix_name"]
    output = ensure_results_dirs(config)

    rows = audit_input(input_path, matrix_name, gene)
    table_path = output["tables"] / "healthy_spatial_dkk3_audit.csv"
    write_csv(table_path, rows)

    summary = {
        "input": str(input_path),
        "gene": gene,
        "matrices": len(rows),
        "matrices_with_gene": sum(row["feature_hits"] == 1 for row in rows),
        "matrices_with_detected_gene": sum(
            row["positive_spots"] > 0 for row in rows
        ),
        "total_positive_spots": sum(row["positive_spots"] for row in rows),
        "total_counts": sum(row["total_counts"] for row in rows),
        "table": str(table_path),
    }
    print(json.dumps(summary, indent=2))
    return 0 if rows and summary["matrices_with_gene"] == len(rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
