#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import sys
import zipfile
from collections import Counter
from pathlib import Path

from common import (
    configured_path,
    ensure_results_dirs,
    load_config,
    normalize_condition,
    resolve_column,
    scrinshot_cell_map_members,
)


SCRINSHOT_METADATA_COLUMNS = {
    "donor",
    "sample",
    "disease",
    "smoking",
    "gender",
    "roi",
    "x",
    "y",
    "cellclass",
    "celltype",
    "cellsubtype",
}


def audit_scrna(path: Path, config: dict) -> dict:
    result: dict = {
        "path": str(path),
        "exists": path.exists(),
        "gene": config["project"]["gene"],
    }
    if not path.exists():
        result["status"] = "missing"
        return result

    try:
        import anndata as ad
    except ImportError as exc:
        result.update(
            status="dependency_missing",
            error="Install anndata to inspect the configured h5ad.",
        )
        return result

    adata = ad.read_h5ad(path, backed="r")
    try:
        gene = config["project"]["gene"]
        obs_columns = list(map(str, adata.obs.columns))
        sample_col = resolve_column(
            obs_columns,
            config["metadata"]["sample_candidates"],
            label="scRNA sample/donor column",
            required=False,
        )
        condition_col = resolve_column(
            obs_columns,
            config["metadata"]["condition_candidates"],
            label="scRNA condition column",
            required=False,
        )
        study_col = resolve_column(
            obs_columns,
            config["metadata"]["study_candidates"],
            label="scRNA study column",
            required=False,
        )
        cell_type_col = resolve_column(
            obs_columns,
            config["metadata"]["cell_type_candidates"],
            label="scRNA cell-type column",
            required=False,
        )

        result.update(
            status="ok",
            n_obs=int(adata.n_obs),
            n_vars=int(adata.n_vars),
            gene_present=bool(gene in adata.var_names),
            layers=sorted(map(str, adata.layers.keys())),
            resolved_columns={
                "sample": sample_col,
                "condition": condition_col,
                "study": study_col,
                "cell_type": cell_type_col,
            },
            counts_layer_present=config["scrna"]["counts_layer"] in adata.layers,
            expression_layer_present=(
                config["scrna"]["expression_layer"] in adata.layers
            ),
        )

        if condition_col:
            normalized = adata.obs[condition_col].map(
                lambda value: normalize_condition(value, config)
            )
            result["cells_by_condition"] = {
                str(key): int(value)
                for key, value in normalized.value_counts(dropna=False).items()
            }

        if sample_col and condition_col:
            sample_frame = adata.obs[[sample_col, condition_col]].copy()
            sample_frame["normalized_condition"] = sample_frame[condition_col].map(
                lambda value: normalize_condition(value, config)
            )
            pairs = sample_frame[
                [sample_col, "normalized_condition"]
            ].drop_duplicates()
            result["samples_by_condition"] = {
                str(key): int(value)
                for key, value in pairs["normalized_condition"]
                .value_counts(dropna=False)
                .items()
            }

        if sample_col and condition_col and study_col:
            study_frame = adata.obs[
                [sample_col, condition_col, study_col]
            ].drop_duplicates()
            study_frame["normalized_condition"] = study_frame[condition_col].map(
                lambda value: normalize_condition(value, config)
            )
            counts = (
                study_frame.groupby(
                    [study_col, "normalized_condition"], observed=True
                )[sample_col]
                .nunique()
                .reset_index(name="n_samples")
            )
            result["samples_by_study_condition"] = counts.to_dict("records")
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()

    return result


def cell_map_members(archive: zipfile.ZipFile) -> list[str]:
    return scrinshot_cell_map_members(archive.namelist())


def read_csv_header(archive: zipfile.ZipFile, member: str) -> list[str]:
    with archive.open(member) as raw:
        line = raw.readline().decode("utf-8-sig", errors="replace")
    return next(csv.reader([line]))


def read_csv_first_row(
    archive: zipfile.ZipFile, member: str
) -> tuple[list[str], list[str] | None]:
    with archive.open(member) as raw:
        header_line = raw.readline().decode("utf-8-sig", errors="replace")
        first_line = raw.readline().decode("utf-8-sig", errors="replace")
    header = next(csv.reader([header_line]))
    first = next(csv.reader([first_line])) if first_line else None
    return header, first


def audit_scrinshot(
    path: Path,
    config: dict,
) -> tuple[dict, list[dict], list[str]]:
    result: dict = {
        "path": str(path),
        "exists": path.exists(),
        "gene": config["project"]["gene"],
    }
    sections: list[dict] = []
    gene_panel: set[str] = set()

    if not path.exists():
        result["status"] = "missing"
        return result, sections, []

    if not zipfile.is_zipfile(path):
        result.update(status="invalid", error="Configured file is not a ZIP archive.")
        return result, sections, []

    with zipfile.ZipFile(path) as archive:
        members = cell_map_members(archive)
        result["archive_members"] = len(archive.infolist())
        result["cell_map_csvs"] = len(members)

        donor_counter: Counter[str] = Counter()
        condition_counter: Counter[str] = Counter()

        for member in members:
            header, first = read_csv_first_row(archive, member)
            gene_columns = [
                column
                for column in header
                if column.strip().casefold() not in SCRINSHOT_METADATA_COLUMNS
                and not column.strip().casefold().startswith("unnamed:")
            ]
            gene_panel.update(gene_columns)

            row = dict(zip(header, first or []))
            donor = row.get("Donor", "")
            condition = normalize_condition(row.get("Disease", ""), config)
            sample = row.get("Sample", Path(member).stem)
            donor_counter[donor] += 1
            condition_counter[condition] += 1
            sections.append(
                {
                    "archive_member": member,
                    "donor": donor,
                    "sample": sample,
                    "condition": condition,
                    "n_expression_columns": len(gene_columns),
                }
            )

        dot_genes = sorted(
            {
                Path(name).stem
                for name in archive.namelist()
                if "/Dot-coordinates/" in name
                and name.casefold().endswith(".csv")
            }
        )

    gene = config["project"]["gene"]
    result.update(
        status="ok",
        processed_expression_columns=len(gene_panel),
        dot_coordinate_genes=len(dot_genes),
        gene_present_in_processed_tables=gene in gene_panel,
        gene_present_in_dot_coordinates=gene in dot_genes,
        direct_spatial_gene_analysis_supported=(
            gene in gene_panel or gene in dot_genes
        ),
        sections_by_donor=dict(sorted(donor_counter.items())),
        sections_by_condition=dict(sorted(condition_counter.items())),
    )
    return result, sections, sorted(gene_panel | set(dot_genes))


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit COPD scRNA and Firsova SCRINSHOT inputs."
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit non-zero when required inputs or DKK3 scRNA coverage are missing.",
    )
    args = parser.parse_args()

    config = load_config()
    output = ensure_results_dirs(config)
    scrna_path = configured_path(
        config, "scrna_h5ad", environment_variable="COPD_SCRNA_H5AD"
    )
    scrinshot_path = configured_path(
        config, "scrinshot_zip", environment_variable="COPD_SCRINSHOT_ZIP"
    )

    scrna = audit_scrna(scrna_path, config)
    scrinshot, sections, genes = audit_scrinshot(scrinshot_path, config)
    audit = {
        "project": config["project"],
        "scrna": scrna,
        "scrinshot": scrinshot,
        "interpretation": {
            "direct_dkk3_source": "scRNA" if scrna.get("gene_present") else None,
            "spatial_dkk3_source": (
                "Firsova SCRINSHOT"
                if scrinshot.get("direct_spatial_gene_analysis_supported")
                else None
            ),
            "scrinshot_role": (
                "direct spatial DKK3"
                if scrinshot.get("direct_spatial_gene_analysis_supported")
                else "cellular-neighborhood context only"
            ),
        },
    }

    audit_path = output["meta"] / "input_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2, default=str), encoding="utf-8")

    if sections:
        write_csv(
            output["meta"] / "scrinshot_sections.csv",
            sections,
            [
                "archive_member",
                "donor",
                "sample",
                "condition",
                "n_expression_columns",
            ],
        )
    if genes:
        write_csv(
            output["meta"] / "scrinshot_gene_panel.csv",
            [{"gene": gene} for gene in genes],
            ["gene"],
        )

    print(json.dumps(audit, indent=2, default=str))
    print(f"\nWrote {audit_path}")

    failures: list[str] = []
    if not scrna.get("exists"):
        failures.append("scRNA input is missing")
    elif not scrna.get("gene_present"):
        failures.append(f"{config['project']['gene']} is absent from scRNA")
    elif not scrna.get("counts_layer_present"):
        failures.append("configured scRNA raw-count layer is absent")
    if not scrinshot.get("exists"):
        failures.append("SCRINSHOT archive is missing")

    if args.strict and failures:
        print("Strict audit failed: " + "; ".join(failures), file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
