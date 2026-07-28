#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import shutil
import sys
import tarfile
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Iterable

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from common import configured_path, ensure_results_dirs, load_config, project_path


DEFAULT_GSE = "GSE292993"


@dataclass(frozen=True)
class DccPayload:
    name: str
    size_bytes: int
    stream: BinaryIO


def accession_download_dir(accession: str) -> Path:
    return project_path(f"data/raw/downloads/geo/{accession}/supplementary")


def find_first_existing(candidates: Iterable[Path]) -> Path | None:
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def open_text(path: Path):
    if path.name.endswith(".gz"):
        return gzip.open(path, mode="rt", encoding="utf-8", errors="replace")
    return path.open("rt", encoding="utf-8", errors="replace")


def text_has_gene(line: str, gene: str) -> bool:
    return re.search(rf"(^|[^A-Za-z0-9_]){re.escape(gene)}([^A-Za-z0-9_]|$)", line, re.I) is not None


def inspect_pkc(path: Path | None, gene: str) -> dict:
    result = {
        "path": str(path) if path else None,
        "exists": bool(path and path.exists()),
        "gene": gene,
    }
    if not path or not path.exists():
        result["status"] = "missing"
        return result

    line_count = 0
    gene_line_count = 0
    first_gene_line = None
    first_nonempty_line = None
    with open_text(path) as handle:
        for line in handle:
            line_count += 1
            stripped = line.strip()
            if stripped and first_nonempty_line is None:
                first_nonempty_line = stripped[:250]
            if text_has_gene(line, gene):
                gene_line_count += 1
                if first_gene_line is None:
                    first_gene_line = stripped[:500]

    result.update(
        status="ok",
        size_bytes=path.stat().st_size,
        line_count=line_count,
        gene_present=gene_line_count > 0,
        gene_line_count=gene_line_count,
        first_gene_line=first_gene_line,
        first_nonempty_line=first_nonempty_line,
    )
    return result


def tar_dcc_members(path: Path) -> list[tarfile.TarInfo]:
    with tarfile.open(path) as archive:
        return [
            member
            for member in archive.getmembers()
            if member.isfile() and member.name.casefold().endswith(".dcc")
        ]


def iter_dcc_from_tar(path: Path) -> Iterable[DccPayload]:
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.isfile() or not member.name.casefold().endswith(".dcc"):
                continue
            stream = archive.extractfile(member)
            if stream is None:
                continue
            with stream:
                yield DccPayload(member.name, int(member.size), stream)


def iter_dcc_from_dir(path: Path) -> Iterable[DccPayload]:
    for dcc_path in sorted(path.rglob("*.dcc")):
        with dcc_path.open("rb") as stream:
            yield DccPayload(str(dcc_path.relative_to(path)), dcc_path.stat().st_size, stream)


def inspect_dcc_payload(payload: DccPayload, gene: str) -> dict:
    line_count = 0
    gene_line_count = 0
    first_gene_line = None
    first_nonempty_line = None
    contains_negative = False
    contains_no_template = False

    wrapper = io.TextIOWrapper(payload.stream, encoding="utf-8", errors="replace")
    for line in wrapper:
        line_count += 1
        stripped = line.strip()
        if stripped and first_nonempty_line is None:
            first_nonempty_line = stripped[:250]
        folded = stripped.casefold()
        contains_negative = contains_negative or "negative" in folded
        contains_no_template = contains_no_template or "no template" in folded
        if text_has_gene(line, gene):
            gene_line_count += 1
            if first_gene_line is None:
                first_gene_line = stripped[:500]

    filename = Path(payload.name).name
    return {
        "dcc_member": payload.name,
        "dcc_filename": filename,
        "roi_id_guess": Path(filename).stem,
        "size_bytes": payload.size_bytes,
        "line_count": line_count,
        "gene_present": gene_line_count > 0,
        "gene_line_count": gene_line_count,
        "first_gene_line": first_gene_line,
        "first_nonempty_line": first_nonempty_line,
        "contains_negative_probe_text": contains_negative,
        "contains_no_template_control_text": contains_no_template,
    }


def inspect_dcc_inputs(raw_tar: Path | None, dcc_dir: Path, gene: str) -> tuple[dict, list[dict]]:
    result = {
        "raw_tar": str(raw_tar) if raw_tar else None,
        "dcc_dir": str(dcc_dir),
        "gene": gene,
    }
    rows: list[dict] = []

    if raw_tar and raw_tar.exists():
        members = tar_dcc_members(raw_tar)
        result.update(
            source="tar",
            status="ok",
            raw_tar_size_bytes=raw_tar.stat().st_size,
            dcc_count=len(members),
            first_dcc_members=[member.name for member in members[:10]],
        )
        for payload in iter_dcc_from_tar(raw_tar):
            rows.append(inspect_dcc_payload(payload, gene))
    elif dcc_dir.exists():
        dcc_files = sorted(dcc_dir.rglob("*.dcc"))
        result.update(
            source="directory",
            status="ok" if dcc_files else "missing",
            dcc_count=len(dcc_files),
            first_dcc_members=[str(path.relative_to(dcc_dir)) for path in dcc_files[:10]],
        )
        for payload in iter_dcc_from_dir(dcc_dir):
            rows.append(inspect_dcc_payload(payload, gene))
    else:
        result.update(source=None, status="missing", dcc_count=0, first_dcc_members=[])

    result["dccs_with_gene"] = sum(1 for row in rows if row["gene_present"])
    result["dccs_with_negative_probe_text"] = sum(
        1 for row in rows if row["contains_negative_probe_text"]
    )
    result["dccs_with_no_template_control_text"] = sum(
        1 for row in rows if row["contains_no_template_control_text"]
    )
    return result, rows


def copy_geo_inputs(raw_tar: Path | None, pkc_path: Path | None, config: dict) -> dict:
    dcc_dir = configured_path(config, "geomx_dcc_dir")
    pkc_dir = configured_path(config, "geomx_pkc_dir")
    metadata_dir = configured_path(config, "geomx_metadata_dir")
    copied = {"dcc_files": 0, "pkc_files": 0, "metadata_files": 0}

    if raw_tar and raw_tar.exists():
        dcc_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(raw_tar) as archive:
            for member in archive.getmembers():
                if not member.isfile() or not member.name.casefold().endswith(".dcc"):
                    continue
                target = dcc_dir / Path(member.name).name
                stream = archive.extractfile(member)
                if stream is None:
                    continue
                with stream, target.open("wb") as output:
                    shutil.copyfileobj(stream, output)
                copied["dcc_files"] += 1

    if pkc_path and pkc_path.exists():
        pkc_dir.mkdir(parents=True, exist_ok=True)
        target = pkc_dir / pkc_path.name
        if pkc_path.resolve() != target.resolve():
            shutil.copy2(pkc_path, target)
            copied["pkc_files"] += 1

    download_root = accession_download_dir(DEFAULT_GSE).parent
    if download_root.exists():
        metadata_dir.mkdir(parents=True, exist_ok=True)
        for metadata_path in sorted(download_root.rglob(f"{DEFAULT_GSE}_family.*")):
            target = metadata_dir / metadata_path.name
            if metadata_path.resolve() != target.resolve():
                shutil.copy2(metadata_path, target)
                copied["metadata_files"] += 1

    return copied


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Audit GSE292993 GeoMx PKC and DCC inputs before downstream QC."
    )
    parser.add_argument(
        "--copy-geo-inputs",
        action="store_true",
        help="Copy DCC, PKC, and GEO metadata from data/raw/downloads into the structured raw folders.",
    )
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    config = load_config()
    gene = config["project"]["gene"]
    output = ensure_results_dirs(config)
    dcc_dir = configured_path(config, "geomx_dcc_dir")
    pkc_dir = configured_path(config, "geomx_pkc_dir")

    raw_tar = find_first_existing(
        [
            accession_download_dir(DEFAULT_GSE) / f"{DEFAULT_GSE}_RAW.tar",
            project_path(f"data/raw/{DEFAULT_GSE}_RAW.tar"),
            project_path(f"data/raw/gse292993/{DEFAULT_GSE}_RAW.tar"),
        ]
    )
    pkc_path = find_first_existing(
        [
            accession_download_dir(DEFAULT_GSE) / f"{DEFAULT_GSE}_Hs_R_NGS_WTA_v1.0.pkc.gz",
            *sorted(pkc_dir.glob("*.pkc.gz")),
            *sorted(pkc_dir.glob("*.pkc")),
        ]
    )

    copy_summary = None
    if args.copy_geo_inputs:
        copy_summary = copy_geo_inputs(raw_tar, pkc_path, config)

    pkc = inspect_pkc(pkc_path, gene)
    dcc, dcc_rows = inspect_dcc_inputs(raw_tar, dcc_dir, gene)

    audit = {
        "project": config["project"],
        "dataset": DEFAULT_GSE,
        "pkc": pkc,
        "dcc": dcc,
        "copy_summary": copy_summary,
        "interpretation": {
            "panel_supports_direct_dkk3": bool(pkc.get("gene_present")),
            "dcc_supports_direct_dkk3": dcc.get("dccs_with_gene", 0) > 0,
            "ready_for_geomx_qc": bool(pkc.get("gene_present"))
            and dcc.get("dcc_count", 0) > 0,
        },
    }

    audit_path = output["meta"] / "gse292993_geomx_input_audit.json"
    dcc_table = output["tables"] / "gse292993_dcc_input_manifest.csv"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")
    write_csv(dcc_table, dcc_rows)

    print(json.dumps(audit, indent=2))
    print(f"\nWrote {audit_path}")
    print(f"Wrote {dcc_table}")

    failures = []
    if pkc.get("status") != "ok":
        failures.append("PKC file missing or unreadable")
    elif not pkc.get("gene_present"):
        failures.append(f"{gene} was not found in the PKC panel file")
    if dcc.get("status") != "ok" or dcc.get("dcc_count", 0) == 0:
        failures.append("No DCC files found in GSE292993_RAW.tar or configured dcc_dir")

    if args.strict and failures:
        print("Strict audit failed: " + "; ".join(failures))
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
