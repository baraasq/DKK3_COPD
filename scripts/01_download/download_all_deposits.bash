#!/usr/bin/env bash
set -Eeuo pipefail

# Download the public deposits used by the COPD project:
#   PRJNA1282758  raw SRA reads (FASTQ)
#   GSE302339     Cell Ranger processed count matrices
#   GSE237120     GeoMx data
#   GSE292993     GeoMx data
#   Zenodo 16341197 analysis code
#
# The script is restartable. curl resumes partial files, and completed SRA runs
# receive marker files. Raw FASTQ conversion is intentionally gated because it
# can require substantially more disk than the compressed SRA archives.

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="${COPD_PUBLIC_ROOT:-$(cd -- "${SCRIPT_DIR}/../.." && pwd)}"
THREADS="${THREADS:-8}"
SRA_PROJECT="PRJNA1282758"
ZENODO_RECORD="16341197"
GEO_ACCESSIONS=(GSE302339 GSE237120 GSE292993)

DOWNLOAD_ROOT="${PROJECT_ROOT}/data/raw/downloads"
GEO_ROOT="${DOWNLOAD_ROOT}/geo"
SRA_ROOT="${DOWNLOAD_ROOT}/sra/${SRA_PROJECT}"
ZENODO_ROOT="${DOWNLOAD_ROOT}/zenodo/${ZENODO_RECORD}"

usage() {
    cat <<'EOF'
Usage:
  bash scripts/01_download/download_all_deposits.bash TARGET

TARGET:
  metadata  Download GEO metadata/supplementary files, Zenodo code, and the
            PRJNA1282758 SRA run table (no FASTQ conversion).
  geo       Download every supplementary file and family metadata for
            GSE302339, GSE237120, and GSE292993.
  zenodo    Download every file from Zenodo record 16341197 and verify checksums.
  sra-info  Download the PRJNA1282758 SRA run table and print its size summary.
  sra       Download all PRJNA1282758 SRA archives and convert every run to
            gzipped FASTQ.
  all       Run geo, zenodo, and sra.

Large-download safety:
  The "sra" and "all" targets require CONFIRM_SRA_FASTQ=YES.

Examples:
  bash scripts/01_download/download_all_deposits.bash metadata
  THREADS=16 CONFIRM_SRA_FASTQ=YES \
    bash scripts/01_download/download_all_deposits.bash sra
EOF
}

log() {
    printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

download_file() {
    local url="$1"
    local destination="$2"
    mkdir -p -- "$(dirname -- "${destination}")"
    log "Downloading ${url}"
    curl --fail --location --retry 5 --retry-delay 5 \
        --continue-at - --output "${destination}" "${url}"
}

geo_bucket() {
    local accession="$1"
    printf '%snnn\n' "${accession:0:${#accession}-3}"
}

download_geo_supplementary() {
    local accession="$1"
    local bucket
    local base_url
    local destination
    local listing

    bucket="$(geo_bucket "${accession}")"
    base_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${bucket}/${accession}"
    destination="${GEO_ROOT}/${accession}"
    listing="${destination}/supplementary_files.tsv"
    mkdir -p -- "${destination}/supplementary" "${destination}/metadata"

    log "Discovering supplementary files for ${accession}"
    python3 - "${base_url}/suppl/" "${listing}" <<'PY'
import html.parser
import pathlib
import sys
import urllib.parse
import urllib.request

base_url, output_path = sys.argv[1:]

class Links(html.parser.HTMLParser):
    def __init__(self):
        super().__init__()
        self.hrefs = []

    def handle_starttag(self, tag, attrs):
        if tag.lower() != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.hrefs.append(href)

request = urllib.request.Request(
    base_url, headers={"User-Agent": "COPD_public reproducible downloader"}
)
with urllib.request.urlopen(request) as response:
    page = response.read().decode("utf-8", errors="replace")

parser = Links()
parser.feed(page)
rows = []
for href in parser.hrefs:
    if href.startswith(("?", "/", "#")) or href in ("../", "./") or href.endswith("/"):
        continue
    url = urllib.parse.urljoin(base_url, href)
    name = pathlib.PurePosixPath(urllib.parse.unquote(urllib.parse.urlparse(url).path)).name
    if name:
        rows.append((url, name))

rows = sorted(set(rows))
if not rows:
    raise SystemExit(f"No supplementary files discovered at {base_url}")

with open(output_path, "w", encoding="utf-8", newline="\n") as handle:
    handle.write("url\tfilename\n")
    for url, name in rows:
        handle.write(f"{url}\t{name}\n")
PY

    while IFS=$'\t' read -r url filename; do
        [[ "${url}" == "url" ]] && continue
        download_file "${url}" "${destination}/supplementary/${filename}"
    done < "${listing}"

    download_file \
        "${base_url}/soft/${accession}_family.soft.gz" \
        "${destination}/metadata/${accession}_family.soft.gz"
    download_file \
        "${base_url}/miniml/${accession}_family.xml.tgz" \
        "${destination}/metadata/${accession}_family.xml.tgz"
}

download_all_geo() {
    require_command curl
    require_command python3
    local accession
    for accession in "${GEO_ACCESSIONS[@]}"; do
        download_geo_supplementary "${accession}"
    done
}

download_zenodo() {
    require_command curl
    require_command python3
    local record_json="${ZENODO_ROOT}/record.json"
    local manifest="${ZENODO_ROOT}/files.tsv"

    mkdir -p -- "${ZENODO_ROOT}/files"
    download_file \
        "https://zenodo.org/api/records/${ZENODO_RECORD}" \
        "${record_json}"

    python3 - "${record_json}" "${manifest}" <<'PY'
import json
import sys

record_path, manifest_path = sys.argv[1:]
with open(record_path, encoding="utf-8") as handle:
    record = json.load(handle)

files = record.get("files", [])
if not files:
    raise SystemExit("Zenodo record contains no downloadable files")

with open(manifest_path, "w", encoding="utf-8", newline="\n") as handle:
    handle.write("url\tfilename\tsize_bytes\tchecksum\n")
    for item in files:
        links = item.get("links", {})
        url = links.get("content") or links.get("download") or links.get("self")
        if not url:
            raise SystemExit(f"No download link for {item.get('key')}")
        filename = item["key"].replace("\\", "/").split("/")[-1]
        handle.write(
            f"{url}\t{filename}\t{item.get('size', '')}\t{item.get('checksum', '')}\n"
        )
PY

    while IFS=$'\t' read -r url filename size checksum; do
        [[ "${url}" == "url" ]] && continue
        download_file "${url}" "${ZENODO_ROOT}/files/${filename}"
    done < "${manifest}"

    python3 - "${manifest}" "${ZENODO_ROOT}/files" <<'PY'
import hashlib
import pathlib
import sys

manifest_path, files_dir = sys.argv[1:]
failures = []
with open(manifest_path, encoding="utf-8") as handle:
    next(handle)
    for line in handle:
        _, filename, expected_size, checksum = line.rstrip("\n").split("\t")
        path = pathlib.Path(files_dir, filename)
        if not path.is_file():
            failures.append(f"{filename}: missing")
            continue
        if expected_size and path.stat().st_size != int(expected_size):
            failures.append(
                f"{filename}: size {path.stat().st_size}, expected {expected_size}"
            )
        if checksum:
            algorithm, expected = checksum.split(":", 1)
            digest = hashlib.new(algorithm)
            with path.open("rb") as stream:
                for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
                    digest.update(block)
            if digest.hexdigest().lower() != expected.lower():
                failures.append(f"{filename}: {algorithm} checksum mismatch")

if failures:
    raise SystemExit("Zenodo verification failed:\n" + "\n".join(failures))
print("Zenodo files verified successfully.")
PY
}

download_sra_runinfo() {
    require_command python3
    local runinfo="${SRA_ROOT}/metadata/SraRunInfo.csv"
    local runlist="${SRA_ROOT}/metadata/run_accessions.txt"
    local manifest="${SRA_ROOT}/metadata/run_manifest.tsv"
    mkdir -p -- "${SRA_ROOT}/metadata"

    log "Querying NCBI SRA for ${SRA_PROJECT}"
    python3 - "${SRA_PROJECT}" "${runinfo}" <<'PY'
import shutil
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

project, output_path = sys.argv[1:]
base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
agent = "COPD_public reproducible downloader"

search_params = urllib.parse.urlencode(
    {
        "db": "sra",
        "term": f"{project}[BioProject]",
        "usehistory": "y",
        "retmax": "0",
    }
)
request = urllib.request.Request(
    f"{base}esearch.fcgi?{search_params}", headers={"User-Agent": agent}
)
with urllib.request.urlopen(request) as response:
    root = ET.fromstring(response.read())

count = int(root.findtext("Count", "0"))
query_key = root.findtext("QueryKey")
webenv = root.findtext("WebEnv")
if count == 0 or not query_key or not webenv:
    raise SystemExit(f"No public SRA runs found for {project}")

fetch_params = urllib.parse.urlencode(
    {
        "db": "sra",
        "query_key": query_key,
        "WebEnv": webenv,
        "retmax": str(count),
        "rettype": "runinfo",
        "retmode": "text",
    }
)
request = urllib.request.Request(
    f"{base}efetch.fcgi?{fetch_params}", headers={"User-Agent": agent}
)
with urllib.request.urlopen(request) as response, open(output_path, "wb") as output:
    shutil.copyfileobj(response, output)
print(f"NCBI returned {count} SRA records.")
PY

    python3 - "${runinfo}" "${runlist}" "${manifest}" <<'PY'
import csv
import sys

runinfo_path, runlist_path, manifest_path = sys.argv[1:]
with open(runinfo_path, encoding="utf-8-sig", newline="") as handle:
    rows = list(csv.DictReader(handle))

rows = [row for row in rows if row.get("Run", "").startswith(("SRR", "ERR", "DRR"))]
if not rows:
    raise SystemExit("SraRunInfo.csv did not contain any run accessions")

with open(runlist_path, "w", encoding="utf-8", newline="\n") as handle:
    for row in rows:
        handle.write(row["Run"] + "\n")

columns = [
    "Run",
    "BioSample",
    "SampleName",
    "LibraryName",
    "LibraryLayout",
    "spots",
    "bases",
    "size_MB",
]
with open(manifest_path, "w", encoding="utf-8", newline="\n") as handle:
    handle.write("\t".join(columns) + "\n")
    for row in rows:
        handle.write("\t".join(row.get(column, "") for column in columns) + "\n")

total_mb = sum(float(row.get("size_MB") or 0) for row in rows)
total_bases = sum(int(float(row.get("bases") or 0)) for row in rows)
print(f"Runs: {len(rows)}")
print(f"Compressed SRA size reported by NCBI: {total_mb / 1024:.1f} GiB")
print(f"Total sequenced bases: {total_bases:,}")
print("Allow substantially more free space for SRA archives, temporary files, and FASTQ.")
PY
}

download_sra_fastq() {
    [[ "${CONFIRM_SRA_FASTQ:-NO}" == "YES" ]] || die \
        "Set CONFIRM_SRA_FASTQ=YES after reviewing the sra-info size summary."

    require_command prefetch
    require_command fasterq-dump
    require_command vdb-validate
    require_command python3

    download_sra_runinfo

    local runlist="${SRA_ROOT}/metadata/run_accessions.txt"
    local archive_root="${SRA_ROOT}/archive"
    local fastq_root="${SRA_ROOT}/fastq"
    local temp_root="${SRA_ROOT}/tmp"
    local done_root="${SRA_ROOT}/done"
    local run
    local archive
    local fastq
    local -a fastqs

    mkdir -p -- "${archive_root}" "${fastq_root}" "${temp_root}" "${done_root}"

    while IFS= read -r run; do
        [[ -n "${run}" ]] || continue
        if [[ -f "${done_root}/${run}.complete" ]]; then
            log "Skipping completed run ${run}"
            continue
        fi

        log "Downloading SRA archive ${run}"
        prefetch "${run}" --max-size u --output-directory "${archive_root}"
        archive="${archive_root}/${run}/${run}.sra"
        [[ -f "${archive}" ]] || die "Expected archive was not created: ${archive}"
        vdb-validate "${archive}"

        log "Converting ${run} to FASTQ with ${THREADS} threads"
        mkdir -p -- "${temp_root}/${run}"
        fasterq-dump "${archive}" \
            --split-files \
            --force \
            --threads "${THREADS}" \
            --temp "${temp_root}/${run}" \
            --outdir "${fastq_root}"

        shopt -s nullglob
        fastqs=("${fastq_root}/${run}"*.fastq)
        shopt -u nullglob
        (( ${#fastqs[@]} > 0 )) || die "No FASTQ files produced for ${run}"

        if command -v pigz >/dev/null 2>&1; then
            pigz --force --processes "${THREADS}" -- "${fastqs[@]}"
        else
            gzip --force -- "${fastqs[@]}"
        fi
        printf 'completed_utc\t%s\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" \
            > "${done_root}/${run}.complete"
        log "Completed ${run}"
    done < "${runlist}"
}

main() {
    local target="${1:-}"

    case "${target}" in
        metadata)
            download_all_geo
            download_zenodo
            download_sra_runinfo
            ;;
        geo)
            download_all_geo
            ;;
        zenodo)
            download_zenodo
            ;;
        sra-info)
            download_sra_runinfo
            ;;
        sra)
            download_sra_fastq
            ;;
        all)
            download_all_geo
            download_zenodo
            download_sra_fastq
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            usage >&2
            exit 2
            ;;
    esac
}

main "$@"
