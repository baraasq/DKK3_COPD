# Pipeline scripts

Scripts are grouped by execution stage:

1. `00_setup`
2. `01_download`
3. `02_metadata`
4. `03_geomx_qc`
5. `04_normalization`
6. `05_dkk3`
7. `06_celltype`
8. `07_communication`
9. `08_validation`
10. `09_figures`

Shared code belongs in `utils`. Existing top-level scripts are the earlier
Firsova/SCRINSHOT validation prototypes and will remain separate from the
GSE292993 primary workflow.

## Public-data download

From the repository root on the Linux server:

```bash
# Smaller GEO/Zenodo deposits plus the SRA run table; no FASTQ conversion
bash scripts/01_download/download_all_deposits.bash metadata

# Review the current number of runs and compressed SRA size
bash scripts/01_download/download_all_deposits.bash sra-info

# Download and convert every PRJNA1282758 run after checking free disk space
THREADS=16 CONFIRM_SRA_FASTQ=YES \
  bash scripts/01_download/download_all_deposits.bash sra
```

For quota-limited filesystems, put the `fasterq-dump` temporary directory on a
larger scratch volume and remove each `.sra` archive after its gzipped FASTQ
files have been written:

```bash
THREADS=16 CONFIRM_SRA_FASTQ=YES \
  SRA_TEMP_ROOT=/scratch/baraa/COPD_public_sra_tmp \
  KEEP_SRA_ARCHIVE=NO \
  bash scripts/01_download/download_all_deposits.bash sra
```

If `fasterq-dump` fails with `disk-limit exceeded`, remove the failed run's
temporary files and partial FASTQs before rerunning:

```bash
rm -rf data/raw/downloads/sra/PRJNA1282758/tmp/SRR34233583
rm -f data/raw/downloads/sra/PRJNA1282758/fastq/SRR34233583*.fastq \
      data/raw/downloads/sra/PRJNA1282758/fastq/SRR34233583*.fastq.gz
```

The downloader resumes partial HTTP files, skips completed SRA runs, retains
all read files emitted by `fasterq-dump --split-files`, and verifies the
checksums published by Zenodo.

## GSE292993 GeoMx input audit

After the GEO files are downloaded, audit the PKC panel and DCC files before
building downstream GeoMx objects:

```bash
python scripts/03_geomx_qc/00_audit_gse292993_inputs.py --copy-geo-inputs --strict
```

This writes:

- `results/meta/gse292993_geomx_input_audit.json`
- `results/tables/gse292993_dcc_input_manifest.csv`

The audit checks that the WTA PKC file is readable, `DKK3` is present in the
panel, DCC files are discoverable from `GSE292993_RAW.tar` or
`data/raw/gse292993/dcc`, and each DCC can be scanned for gene/control text.

Then profile the copied DCC files and GEO sample metadata:

```bash
python scripts/03_geomx_qc/01_profile_gse292993_objects.py --strict
```

This writes:

- `results/meta/gse292993_object_profile.json`
- `results/tables/gse292993_geo_sample_metadata.csv`
- `results/tables/gse292993_roi_metadata_initial.csv`

The profile checks that all DCC filenames can be matched back to GEO sample
records and captures a preview of the DCC structure before downstream GeoMx
object construction.
