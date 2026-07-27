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

The downloader resumes partial HTTP files, skips completed SRA runs, retains
all read files emitted by `fasterq-dump --split-files`, and verifies the
checksums published by Zenodo.
