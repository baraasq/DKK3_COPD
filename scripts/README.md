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

Finally extract per-ROI DCC QC metrics and PKC code mappings:

```bash
python scripts/03_geomx_qc/02_extract_gse292993_dcc_qc.py --strict
```

This writes:

- `results/meta/gse292993_dcc_qc_summary.json`
- `results/tables/gse292993_pkc_code_map.csv`
- `results/tables/gse292993_dcc_qc_metrics.csv`

The DCC QC extraction parses `Scan_Attributes`, `NGS_Processing_Attributes`,
and `<Code_Summary>` from each `.dcc.gz` file. It reports raw, trimmed,
stitched, and aligned reads; trimmed, stitched, and aligned fractions; UMI/RTS
Q30; total counts; detected code count; `DKK3` counts; and negative-probe
summaries when the PKC labels negative-control RTS codes. `DKK3` counts are
reported only when the PKC file resolves one or more RTS code IDs for `DKK3`.

After extracting DCC QC metrics, merge them with GEO sample metadata and create
conservative ROI-level QC flags:

```bash
python scripts/03_geomx_qc/03_flag_gse292993_roi_qc.py --strict
```

This writes:

- `results/meta/gse292993_roi_qc_flag_summary.json`
- `results/tables/gse292993_roi_qc_flags.csv`

The default thresholds are intentionally conservative and mainly remove obvious
technical failures: at least 100,000 aligned reads, 10,000 detected code counts,
10,000 total code counts, trimmed fraction of at least 0.90, stitched fraction
of at least 0.80, aligned/stitched fraction of at least 0.80, and UMI/RTS Q30 of
at least 0.98.

After ROI QC, compute DKK3 detectability relative to negative-probe background:

```bash
python scripts/03_geomx_qc/04_compute_gse292993_dkk3_loq.py --strict
```

This writes:

- `results/meta/gse292993_dkk3_loq_summary.json`
- `results/tables/gse292993_dkk3_loq_flags.csv`

The LOQ step keeps two background calls side by side: a geometric negative-probe
LOQ using `exp(mean(log(negative + 1)) + 2*sd(log(negative + 1))) - 1`, and an
arithmetic `mean + 2*sd` call. Do not discard low-DKK3 ROIs automatically from
expression modeling; use these flags to interpret detectability and sensitivity.

## GSE292993 DKK3 summaries

After QC and LOQ calling, summarize DKK3 at ROI and donor-compartment levels:

```bash
python scripts/05_dkk3/00_summarize_gse292993_dkk3.py --strict
```

This writes:

- `results/meta/gse292993_dkk3_signal_summary.json`
- `results/tables/gse292993_dkk3_roi_signal.csv`
- `results/tables/gse292993_dkk3_primary_roi_by_group.csv`
- `results/tables/gse292993_dkk3_donor_compartment_summary.csv`
- `results/tables/gse292993_dkk3_donor_diagnosis_compartment_summary.csv`
- `results/tables/gse292993_dkk3_donor_overall_summary.csv`
- `results/tables/gse292993_dkk3_donor_diagnosis_overall_summary.csv`

The primary summary keeps QC-passing ROIs with known COPD/control labels,
airway/parenchyma/vessel compartments, and known donor IDs. ROI summaries are
descriptive; donor-compartment summaries are the unit to carry forward into
inference.

Then run donor-aware COPD vs control comparisons within each compartment:

```bash
python scripts/05_dkk3/01_test_gse292993_dkk3_donor_effects.py --strict
```

This writes:

- `results/meta/gse292993_dkk3_donor_effect_summary.json`
- `results/tables/gse292993_dkk3_donor_effect_tests.csv`

The primary metric is donor-level median `log1p(DKK3 CPM)` within each
compartment. The script also tests donor-level median CPM, median raw count, and
above-LOQ fractions using permutation p-values and bootstrap confidence
intervals. These are still unadjusted exploratory tests; treat them as the
first donor-aware screen before fuller mixed or multivariable models.

If COPD/control grouping disagrees with prior knowledge or appears sensitive to
smoker controls, run the stratified donor-aware comparisons:

```bash
python scripts/05_dkk3/02_test_gse292993_dkk3_smoking_strata.py --strict
```

This writes:

- `results/meta/gse292993_dkk3_smoking_strata_effect_summary.json`
- `results/tables/gse292993_dkk3_smoking_strata_effect_tests.csv`

The default pairwise contrasts are `COPD` vs `Non Smoker`, `COPD` vs `Smoker`,
and `Smoker` vs `Non Smoker` within each compartment, using the donor-level
tables that preserve the raw GEO `characteristics_condition` labels.

## DKK3 Figures

Create the donor-level parenchymal DKK3 summary figure after the stratified
tests are available:

```bash
python scripts/09_figures/00_plot_gse292993_dkk3_parenchyma.py
```

This writes PNG, SVG, and PDF versions of:

- `results/figures/gse292993_dkk3/gse292993_dkk3_parenchyma_donor_signal`

The figure shows donor-level `Non Smoker`, `Smoker`, and `COPD` distributions
for median `log1p(DKK3 CPM)`, above-LOQ fraction, and median raw DKK3 count.
