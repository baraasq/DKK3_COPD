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

By default this exports the annotation, meta-merge, and final merged cell-type
notebooks when available (`2_celltype_annotation`, `8_meta_merge`, and
`10_pertpy_celltype_merged`). This writes:

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

To plot GeoMx-native ROI QC metrics and sample balance:

```bash
python scripts/03_geomx_qc/05_plot_gse292993_roi_qc.py --strict
```

This writes PNG, SVG, and PDF versions of:

- `results/figures/gse292993_qc/gse292993_roi_qc_metric_distributions`
- `results/figures/gse292993_qc/gse292993_roi_qc_sample_balance`
- `results/meta/gse292993_roi_qc_plot_summary.json`

The closest GeoMx equivalents to standard spatial/scRNA QC are
`total_code_counts` for total recovered signal, `n_code_counts` for detected
targets, `umi_q30` for UMI quality, and negative-probe summaries for background.
GeoMx DCC files do not directly provide a `% mitochondrial` metric, and nuclei
count / ROI area require image or ROI annotation metadata if available. Rows
labeled `unknown diagnosis` or `unknown compartment` are retained here only for
QC auditing; they indicate unresolved GEO metadata labels, not a biological
group to interpret.

To audit whether the GeoMx ROI metadata has direct emphysema labels or only
proxy fields such as LAA950/GOLD/FEV1:

```bash
python scripts/03_geomx_qc/06_audit_gse292993_phenotype_labels.py \
  --include-qc-only \
  --strict
```

This writes:

- `results/meta/gse292993_phenotype_label_audit_summary.json`
- `results/tables/gse292993_phenotype_label_column_summary.csv`
- `results/tables/gse292993_phenotype_label_group_counts.csv`
- `results/tables/gse292993_laa_threshold_sensitivity_counts.csv`

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

## GSE292993 WTA matrix and scRNA reference

Before cell-type composition analysis, build the full GeoMx WTA target matrix
from the DCC files:

```bash
python scripts/04_normalization/00_build_gse292993_geomx_wta_matrix.py --strict
```

This writes:

- `data/processed/gse292993/gse292993_geomx_counts_by_roi.tsv.gz`
- `data/processed/gse292993/gse292993_geomx_logcpm_by_roi.tsv.gz`
- `results/tables/gse292993_geomx_feature_manifest.csv`
- `results/tables/gse292993_geomx_matrix_roi_manifest.csv`
- `results/meta/gse292993_geomx_matrix_summary.json`

By default, negative/control PKC features are excluded from the expression
matrix. Counts from multiple RTS codes mapping to the same target are summed.
PKC records mapping implausibly many RTS codes to one target are dropped by
default because they usually represent broad panel metadata rather than a
single gene-like WTA target. Dropped targets are written to
`results/tables/gse292993_geomx_dropped_feature_manifest.csv`.

Then audit the scRNA-seq reference and, when an annotated `.h5ad` is available,
build logCPM cell-type signatures:

```bash
python scripts/06_celltype/00_audit_scrna_reference.py --strict
```

This writes:

- `results/meta/gse302339_scrna_reference_audit_summary.json`
- `results/tables/gse302339_scrna_reference_file_manifest.csv`
- `results/tables/gse302339_scrna_celltype_counts.csv`
- `results/tables/gse302339_scrna_signature_manifest.csv`
- `data/processed/gse292993/gse302339_scrna_reference_signatures_logcpm.csv`

If the downloaded scRNA deposit only contains Cell Ranger count matrices, the
audit will report that it is counts-only and not ready for deconvolution. In
that case, provide an annotated `.h5ad` with a cell-type column via
`--reference-h5ad` and, if needed, `--cell-type-column`.

If the Zenodo Scanpy workflow is present but no final annotated `.h5ad` was
deposited, inspect the authors' notebooks for read/write paths and annotation
clues:

```bash
python scripts/06_celltype/02_inspect_gse302339_scanpy_notebooks.py --strict
```

This writes:

- `results/meta/gse302339_scanpy_notebook_inspection_summary.json`
- `results/tables/gse302339_scanpy_notebook_summary.csv`
- `results/tables/gse302339_scanpy_notebook_code_matches.csv`
- `results/tables/gse302339_scanpy_notebook_artifact_paths.csv`

Use this after the scRNA reference audit fails with a missing `.h5ad`. It scans
only notebook code cells, so embedded UMAP/dotplot image payloads are ignored.
Notebooks mentioning `write_h5ad` or `.write` are the best candidates for
reconstructing an annotated reference; notebooks with only `read_10x_h5`,
`leiden`, and marker-gene terms are likely annotation workflow notebooks rather
than deposited reference objects. The script also reports no-extension
`output/...` paths opened in `wb`/`rb` mode, because the authors' notebooks may
pickle intermediate annotated AnnData objects instead of writing `.h5ad` files.

When the notebook inspection shows annotation dictionaries or no-extension
`output/...` artifacts, export the full code cells from the key notebooks:

```bash
python scripts/06_celltype/03_extract_gse302339_annotation_code.py --strict
```

This writes:

- `results/meta/gse302339_annotation_code_export_summary.json`
- `results/tables/gse302339_annotation_code_cells.csv`
- `intermediate/gse302339_scanpy_workflow_code/2_celltype_annotation.py`
- `intermediate/gse302339_scanpy_workflow_code/8_meta_merge.py`
- `intermediate/gse302339_scanpy_workflow_code/10_pertpy_celltype_merged.py`

Those exported `.py` files are audit artifacts, not cleaned runnable scripts.
Use them to recover the authors' cluster-to-cell-type dictionaries and the
pickle-style intermediate object names before deciding whether to reconstruct
their annotated scRNA reference or switch to an external annotated lung atlas.
Notebook-only shell/magic lines such as `!rm ...` or `%matplotlib` are replaced
with no-op `pass` lines during export so the files can be syntax-checked with
`python -m py_compile`. The unused optional `scvi` import/seed and `mudata`
import lines in the preprocessing notebook are also skipped during export
because they can conflict with the active environment and are not used elsewhere
in the exported preprocessing code.

To turn the exported author `celldict_*` assignments into tidy cluster-to-cell
type mapping tables:

```bash
python scripts/06_celltype/04_extract_gse302339_author_celltype_maps.py --strict
```

This writes:

- `results/meta/gse302339_author_celltype_map_summary.json`
- `results/tables/gse302339_author_celltype_cluster_maps.csv`
- `results/tables/gse302339_author_celltype_dictionary_blocks.csv`

The parser records the downstream `.obs` column assigned after each dictionary.
With no `--code-file` arguments, it scans all exported notebook `.py` files in
`intermediate/gse302339_scanpy_workflow_code`. For the pneumocyte/fibroblast
question, use a map only when the saved-object audit confirms that the map and
object use the same clustering.

To trace the exact object-local clustering steps, including subset assignments,
Scanpy/Harmony/Leiden calls, `.obs` assignments, cell-type dictionaries, and
pickle read/write artifacts:

```bash
python scripts/06_celltype/13_trace_gse302339_scanpy_clustering_workflow.py --strict
```

This writes:

- `results/meta/gse302339_scanpy_clustering_workflow_summary.json`
- `results/tables/gse302339_scanpy_clustering_workflow_steps.csv`

Use this table to confirm which object variable each Leiden clustering belongs
to before applying an extracted cluster-to-cell-type dictionary.

Before rerunning the authors' preprocessing notebook export, prepare the small
`input/` helper files it expects:

```bash
python scripts/06_celltype/06_prepare_gse302339_author_input_helpers.py --strict
```

This legacy helper reconstructs `input/GOCC_RIBOSOMAL_SUBUNIT.v2023.1.Hs.csv`
from broad ribosomal prefixes. That approximation is adequate for dependency
discovery, but it is **not** suitable for an exact author-annotation replay.
Script 15 instead extracts the authors' exact 185-gene tuple from the deposited
code. This helper also extracts
`input/meta_cr8.csv` from the deposited Zenodo `scanpy_workflow.zip` if that
sidecar is present. The script scans the notebooks for literal `input/...`
file reads first, so unresolved sidecars are reported before rerunning a long
preprocessing step. It writes:

- `results/meta/gse302339_author_input_helper_summary.json`
- `results/tables/gse302339_ribosomal_gene_helper_manifest.csv`
- `results/tables/gse302339_author_notebook_sidecar_manifest.csv`

If the sidecar audit reports that `input/meta_cr8.csv` is missing from the
Zenodo notebook archive, reconstruct it from the Cell Ranger filenames and GEO
sample metadata:

```bash
python scripts/06_celltype/07_reconstruct_gse302339_meta_cr8.py \
  --download-soft \
  --strict
```

This downloads/uses the GSE302339 family SOFT file, merges GEO sample metadata
onto the 65 extracted Cell Ranger `.h5` files, and checks whether the exported
author code references any `meta[...]` columns missing from the reconstructed
table. It writes:

- `input/meta_cr8.csv`
- `results/meta/gse302339_meta_cr8_reconstruction_summary.json`
- `results/tables/gse302339_meta_cr8_reconstruction_manifest.csv`

Before rerunning long author notebook exports, audit the active Python
environment for imported and indirect runtime dependencies:

```bash
python scripts/06_celltype/08_audit_gse302339_author_dependencies.py --strict
```

The dependency audit scans the exported notebook `.py` files, adds known
runtime-only extras such as `harmonypy` when `harmony_integrate` is present,
and prints a `python -m pip install ...` command for any required modules that
are missing. It writes:

- `results/meta/gse302339_author_dependency_audit_summary.json`
- `results/tables/gse302339_author_dependency_manifest.csv`
- `results/tables/gse302339_author_dependency_source_manifest.csv`

If the installed Scanpy/HarmonyPy combination raises an `X_pca_harmony` shape
mismatch, patch the exported preprocessing code to call HarmonyPy directly and
validate the output orientation:

```bash
python scripts/06_celltype/09_patch_gse302339_author_preprocessing_runtime.py --strict
```

This modifies only
`intermediate/gse302339_scanpy_workflow_code/1_preprocessing_doublet_detection.py`
and writes:

- `results/meta/gse302339_author_preprocessing_patch_summary.json`

If the cell-type annotation notebook reaches the main annotation outputs but
then tries to read the optional downstream ABT/meta-merge artifact
`output/adata_mergedmeta_abt_cr8`, stop it cleanly before that section for the
parenchymal-signature workflow:

```bash
python scripts/06_celltype/10_patch_gse302339_celltype_annotation_for_signature_run.py --strict
```

This modifies only
`intermediate/gse302339_scanpy_workflow_code/2_celltype_annotation.py` and
writes:

- `results/meta/gse302339_celltype_annotation_signature_patch_summary.json`

After the authors' annotation workflow has been rerun, or if the pickle/H5AD
artifact is otherwise available, first audit the saved author objects against
the extracted cluster-to-cell-type dictionaries:

```bash
python scripts/06_celltype/11_audit_gse302339_saved_author_objects.py --strict
```

This is a guardrail before deconvolution. It checks whether the saved object
naturally contains the required labels (`AT1`, `AT2`, and `Fibroblast`) in its
own `.obs` annotation column, and separately reports whether cluster IDs from
the extracted notebook dictionaries overlap the saved object's Leiden clusters.
Do not force a cluster map onto a saved object unless the audit and notebook
context show that both came from the same clustering.

It writes:

- `results/meta/gse302339_saved_author_object_audit_summary.json`
- `results/tables/gse302339_saved_author_object_summary.csv`
- `results/tables/gse302339_saved_object_cluster_map_compatibility.csv`
- `results/tables/gse302339_saved_object_required_celltype_compatibility.csv`

Before building signatures, also audit whether the saved objects contain a
usable expression source. The reconstructed author objects may have scaled or
regressed values in `.X`, which are not safe for logCPM signature construction:

```bash
python scripts/06_celltype/12_audit_gse302339_author_expression_slots.py
```

This writes:

- `results/meta/gse302339_author_expression_slot_audit_summary.json`
- `results/tables/gse302339_author_expression_slot_audit.csv`

To compare the author notebook's recorded clustering outputs against your
locally reconstructed saved objects, run:

```bash
python scripts/06_celltype/14_compare_gse302339_notebook_expected_clusters.py --strict
```

This is a hard guardrail for the GSE302339 reference. The deposited notebook
recorded 41 parenchyma Leiden clusters before applying the
`parenchyma_celltype_level1` dictionary. If your reconstructed
`output/parenchyma_harmony_annotated_cr8` has fewer clusters, the author
cluster IDs are not compatible with that object and should not be used to build
AT1/AT2/Fibroblast signatures. It writes:

- `results/meta/gse302339_notebook_expected_cluster_comparison_summary.json`
- `results/tables/gse302339_notebook_expected_cluster_comparison.csv`

### Checkpoint-matched author annotation replay

The public deposits do not contain an annotated H5AD/pickle or a barcode-level
cell-type table. To recover the deposited notebook annotations without manual
marker reannotation or cross-object cluster-map transfer, create the isolated
author runtime and prepare a deterministic replay:

```bash
conda env create -f environment.gse302339-author.yml
conda activate gse302339_author

python scripts/06_celltype/15_prepare_gse302339_author_exact_replay.py \
  --prepare-code \
  --strict
```

If the first environment creation failed during pip installation, remove that
partial env and recreate it from the current file:

```bash
conda deactivate
conda env remove -n gse302339_author
conda env create -f environment.gse302339-author.yml
conda activate gse302339_author
```

The preflight reconstructs the authors' 65-file processing order from unique
raw cell counts, restores the global seed, extracts the exact 185-gene
ribosomal tuple, checks the metadata join, and generates guarded notebook code.
It also requires at least 25 GiB free and skips serialization of the unused
2.3-GiB concatenated checkpoint. It must report
`ready_to_run_exact_replay: true` before the long run:

```bash
set -o pipefail
python intermediate/gse302339_author_exact_code/1_preprocessing_doublet_detection.py \
  2>&1 | tee logs/gse302339_author_exact_preprocessing.log
python intermediate/gse302339_author_exact_code/2_celltype_annotation.py \
  2>&1 | tee logs/gse302339_author_exact_annotation.log
```

Do not run notebook 2 if notebook 1 exits non-zero. The annotation step depends
on the integrated object written by the preprocessing replay after all gates pass.

The generated scripts check every sample plus the 160,620-by-18,941 concat,
160,620-by-2,323 HVG, full 62-cluster, parenchyma 41-cluster, and immune
38-cluster checkpoints before assigning labels. Then export a stable
sample-plus-barcode annotation table:

```bash
python scripts/06_celltype/16_export_gse302339_author_exact_annotations.py --strict
```

The export also requires the 26-iteration Harmony checkpoint and guarded
annotation-completion message in the two replay logs.

Full rationale and commands are in
`docs/gse302339_author_annotation_replay.md`. This reconstructs the deposited
notebook annotation (160,620 cells), not the unexplained paper-final 128,433-cell
subset.

If the exact-replay object passes that audit, build a GeoMx-compatible parenchymal
signature matrix from the reconstructed author object:

```bash
python scripts/06_celltype/05_build_gse302339_author_signatures.py \
  --input-object output/gse302339_author_exact/parenchyma_harmony_annotated_cr8 \
  --cell-type-column parenchyma_celltype_level1 \
  --expression-source raw \
  --signature-transform logcpm \
  --signature-output data/processed/gse292993/gse302339_author_exact_parenchyma_signatures_rawX_logcpm.csv \
  --no-write-h5ad \
  --strict
```

By default this expects `output/parenchyma_harmony_annotated_cr8` and the
`parenchyma_celltype_level1` column, then requires `AT1`, `AT2`, and
`Fibroblast` to be present in the emitted signature. It writes:

- `results/meta/gse302339_author_signature_reference_summary.json`
- `results/tables/gse302339_author_signature_celltype_counts.csv`
- `results/tables/gse302339_author_signature_manifest.csv`
- `data/processed/gse292993/gse302339_author_exact_parenchyma_signatures_rawX_logcpm.csv`

Do not apply a cluster map across different saved objects. For example, the
`parenchyma_celltype_level1` map extracted from `2_celltype_annotation.py` is
assigned to the notebook variable `parenchyma` after its own parenchymal
subclustering. It should not be applied to `output/adata_harmony_annotated_cr8`,
even if some cluster IDs overlap numerically. Cluster IDs are not stable cell
type identifiers across object states, graph construction, package versions, or
Leiden runs. The signature builder refuses this cross-object relabeling by
default.

Only build a pneumocyte/fibroblast signature from an object that naturally
contains the required labels in its own `.obs` column, or from the exact object
whose clustering generated the author map. If the current reconstructed
`parenchyma_harmony_annotated_cr8` object lacks `Fibroblast`, this command
should fail under `--strict`. Do not force the extracted cluster dictionary onto
that divergent object. Use the checkpoint-matched script-15 replay above, or
obtain the authors' saved annotated object.

The signature builder uses `--expression-source auto` by default: it prefers the
requested layer, then `raw.X`, then `.X`. Use the expression-slot audit above to
confirm which source was selected. In notebook 1, `.raw` is assigned before
`normalize_total` and `log1p`, so the checkpoint-matched replay's `raw.X` is the
count matrix. Use `--signature-transform logcpm` for that source.

Use that signature matrix with NNLS deconvolution:

```bash
python scripts/06_celltype/01_deconvolve_gse292993_compartment_nnls.py \
  --signature-csv data/processed/gse292993/gse302339_author_exact_parenchyma_signatures_rawX_logcpm.csv \
  --strict
```

To project the same reference signatures onto all three GeoMx ROI compartments,
run the NNLS step once per compartment:

```bash
for compartment in airway parenchyma vessel; do
  python scripts/06_celltype/01_deconvolve_gse292993_compartment_nnls.py \
    --compartment "$compartment" \
    --signature-csv data/processed/gse292993/gse302339_author_exact_parenchyma_signatures_rawX_logcpm.csv \
    --strict
done
```

Then plot donor-balanced cell composition across airway, parenchyma, and vessel
ROIs:

```bash
python scripts/09_figures/01_plot_gse292993_cell_composition_compartments.py --strict
```

This writes:

- `results/figures/gse292993_cell_composition/gse292993_nnls_cell_composition_all_compartments_stacked`
- `results/figures/gse292993_cell_composition/gse292993_nnls_cell_composition_all_compartments_heatmap`
- `results/tables/gse292993_nnls_cell_composition_all_compartments_summary.csv`
- `results/tables/gse292993_nnls_cell_composition_all_compartments_donor_means.csv`
- `results/meta/gse292993_nnls_cell_composition_all_compartments_plot_summary.json`

The plotting script first averages ROI fractions within each donor, then
averages donors within each diagnosis/compartment. This makes the visual less
sensitive to the larger COPD ROI count.

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
airway/parenchyma/vessel compartment labels, and known donor IDs. Unknown
diagnosis or compartment rows are excluded from the DKK3 biological summaries
and retained only in QC plots/tables. ROI summaries are descriptive;
donor-compartment summaries are the unit to carry forward into inference.

After the WTA matrix, scRNA signatures, and DKK3 ROI summaries exist, run the
baseline NNLS deconvolution for the selected GeoMx compartment. The default is
`parenchyma` because the first biological question is pneumocyte/fibroblast
composition:

```bash
python scripts/06_celltype/01_deconvolve_gse292993_compartment_nnls.py --strict
```

This writes:

- `results/meta/gse292993_parenchyma_nnls_deconvolution_summary.json`
- `results/tables/gse292993_parenchyma_nnls_deconvolution_roi.csv`
- `results/tables/gse292993_parenchyma_nnls_deconvolution_donor.csv`
- `results/tables/gse292993_parenchyma_nnls_deconvolution_dkk3_correlations.csv`
- `results/tables/gse292993_parenchyma_nnls_deconvolution_celltype_manifest.csv`

The NNLS deconvolution is a first-pass composition estimate, not a final cell
assignment. Use it to ask whether parenchymal DKK3 tracks fibroblast/stromal,
pneumocyte, endothelial, or immune composition before moving to
ligand-receptor analysis.

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

Create the donor-level DKK3 summary figure across airway, parenchyma, and
vessel compartments after the stratified tests are available:

```bash
python scripts/09_figures/00_plot_gse292993_dkk3_compartments.py
```

This writes PNG, SVG, and PDF versions of:

- `results/figures/gse292993_dkk3/gse292993_dkk3_all_compartments_donor_signal`

The figure shows donor-level `Non Smoker`, `Smoker`, and `COPD` distributions
for median `log1p(DKK3 CPM)`, above-LOQ fraction, and median raw DKK3 count in
each biological compartment. To recreate a single-compartment figure, pass a
compartment explicitly:

```bash
python scripts/09_figures/00_plot_gse292993_dkk3_compartments.py --compartment parenchyma
```

To deliberately include unresolved compartment labels for a diagnostic plot:

```bash
python scripts/09_figures/00_plot_gse292993_dkk3_compartments.py --compartment diagnostic
```
