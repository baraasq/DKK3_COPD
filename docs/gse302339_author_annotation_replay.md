# GSE302339 author-annotation replay

## Goal

Recover the cell labels assigned by the deposited GSE302339 Scanpy workflow
without manually reannotating markers and without applying cluster IDs to a
different clustering.

No annotated H5AD, pickle, or barcode-to-cell-type table is present in the GEO
or Zenodo deposits. The recoverable public target is therefore the **deposited
notebook annotation**, reconstructed in the notebook's own runtime and accepted
only when its recorded checkpoints match.

## Why the earlier reconstruction diverged

The executed notebook records Python 3.12.8, Scanpy 1.10.4, AnnData 0.11.1,
DoubletDetection 4.2, and an author-era HarmonyPy release. The earlier run used
newer packages, skipped the global seed set by `scvi.settings.seed = 0`, used a
different `os.listdir()` sample order, and approximated the missing ribosomal
sidecar. Those differences changed doublet calls, Harmony, and Leiden clusters.

The exact replay fixes all four sources of divergence:

1. an isolated, pinned environment;
2. `random.seed(0)` and `numpy.random.seed(0)` as the lightweight equivalent of
   the notebook's seed-setting call for the used code paths;
3. the 65-file author order reconstructed from the unique raw-cell totals in
   notebook output;
4. the exact 185-gene ribosomal tuple embedded in the deposited preprocessing
   code.

## Run on m5server

From the project root:

```bash
cd /mnt/flowlib/baraa/projects/COPD_public

conda env create -f environment.gse302339-author.yml
conda activate gse302339_author

python scripts/06_celltype/15_prepare_gse302339_author_exact_replay.py \
  --prepare-code \
  --strict
```

If an earlier environment creation failed partway through pip installation,
remove the partial environment before recreating it:

```bash
conda deactivate
conda env remove -n gse302339_author
conda env create -f environment.gse302339-author.yml
conda activate gse302339_author
```

The preflight must report `ready_to_run_exact_replay: true`. It audits every
installed version, maps all 65 H5 files into the author order, checks the exact
metadata join keys, writes the 185-gene sidecar, and generates guarded copies of
notebooks 1 and 2. It also requires at least 25 GiB free by default. The replay
does not serialize the unused 2.3-GiB concatenated checkpoint, but retains the
integrated and annotated objects needed for validation. Override the threshold
only after checking both filesystem space (`df -h`) and the per-user quota
(`quota -s`, where available).

Run the generated workflow with pipe failures preserved:

```bash
set -o pipefail

python intermediate/gse302339_author_exact_code/1_preprocessing_doublet_detection.py \
  2>&1 | tee logs/gse302339_author_exact_preprocessing.log

python intermediate/gse302339_author_exact_code/2_celltype_annotation.py \
  2>&1 | tee logs/gse302339_author_exact_annotation.log
```

The generated code stops immediately if any of these deposited-notebook
checkpoints diverge:

- all 65 per-sample raw, post-low-quality-filter, and post-doublet dimensions;
- concatenated matrix: 160,620 cells by 18,941 genes;
- highly variable-gene object: 160,620 cells by 2,323 genes;
- full object: 62 Leiden clusters with IDs 0 through 61;
- parenchyma object: 41 Leiden clusters;
- immune object: 38 Leiden clusters.

Only after the corresponding object passes its gates does the generated code
apply the authors' object-local dictionary. This is a replay of the authors'
annotation operation, not cross-object cluster-ID transfer.

## Make labels independent of cluster IDs

After both generated scripts finish:

```bash
python scripts/06_celltype/16_export_gse302339_author_exact_annotations.py --strict
```

This rechecks all object and label gates and writes:

```text
data/external/scrna_reference/gse302339_author_exact_cell_annotations.tsv.gz
```

The table is keyed by `batch + original 10x barcode`, contains broad and
detailed author labels, and is the stable annotation artifact for downstream
work. It also verifies that the preprocessing log contains the deposited
Harmony convergence checkpoint (26 iterations) and that the annotation log
reached its guarded completion point. Once exported, downstream analyses do not
use Leiden IDs.

## Build the parenchymal reference

Build signatures from all author-labeled parenchymal classes, retaining AT1,
AT2, fibroblast, endothelial, smooth-muscle, mesothelial, club, and airway
epithelial backgrounds:

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

Do not use the earlier forced full-object signature or its NNLS outputs.

## Scope limitation

This reconstructs the deposited notebook annotation. The paper reports 128,433
final cells and different QC/version details, whereas the deposited executed
notebook records 160,620 cells and does not include the later paper-final
filter. A paper-final, bit-for-bit label set requires an author-provided object
or a small table containing sample, original barcode, and final cell type.
