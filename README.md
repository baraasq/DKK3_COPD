# DKK3 COPD Spatial Revisit

This repository adapts the DKK3-focused IPF workflow in
[`baraasq/DKK3_IPF`](https://github.com/baraasq/DKK3_IPF) to COPD.

## Primary dataset

The primary disease-spatial dataset is **GSE292993**, a NanoString GeoMx
whole-transcriptome study with COPD and control regions of interest from lung
parenchyma, airways, and vessels.

`DKK3` is included in the whole-transcriptome panel. This makes GSE292993 the
primary dataset for direct spatial tests of DKK3 in COPD. The earlier
Firsova/SCRINSHOT work is retained as a secondary cell-resolved neighborhood
validation track; it does not measure DKK3.

## Analysis units

- **Biological replicate:** donor
- **Repeated observations:** GeoMx ROIs nested within donor
- **Primary compartments:** parenchyma, airway, and vessel
- **Primary question:** whether COPD-associated DKK3 expression differs by
  anatomical compartment and is linked to pneumocyte/fibroblast programs
- **Primary inference:** donor-aware models, not ROI-as-replicate tests

## Project structure

```text
.
├── config
│   ├── gse292993.yaml
│   ├── pipeline.yaml
│   └── project.toml
├── data
│   ├── raw/gse292993
│   │   ├── dcc
│   │   ├── geo_metadata
│   │   └── pkc
│   ├── processed/gse292993
│   ├── external
│   │   ├── scrna_reference
│   │   └── spatial_validation
│   └── spatial/gse292993
│       ├── images
│       └── roi_annotations
├── intermediate
│   ├── qc
│   ├── geomx_objects
│   ├── normalized
│   ├── pseudobulk
│   ├── dkk3
│   ├── deconvolution
│   ├── communication
│   │   ├── cellchat_dkk3
│   │   ├── dialogue_dkk3
│   │   └── misty_dkk3
│   ├── modules
│   └── validation
├── logs
├── meta
│   ├── source_manifest.tsv
│   ├── sample_sheet.csv
│   └── roi_metadata.csv
├── results
│   ├── archived
│   ├── qc
│   ├── dkk3
│   ├── clustering
│   ├── compartments
│   ├── figures
│   ├── niches
│   ├── pathways
│   └── tables
├── scripts
│   ├── 00_setup
│   ├── 01_download
│   ├── 02_metadata
│   ├── 03_geomx_qc
│   ├── 04_normalization
│   ├── 05_dkk3
│   ├── 06_celltype
│   ├── 07_communication
│   ├── 08_validation
│   ├── 09_figures
│   └── utils
├── tests
├── docs
├── dir_engineering.bash
└── environment.yml
```

Raw files are immutable. Derived analysis objects go in `data/processed` or
`intermediate`; only compact tables, figures, and reports belong in `results`.

## Reproduce the structure on m5server

```bash
bash dir_engineering.bash /mnt/flowlib/baraa/projects/COPD_public
```

The bootstrap script is idempotent: it creates missing directories and does
not delete or overwrite existing data.

## Immediate workflow

1. Download the public deposits with
   `scripts/01_download/download_all_deposits.bash`; inspect `sra-info` before
   authorizing the much larger PRJNA1282758 FASTQ conversion.
2. build donor-, slide-, segment-, ROI-, compartment-, and diagnosis-level
   metadata.
3. Audit the DCC files and WTA panel before filtering.
4. Run GeoMx QC and normalization without treating ROIs as independent
   patients.
5. Map DKK3 by compartment and disease at both ROI and donor-summary levels.
6. Reconstruct the deposited GSE302339 author annotations with the
   checkpoint-matched replay in
   [`docs/gse302339_author_annotation_replay.md`](docs/gse302339_author_annotation_replay.md),
   then resolve pneumocyte/fibroblast signals without cross-object cluster-ID
   transfer.
7. Add communication and cross-dataset validation only after the primary
   expression model is stable.

## Existing validation code

The original top-level `scripts/00_*` through `scripts/03_*` audit Firsova
SCRINSHOT, COPD scRNA, and healthy Visium/RRST inputs. They are retained for
the secondary validation track and are not the GSE292993 primary pipeline.
