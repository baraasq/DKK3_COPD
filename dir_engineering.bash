#!/usr/bin/env bash
set -euo pipefail

# Idempotent project scaffold for the DKK3 COPD spatial revisit.
# Usage:
#   bash dir_engineering.bash /mnt/flowlib/baraa/projects/COPD_public

PROJECT_ROOT="${1:-/mnt/flowlib/baraa/projects/COPD_public}"

mkdir -p "${PROJECT_ROOT}"
cd "${PROJECT_ROOT}"

mkdir -p \
  config \
  data/raw/gse292993/dcc \
  data/raw/gse292993/pkc \
  data/raw/gse292993/geo_metadata \
  data/processed/gse292993 \
  data/external/scrna_reference \
  data/external/spatial_validation \
  data/spatial/gse292993/roi_annotations \
  data/spatial/gse292993/images \
  intermediate/qc \
  intermediate/geomx_objects \
  intermediate/normalized \
  intermediate/pseudobulk \
  intermediate/dkk3 \
  intermediate/deconvolution \
  intermediate/communication/cellchat_dkk3 \
  intermediate/communication/dialogue_dkk3 \
  intermediate/communication/misty_dkk3 \
  intermediate/modules \
  intermediate/validation \
  logs \
  meta \
  results/archived \
  results/qc \
  results/dkk3 \
  results/clustering \
  results/compartments \
  results/figures \
  results/niches \
  results/pathways \
  results/tables \
  scripts/00_setup \
  scripts/01_download \
  scripts/02_metadata \
  scripts/03_geomx_qc \
  scripts/04_normalization \
  scripts/05_dkk3 \
  scripts/06_celltype \
  scripts/07_communication \
  scripts/08_validation \
  scripts/09_figures \
  scripts/utils \
  tests \
  docs

printf 'COPD project structure ready at %s\n' "${PROJECT_ROOT}"

