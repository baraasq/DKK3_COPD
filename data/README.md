# Data directories

- `raw/gse292993`: immutable GEO downloads, including DCC and PKC files.
- `processed/gse292993`: analysis-ready expression objects and matrices.
- `spatial/gse292993`: ROI annotations, slide layouts, masks, and images.
- `external/scrna_reference`: external single-cell references.
- `external/spatial_validation`: secondary spatial cohorts.

Large data files are intentionally ignored by Git. Record every download,
checksum, and local path in `meta/source_manifest.tsv`.

