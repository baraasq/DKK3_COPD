# GSE302339 author cell-type marker provenance

This project uses a marker-program fallback when the deposited GSE302339 code
cannot reproduce the authors' exact final cell-type clusters.

The paper's Methods section states that cell-type prediction was first performed
with the Python implementation of decoupleR using PanglaoDB, then manually
examined with a canonical cell-type marker dictionary. The paper's Table 3 and
the deposited `scanpy_workflow/2_celltype_annotation.ipynb` contain the
canonical marker dictionary.

## Markers explicitly listed by the authors

The extracted author notebook marker dictionary includes:

- Alveolar epithelia: `SFTPC`, `SFTPB`, `EPCAM`, `AGER`
- Airway epithelia: `SCGB1A1`, `MUC4`
- Lung endothelia: `PECAM1`, `LYVE1`
- Fibroblast: `ACTA2`, `PDGFRB`, `COL5A1`, `COL3A1`, `COL18A1`, `MPG`
- Mesothelia: `ITLN1`, `UPK3B`
- Basophil/Mast cell: `CPA3`, `TPSAB1`
- B cell: `CD19`, `MS4A1`, `BANK1`, `CD79A`
- Plasma cells: `MZB1`, `HSP90B1`, `FNDC3B`, `PRDM1`, `IGKC`, `JCHAIN`
- pDC: `IL3RA`, `GZMB`, `COBLL1`, `TCF4`
- cDC/cDC1/cDC2/Migratory DC: `BATF3`, `CLEC9A`, `CADM1`, `CST3`, `COTL1`,
  `LYZ`, `CLEC10A`, `FCER1A`, `DMXL2`, `CCR7`, `FSCN1`, `LAMP3`, `CCL22`
- Monocyte/macrophage/neutrophil: `FCN1`, `CD14`, `TCF7L2`, `FCGR3A`, `LYN`,
  `VCAN`, `CCR2`, `CX3CR1`, `CSF1R`, `CD74`, `MSR1`, `MARCO`, `FBP1`, `APOE`,
  `FABP4`, `LYZ`, `PPARG`, `MRC1`, `LTA4H`, `CTSD`, `CTSL`, `AQP9`, `FCGR3B`,
  `CXCR2`, `IL1R2`, `MMP9`, `CSF3R`, `S100A8`, `S100A9`, `CST3`
- T/NK/ILC/MAIT: `CD3D`, `CD3E`, `CD3G`, `TRAC`, `TRBC1`, `TRBC2`, `LCK`,
  `FYN`, `CD4`, `CD40LG`, `FOXP3`, `IL2RA`, `MKI67`, `CD8A`, `CD8B`, `GNLY`,
  `NKG7`, `NCAM1`, `CD247`, `GRIK4`, `FCER1G`, `TYROBP`, `KLRG1`, `TRDV2`,
  `TRGV9`, `TRGV10`, `TRDC`, `CD7`, `IL7R`, `ID2`, `PLCG2`, `SYNE1`,
  `SLC4A10`
- Megakaryocyte: `PPBP`

## How this is used here

For deconvolution, the marker-program script keeps AT1 and AT2 as separate
programs because GeoMx deconvolution benefits from AT1/AT2 distinction. The
authors' `Alveolar epithelia` marker set is therefore split cautiously:

- AT1 keeps AT1-enriched markers such as `AGER`, `HOPX`, `AQP5`, `CLDN18`.
- AT2 keeps AT2-enriched markers such as `SFTPC`, `SFTPB`, `SFTPA1`, `SFTPA2`,
  `ABCA3`, `NAPSA`, `LPCAT1`.
- `EPCAM` is not used as an AT1/AT2 deciding marker because it is broad
  epithelial and can reduce separation.

The author markers were added mainly to rescue under-selected stromal,
endothelial, airway, and immune cells while retaining marker-validation checks.
