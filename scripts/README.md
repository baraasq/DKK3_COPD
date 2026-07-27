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

