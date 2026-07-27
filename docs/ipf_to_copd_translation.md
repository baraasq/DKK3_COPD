# IPF-to-COPD workflow translation

| Mature IPF component | COPD implementation with GSE292993 | Main change |
|---|---|---|
| Direct spatial DKK3 mapping | Map DKK3 across GeoMx ROIs | Use ROI/segment geometry instead of Visium spots |
| DKK3-high and DKK3-low regions | Define within compartment and preferably within donor | Avoid thresholds driven by airway/parenchyma mixing |
| Pneumocyte-fibroblast interface | Score epithelial and fibroblast programs in parenchymal ROIs | GeoMx ROIs are multicellular mixtures |
| Donor-aware disease inference | Model repeated ROIs nested within donor or aggregate donor-by-compartment | Essential |
| Spatial graph analysis | Use ROI layouts only if coordinates or adjacency are deposited | Do not invent spot-level adjacency |
| Cell-type deconvolution | Anchor to a COPD sc/snRNA reference | Validate reference/platform compatibility |
| CellChat/NicheNet | Run in the cell-resolved reference, prioritized by GeoMx DKK3 programs | Not directly identifiable from bulk ROI profiles |
| DIALOGUE programs | Derive multicellular programs in sc/snRNA and test scores in GeoMx | Requires gene-overlap and donor coverage audits |
| MISTy | Use only when valid ROI neighborhoods can be reconstructed | May not be supported by GEO files |
| Pathway analysis | Test DKK3-associated programs separately by compartment | Account for ROI composition |
| Patient reproducibility | Display donor-level estimates and leave-one-donor-out sensitivity | ROIs are not patients |
| Independent spatial validation | Use Xenium/SCRINSHOT for cell neighborhoods | These panels do not measure DKK3 |

## Recommended execution order

1. Download and checksum GSE292993.
2. Build the donor-slide-ROI-segment metadata table.
3. Verify DKK3 in the PKC and raw DCC files.
4. Perform GeoMx QC with donor and compartment balance audits.
5. Normalize and create donor-aware DKK3 summaries.
6. Fit the primary parenchymal COPD-versus-Control model.
7. Run airway and vessel analyses as secondary contrasts.
8. Score pneumocyte, fibroblast, ECM, WNT, and injury programs.
9. Orient DKK3 source/receiver states in external COPD sc/snRNA data.
10. Validate relevant cellular neighborhoods in Xenium/SCRINSHOT.
11. Add communication analyses only after the direct DKK3 result is stable.

