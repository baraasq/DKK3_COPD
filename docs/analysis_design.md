# Analysis design

## Primary question

Does spatial DKK3 expression differ between COPD and control lung, within
parenchymal, airway, and vascular compartments, and is it associated with
pneumocyte/fibroblast remodeling programs?

The primary dataset is GSE292993, a GeoMx whole-transcriptome dataset in which
DKK3 is measured directly in COPD tissue.

## Evidence layer A: GSE292993 data audit

Before biological filtering:

1. Download the GEO metadata, DCC files, PKC panel, and any slide/ROI
   annotations.
2. Generate checksums and record them in `meta/source_manifest.tsv`.
3. Resolve donor, diagnosis, GOLD stage, smoking status, slide, ROI, segment,
   compartment, ROI area, and nuclei count.
4. Confirm DKK3 in the panel and in the raw DCC files.
5. Reconcile the number of donors, slides, ROIs, and segments against the
   publication and GEO record.
6. Preserve raw files unchanged.

No downstream result is considered valid until each ROI maps to a donor and
an anatomical compartment.

## Evidence layer B: GeoMx QC and normalization

QC thresholds must be derived from the study distributions and negative
probes rather than copied from Visium or single-cell pipelines.

Report QC by:

- donor;
- diagnosis and GOLD stage;
- slide;
- compartment;
- ROI and segment;
- sequencing/binding metrics;
- nuclei count and ROI area when available.

The primary normalized representation is Q3. TMM is a sensitivity analysis.
Raw counts, QC-filtered counts, and normalized values remain separate.

## Evidence layer C: direct spatial DKK3 analysis

### Descriptive analysis

- DKK3 count, normalized expression, and detection by ROI.
- DKK3 maps on available slide/ROI layouts.
- Donor-level DKK3 summaries within parenchyma, airway, and vessel.
- COPD severity trends shown as donor-level points.
- DKK3 correlations with epithelial, stromal, immune, and injury modules.

### Primary inference

ROIs are repeated observations nested within donors. The primary disease
contrast must therefore account for donor clustering or first aggregate
appropriate ROIs to donor-by-compartment summaries.

The prespecified primary contrast is:

```text
COPD versus Control within parenchymal ROIs
```

Secondary contrasts evaluate airway and vessel compartments. A pooled model
requires explicit disease-by-compartment terms and donor-aware uncertainty.

Report:

- donor counts;
- ROI counts;
- effect size and confidence interval;
- model specification;
- compartment;
- sensitivity to normalization and QC choices.

An ROI-only p-value that treats every ROI as an independent patient is not a
valid disease-level result.

## Evidence layer D: pneumocyte/fibroblast programs

GeoMx measures multicellular ROIs rather than individual cells. Pneumocyte and
fibroblast biology will therefore be inferred in stages:

1. score prespecified AT1, AT2, transitional epithelial, fibroblast, activated
   fibroblast, ECM, WNT, and injury modules;
2. test their association with DKK3 within each compartment;
3. adjust for broad composition estimates where supported;
4. anchor cell-state interpretation to an external COPD sc/snRNA-seq
   reference;
5. verify that conclusions reproduce across donors.

DKK3-to-module association does not by itself establish which cell secreted
DKK3 or which cell received a signal.

## Evidence layer E: communication analysis

CellChat, NicheNet, DIALOGUE, or related tools should be run on a cell-resolved
COPD reference, not directly on GeoMx ROIs as if each ROI were a cell.

The GeoMx analysis supplies:

- disease- and compartment-specific DKK3 evidence;
- DKK3-associated tissue programs;
- spatially defined contexts to prioritize.

The single-cell reference supplies:

- plausible DKK3-producing states;
- pneumocyte and fibroblast receiver states;
- candidate receptor/pathway programs;
- donor-aware cell-state contrasts.

Communication results remain hypothesis-generating unless supported by direct
receptor/pathway evidence and spatial co-occurrence.

## Secondary spatial validation

### GSE313006 Xenium

Use for cell-resolved pneumocyte/fibroblast neighborhoods and disease-associated
cellular communities. DKK3 is absent from its targeted panel, so it cannot
validate DKK3 expression directly.

### Firsova COPD SCRINSHOT

Use for cell-resolved COPD neighborhood composition. DKK3 is absent from the
COPD SCRINSHOT panel.

### GSE237120 lymphoid follicles

Use only for follicle-specific sensitivity questions. It should not replace
the parenchymal primary analysis because the ROI selection answers a different
biological question.

## Interpretation guardrails

- Donors, not ROIs or cells, are biological replicates.
- Multiple ROIs and slides from one donor remain nested within that donor.
- GeoMx ROIs are multicellular mixtures, not single-cell measurements.
- Compartment differences must not be interpreted as disease differences.
- Detection of DKK3 in a whole-transcriptome panel does not identify its source
  cell without supporting evidence.
- Cross-sectional associations do not establish signaling direction or
  causality.

