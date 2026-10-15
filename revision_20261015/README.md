# Revised analysis code: navigation and execution

This directory records the scripts used for the manuscript revision. Archived source files are unchanged and checksummed. They include historical stages required to understand the completed analysis; do not substitute an early 100/1,000-permutation stage for the final 10,000-permutation outputs.

## Project layout

`original_project/` retains relative paths from the original project. `cellchat_10000/original_project/` holds the final CellChat drivers and computational dependencies. To rerun, copy both trees into a **separate working project**, maintaining their relative paths, and provide the source data and saved inputs. The archive is not an installer or an automated end-to-end workflow. Do not execute scripts in-place expecting input data to be bundled.

1. `Cancers_ATRT_TME_Submission/Revision_20261005/analysis/code/01_extract_inputs.py`, `02_annotation.py`, `03_sensitivity.py`, and `08_collaborator_cohorts.py`: discovery extraction, TME-only identity refinement, program/threshold/donor sensitivities, and cell-ID exports. `analysis_config.json` and `reference/reference_gene_sets.json` supply the marker and reference definitions. Input counts/metadata must be restored under the referenced `NCBIGEO/snRNA/` paths. These scripts also reference existing original scoring inputs; inspect their paths before execution.
2. `ALSF-ATRT-snRNAseq/analysis_20261006/code/01_prepare.py` through `06_export.py`: external count/CNV preparation, scoring, annotation, sensitivities, expression checks, and export. Requires the downloaded ScPCA objects and source inferCNV results. The initial preparation expects the reference/configuration outputs from step 1.
3. `Cancers_ATRT_TME_Submission/Revision_20261007/analysis/code/reviewer1_analysis.py` and `Revision_20261008_AllReviewers/analysis/code/`: donor and external expression analyses, initial model preparation, denominator reconciliation, and metadata/C2S audit. Early CellChat preparation scripts retain their historical permutation counts; final inference uses step 4.
4. [Final CellChat instructions](cellchat_10000/README.txt): pooled subtype models, controlled addition of Unannotated nuclei only, and unified donor/conservative/external fits. Requires the saved input CellChat objects and cell manifests documented in that guide. The controlled comparison keeps original patients, glial cells, per-patient tumor counts and the 3,220-interaction database family fixed; the separate all-malignant comparison includes a broader tumor definition.
5. `Revision_20261010_CellChat10000/s16_microglia_update_20261012/code/original_sensitivity_analysis.py`: archived microglial sensitivity implementation supporting current Table S15 (historical S16).

## Reproducibility boundaries

- Original scripts are preserved byte-for-byte. Absolute paths, local Rcpp cache/compiler configuration, and original output directory assumptions must be adapted on another computer. Original manuscript figure/table numbers inside those scripts remain historical; use the root README numbering map.
- The CellChat package includes archived sessionInfo and equivalence-validation records, not expression data or serialized input models. Its checksum list covers the unmodified package contents. `ARCHIVE_SHA256SUMS.txt` covers the additional archived sources.
- Use source gene panels jointly with lineage exclusions, matched-null evidence, and reference/neighborhood support; the stricter population is not assumed to be ground truth.
- Annotation and scoring sensitivity sets overlap. Cells, donors, specimens, and patients are different denominators.
- CellChat within-model BH support does not establish between-subtype differences or functional signaling.
- The original scripts may output result figures as well as data. Manuscript editing and reviewer-response automation are deliberately excluded from this archive.
