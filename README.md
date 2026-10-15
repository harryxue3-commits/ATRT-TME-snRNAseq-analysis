# ATRT tumor-microenvironment analysis

Analysis code for astrocyte and microglial expression programs and expression-supported glia–tumor interaction candidates in GSE283839, with an external ScPCA comparison. Updated 15 October 2026 to accompany the revised manuscript.

## Start here

- [Revision analyses and run guide](revision_20261015/README.md): TME-restricted annotation, conservative/broader sensitivities, donor-level analyses, and external-cohort methods.
- [10,000-permutation CellChat code and original run instructions](revision_20261015/cellchat_10000/README.txt): discovery subtype, donor-specific, conservative-annotation, and external runs; controlled Unannotated inclusion; validation and session records.
- [Original pipeline commands](revision_20261015/README_AUGUST.md): the original 01–06 analyses. These alone do **not** reproduce all revised results. Paths in this archived guide are relative to the repository root.

## Current analysis conventions

Primary discovery analyses retain deposited cell labels and all malignant assignments. Marker-based refinement is restricted to the 516 source TME nuclei. Conservative and broader populations overlap and are sensitivity definitions, not independent replicates or definitive cell identities. Macrophage-like and BAM-like labels denote expression phenotypes, not proven origin.

Microglial scoring uses 75 retained gene–program entries representing 73 distinct genes (STAT1 and TREM2 each appear twice), with detection in at least nine of 288 nuclei. Continuous program scores are primary; historical `dominant_program` output keys mean highest-relative-score display assignment, not dominant biological function. M1/M2 keys retained in archived code correspond to inflammatory/immunoregulatory marker associations.

Astrocyte pseudobulk module scores average genes actually present in the count matrix. The support-module script now explicitly distinguishes absent genes from measured zeros and exports gene availability: GPX1 is unavailable, leaving 14 of 15 requested antioxidant/glutathione genes. See Table S12 for the requested-panel denominator sensitivity.

Current CellChat fits use **10,000 permutations**, finite P = (b + 1)/10,001, and BH correction across the full tested source–target–interaction family within each model before extracting highlighted candidates. Here b is the number of null scores **strictly greater** than the observed score, matching CellChat 2.1.2; zero observed scores receive P = 1. The native pathway summaries retain native P values; use the explicit finite-P/BH interaction outputs for revised inference. The separate donor-level correlation permutation analysis has its own settings and is not a CellChat test.

The original `06_cellchat_analysis/CellChat_analysis.R` and original pooled probabilities are retained as historical analyses, not relabeled as new 10,000-permutation results. Use the revision drivers linked above for current inference. Receptor-complex checks use complete database definitions, including SPP1–(ITGAV + ITGB1); candidate-panel single-subunit shorthand is not a complete receptor definition.

The external ScPCA cohort comprises eight specimens from four patients. Repeated specimens are linked before patient-level summaries; myeloid evidence is not automatically microglia-specific. Cell2Sentence provides a supplementary comparison, not independent validation. Its archived template settings are distinguished from unavailable historical execution logs, immutable model revision, and package versions (Table S18).

## Current figure/table map

| Analysis | Current manuscript location |
|---|---|
| Overall workflow / cohort landscape | Figures 1 / 2 |
| Astrocyte programs / paired support modules | Figures 3 / 4; Tables S2–S3, S8, S12 |
| Microglial programs | Figure 5; Tables S2, S10–S11, S15 |
| CellChat / summary schematic | Figures 6 / 7; Tables S4–S5, S13, S16–S17 |
| Conservative identity markers | Main Table 1 (formerly Table S15) |
| Annotation sensitivities / donor CellChat / external / C2S reproducibility | Tables S15 / S16 / S17 / S18 (formerly S16 / S17 / S18 / S19) |
| Gene terminology and abbreviations | Table S7 |

Archived filenames retain their original numbering for provenance; consult this map when reading them.

## Inputs, environment, and verification

Raw expression data, individual cell/patient metadata, manuscripts, credentials, and saved model objects are not included. Obtain the source datasets and restore required inputs before running the archive. Some original drivers retain machine-specific paths; the revision guide documents this limitation rather than claiming a turnkey rerun. Use the archived R session records for the CellChat environment; the original environment.yml is not a complete environment lock for the revision.

This repository update was checked for Python/R syntax and archive integrity, and the available-gene correction was tested on a synthetic matrix. The full computational analyses were not rerun during repository preparation. Expression-supported candidates require spatial, protein-level, and functional follow-up to establish biological communication.
