# ATRT tumor-microenvironment analysis code

This repository contains reproducible analysis scripts for the cohort-composition, astrocyte, microglial, Cell2Sentence, and CellChat analyses of the GSE283839 ATRT single-nucleus RNA-sequencing dataset.

Raw GEO data and generated results are not included in this repository.

## Repository layout

```text
01_cohort_composition/
02_astrocyte_state_analysis/
03_astrocyte_support_analysis/
04_microglia_state_analysis/
05_cell2sentence_analysis/
06_cellchat_analysis/
environment/
```

## Environment

Create the Python/R analysis environment with:

```bash
conda env create -f environment/environment.yml
conda activate atrt-tme-analysis
```

The Cell2Sentence 27B inference script is designed for Google Colab and contains its own installation cell. The CellChat script uses additional R dependencies listed in `06_cellchat_analysis/CellChat_analysis.R`.

## Input data

Download the GSE283839 raw count matrix and metadata from GEO:

- `GSE283839_ATRT_RNA_V3_counts_raw.tsv`
- `GSE283839_ATRT_RNA_V3_metadata.tsv`

Commands below use illustrative paths under `data/`; users may place files elsewhere and supply the corresponding paths.

## 1. Cohort composition

```bash
python 01_cohort_composition/calculate_cohort_composition.py \
  --metadata data/GSE283839_ATRT_RNA_V3_metadata.tsv \
  --counts data/GSE283839_ATRT_RNA_V3_counts_raw.tsv.gz \
  --output-dir results/cohort_composition
```

## 2. Astrocyte state programs

```bash
python 05_cell2sentence_analysis/prepare_tme_cell2sentence_anndata.py \
  --counts data/GSE283839_ATRT_RNA_V3_counts_raw.tsv.gz \
  --metadata data/GSE283839_ATRT_RNA_V3_metadata.tsv.gz \
  --output results/cell2sentence/GSE283839_TME_C2S_ready.h5ad

python 02_astrocyte_state_analysis/score_astrocyte_state_programs.py \
  --input results/cell2sentence/GSE283839_TME_C2S_ready.h5ad \
  --output-dir results/astrocyte_state_programs
```

## 3. Astrocyte support modules

```bash
python 03_astrocyte_support_analysis/analyze_astrocyte_support_modules.py \
  --counts data/GSE283839_ATRT_RNA_V3_counts_raw.tsv.gz \
  --metadata data/GSE283839_ATRT_RNA_V3_metadata.tsv.gz \
  --output-dir results/astrocyte_support_modules
```

## 4. Microglial state programs

```bash
python 04_microglia_state_analysis/extract_microglia_analysis_matrices.py \
  --counts data/GSE283839_ATRT_RNA_V3_counts_raw.tsv.gz \
  --metadata data/GSE283839_ATRT_RNA_V3_metadata.tsv \
  --output-dir results/microglia_inputs

Rscript 04_microglia_state_analysis/score_microglia_state_programs.R \
  data/GSE283839_ATRT_RNA_V3_metadata.tsv \
  results/microglia_inputs/microglia_raw_counts.tsv.gz \
  results/microglia_inputs/microglia_state_program_gene_sets.tsv \
  results/microglia_state_programs
```

## 5. Cell2Sentence analysis

```bash
python 05_cell2sentence_analysis/generate_cell2sentence_ranked_gene_sentences.py \
  --input results/cell2sentence/GSE283839_TME_C2S_ready.h5ad \
  --output results/cell2sentence/GSE283839_TME_cell_sentences_top1000.csv
```

Run `05_cell2sentence_analysis/run_cell2sentence_27b_cell_type_prediction_colab.py` in Google Colab. Change `INPUT_CSV`, `OUTPUT_CSV`, and `CHECKPOINT_CSV` to the corresponding Google Drive locations.

Then calculate supplementary rank-weighted state-program scores:

```bash
python 05_cell2sentence_analysis/score_cell2sentence_state_programs.py \
  --metadata data/GSE283839_ATRT_RNA_V3_metadata.tsv \
  --predictions results/cell2sentence/GSE283839_TME_27B_predictions.csv \
  --output-dir results/cell2sentence/state_programs
```

## 6. CellChat analysis

The CellChat v2 analysis is in:

```bash
06_cellchat_analysis/CellChat_analysis.R
```

Place the raw count matrix and metadata in the working directory expected by the script, then run:

```r
source("06_cellchat_analysis/CellChat_analysis.R")
```

CellChat outputs are written to `cellchat_results/` with subfolders for figures, serialized objects, and tables.

## Reproducibility notes

- Required input and output paths are supplied explicitly where possible.
- Raw data, generated result directories, large GEO files, and local system files are intentionally excluded.
- Random seeds and analysis thresholds are defined in the corresponding scripts.
- Communication results represent transcriptional compatibility or computational predictions, not demonstrated signaling.
