# ATRT Tumor Microenvironment Cell-Cell Communication Analysis using snRNAseq data

Reproducible analysis code for astrocyte, microglial, and cell–cell communication analyses of the GSE283839 ATRT single-nucleus RNA-sequencing dataset.


CellChat v2 analysis of ligand-receptor interactions between tumor cells, astrocytes and microglia in atypical teratoid/rhabdoid tumor (ATRT) snRNA-seq data (GSE283839).

## Data

Download the following two files from [GEO accession GSE283839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283839) and place them in the same directory as the script:

| File | Description |
|---|---|
| `GSE283839_ATRT_RNA_V3_counts_raw.tsv` | Raw gene-by-nucleus count matrix |
| `GSE283839_ATRT_RNA_V3_metadata.tsv` | Cell-level annotations and sample metadata |

## Requirements

- R >= 4.3
- Bioconductor >= 3.18

### Install dependencies

Run the following once before executing the pipeline:

```r
# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"))

# CellChat v2
if (!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("jinworks/CellChat")

# CRAN packages
install.packages(c(
  "Seurat", "data.table", "future",
  "ggplot2", "patchwork", "circlize",
  "scales", "RColorBrewer", "viridis",
  "dplyr", "tidyr"
))
```

## Usage

```r
source("CellChat_analysis.R")
```

All outputs are written to `cellchat_results/` in the working directory:

```
cellchat_results/
  figures/    # PDF and PNG publication figures
  rds/        # Serialised CellChat objects
  tables/     # CSV interaction tables
```

## Citation

If you use this code, please cite the source dataset:

> Blanco-Carmona E et al. A cycling, progenitor-like cell population at the base of atypical teratoid rhabdoid tumor subtype differentiation trajectories. *Neuro Oncol*, 2025.

and CellChat v2:

> Jin S et al. Inference and analysis of cell-cell communication using CellChat. *Nature Communications*, 2021.
