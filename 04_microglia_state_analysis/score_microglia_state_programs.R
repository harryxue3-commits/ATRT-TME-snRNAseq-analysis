#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop(paste(
    "Usage: Rscript score_microglia_state_programs.R",
    "<metadata.tsv> <microglia_raw_counts.tsv.gz>",
    "<microglia_state_program_gene_sets.tsv> <output_dir>"
  ))
}

metadata_path <- normalizePath(args[[1]], mustWork = TRUE)
counts_path <- normalizePath(args[[2]], mustWork = TRUE)
gene_sets_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- normalizePath(args[[4]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_gzip_tsv <- function(path) {
  fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE)
}

metadata <- fread(metadata_path, sep = "\t")
microglia_metadata <- metadata[
  Final_Annotation == "Microglia" & Final_Annotation_Focus == "TME"
]
stopifnot(nrow(microglia_metadata) == 288L)

count_table <- read_gzip_tsv(counts_path)
genes <- count_table[[1]]
counts <- as.matrix(count_table[, -1])
storage.mode(counts) <- "numeric"
rownames(counts) <- genes
stopifnot(setequal(colnames(counts), microglia_metadata$Cell))
counts <- counts[, microglia_metadata$Cell, drop = FALSE]

library_sizes <- colSums(counts)
log_cp10k <- log1p(sweep(counts, 2, pmax(library_sizes, 1), "/") * 1e4)

gene_sets <- fread(gene_sets_path)
program_order <- unique(gene_sets$module)
program_scores <- list()
gene_usage <- list()

for (program in program_order) {
  requested <- gene_sets[module == program, gene]
  available <- intersect(requested, rownames(log_cp10k))
  retained <- available[rowSums(counts[available, , drop = FALSE] > 0) >= 9L]
  if (length(retained) == 0L) stop(sprintf("No genes retained for program: %s", program))
  z_scores <- t(scale(t(log_cp10k[retained, , drop = FALSE])))
  z_scores[!is.finite(z_scores)] <- 0
  program_scores[[program]] <- colMeans(z_scores)
  gene_usage[[program]] <- retained
}

cell_scores <- data.table(Cell = colnames(log_cp10k))
for (program in names(program_scores)) cell_scores[[program]] <- program_scores[[program]]
cell_scores <- merge(
  cell_scores,
  microglia_metadata[, .(Cell, ID, subtype)],
  by = "Cell"
)

cell_scores[, m1_m2_category := fifelse(
  m1_inflammatory > 0 & m2_immunoregulatory > 0, "Both above zero",
  fifelse(
    m1_inflammatory > 0, "M1 only",
    fifelse(m2_immunoregulatory > 0, "M2 only", "Neither above zero")
  )
)]

score_matrix <- as.matrix(cell_scores[, ..program_order])
cell_scores[, dominant_program := program_order[max.col(score_matrix, ties.method = "first")]]
fwrite(cell_scores, file.path(output_dir, "microglia_cell_state_scores.tsv"), sep = "\t")

usage <- rbindlist(lapply(names(gene_usage), function(program) {
  data.table(
    module = program,
    n_genes = length(gene_usage[[program]]),
    genes_used = paste(gene_usage[[program]], collapse = ",")
  )
}))
fwrite(usage, file.path(output_dir, "microglia_program_genes_retained.tsv"), sep = "\t")

gene_expression <- rbindlist(lapply(program_order, function(program) {
  program_genes <- intersect(gene_sets[module == program, gene], rownames(counts))
  rbindlist(lapply(program_genes, function(gene) {
    values <- counts[gene, ]
    data.table(
      module = program,
      gene = gene,
      retained_for_scoring = gene %in% gene_usage[[program]],
      cells_detected = sum(values > 0),
      percent_detected = 100 * mean(values > 0),
      mean_raw_count = mean(values),
      mean_log_cp10k = mean(log_cp10k[gene, ])
    )
  }))
}))
fwrite(gene_expression, file.path(output_dir, "microglia_program_gene_expression.tsv"), sep = "\t")

category_summary <- cell_scores[, .(n = .N), by = m1_m2_category]
category_summary[, percent := 100 * n / sum(n)]
fwrite(category_summary, file.path(output_dir, "microglia_m1_m2_category_summary.tsv"), sep = "\t")

dominant_summary <- cell_scores[, .(n = .N), by = dominant_program]
dominant_summary[, percent := 100 * n / sum(n)]
setorder(dominant_summary, -n, dominant_program)
fwrite(dominant_summary, file.path(output_dir, "microglia_dominant_program_summary.tsv"), sep = "\t")

long_scores <- melt(
  cell_scores,
  id.vars = c("Cell", "ID", "subtype", "m1_m2_category", "dominant_program"),
  measure.vars = program_order,
  variable.name = "module",
  value.name = "score"
)
sample_scores <- long_scores[, .(mean_score = mean(score)), by = .(ID, subtype, module)]
fwrite(sample_scores, file.path(output_dir, "microglia_sample_level_program_scores.tsv"), sep = "\t")

subtype_tests <- sample_scores[, {
  test <- tryCatch(kruskal.test(mean_score ~ subtype), error = function(error) NULL)
  .(p_value = if (is.null(test)) NA_real_ else test$p.value)
}, by = module]
subtype_tests[, adjusted_p_value := p.adjust(p_value, method = "BH")]
fwrite(subtype_tests, file.path(output_dir, "microglia_subtype_program_tests.tsv"), sep = "\t")

correlation <- data.table(
  comparison = "M1/inflammatory versus M2/immunoregulatory",
  method = "Spearman rank correlation",
  n_microglia = nrow(cell_scores),
  rho = cor(cell_scores$m1_inflammatory, cell_scores$m2_immunoregulatory, method = "spearman")
)
fwrite(correlation, file.path(output_dir, "microglia_m1_m2_correlation.tsv"), sep = "\t")

cat(sprintf("Scored %d microglial nuclei across %d programs.\n", nrow(cell_scores), length(program_order)))
