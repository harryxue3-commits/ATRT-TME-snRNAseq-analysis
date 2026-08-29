# =============================================================================
# CellChat v2 — ATRT snRNA-seq (GSE283839)
# Cell-cell communication: Tumor × Astrocytes × Microglia
#
# Dataset  : GSE283839 (Jin et al.) — 17,564 nuclei, 12 ATRT samples
# Subtypes : ATRT-SHH (n=5), ATRT-TYR (n=4), ATRT-MYC (n=2, after QC)
#
# Analysis overview
#   1. Load raw counts + metadata, log-normalise, define focal populations
#   2. Sample balancing: per-group threshold exclusion + downsampling
#   3. CellChat objects: pan-cancer (cc_pan), per-subtype (cc_SHH/TYR/MYC),
#      and Microglia × tumour-state (cc_states)
#   4. Publication figures and interaction tables
#
# Output: cellchat_results/  figures/  rds/  tables/
# =============================================================================


# ── 0.  INSTALL (run once) ────────────────────────────────────────────────────
if (FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("ComplexHeatmap", "BiocNeighbors"))
  if (!requireNamespace("remotes", quietly = TRUE))
    install.packages("remotes")
  remotes::install_github("jinworks/CellChat")
  install.packages(c("Seurat", "data.table", "future",
                     "ggplot2", "patchwork", "circlize",
                     "scales", "RColorBrewer", "viridis",
                     "dplyr", "tidyr"))
}


# ── 1.  PACKAGES ─────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(CellChat); library(Seurat); library(SeuratObject)
  library(data.table); library(future)
  library(ggplot2); library(patchwork); library(ComplexHeatmap)
  library(circlize); library(scales); library(RColorBrewer)
  library(viridis); library(dplyr); library(tidyr)
})

# Seurat v4/v5 compatibility
SEURAT_V5 <- packageVersion("Seurat") >= "5.0.0"
get_norm <- function(obj) {
  if (SEURAT_V5) GetAssayData(obj, assay = "RNA", layer = "data")
  else           GetAssayData(obj, assay = "RNA", slot  = "data")
}

plan("sequential")
options(future.globals.maxSize = 8e9)

# ── 2.  PATHS & CONSTANTS ────────────────────────────────────────────────────
DATA_DIR <- "."
OUT_DIR  <- file.path(DATA_DIR, "cellchat_results")
FIG_DIR  <- file.path(OUT_DIR, "figures")
RDS_DIR  <- file.path(OUT_DIR, "rds")
TAB_DIR  <- file.path(OUT_DIR, "tables")
for (d in c(FIG_DIR, RDS_DIR, TAB_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Tumour progenitor state labels from the original annotation
TUMOR_LABELS <- c("NPC-like","IPC-like","OPC-like","RG-like",
                  "CP-like","Hypoxic","Cilia-like","Mesenchymal-like")

SUBTYPES <- c("ATRT-SHH","ATRT-TYR","ATRT-MYC")
ST_SHORT <- c("SHH","TYR","MYC")

# Minimum nuclei per CellChat group
MIN_CELLS <- 10

# Sample excluded globally: critically low tumour cell recovery
EXCLUDE_SAMPLES <- c("ATRT-MYC-2")

# Samples with >= 5 astrocyte nuclei — used for all Astrocyte group analyses
ASTRO_SAMPLES <- c("ATRT-SHH-4", "ATRT-TYR-2", "ATRT-MYC-3")

# Per-sample minimum nuclei for inclusion in each group (downsampling step)
MIN_TUMOR_PER_SAMPLE     <- 150
MIN_MICROGLIA_PER_SAMPLE <- 10

# Colour palettes
POP_COLORS <- c(Tumor = "#C0392B", Astrocytes = "#2471A3", Microglia = "#1E8449")
SUBTYPE_COLORS <- c("ATRT-SHH" = "#8E44AD", "ATRT-TYR" = "#E67E22", "ATRT-MYC" = "#2E86C1")

# Publication ggplot2 theme
pub_theme <- theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size =  9, colour = "grey40", hjust = 0.5),
    axis.text        = element_text(size =  9),
    axis.title       = element_text(size = 10),
    legend.text      = element_text(size =  8),
    legend.title     = element_text(size =  9, face = "bold"),
    strip.text       = element_text(size =  9, face = "bold"),
    strip.background = element_rect(fill = "grey95", colour = "grey70"),
    panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )

# Save figure as PDF + PNG
save_fig <- function(p, stem, w = 8, h = 6) {
  ggsave(file.path(FIG_DIR, paste0(stem, ".pdf")), p,
         width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG_DIR, paste0(stem, ".png")), p,
         width = w, height = h, dpi = 300)
  invisible(NULL)
}


# ── 3.  LOAD DATA ────────────────────────────────────────────────────────────
cat("Loading metadata...\n")
meta <- fread(file.path(DATA_DIR, "GSE283839_ATRT_RNA_V3_metadata.tsv"),
              sep = "\t", data.table = FALSE)
rownames(meta) <- meta$Cell

cat("Loading raw counts...\n")
dt         <- fread(file.path(DATA_DIR, "GSE283839_ATRT_RNA_V3_counts_raw.tsv"),
                    sep = "\t", data.table = FALSE)
gene_names <- dt[[1]]
dt         <- dt[, -1]
# Normalise barcode separator (counts file uses '.', metadata uses '-')
barcodes   <- sub("([ACGT]{16})\\.", "\\1-", colnames(dt))
mat        <- as.matrix(dt)
rm(dt); gc()
rownames(mat) <- gene_names
colnames(mat) <- barcodes
cat(sprintf("  Matrix: %d genes x %d cells\n", nrow(mat), ncol(mat)))


# ── 4.  BUILD SEURAT OBJECT ──────────────────────────────────────────────────
shared <- intersect(colnames(mat), rownames(meta))
cat(sprintf("  Matched barcodes: %d\n", length(shared)))

so <- CreateSeuratObject(counts = mat[, shared], meta.data = meta[shared, ])
rm(mat); gc()

so <- NormalizeData(so, normalization.method = "LogNormalize",
                    scale.factor = 10000, verbose = FALSE)

# Assign focal cell group labels
so$cell_group <- dplyr::case_when(
  so$Final_Annotation %in% TUMOR_LABELS ~ "Tumor",
  so$Final_Annotation == "Astrocytes"   ~ "Astrocytes",
  so$Final_Annotation == "Microglia"    ~ "Microglia",
  TRUE                                   ~ NA_character_
)

# Subset to focal populations, removing the excluded sample
so_sub <- subset(so,
                 cells = rownames(so@meta.data)[
                   !is.na(so$cell_group) &
                     !so@meta.data$ID %in% EXCLUDE_SAMPLES])
rm(so); gc()
cat(sprintf("Focal cells after exclusion: %d\n", ncol(so_sub)))

# Astrocyte-restricted subset: samples with >= 5 annotated astrocytes
so_sub_astro <- subset(so_sub,
                       cells = rownames(so_sub@meta.data)[
                         so_sub@meta.data$ID %in% ASTRO_SAMPLES])

cat("\nCell counts per sample:\n")
so_sub@meta.data %>%
  filter(cell_group %in% c("Astrocytes","Microglia","Tumor")) %>%
  count(ID, subtype, cell_group) %>%
  tidyr::pivot_wider(names_from = cell_group, values_from = n,
                     values_fill = 0L) %>%
  arrange(subtype, ID) %>%
  print(n = Inf)


# ── 5.  HELPER FUNCTIONS ─────────────────────────────────────────────────────

#' Exclude low-cell samples and downsample remaining samples per group.
#'
#' For each cell group independently:
#'   - Samples below the group-specific threshold are excluded.
#'   - Remaining samples are downsampled to the minimum per-sample count,
#'     ensuring no single donor dominates the expression averages.
#'
#' Thresholds are defined by MIN_TUMOR_PER_SAMPLE and MIN_MICROGLIA_PER_SAMPLE.
#' Astrocytes are passed through unchanged (threshold = 0).
downsample_per_sample <- function(obj,
                                  group_col  = "cell_group",
                                  sample_col = "ID",
                                  seed       = 42) {
  set.seed(seed)
  meta       <- obj@meta.data
  meta$.cell <- rownames(meta)
  
  thresholds <- c(Tumor      = MIN_TUMOR_PER_SAMPLE,
                  Microglia  = MIN_MICROGLIA_PER_SAMPLE,
                  Astrocytes = 0L)
  
  keep_all <- c()
  
  for (grp in sort(unique(meta[[group_col]]))) {
    grp_meta  <- meta[meta[[group_col]] == grp, ]
    threshold <- thresholds[grp]
    if (is.na(threshold)) threshold <- 0L
    
    samp_counts <- tapply(grp_meta$.cell, grp_meta[[sample_col]], length)
    valid_samps <- names(samp_counts[samp_counts >= threshold])
    excluded    <- names(samp_counts[samp_counts <  threshold])
    
    if (length(excluded) > 0)
      cat(sprintf("  [%s] Excluding samples (< %d cells): %s\n",
                  grp, threshold, paste(excluded, collapse = ", ")))
    
    if (length(valid_samps) == 0) {
      cat(sprintf("  [%s] No qualifying samples — group dropped\n", grp))
      next
    }
    
    target <- min(samp_counts[valid_samps])
    cat(sprintf("  [%s] %d samples, %d cells each -> %d total\n",
                grp, length(valid_samps), target,
                target * length(valid_samps)))
    
    kept <- unlist(lapply(valid_samps, function(s) {
      cells <- grp_meta$.cell[grp_meta[[sample_col]] == s]
      if (length(cells) <= target) cells else sample(cells, target)
    }))
    keep_all <- c(keep_all, kept)
  }
  
  cat(sprintf("  Downsampled: %d -> %d cells\n", nrow(meta), length(keep_all)))
  subset(obj, cells = keep_all)
}


#' Run the full CellChat v2 pipeline on a Seurat object.
#'
#' Groups with fewer than min.cells nuclei are dropped before inference.
#' Returns NULL if fewer than two valid groups remain.
#'
#' @param pop.size Passed to computeCommunProb. Use TRUE for single-condition
#'   analyses (corrects for group abundance); FALSE for cross-condition
#'   comparisons (removes abundance as a confound).
run_cellchat <- function(obj, label = "", min.cells = MIN_CELLS,
                         pop.size = FALSE) {
  gc_tab <- table(obj$cell_group)
  valid  <- names(gc_tab[gc_tab >= min.cells])
  
  cat(sprintf("\n[CellChat] %s | %d cells | groups: %s\n",
              label, ncol(obj), paste(valid, collapse = ", ")))
  
  if (length(valid) < 2) { cat("  Skipping: fewer than 2 valid groups.\n"); return(NULL) }
  
  dropped <- setdiff(names(gc_tab), valid)
  if (length(dropped))
    cat(sprintf("  Dropping (< %d cells): %s\n",
                min.cells, paste(dropped, collapse = ", ")))
  
  obj <- subset(obj,
                cells = rownames(obj@meta.data)[obj$cell_group %in% valid])
  
  cc <- createCellChat(object   = get_norm(obj),
                       meta     = obj@meta.data,
                       group.by = "cell_group")
  cc@DB <- CellChatDB.human
  
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  
  cc <- computeCommunProb(cc,
                          type            = "truncatedMean",
                          trim            = 0.1,
                          nboot           = 100,
                          seed.use        = 42,
                          population.size = pop.size)
  cc <- filterCommunication(cc, min.cells = min.cells)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
  
  cat(sprintf("  -> %d interactions | %d pathways\n",
              sum(cc@net$count, na.rm = TRUE),
              length(cc@netP$pathways)))
  cc
}


# ── 6.  CELLCHAT RUNS ────────────────────────────────────────────────────────

# ── 6a. Pan-cancer ───────────────────────────────────────────────────────────
# Restricted to astrocyte-qualifying samples; no downsampling applied
# (three samples represent distinct subtypes, not biological replicates).
cat("\n=== Pan-cancer CellChat ===\n")
cat("Samples:", paste(ASTRO_SAMPLES, collapse = ", "), "\n")

cc_pan <- run_cellchat(so_sub_astro, label = "Pan-cancer", pop.size = TRUE)
stopifnot(!is.null(cc_pan))
saveRDS(cc_pan, file.path(RDS_DIR, "cellchat_pan.rds"))

pan_comm <- subsetCommunication(cc_pan)
write.csv(pan_comm, file.path(TAB_DIR, "pan_interactions.csv"), row.names = FALSE)


# ── 6b. Per-subtype ──────────────────────────────────────────────────────────
# Astrocytes: qualifying sample only (one per subtype).
# Tumor + Microglia: all samples within the subtype, threshold-filtered
# and downsampled to equalise donor contributions.
cat("\n=== Per-subtype CellChat ===\n")
cc_sub_list <- list()

for (st in SUBTYPES) {
  st_short <- sub("ATRT-", "", st)
  
  cells_tm <- rownames(so_sub@meta.data)[
    so_sub@meta.data$subtype == st &
      so_sub@meta.data$cell_group %in% c("Tumor","Microglia")]
  
  cells_astro <- rownames(so_sub@meta.data)[
    so_sub@meta.data$subtype == st &
      so_sub@meta.data$cell_group == "Astrocytes" &
      so_sub@meta.data$ID %in% ASTRO_SAMPLES]
  
  if (length(cells_astro) == 0) {
    cat(sprintf("  Skipping %s: no qualifying Astrocyte sample\n", st)); next
  }
  
  obj_st    <- subset(so_sub, cells = unique(c(cells_tm, cells_astro)))
  cat(sprintf("\n--- %s ---\n", st))
  cat(sprintf("  Astrocytes from: %s\n",
              paste(unique(so_sub@meta.data[cells_astro, "ID"]), collapse = ", ")))
  
  obj_st_ds <- downsample_per_sample(obj_st)
  rm(obj_st)
  
  cc_st <- run_cellchat(obj_st_ds, label = st, min.cells = MIN_CELLS,
                        pop.size = FALSE)
  rm(obj_st_ds); gc()
  if (is.null(cc_st)) next
  
  saveRDS(cc_st, file.path(RDS_DIR, paste0("cellchat_", st_short, ".rds")))
  comm_st         <- subsetCommunication(cc_st)
  comm_st$subtype <- st_short
  write.csv(comm_st,
            file.path(TAB_DIR, paste0("interactions_", st_short, ".csv")),
            row.names = FALSE)
  cc_sub_list[[st_short]] <- cc_st
}
cat(sprintf("\nSubtypes completed: %s\n", paste(names(cc_sub_list), collapse = ", ")))


# ── 6c. Merge subtype objects ─────────────────────────────────────────────────
cc_merged <- NULL
if (length(cc_sub_list) >= 2) {
  ref_groups <- levels(cc_pan@idents)
  lifted <- lapply(names(cc_sub_list), function(st) {
    cc      <- cc_sub_list[[st]]
    missing <- setdiff(ref_groups, levels(cc@idents))
    if (length(missing))
      cat(sprintf("  liftCellChat %s — adding: %s\n",
                  st, paste(missing, collapse = ", ")))
    if (length(missing)) liftCellChat(cc, group.new = ref_groups) else cc
  })
  names(lifted) <- names(cc_sub_list)
  cc_merged <- mergeCellChat(lifted, add.names = names(lifted))
  saveRDS(cc_merged, file.path(RDS_DIR, "cellchat_merged.rds"))
  cat("Merged CellChat object saved.\n")
}


# ══════════════════════════════════════════════════════════════════════════════
# ── 7.  PUBLICATION FIGURES ───────────────────────────────────────────────────
# ══════════════════════════════════════════════════════════════════════════════
cat("\n=== Generating figures ===\n")


# ── Fig 1: Cell composition ───────────────────────────────────────────────────
comp <- so_sub@meta.data %>%
  count(cell_group, subtype) %>%
  group_by(subtype) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

p1a <- ggplot(comp, aes(x = factor(sub("ATRT-","", subtype), ST_SHORT),
                        y = pct, fill = cell_group)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = POP_COLORS, name = "Population") +
  scale_y_continuous(labels = \(x) paste0(x, "%"), expand = c(0, 0)) +
  labs(title    = "Cell composition by ATRT subtype",
       subtitle = paste("Excluded:", paste(EXCLUDE_SAMPLES, collapse = ", ")),
       x = NULL, y = "% cells") +
  pub_theme
save_fig(p1a, "fig_01a_composition", w = 5, h = 5)

heat_cnt <- so_sub@meta.data %>%
  count(ID, subtype, cell_group) %>%
  mutate(ID = factor(ID, levels = unique(ID[order(subtype, ID)])))

p1b <- ggplot(heat_cnt, aes(x = cell_group, y = ID, fill = n)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = n, colour = n > 300), size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "black","TRUE" = "white")) +
  scale_fill_viridis_c(option = "D", trans = "sqrt", name = "n (sq rt)") +
  facet_grid(subtype ~ ., scales = "free_y", space = "free_y") +
  labs(title    = "Cells per sample",
       subtitle = "Bold = qualifying sample for Astrocyte analyses",
       x = NULL, y = NULL) +
  pub_theme +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x  = element_text(angle = 30, hjust = 1),
        axis.text.y  = element_text(
          face = ifelse(levels(heat_cnt$ID) %in% ASTRO_SAMPLES,
                        "bold", "plain")))
save_fig(p1b, "fig_01b_counts_heatmap", w = 5, h = 8)
cat("  Fig 1 done\n")


# ── Fig 2: Pan-cancer interaction network ────────────────────────────────────
grp_sz <- as.numeric(table(cc_pan@idents))

pdf(file.path(FIG_DIR, "fig_02a_pan_circle_count.pdf"), w = 6, h = 6)
par(mar = c(0, 0, 2.5, 0))
netVisual_circle(cc_pan@net$count, vertex.weight = grp_sz,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = POP_COLORS,
                 title.name = "Pan-cancer — number of interactions")
dev.off()

pdf(file.path(FIG_DIR, "fig_02b_pan_circle_weight.pdf"), w = 6, h = 6)
par(mar = c(0, 0, 2.5, 0))
netVisual_circle(cc_pan@net$weight, vertex.weight = grp_sz,
                 weight.scale = TRUE, label.edge = FALSE,
                 color.use = POP_COLORS,
                 title.name = "Pan-cancer — interaction strength")
dev.off()

p2c <- netVisual_heatmap(cc_pan, measure = "count",
                         color.use = POP_COLORS, title.name = "# interactions")
p2d <- netVisual_heatmap(cc_pan, measure = "weight",
                         color.use = POP_COLORS, title.name = "Interaction strength")
pdf(file.path(FIG_DIR, "fig_02c_pan_heatmap_count.pdf"),  w = 4.5, h = 4)
print(p2c); dev.off()
pdf(file.path(FIG_DIR, "fig_02d_pan_heatmap_weight.pdf"), w = 4.5, h = 4)
print(p2d); dev.off()
cat("  Fig 2 done\n")


# ── Fig 3: Pan-cancer L-R bubble plots ───────────────────────────────────────
pan_pairs <- list(
  list(src = "Tumor",      tgt = "Astrocytes", nm = "Tumor_Astrocytes"),
  list(src = "Tumor",      tgt = "Microglia",  nm = "Tumor_Microglia"),
  list(src = "Astrocytes", tgt = "Microglia",  nm = "Astrocytes_Microglia")
)
for (pp in pan_pairs) {
  p <- tryCatch(
    netVisual_bubble(cc_pan,
                     sources.use    = pp$src,
                     targets.use    = pp$tgt,
                     bidirection    = TRUE,
                     remove.isolate = TRUE,
                     font.size      = 9,
                     title.name     = paste0("Pan-cancer | ",
                                             gsub("_", " <-> ", pp$nm))) +
      pub_theme +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)),
    error = \(e) { cat("  [warn]", pp$nm, ":", e$message, "\n"); NULL }
  )
  if (!is.null(p))
    save_fig(p, paste0("fig_03_pan_bubble_", pp$nm), w = 11, h = 9)
}
cat("  Fig 3 done\n")


# ── Fig 4: Signaling roles and top pathways ───────────────────────────────────
p4a <- tryCatch(
  netAnalysis_signalingRole_scatter(cc_pan, color.use = POP_COLORS,
                                    title = "Signaling roles — Pan-cancer") +
    pub_theme,
  error = \(e) NULL)
if (!is.null(p4a)) save_fig(p4a, "fig_04a_pan_role_scatter", w = 6, h = 5)

pdf(file.path(FIG_DIR, "fig_04b_pan_role_heatmap_out.pdf"), w = 8, h = 5)
tryCatch(print(netAnalysis_signalingRole_heatmap(
  cc_pan, pattern = "outgoing", color.use = POP_COLORS,
  title = "Outgoing signaling")), error = \(e) NULL)
dev.off()

pdf(file.path(FIG_DIR, "fig_04c_pan_role_heatmap_in.pdf"), w = 8, h = 5)
tryCatch(print(netAnalysis_signalingRole_heatmap(
  cc_pan, pattern = "incoming", color.use = POP_COLORS,
  title = "Incoming signaling")), error = \(e) NULL)
dev.off()

# Focal interactions: all pairs between the three populations
focal_pan <- pan_comm %>%
  filter(
    (source == "Tumor"      & target %in% c("Astrocytes","Microglia")) |
      (source == "Astrocytes" & target %in% c("Tumor","Microglia"))      |
      (source == "Microglia"  & target %in% c("Tumor","Astrocytes"))
  )
write.csv(focal_pan, file.path(TAB_DIR, "pan_focal_interactions.csv"),
          row.names = FALSE)

top_pw_df <- focal_pan %>%
  mutate(pair = paste0(source, " -> ", target)) %>%
  group_by(pair, pathway_name) %>%
  summarise(total = sum(prob), .groups = "drop") %>%
  group_by(pair) %>%
  slice_max(total, n = 8, with_ties = FALSE) %>% ungroup()

p4d <- ggplot(top_pw_df, aes(x = reorder(pathway_name, total),
                             y = total, fill = pair)) +
  geom_col(width = 0.7, show.legend = FALSE) +
  facet_wrap(~pair, scales = "free", ncol = 2) +
  coord_flip() +
  scale_fill_brewer(palette = "Set1") +
  labs(title = "Top pathways per focal pair (pan-cancer)",
       x = NULL, y = "Summed probability") +
  pub_theme
save_fig(p4d, "fig_04d_pan_top_pathways", w = 11, h = 8)

top6_pw <- focal_pan %>%
  group_by(pathway_name) %>% summarise(t = sum(prob)) %>%
  slice_max(t, n = 6, with_ties = FALSE) %>% pull(pathway_name)

pdf(file.path(FIG_DIR, "fig_04e_pan_top_pathway_chords.pdf"), w = 14, h = 9)
par(mfrow = c(2, 3))
for (pw in top6_pw) {
  par(mar = c(0, 0, 3, 0))
  tryCatch(netVisual_aggregate(cc_pan, signaling = pw, layout = "chord",
                               color.use = POP_COLORS,
                               vertex.label.cex = 1, title.name = pw),
           error = \(e) NULL)
}
dev.off()
cat("  Fig 4 done\n")


# ── Fig 5: Per-subtype circle plots ──────────────────────────────────────────
pdf(file.path(FIG_DIR, "fig_05_subtype_circles.pdf"),
    w = 5.5 * length(cc_sub_list), h = 5.5)
par(mfrow = c(1, length(cc_sub_list)))
for (st in names(cc_sub_list)) {
  cc_st <- cc_sub_list[[st]]
  grp_s <- as.numeric(table(cc_st@idents))
  par(mar = c(0, 0, 3, 0))
  netVisual_circle(cc_st@net$count,
                   vertex.weight = grp_s, weight.scale = TRUE,
                   label.edge    = FALSE,
                   color.use     = POP_COLORS[levels(cc_st@idents)],
                   title.name    = paste0("ATRT-", st))
}
dev.off()
cat("  Fig 5 done\n")


# ── Fig 6: Subtype comparison ─────────────────────────────────────────────────
if (!is.null(cc_merged)) {
  n_sub   <- length(cc_sub_list)
  sub_col <- unname(SUBTYPE_COLORS[paste0("ATRT-", names(cc_sub_list))])
  
  p6a <- tryCatch(
    compareInteractions(cc_merged, show.legend = FALSE,
                        group = seq_len(n_sub), color.use = sub_col) +
      pub_theme + labs(title = "# interactions by subtype"),
    error = \(e) NULL)
  if (!is.null(p6a)) save_fig(p6a, "fig_06a_compare_count", w = 5, h = 4)
  
  p6b <- tryCatch(
    compareInteractions(cc_merged, show.legend = FALSE, measure = "weight",
                        group = seq_len(n_sub), color.use = sub_col) +
      pub_theme + labs(title = "Interaction strength by subtype"),
    error = \(e) NULL)
  if (!is.null(p6b)) save_fig(p6b, "fig_06b_compare_weight", w = 5, h = 4)
  
  sub_pairs <- combn(seq_len(n_sub), 2, simplify = FALSE)
  sub_names <- names(cc_sub_list)
  for (sp in sub_pairs) {
    nm_pair <- paste0(sub_names[sp[1]], "_vs_", sub_names[sp[2]])
    pdf(file.path(FIG_DIR, paste0("fig_06c_diff_", nm_pair, ".pdf")), w = 5, h = 4)
    tryCatch(print(netVisual_diffInteraction(
      cc_merged, comparison = sp, color.use = POP_COLORS,
      title.name = nm_pair)), error = \(e) NULL)
    dev.off()
  }
  
  pdf(file.path(FIG_DIR, "fig_06d_pathway_rank.pdf"), w = 10, h = 6)
  tryCatch(rankNet(cc_merged, mode = "comparison", stacked = TRUE,
                   do.stat = TRUE, color.use = sub_col),
           error = \(e) NULL)
  dev.off()
  
  make_sub_bubble <- function(st, src, tgt) {
    cc <- cc_sub_list[[st]]
    if (!all(c(src, tgt) %in% levels(cc@idents))) return(NULL)
    tryCatch(
      netVisual_bubble(cc, sources.use = src, targets.use = tgt,
                       bidirection = TRUE, remove.isolate = TRUE,
                       font.size = 8, title.name = paste0("ATRT-", st)) +
        pub_theme +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)),
      error = \(e) NULL)
  }
  
  for (pair in list(
    list(src = "Tumor",      tgt = "Astrocytes", nm = "Tumor_Astrocytes"),
    list(src = "Tumor",      tgt = "Microglia",  nm = "Tumor_Microglia"),
    list(src = "Astrocytes", tgt = "Microglia",  nm = "Astrocytes_Microglia")
  )) {
    panels <- Filter(Negate(is.null),
                     lapply(names(cc_sub_list),
                            \(st) make_sub_bubble(st, pair$src, pair$tgt)))
    if (length(panels) >= 2) {
      p <- wrap_plots(panels, nrow = 1) +
        plot_annotation(
          title    = paste0(gsub("_", " <-> ", pair$nm), " — per subtype"),
          subtitle = "Astrocytes: qualifying sample only | Tumor + Microglia: all samples (balanced)",
          theme    = theme(
            plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
            plot.subtitle = element_text(hjust = 0.5, size = 9, colour = "grey40")))
      save_fig(p, paste0("fig_06e_sub_bubble_", pair$nm),
               w = 7 * length(panels), h = 9)
    }
  }
  cat("  Fig 6 done\n")
}


# ── Fig 7: Interaction presence dotplot ──────────────────────────────────────
all_sub_comm <- bind_rows(lapply(names(cc_sub_list), function(st) {
  df <- subsetCommunication(cc_sub_list[[st]])
  df$subtype <- st; df
}))

top50_pan <- focal_pan %>%
  slice_max(prob, n = 50, with_ties = FALSE) %>%
  pull(interaction_name) %>% unique()

presence <- all_sub_comm %>%
  filter(
    interaction_name %in% top50_pan,
    (source == "Tumor"      & target %in% c("Astrocytes","Microglia")) |
      (source == "Astrocytes" & target %in% c("Tumor","Microglia"))      |
      (source == "Microglia"  & target %in% c("Tumor","Astrocytes"))
  ) %>%
  mutate(direction = paste0(source, " -> ", target),
         lr_label  = paste0(interaction_name, "  [", pathway_name, "]"))

p7 <- ggplot(presence,
             aes(x = factor(subtype, ST_SHORT),
                 y = reorder(lr_label, prob),
                 colour = direction, size = prob)) +
  geom_point(alpha = 0.85) +
  scale_colour_brewer(palette = "Dark2", name = "Direction") +
  scale_size_continuous(range = c(1, 6), name = "Prob.") +
  facet_wrap(~direction, scales = "free_y", ncol = 3) +
  labs(title    = "Interaction presence across ATRT subtypes",
       subtitle = "Top-50 pan-cancer LR pairs | absent = not detected in that subtype",
       x = "Subtype", y = NULL) +
  pub_theme + theme(axis.text.y = element_text(size = 7))
save_fig(p7, "fig_07_presence_dotplot", w = 16, h = 14)
cat("  Fig 7 done\n")


# ── Fig 8: Signaling role scatter per subtype ─────────────────────────────────
role_panels <- Filter(Negate(is.null), lapply(names(cc_sub_list), function(st) {
  cc <- cc_sub_list[[st]]
  tryCatch(
    netAnalysis_signalingRole_scatter(
      cc, color.use = POP_COLORS[levels(cc@idents)],
      title = paste0("ATRT-", st)) +
      pub_theme + theme(legend.position = "none"),
    error = \(e) NULL)
}))
if (length(role_panels) >= 2) {
  p8 <- wrap_plots(role_panels, nrow = 1) +
    plot_annotation(
      title = "Signaling roles by subtype",
      theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13)))
  save_fig(p8, "fig_08_role_scatter_subtype",
           w = 5 * length(role_panels), h = 5)
}
cat("  Fig 8 done\n")


# ── Fig 9: Sample audit ───────────────────────────────────────────────────────
sample_audit <- so_sub@meta.data %>%
  count(ID, subtype, cell_group) %>%
  group_by(cell_group) %>%
  mutate(pct = n / sum(n) * 100) %>% ungroup()

p9 <- ggplot(mutate(sample_audit,
                    ID = factor(ID, levels = rev(sort(unique(ID))))),
             aes(x = cell_group, y = ID, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = n, colour = n > 300), size = 2.7, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "black","TRUE" = "white")) +
  scale_fill_viridis_c(option = "D", trans = "sqrt", name = "n (sq rt)") +
  facet_grid(subtype ~ ., scales = "free_y", space = "free_y") +
  labs(title    = "Cell counts per sample",
       subtitle = paste("Excluded:", paste(EXCLUDE_SAMPLES, collapse = ", ")),
       x = NULL, y = NULL) +
  pub_theme +
  theme(strip.text.y = element_text(angle = 0),
        axis.text.x  = element_text(angle = 30, hjust = 1))
save_fig(p9, "fig_09_sample_audit", w = 5.5, h = 8)
cat("  Fig 9 done\n")


# ── Fig 10: Fine-grained tumour cell states ───────────────────────────────────
so_sub_astro@meta.data$cell_group_fine <- dplyr::case_when(
  so_sub_astro@meta.data$Final_Annotation == "Astrocytes" ~ "Astrocytes",
  so_sub_astro@meta.data$Final_Annotation == "Microglia"  ~ "Microglia",
  so_sub_astro@meta.data$Final_Annotation %in% TUMOR_LABELS ~
    so_sub_astro@meta.data$Final_Annotation,
  TRUE ~ NA_character_
)
fg_keep  <- names(Filter(\(x) x >= MIN_CELLS * 2,
                         table(so_sub_astro@meta.data$cell_group_fine)))
fg_cells <- rownames(so_sub_astro@meta.data)[
  !is.na(so_sub_astro@meta.data$cell_group_fine) &
    so_sub_astro@meta.data$cell_group_fine %in% fg_keep]
so_fine  <- subset(so_sub_astro, cells = fg_cells)
so_fine@meta.data$cell_group <- so_fine@meta.data$cell_group_fine

cc_fine <- run_cellchat(so_fine, label = "Fine-grained tumour states",
                        min.cells = MIN_CELLS * 2)
rm(so_fine); gc()

if (!is.null(cc_fine)) {
  saveRDS(cc_fine, file.path(RDS_DIR, "cellchat_fine.rds"))
  t_grps <- sort(fg_keep[!fg_keep %in% c("Astrocytes","Microglia")])
  for (tgt in c("Astrocytes","Microglia")) {
    if (!tgt %in% levels(cc_fine@idents)) next
    p <- tryCatch(
      netVisual_bubble(cc_fine, sources.use = t_grps, targets.use = tgt,
                       remove.isolate = TRUE, font.size = 8,
                       title.name = paste0("Tumour states -> ", tgt)) +
        pub_theme,
      error = \(e) NULL)
    if (!is.null(p))
      save_fig(p, paste0("fig_10_fine_bubble_to_", tgt), w = 12, h = 9)
  }
  cat("  Fig 10 done\n")
}


# ── Top-10 L-R pairs per subtype (sender -> Tumor) ────────────────────────────
# Minimum group sizes for figure inclusion (applied post-CellChat)
MIN_SENDER_CELLS <- 5
MIN_TUMOR_CELLS  <- 20

sub_grp_n <- lapply(cc_sub_list, function(cc) {
  n <- as.integer(table(cc@idents)); names(n) <- levels(cc@idents); n
})

# Build a consistent pathway colour map across all figure pages
all_pathways <- sort(unique(unlist(lapply(names(cc_sub_list), function(st) {
  subsetCommunication(cc_sub_list[[st]]) %>%
    filter(ligand != receptor,
           (source %in% c("Astrocytes","Microglia") & target == "Tumor") |
             (source == "Tumor" & target %in% c("Astrocytes","Microglia"))) %>%
    pull(pathway_name) %>% unique()
}))))

kelly_cols <- c(
  "#F3C300","#875692","#F38400","#A1CAF1","#BE0032",
  "#C2B280","#848482","#008856","#E68FAC","#0067A5",
  "#F99379","#604E97","#F6A600","#B3446C","#DCD300",
  "#882D17","#8DB600","#654522","#E25822","#2B3D26",
  "#222222","#F2F3F4"
)
PATHWAY_COLORS <- if (length(all_pathways) <= length(kelly_cols)) {
  setNames(kelly_cols[seq_along(all_pathways)], all_pathways)
} else {
  setNames(colorRampPalette(brewer.pal(8,"Set1"))(length(all_pathways)),
           all_pathways)
}

for (sender in c("Astrocytes","Microglia")) {
  pdf_path <- file.path(FIG_DIR,
                        paste0("fig_top10_subtype_", sender, "_to_Tumor.pdf"))
  pdf(pdf_path, width = 7, height = 6)
  
  for (st in names(cc_sub_list)) {
    grp_n   <- sub_grp_n[[st]]
    n_send  <- grp_n[sender]
    n_tumor <- grp_n["Tumor"]
    
    if (is.na(n_send)  || n_send  < MIN_SENDER_CELLS) next
    if (is.na(n_tumor) || n_tumor < MIN_TUMOR_CELLS)  next
    
    df <- subsetCommunication(cc_sub_list[[st]]) %>%
      filter(source == sender, target == "Tumor", ligand != receptor) %>%
      slice_max(prob, n = 10, with_ties = FALSE) %>%
      mutate(lr_label = paste0(ligand, " -> ", receptor))
    
    if (nrow(df) == 0) next
    
    p <- ggplot(df, aes(x = prob, y = reorder(lr_label, prob),
                        fill = pathway_name)) +
      geom_segment(aes(x = 0, xend = prob,
                       y = reorder(lr_label, prob),
                       yend = reorder(lr_label, prob)),
                   colour = "grey80", linewidth = 0.4) +
      geom_point(shape = 21, colour = "white", stroke = 0.3, size = 5) +
      scale_fill_manual(values = PATHWAY_COLORS, name = "Pathway", drop = TRUE) +
      scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
      labs(title    = paste0("ATRT-", st, "  |  ", sender, " -> Tumor"),
           subtitle = sprintf(
             "Top 10 L-R pairs | homoprotein excluded | n(%s) = %d  n(Tumor) = %d",
             sender, n_send, n_tumor),
           x = "Communication probability", y = NULL) +
      pub_theme +
      theme(axis.text.y        = element_text(size = 10),
            legend.position    = "right",
            panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3))
    print(p)
  }
  dev.off()
  cat(sprintf("  Saved: %s\n", pdf_path))
}


# ── Microglia x tumour cell states — bidirectional heatmap ───────────────────
# All samples and all available cells are used here to maximise statistical
# power for the Wilcoxon-based over-expression test across fine-grained states.
cat("\n=== Microglia x Tumour states ===\n")

so_sub@meta.data$cell_group_states <- dplyr::case_when(
  so_sub@meta.data$Final_Annotation == "Microglia"    ~ "Microglia",
  so_sub@meta.data$Final_Annotation %in% TUMOR_LABELS ~ so_sub@meta.data$Final_Annotation,
  TRUE ~ NA_character_
)

so_states <- subset(so_sub,
                    cells = rownames(so_sub@meta.data)[
                      !is.na(so_sub@meta.data$cell_group_states)])
so_states@meta.data$cell_group <- so_states@meta.data$cell_group_states

cat("Cell counts per group:\n")
print(sort(table(so_states@meta.data$cell_group), decreasing = TRUE))

cc_states <- run_cellchat(so_states, label = "Microglia vs tumour states",
                          min.cells = MIN_CELLS)
rm(so_states); gc()
saveRDS(cc_states, file.path(RDS_DIR, "cellchat_microglia_tumor_states.rds"))

comm_states <- subsetCommunication(cc_states) %>%
  filter(source == "Microglia", target %in% TUMOR_LABELS, ligand != receptor) %>%
  mutate(lr_label = paste0(ligand, " -> ", receptor))
write.csv(comm_states,
          file.path(TAB_DIR, "microglia_tumor_states_interactions.csv"),
          row.names = FALSE)
cat(sprintf("Microglia -> tumour state interactions: %d\n", nrow(comm_states)))

state_order <- intersect(
  c("Cilia-like","CP-like","Hypoxic","IPC-like",
    "Mesenchymal-like","NPC-like","OPC-like","RG-like"),
  levels(cc_states@idents)
)

# Display pairs from the per-subtype top-10 Microglia->Tumor runs
top10_pairs <- unique(unlist(lapply(names(cc_sub_list), function(st) {
  subsetCommunication(cc_sub_list[[st]]) %>%
    filter(source == "Microglia", target == "Tumor", ligand != receptor) %>%
    slice_max(prob, n = 10, with_ties = FALSE) %>%
    mutate(lr_label = paste0(ligand, " -> ", receptor)) %>%
    pull(lr_label)
})))

comm_filtered        <- filter(comm_states, lr_label %in% top10_pairs)
tumor_states_present <- sort(unique(comm_states$target))

lr_order  <- comm_filtered %>%
  group_by(lr_label) %>%
  summarise(max_prob = max(prob), .groups = "drop") %>%
  arrange(desc(max_prob)) %>% pull(lr_label)
all_pairs <- union(lr_order, setdiff(top10_pairs, lr_order))

grid_df <- expand.grid(lr_label    = all_pairs,
                       tumor_state = tumor_states_present,
                       stringsAsFactors = FALSE) %>%
  left_join(select(comm_filtered, lr_label, target, prob, pathway_name),
            by = c("lr_label","tumor_state" = "target")) %>%
  replace_na(list(prob = 0)) %>%
  mutate(lr_label    = factor(lr_label,    levels = all_pairs),
         tumor_state = factor(tumor_state, levels = tumor_states_present))

p_heat <- ggplot(grid_df, aes(x = tumor_state, y = lr_label, fill = prob)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(data = filter(grid_df, prob > 0),
            aes(label = sprintf("%.2f", prob)),
            size = 2.5, colour = "white", fontface = "bold") +
  scale_fill_distiller(palette = "YlOrRd", direction = 1,
                       name = "Communication\nprobability", na.value = "grey97") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(title    = "Microglia -> Tumour cell states",
       subtitle = sprintf(
         "Per-subtype top-10 L-R pairs (%d shown) | grey = not detected",
         length(all_pairs)),
       x = NULL, y = NULL) +
  pub_theme +
  theme(axis.text.x  = element_text(angle = 35, hjust = 1, size = 9),
        axis.text.y  = element_text(size = 9),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5))

ggsave(file.path(FIG_DIR, "fig_heatmap_microglia_tumor_states.pdf"),
       p_heat,
       width     = 3 + length(tumor_states_present) * 1.1,
       height    = 3 + length(all_pairs) * 0.4,
       device    = cairo_pdf, limitsize = FALSE)
cat("  Saved: fig_heatmap_microglia_tumor_states.pdf\n")


# ── Master interaction table ──────────────────────────────────────────────────
n_sub_per_lr <- all_sub_comm %>%
  filter(
    (source == "Tumor"      & target %in% c("Astrocytes","Microglia")) |
      (source == "Astrocytes" & target %in% c("Tumor","Microglia"))      |
      (source == "Microglia"  & target %in% c("Tumor","Astrocytes"))
  ) %>%
  group_by(interaction_name, source, target, pathway_name) %>%
  summarise(n_subtypes = n_distinct(subtype),
            subtypes   = paste(sort(unique(subtype)), collapse = "/"),
            mean_prob  = mean(prob),
            max_prob   = max(prob),
            .groups    = "drop") %>%
  arrange(desc(n_subtypes), desc(mean_prob))
write.csv(n_sub_per_lr,
          file.path(TAB_DIR, "master_interactions_subtype_presence.csv"),
          row.names = FALSE)


# ── Session summary ───────────────────────────────────────────────────────────
cat("\n", strrep("=", 55), "\n", sep = "")
cat("COMPLETE\n")
cat(strrep("=", 55), "\n", sep = "")
cat("Figures : ", FIG_DIR, "\n")
cat("RDS     : ", RDS_DIR, "\n")
cat("Tables  : ", TAB_DIR, "\n")
cat(sprintf("\nPan-cancer : %d interactions | %d pathways\n",
            sum(cc_pan@net$count, na.rm = TRUE), length(cc_pan@netP$pathways)))
for (st in names(cc_sub_list))
  cat(sprintf("%-10s : %d interactions | %d pathways\n", st,
              sum(cc_sub_list[[st]]@net$count, na.rm = TRUE),
              length(cc_sub_list[[st]]@netP$pathways)))
cat(sprintf("\nPan-cancer (all 3 subtypes): %d interactions\n",
            sum(n_sub_per_lr$n_subtypes == 3)))
cat(sprintf("Subtype-specific (1 subtype): %d interactions\n",
            sum(n_sub_per_lr$n_subtypes == 1)))
cat("\nSession info:\n")
sessionInfo()