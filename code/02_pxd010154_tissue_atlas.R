#!/usr/bin/env Rscript
# =============================================================================
# R for Mass Spectrometry — Chapter: Human Tissue Proteome Atlas
# Dataset: PXD010154 (Human Proteome Atlas)
# 30 healthy human tissues, 17 individual organ proteomes, ultra-deep tonsil
# =============================================================================
#
# Pipeline sections:
#   1. Extract & load organ-specific proteome matrices
#   2. Load aggregate 30-tissue proteome
#   3. Tissue comparison — protein abundance heatmap
#   4. Organ similarity — PCA + hierarchical clustering
#   5. Tissue-enriched protein detection
#   6. Full proteome vs RNAseq-based proteome comparison (tonsil)
#   7. Summary tables
#
# Usage from cluster:
#   cd /cluster/data/ck_care/ltran/r4ms_book/analysis
#   Rscript tissue_atlas_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
})

# ── Paths ──────────────────────────────────────────────────────────────────
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  normalizePath(file.path(getwd(), "tissue_atlas_analysis.R"), mustWork = FALSE)
}

analysis_dir <- dirname(get_script_path())
project_dir  <- dirname(analysis_dir)
raw_dir       <- file.path(project_dir, "raw", "PXD010154")
extract_dir   <- file.path(raw_dir, "extracted")
out_fig       <- file.path(analysis_dir, "results", "figures", "tissue_atlas")
out_tbl       <- file.path(analysis_dir, "results", "tables", "tissue_atlas")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== Human Tissue Proteome Atlas — Analysis Pipeline ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ── 1. Extract organ-specific proteome matrices ────────────────────────────
cat("\n── 1. Extracting organ proteome matrices ──\n")

organ_zips <- list.files(raw_dir, pattern = "^P01[0-9]+_.*_fullproteome_RNAseq_txt\\.zip$",
                         full.names = TRUE)
organ_info <- tibble(
  zip_path = organ_zips,
  accession = str_extract(basename(organ_zips), "^P[0-9]+"),
  organ = str_replace(basename(organ_zips), "^P[0-9]+_(.+)_fullproteome_RNAseq_txt\\.zip$", "\\1")
)

cat(sprintf("  Found %d organ-specific ZIP files\n", nrow(organ_info)))

organ_matrices <- list()

for (i in seq_len(nrow(organ_info))) {
  acc   <- organ_info$accession[i]
  organ <- organ_info$organ[i]
  zipf  <- organ_info$zip_path[i]

  cat(sprintf("  [%d/%d] %s (%s)... ", i, nrow(organ_info), organ, acc))

  # Validate ZIP first
  contents <- tryCatch(unzip(zipf, list = TRUE), error = function(e) NULL)
  if (is.null(contents) || nrow(contents) == 0) {
    cat("corrupt ZIP, skipping\n")
    next
  }

  # Find the MaxQuant proteinGroups file
  pg_file <- grep("proteinGroups\\.txt$", contents$Name, value = TRUE)
  if (length(pg_file) == 0) {
    cat("no proteinGroups.txt found, skipping\n")
    next
  }

  # Extract and read
  tmp <- file.path(extract_dir, paste0(organ, "_", acc))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  unzip(zipf, files = pg_file[1], exdir = tmp, overwrite = FALSE)

  pg <- tryCatch({
    read.delim(file.path(tmp, pg_file[1]), check.names = FALSE, stringsAsFactors = FALSE)
  }, error = function(e) NULL)

  if (is.null(pg)) {
    cat("failed to read proteinGroups\n")
    next
  }

  # MaxQuant proteinGroups: columns include "Protein IDs", "Gene names",
  # "iBAQ", "iBAQ log10", "LFQ intensity", "Intensity", etc.
  # Use iBAQ if available, otherwise Intensity
  intensity_col <- grep("^iBAQ$", colnames(pg), value = TRUE)
  if (length(intensity_col) == 0) {
    intensity_col <- grep("^Intensity$", colnames(pg), value = TRUE)
  }
  if (length(intensity_col) == 0) {
    # Try LFQ intensity
    lfq_cols <- grep("^LFQ intensity", colnames(pg), value = TRUE)
    if (length(lfq_cols) > 0) {
      intensity_col <- lfq_cols[1]
    }
  }
  if (length(intensity_col) == 0) {
    cat("no intensity column found\n")
    next
  }

  # Build expression vector: protein → intensity
  prot_ids <- pg[["Protein IDs"]] %||% pg[["Majority protein IDs"]] %||% rownames(pg)
  gene_names <- pg[["Gene names"]] %||% pg[["Gene names"]]

  expr_vec <- as.numeric(pg[[intensity_col[1]]])
  names(expr_vec) <- prot_ids

  # Remove contaminants, reverse hits, and zero/NA intensities
  valid_rows <- !grepl("^CON__|^REV__", prot_ids) &
                !is.na(expr_vec) & expr_vec > 0
  expr_vec <- expr_vec[valid_rows]

  # Log2 transform
  expr_vec <- log2(expr_vec)

  organ_matrices[[organ]] <- expr_vec
  cat(sprintf("%d proteins quantified (%s)\n", length(expr_vec), intensity_col[1]))
}

# ── 2. Build combined tissue expression matrix ─────────────────────────────
cat("\n── 2. Building combined tissue expression matrix ──\n")

# Each organ_matrices entry is a named vector: protein_id → log2(iBAQ)
# Get union of all protein IDs across organs
all_features <- unique(unlist(lapply(organ_matrices, names)))
cat(sprintf("  Union of protein IDs across all organs: %d\n", length(all_features)))

# Build matrix: rows = proteins, cols = organs
combined_matrix <- matrix(NA_real_, nrow = length(all_features),
                          ncol = length(organ_matrices))
rownames(combined_matrix) <- all_features
colnames(combined_matrix) <- names(organ_matrices)

for (organ in names(organ_matrices)) {
  vec <- organ_matrices[[organ]]
  combined_matrix[names(vec), organ] <- vec
}
cat(sprintf("  Combined matrix: %d proteins x %d organs\n",
            nrow(combined_matrix), ncol(combined_matrix)))

# ── 3. QC — feature coverage across organs ─────────────────────────────────
cat("\n── 3. Feature coverage across organs ──\n")
coverage <- rowSums(!is.na(combined_matrix))
cat(sprintf("  Features in ≥10 organs: %d\n", sum(coverage >= 10)))
cat(sprintf("  Features in all %d organs: %d\n", ncol(combined_matrix),
            sum(coverage == ncol(combined_matrix))))

# Filter to features present in at least 3 organs
combined_filt <- combined_matrix[coverage >= 3, ]
cat(sprintf("  Filtered matrix: %d features x %d organs\n",
            nrow(combined_filt), ncol(combined_filt)))

# ── 4. Tissue correlation heatmap ──────────────────────────────────────────
cat("\n── 4. Tissue correlation heatmap ──\n")

cor_matrix <- cor(combined_filt, use = "pairwise.complete.obs", method = "spearman")

pdf(file.path(out_fig, "01_tissue_correlation_heatmap.pdf"), width = 12, height = 10)
pheatmap(cor_matrix,
         main = "Tissue Proteome Correlation (Spearman)",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize_number = 7,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
         breaks = seq(-1, 1, length.out = 101),
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
dev.off()
cat("  ✓ Saved: 01_tissue_correlation_heatmap.pdf\n")

# ── 5. PCA of tissue proteomes ─────────────────────────────────────────────
cat("\n── 5. PCA of tissue proteomes ──\n")

pca_input <- t(combined_filt)
pca_input <- pca_input[, colSums(is.na(pca_input)) == 0]  # complete cases only
pca_res <- prcomp(pca_input, scale. = TRUE, center = TRUE)
pca_scores <- as_tibble(pca_res$x, rownames = "organ")
pca_var <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

pdf(file.path(out_fig, "02_tissue_pca.pdf"), width = 10, height = 8)
p <- ggplot(pca_scores, aes(x = PC1, y = PC2, label = organ)) +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(size = 3, max.overlaps = 20) +
  labs(title = "PCA of Human Tissue Proteomes",
       x = sprintf("PC1 (%s%%)", pca_var[1]),
       y = sprintf("PC2 (%s%%)", pca_var[2])) +
  theme_minimal(base_size = 12)
print(p)
dev.off()
cat("  ✓ Saved: 02_tissue_pca.pdf\n")

# ── 6. Tissue-enriched protein detection ───────────────────────────────────
cat("\n── 6. Tissue-enriched protein detection ──\n")

# For each organ, find proteins with highest relative expression
enrichment_results <- list()

for (organ in colnames(combined_filt)) {
  organ_expr <- combined_filt[, organ]
  other_means <- rowMeans(combined_filt[, setdiff(colnames(combined_filt), organ),
                                        drop = FALSE], na.rm = TRUE)

  # Fold change vs mean of other tissues
  fc <- organ_expr - other_means  # log2 scale

  enrich <- tibble(
    feature_id = names(fc),
    organ = organ,
    log2_expr = organ_expr,
    log2_other_mean = other_means,
    log2FC = fc
  ) %>%
    filter(!is.na(log2FC), is.finite(log2FC)) %>%
    arrange(desc(log2FC))

  enrichment_results[[organ]] <- enrich
}

all_enrichment <- bind_rows(enrichment_results)

# Top 5 enriched proteins per organ
top_per_organ <- all_enrichment %>%
  group_by(organ) %>%
  slice_max(log2FC, n = 5) %>%
  ungroup()

write_csv(top_per_organ, file.path(out_tbl, "tissue_enriched_proteins_top5.csv"))
cat(sprintf("  ✓ Saved: tissue_enriched_proteins_top5.csv (%d entries)\n",
            nrow(top_per_organ)))

# ── 7. Tissue-enriched protein heatmap (top 5 per organ) ───────────────────
cat("\n── 7. Tissue-enriched protein heatmap ──\n")

# Use top 5 enriched per organ, take union, keep proteins in ≥6 organs
top_features <- unique(top_per_organ$feature_id)
top_subset <- combined_filt[intersect(top_features, rownames(combined_filt)), ]
top_complete <- top_subset[rowSums(!is.na(top_subset)) >= ncol(top_subset)/2, , drop = FALSE]

if (nrow(top_complete) >= 5) {
  # Impute remaining NAs for clustering
  for (r in seq_len(nrow(top_complete))) {
    row_vals <- top_complete[r, ]
    nas <- is.na(row_vals)
    if (any(nas) && !all(nas)) row_vals[nas] <- min(row_vals, na.rm = TRUE)
    top_complete[r, ] <- row_vals
  }

  pdf(file.path(out_fig, "03_tissue_enriched_heatmap.pdf"), width = 14, height = 10)
  pheatmap(top_complete,
           main = sprintf("Tissue-Enriched Proteins (%d)", nrow(top_complete)),
           scale = "row",
           fontsize_row = 6,
           fontsize_col = 10,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation")
  dev.off()
  cat(sprintf("  ✓ Saved: 03_tissue_enriched_heatmap.pdf (%d proteins)\n",
              nrow(top_complete)))
} else {
  cat("  Too few complete proteins for heatmap, skipping\n")
}

# ── 8. Organ feature count summary ─────────────────────────────────────────
cat("\n── 8. Organ summary ──\n")

organ_summary <- tibble(
  organ = colnames(combined_filt),
  n_proteins = colSums(!is.na(combined_filt)),
  mean_expr = colMeans(combined_filt, na.rm = TRUE),
  median_expr = apply(combined_filt, 2, median, na.rm = TRUE)
) %>% arrange(desc(n_proteins))

write_csv(organ_summary, file.path(out_tbl, "organ_summary.csv"))
cat(sprintf("  ✓ Saved: organ_summary.csv (%d organs)\n", nrow(organ_summary)))
print(organ_summary, n = Inf)

# ── Done ────────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Tissue Atlas Analysis complete: %s ===\n", Sys.time()))
cat(sprintf("  Figures: %s/\n", out_fig))
cat(sprintf("  Tables:  %s/\n", out_tbl))
