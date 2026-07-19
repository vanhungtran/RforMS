#!/usr/bin/env Rscript
# =============================================================================
# R for Mass Spectrometry — Chapter: Clinical DIA Proteomics with Spectronaut
# Dataset: PXD000547 (DIA/SWATH, paired clinical samples)
# =============================================================================
#
# Pipeline sections:
#   1. Import Spectronaut .sepr export files
#   2. QC — density plots, missing values, sample correlation
#   3. Paired differential expression analysis (limma)
#   4. Volcano plot + protein heatmap
#   5. Clinical biomarker prioritization
#
# Usage from cluster:
#   cd /cluster/data/ck_care/ltran/r4ms_book/analysis
#   Rscript clinical_spectronaut_analysis.R
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
  library(limma)
})

# ── Paths ──────────────────────────────────────────────────────────────────
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) return(normalizePath(sub("^--file=", "", file_arg[1])))
  normalizePath(file.path(getwd(), "clinical_spectronaut_analysis.R"), mustWork = FALSE)
}

analysis_dir <- dirname(get_script_path())
project_dir  <- dirname(analysis_dir)
raw_dir       <- file.path(project_dir, "raw", "PXD000547")
out_fig       <- file.path(analysis_dir, "results", "figures", "clinical")
out_tbl       <- file.path(analysis_dir, "results", "tables", "clinical")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_tbl, recursive = TRUE, showWarnings = FALSE)

cat("=== Clinical DIA Proteomics — PXD000547 ===\n")
cat(sprintf("Started: %s\n", Sys.time()))

# ── 1. Import Spectronaut .sepr files ──────────────────────────────────────
cat("\n── 1. Importing Spectronaut export files ──\n")

sepr_files <- list.files(raw_dir, pattern = "\\.sepr$", full.names = TRUE)
cat(sprintf("  Found %d .sepr files\n", length(sepr_files)))

# Spectronaut .sepr format: tab-delimited, rows = peptides, columns = metadata + samples
load_sepr <- function(path) {
  cat(sprintf("  Loading: %s ...\n", basename(path)))

  # Read header to understand structure
  hdr <- read.delim(path, nrows = 1, header = FALSE, stringsAsFactors = FALSE,
                    check.names = FALSE)
  # Read full file
  raw <- read.delim(path, header = TRUE, stringsAsFactors = FALSE,
                    check.names = FALSE, fill = TRUE)

  cat(sprintf("    %d rows x %d cols\n", nrow(raw), ncol(raw)))

  # Identify key columns
  cols <- colnames(raw)

  # PG.ProteinAccessions = protein ID
  # EG.Qvalue = q-value at protein group level
  # EG.PGQuantity = protein group quantity
  # Any column with "PGQuantity" and sample name = abundance

  protein_col <- grep("PG\\.ProteinAccessions", cols, value = TRUE, ignore.case = TRUE)[1]
  qval_col    <- grep("EG\\.Qvalue", cols, value = TRUE, ignore.case = TRUE)[1]
  quantity_cols <- grep("PG\\.Quantity", cols, value = TRUE, ignore.case = TRUE)

  if (is.na(protein_col) || length(quantity_cols) == 0) {
    cat("    WARNING: Could not find expected Spectronaut columns\n")
    return(NULL)
  }

  # Build expression matrix from quantity columns
  expr_mat <- raw %>%
    select(all_of(quantity_cols)) %>%
    mutate(across(everything(), ~ suppressWarnings(as.numeric(.)))) %>%
    as.matrix()
  rownames(expr_mat) <- raw[[protein_col]]

  # Clean sample names (remove "PG.Quantity." or similar prefixes)
  colnames(expr_mat) <- str_replace(colnames(expr_mat),
                                     "^.*PG\\.Quantity\\.?", "")

  # Log2 transform
  expr_mat <- log2(expr_mat + 1)

  cat(sprintf("    Expression matrix: %d proteins x %d samples\n",
              nrow(expr_mat), ncol(expr_mat)))

  # Return list with expression + metadata
  list(
    expression = expr_mat,
    proteins   = raw[[protein_col]],
    qvalue     = if (!is.na(qval_col)) raw[[qval_col]] else rep(NA_real_, nrow(raw)),
    raw        = raw
  )
}

sepr_data <- lapply(sepr_files, load_sepr)
names(sepr_data) <- basename(sepr_files)
sepr_data <- sepr_data[!vapply(sepr_data, is.null, logical(1))]

if (length(sepr_data) == 0) {
  cat("ERROR: No valid Spectronaut files loaded. Check file format.\n")
  cat("Trying alternative: reading as generic tab-delimited...\n")

  # Fallback: try reading .sepr as generic tab-delimited with auto-detect
  for (f in sepr_files) {
    raw <- read.delim(f, header = TRUE, stringsAsFactors = FALSE,
                      check.names = FALSE, fill = TRUE, sep = "\t")
    cat(sprintf("  %s: %d cols\n", basename(f), ncol(raw)))
    cat(sprintf("  Columns: %s\n", paste(head(colnames(raw), 15), collapse = ", ")))
  }
  quit(status = 1)
}

# ── 2. Merge across files (if multiple) + QC ───────────────────────────────
cat("\n── 2. QC — Expression distributions ──\n")

# Combine if multiple files
if (length(sepr_data) > 1) {
  # Use common proteins
  common_prots <- Reduce(intersect, lapply(sepr_data, function(x) rownames(x$expression)))
  cat(sprintf("  Common proteins across %d files: %d\n",
              length(sepr_data), length(common_prots)))

  combined_expr <- do.call(cbind, lapply(sepr_data, function(x) {
    x$expression[common_prots, , drop = FALSE]
  }))
} else {
  combined_expr <- sepr_data[[1]]$expression
}

cat(sprintf("  Combined: %d proteins x %d samples\n",
            nrow(combined_expr), ncol(combined_expr)))

# Density before normalization
pdf(file.path(out_fig, "01_density_before_norm.pdf"), width = 8, height = 6)
plot(density(combined_expr[, 1], na.rm = TRUE), col = "grey80",
     main = "Expression Density — Raw", xlab = "log2(expression)", ylim = c(0, 0.5))
for (i in 2:ncol(combined_expr)) {
  lines(density(combined_expr[, i], na.rm = TRUE), col = "grey80")
}
lines(density(rowMeans(combined_expr, na.rm = TRUE), na.rm = TRUE), col = "red", lwd = 2)
legend("topright", c("Individual samples", "Mean"), col = c("grey80", "red"), lwd = c(1, 2))
dev.off()
cat("  ✓ Saved: 01_density_before_norm.pdf\n")

# ── 3. Median normalization ────────────────────────────────────────────────
cat("\n── 3. Median normalization ──\n")

medians <- apply(combined_expr, 2, median, na.rm = TRUE)
grand_median <- median(medians, na.rm = TRUE)
combined_norm <- sweep(combined_expr, 2, medians - grand_median)

# Density after normalization
pdf(file.path(out_fig, "02_density_after_norm.pdf"), width = 8, height = 6)
plot(density(combined_norm[, 1], na.rm = TRUE), col = "grey80",
     main = "Expression Density — After Median Normalization",
     xlab = "log2(expression)", ylim = c(0, 0.5))
for (i in 2:ncol(combined_norm)) {
  lines(density(combined_norm[, i], na.rm = TRUE), col = "grey80")
}
lines(density(rowMeans(combined_norm, na.rm = TRUE), na.rm = TRUE), col = "steelblue", lwd = 2)
legend("topright", c("Individual samples", "Mean"), col = c("grey80", "steelblue"), lwd = c(1, 2))
dev.off()
cat("  ✓ Saved: 02_density_after_norm.pdf\n")

# ── 4. Sample correlation ──────────────────────────────────────────────────
cat("\n── 4. Sample correlation ──\n")

# Filter to proteins with <50% missing
missing_pct <- rowMeans(is.na(combined_norm))
expr_filt <- combined_norm[missing_pct < 0.5, ]

cor_samples <- cor(expr_filt, use = "pairwise.complete.obs", method = "pearson")

pdf(file.path(out_fig, "03_sample_correlation.pdf"), width = 10, height = 8)
pheatmap(cor_samples,
         main = "Sample-Sample Correlation",
         display_numbers = TRUE,
         number_format = "%.3f",
         fontsize_number = 7,
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
dev.off()
cat("  ✓ Saved: 03_sample_correlation.pdf\n")

# ── 5. Missing value summary ───────────────────────────────────────────────
cat("\n── 5. Missing values ──\n")

na_counts <- colSums(is.na(combined_norm))
na_pct <- na_counts / nrow(combined_norm) * 100
cat(sprintf("  Missing per sample: %.1f–%.1f%%\n", min(na_pct), max(na_pct)))

pdf(file.path(out_fig, "04_missing_values.pdf"), width = 8, height = 6)
barplot(na_pct, names.arg = names(na_pct), las = 2, cex.names = 0.6,
        main = "Missing Values per Sample (%)", ylab = "% Missing",
        col = "steelblue")
abline(h = 50, lty = 2, col = "red")
dev.off()
cat("  ✓ Saved: 04_missing_values.pdf\n")

# ── 6. Paired differential analysis (if sample names indicate pairs) ───────
cat("\n── 6. Paired differential analysis ──\n")

# Detect paired design from sample names
sample_names <- colnames(combined_norm)
cat(sprintf("  Sample names: %s\n", paste(sample_names, collapse = ", ")))

# Try to infer grouping from sample name patterns
# Common patterns: CC_1, CC_2 (condition pairs); Case_1, Control_1, etc.
group_pattern <- str_extract(sample_names, "^[A-Za-z]+")
cat(sprintf("  Inferred groups: %s\n", paste(unique(group_pattern), collapse = ", ")))

# For now: if we have "CC_1" and "CC_2" style, treat CC as condition prefix
# and numeric suffix as pair ID
pair_id <- str_extract(sample_names, "[0-9]+$")

if (length(unique(group_pattern)) >= 2 && !any(is.na(pair_id))) {
  cat(sprintf("  Paired design detected: %d pairs\n", length(unique(pair_id))))

  # Build design matrix
  group <- factor(group_pattern)
  pair  <- factor(pair_id)

  design <- model.matrix(~ 0 + group + pair)
  colnames(design) <- make.names(colnames(design))

  # Impute missing values for DE (MinProb-like: replace NA with row minimum - 1)
  expr_imp <- combined_norm
  for (i in seq_len(nrow(expr_imp))) {
    row_vals <- expr_imp[i, ]
    if (any(is.na(row_vals))) {
      row_min <- min(row_vals, na.rm = TRUE)
      expr_imp[i, is.na(row_vals)] <- row_min - 1
    }
  }

  # limma paired analysis
  fit <- lmFit(expr_imp, design)
  contrast_formula <- paste(colnames(design)[1], colnames(design)[2], sep = " - ")
  contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2, trend = TRUE)

  de_results <- topTable(fit2, number = Inf, sort.by = "P")
  de_results$protein <- rownames(de_results)

  # Add significance flags
  de_results <- de_results %>%
    mutate(
      significance = case_when(
        adj.P.Val < 0.05 & logFC > 1  ~ "Up in GroupB",
        adj.P.Val < 0.05 & logFC < -1 ~ "Up in GroupA",
        TRUE ~ "NS"
      )
    )

  n_sig <- sum(de_results$adj.P.Val < 0.05, na.rm = TRUE)
  cat(sprintf("  Significant proteins (adj.P.Val < 0.05): %d\n", n_sig))
  cat(sprintf("  |log2FC| > 1 + adj.P.Val < 0.05: %d\n",
              sum(de_results$significance != "NS", na.rm = TRUE)))

  write_csv(de_results, file.path(out_tbl, "DE_paired_results.csv"))

  # ── Volcano plot ─────────────────────────────────────────────────────────
  pdf(file.path(out_fig, "05_volcano.pdf"), width = 10, height = 8)

  de_plot <- de_results %>%
    mutate(
      neg_log10_padj = -log10(adj.P.Val),
      label = if_else(adj.P.Val < 0.05 & abs(logFC) > 1.5, protein, "")
    )

  p <- ggplot(de_plot, aes(x = logFC, y = neg_log10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(aes(label = label), size = 3, max.overlaps = 30) +
    scale_color_manual(values = c("Up in GroupA" = "#2166AC", "Up in GroupB" = "#B2182B",
                                   "NS" = "grey70")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    labs(title = "Paired Differential Expression — PXD000547",
         x = "log2 Fold Change", y = "-log10(adjusted P-value)") +
    theme_minimal(base_size = 12)
  print(p)
  dev.off()
  cat("  ✓ Saved: 05_volcano.pdf\n")

  # ── Heatmap of top DE proteins ───────────────────────────────────────────
  sig_proteins <- de_results %>%
    filter(adj.P.Val < 0.05) %>%
    arrange(adj.P.Val)

  n_top <- min(50, nrow(sig_proteins))
  if (n_top >= 5) {
    top_expr <- expr_imp[sig_proteins$protein[1:n_top], ]
    top_expr <- top_expr - rowMeans(top_expr, na.rm = TRUE)  # center rows

    pdf(file.path(out_fig, "06_heatmap_top_de.pdf"), width = 10, height = 12)
    pheatmap(top_expr,
             main = sprintf("Top %d DE Proteins — PXD000547", n_top),
             scale = "none",
             fontsize_row = 6,
             color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
    dev.off()
    cat(sprintf("  ✓ Saved: 06_heatmap_top_de.pdf (%d proteins)\n", n_top))
  }

} else {
  cat("  Could not detect paired design. Computing basic group comparison...\n")

  # Fallback: simple two-group comparison if exactly 2 groups detected
  if (length(unique(group_pattern)) == 2) {
    groups <- unique(group_pattern)
    g1_cols <- which(group_pattern == groups[1])
    g2_cols <- which(group_pattern == groups[2])

    # Simple t-test per protein
    simple_de <- tibble(
      protein = rownames(combined_norm),
      mean_g1 = rowMeans(combined_norm[, g1_cols, drop = FALSE], na.rm = TRUE),
      mean_g2 = rowMeans(combined_norm[, g2_cols, drop = FALSE], na.rm = TRUE),
      log2FC = mean_g2 - mean_g1
    )

    write_csv(simple_de, file.path(out_tbl, "DE_simple_results.csv"))
    cat(sprintf("  ✓ Saved simple DE: %d proteins\n", nrow(simple_de)))
  }
}

# ── Done ────────────────────────────────────────────────────────────────────
cat(sprintf("\n=== Clinical Analysis complete: %s ===\n", Sys.time()))
cat(sprintf("  Figures: %s/\n", out_fig))
cat(sprintf("  Tables:  %s/\n", out_tbl))
