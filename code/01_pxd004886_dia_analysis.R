#!/usr/bin/env Rscript
# =============================================================================
# R for Mass Spectrometry — Chapter: OpenSWATH DIA Pipeline
# Dataset: PXD004886 (Bruderer et al., Nat Commun 2017)
# 11-site inter-lab SWATH-MS reproducibility study
# =============================================================================
#
# Pipeline sections:
#   1. Load OpenSWATH peptide-level data
#   2. Filter (remove decoys, m_score threshold)
#   3. Log2 transform + median normalisation
#   4. Aggregate peptide → protein (median polish)
#   5. Missing value handling (MinProb imputation)
#   6. QC plots (density, PCA, correlation heatmap, missing values)
#   7. Cross-site differential abundance (limma)
#   8. Volcano plot + heatmap
#
# Usage from cluster:
#   cd /cluster/data/ck_care/ltran/r4ms_book/analysis
#   Rscript proteomics_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
  library(RColorBrewer)
  library(limma)
  library(matrixStats)
})

has_cowplot <- requireNamespace("cowplot", quietly = TRUE)
if (has_cowplot) library(cowplot)

# ── Path resolution (portable) ────────────────────────────────────────────────
get_script_path <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  if (length(file_arg) > 0) {
    return(normalizePath(sub("^--file=", "", file_arg[1])))
  }
  frame_files <- vapply(sys.frames(), function(env) {
    if (!is.null(env$ofile)) env$ofile else ""
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0) return(normalizePath(frame_files[1]))
  normalizePath(file.path(getwd(), "proteomics_analysis.R"), mustWork = FALSE)
}

analysis_dir <- dirname(get_script_path())
project_dir <- dirname(analysis_dir)
project_path <- function(...) file.path(project_dir, ...)
analysis_path <- function(...) file.path(analysis_dir, ...)

# ── Output directories ────────────────────────────────────────────────────────
dir.create(analysis_path("results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(analysis_path("results", "tables"),  recursive = TRUE, showWarnings = FALSE)

theme_base <- if (has_cowplot) theme_cowplot(12) else theme_bw(base_size = 12)

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 1: Load OpenSWATH peptide-level data
# ═══════════════════════════════════════════════════════════════════════════════
cat("=== 1. Loading data ===\n")

# Use a single site for demonstration (site01, 659 MB)
# For the full analysis, switch to all_sites_global_q_0.01_applied_to_local_global.txt
data_file <- project_path("raw", "PXD004886", "site01_global_q_0.01_applied_to_local_global.txt")

cat(sprintf("  Reading: %s\n", basename(data_file)))
raw <- read.delim(data_file, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, check.names = FALSE)

cat(sprintf("  Rows: %s  |  Columns: %d\n", format(nrow(raw), big.mark = ","), ncol(raw)))
cat(sprintf("  Unique runs: %d\n", length(unique(raw$filename))))
cat(sprintf("  Unique proteins: %d\n", length(unique(raw$ProteinName))))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 2: Filter
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 2. Filtering ===\n")

filtered <- raw %>%
  filter(
    decoy == 0,                                   # Remove decoy entries
    !is.na(Intensity), Intensity > 0,              # Remove missing/zero intensities
    transition_group_id_m_score < 0.01             # m_score threshold (FDR 1%)
  )

cat(sprintf("  After filter: %s rows (%d%% retained)\n",
    format(nrow(filtered), big.mark = ","),
    round(nrow(filtered) / nrow(raw) * 100)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 3: Aggregate peptide → protein (median per protein per run)
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 3. Peptide → protein aggregation ===\n")

# Simplify run names: extract core identifier from full path
filtered <- filtered %>%
  mutate(
    Run = basename(filename) %>% str_remove("\\.mzXML\\.gz$"),
    # Parse site info from the run name for metadata
    Site = str_extract(Run, "Site\\d+")
  )

# Median summarisation: for each protein, take median intensity across its peptides
prot_long <- filtered %>%
  group_by(ProteinName, Run, Site) %>%
  summarise(
    Intensity = median(Intensity, na.rm = TRUE),
    n_peptides = n(),
    .groups = "drop"
  )

# Pivot to wide matrix: proteins × runs
mat_wide <- prot_long %>%
  select(ProteinName, Run, Intensity) %>%
  pivot_wider(names_from = Run, values_from = Intensity) %>%
  column_to_rownames("ProteinName")

cat(sprintf("  Protein matrix: %d proteins × %d runs\n", nrow(mat_wide), ncol(mat_wide)))

# Build sample metadata
meta <- prot_long %>%
  distinct(Run, Site) %>%
  column_to_rownames("Run")

meta <- meta[colnames(mat_wide), , drop = FALSE]

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 4: Missing value filtering
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 4. Missing value filtering (≥50% valid per condition) ===\n")

# For demo: treat each dynamic range level as a pseudo-condition
# In real analysis, group by site or biological condition
conditions <- unique(meta$Site)
keep <- rep(FALSE, nrow(mat_wide))

for (cond in conditions) {
  cond_runs <- rownames(meta)[meta$Site == cond]
  if (length(cond_runs) > 0) {
    valid_frac <- rowSums(!is.na(mat_wide[, cond_runs, drop = FALSE])) / length(cond_runs)
    keep <- keep | (valid_frac >= 0.5)
  }
}

mat_filt <- mat_wide[keep, ]
cat(sprintf("  Proteins kept: %d / %d\n", nrow(mat_filt), nrow(mat_wide)))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 5: Log2 transform + median normalisation
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 5. Log2 transform + median normalisation ===\n")

mat_log <- log2(mat_filt)

# Median normalisation: center each column to the global median
sample_medians <- apply(mat_log, 2, median, na.rm = TRUE)
global_median  <- median(sample_medians)
mat_norm <- sweep(mat_log, 2, sample_medians - global_median, "-")

cat(sprintf("  Global median: %.2f\n", global_median))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 6: MinProb imputation
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 6. MinProb imputation ===\n")

set.seed(42)
mat_imp <- mat_norm
for (j in seq_len(ncol(mat_imp))) {
  col <- mat_imp[, j]
  miss <- is.na(col)
  if (any(miss)) {
    col_min <- quantile(col, 0.01, na.rm = TRUE)
    col_sd  <- sd(col, na.rm = TRUE) * 0.3
    if (is.na(col_sd) || col_sd == 0) col_sd <- 0.01
    mat_imp[miss, j] <- rnorm(sum(miss), mean = col_min, sd = col_sd)
  }
}

n_imputed <- sum(is.na(mat_norm))
cat(sprintf("  Values imputed: %s (%.1f%%)\n",
    format(n_imputed, big.mark = ","),
    n_imputed / (nrow(mat_norm) * ncol(mat_norm)) * 100))

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 7: QC Plots
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 7. QC Plots ===\n")

# 7a. Density: before vs after normalisation
plot_density <- function(mat, title_str) {
  df <- as.data.frame(mat) %>%
    mutate(Protein = rownames(.)) %>%
    pivot_longer(-Protein, names_to = "Run", values_to = "Intensity") %>%
    filter(!is.na(Intensity))
  ggplot(df, aes(x = Intensity, group = Run)) +
    geom_density(alpha = 0.6, linewidth = 0.3, colour = "steelblue") +
    labs(title = title_str, x = "log2(Intensity)", y = "Density") +
    theme_base
}

p1 <- plot_density(mat_log, "Before Normalisation")
p2 <- plot_density(mat_norm, "After Median Normalisation")

pdf(analysis_path("results", "figures", "01_density_before_norm.pdf"), 8, 5); print(p1); dev.off()
pdf(analysis_path("results", "figures", "01_density_after_norm.pdf"),  8, 5); print(p2); dev.off()
cat("  ✓ Density plots saved\n")

# 7b. Sample correlation heatmap
cor_mat <- cor(mat_imp, use = "pairwise.complete.obs", method = "pearson")
pdf(analysis_path("results", "figures", "02_sample_correlation_heatmap.pdf"), 10, 8)
pheatmap(cor_mat,
         clustering_method = "ward.D2",
         color = colorRampPalette(brewer.pal(9, "Blues"))(100),
         main = "Sample Correlation (Pearson)",
         fontsize = 7)
dev.off()
cat("  ✓ Correlation heatmap saved\n")

# 7c. PCA
pca <- prcomp(t(mat_imp), center = TRUE, scale. = TRUE)
pca_df <- as.data.frame(pca$x) %>%
  rownames_to_column("Run") %>%
  left_join(meta %>% rownames_to_column("Run"), by = "Run")

var_exp <- round(summary(pca)$importance[2, 1:2] * 100, 1)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, label = Run)) +
  geom_point(size = 2, colour = "steelblue") +
  geom_text_repel(size = 2.5, max.overlaps = 15) +
  labs(
    title = "PCA — Protein Expression",
    x = sprintf("PC1 (%.1f%%)", var_exp[1]),
    y = sprintf("PC2 (%.1f%%)", var_exp[2])
  ) + theme_base

pdf(analysis_path("results", "figures", "03_PCA.pdf"), 8, 6); print(p_pca); dev.off()
cat("  ✓ PCA saved\n")

# 7d. Missing values per sample
missing_df <- as.data.frame(mat_norm) %>%
  summarise(across(everything(), ~ sum(is.na(.)) / n() * 100)) %>%
  pivot_longer(everything(), names_to = "Run", values_to = "PctMissing")

p_miss <- ggplot(missing_df, aes(x = reorder(Run, PctMissing), y = PctMissing)) +
  geom_col(fill = "steelblue") +
  labs(title = "Missing Values per Run", x = "Run", y = "% Missing") +
  theme_base + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 6))

pdf(analysis_path("results", "figures", "04_missing_values.pdf"), 12, 5); print(p_miss); dev.off()
cat("  ✓ Missing value histogram saved\n")

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 8: Differential Abundance — Site Comparison
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 8. Cross-site differential abundance (limma) ===\n")

# For this demo, compare first half of runs vs second half as a proxy
# In a real analysis, you'd have biological conditions (e.g., treated vs control)
runs <- colnames(mat_imp)
n_runs <- length(runs)
group <- factor(c(rep("GroupA", n_runs %/% 2), rep("GroupB", n_runs - n_runs %/% 2)))

cat(sprintf("  Design: %d runs in GroupA, %d in GroupB\n", sum(group == "GroupA"), sum(group == "GroupB")))

design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "GroupB_vs_GroupA")

fit <- eBayes(lmFit(mat_imp, design), trend = TRUE, robust = TRUE)

de_results <- topTable(fit, coef = "GroupB_vs_GroupA", number = Inf, sort.by = "P") %>%
  rownames_to_column("Protein") %>%
  rename(log2FC = logFC, pval = P.Value, padj = adj.P.Val) %>%
  mutate(
    significance = case_when(
      padj < 0.05 & log2FC >  1 ~ "Up in GroupB",
      padj < 0.05 & log2FC < -1 ~ "Up in GroupA",
      TRUE ~ "NS"
    ),
    neg_log10_pval = -log10(pval)
  )

cat(sprintf("  DE proteins (padj<0.05, |log2FC|>1): %d\n",
    sum(de_results$significance != "NS")))

# Save DE table
write.csv(de_results, analysis_path("results", "tables", "DE_results.csv"), row.names = FALSE)
cat("  ✓ DE table saved\n")

# Volcano plot
p_volcano <- ggplot(de_results, aes(x = log2FC, y = neg_log10_pval, colour = significance)) +
  geom_point(alpha = 0.5, size = 0.8) +
  scale_colour_manual(values = c("Up in GroupA" = "#4DBBD5", "Up in GroupB" = "#E64B35", "NS" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.3) +
  labs(title = "Volcano Plot — GroupB vs GroupA",
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_base + theme(legend.position = "bottom")

pdf(analysis_path("results", "figures", "05_volcano.pdf"), 8, 6); print(p_volcano); dev.off()
cat("  ✓ Volcano plot saved\n")

# Heatmap: top 50 most variable proteins
top50 <- de_results %>%
  filter(significance != "NS") %>%
  arrange(pval) %>%
  head(50) %>%
  pull(Protein)

if (length(top50) >= 2) {
  mat_top <- mat_imp[top50, , drop = FALSE]
  # z-score per row
  mat_top_z <- t(scale(t(mat_top)))

  pdf(analysis_path("results", "figures", "06_heatmap_top50.pdf"), 10, 12)
  pheatmap(mat_top_z,
           clustering_method = "ward.D2",
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           main = "Top 50 DE Proteins (z-score)",
           fontsize_row = 6, fontsize_col = 6)
  dev.off()
  cat("  ✓ Heatmap saved\n")
}

# ═══════════════════════════════════════════════════════════════════════════════
# SECTION 9: Protein annotation
# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== 9. Protein annotation ===\n")

source(analysis_path("annotation_utils.R"))

ann <- annotate_proteins(de_results$Protein)

# Merge annotations into DE results
de_annotated <- de_results %>%
  left_join(ann %>% select(Protein, Genename, Protein.names), by = "Protein") %>%
  select(Protein, Genename, Protein.names, log2FC, pval, padj, significance, everything())

# Save annotated DE table
write.csv(de_annotated, analysis_path("results", "tables", "DE_results_annotated.csv"), row.names = FALSE)
cat("  ✓ Annotated DE table saved\n")

# Show top annotated hits
cat("\n  Top DE proteins with gene symbols:\n")
de_annotated %>%
  filter(significance != "NS") %>%
  arrange(pval) %>%
  select(Genename, Protein.names, log2FC, padj, significance) %>%
  head(10) %>%
  as.data.frame() %>%
  print()

# ═══════════════════════════════════════════════════════════════════════════════
cat("\n=== Pipeline complete ===\n")
cat(sprintf("  Figures: %s\n", analysis_path("results", "figures")))
cat(sprintf("  Tables:  %s\n", analysis_path("results", "tables")))
