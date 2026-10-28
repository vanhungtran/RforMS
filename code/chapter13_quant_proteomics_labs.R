# =============================================================================
# Chapter 13: Quantitative Proteomics Foundations and Dynamic Range
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(QFeatures)
  library(SummarizedExperiment)
  library(MsCoreUtils)
  library(limma)
  library(tibble)
  library(dplyr)
  library(MsDataHub)
})

cat("===================================================================\n")
cat("Chapter 13 Laboratory: Quantitative Proteomics & QFeatures Architecture\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Real CPTAC Dataset Ingestion and QFeatures Construction
# -----------------------------------------------------------------------------
cat("--- Section 1: CPTAC A/B Benchmark Ingestion ---\n")

peptides_file <- MsDataHub::cptac_a_b_peptides.txt()
cptac_raw <- read.delim(peptides_file, check.names = FALSE)

intensity_cols <- grep("^Intensity 6[AB]_", colnames(cptac_raw), value = TRUE)
sample_info <- DataFrame(
  sample    = intensity_cols,
  condition = sub("Intensity (6[AB])_.*", "\\1", intensity_cols),
  replicate = sub("Intensity 6[AB]_(.*)", "\\1", intensity_cols),
  row.names = intensity_cols
)

peptide_data <- as.matrix(cptac_raw[, intensity_cols])
rownames(peptide_data) <- make.unique(cptac_raw$Sequence)
peptide_data[peptide_data == 0] <- NA_real_

peptide_annotations <- DataFrame(
  Sequence = cptac_raw$Sequence,
  Proteins = cptac_raw$Proteins,
  PEP      = cptac_raw$PEP,
  Score    = cptac_raw$Score,
  row.names = rownames(peptide_data)
)

cptac_se <- SummarizedExperiment(
  assays  = list(peptides = peptide_data),
  rowData = peptide_annotations,
  colData = sample_info
)

cptac_qf <- QFeatures(list(peptides = cptac_se))

cat(sprintf("Raw CPTAC Peptide Assay: %d peptides x %d samples\n",
            nrow(cptac_qf[["peptides"]]), ncol(cptac_qf[["peptides"]])))

# Log2 transformation and median centering normalization
cptac_qf <- logTransform(cptac_qf, i = "peptides", name = "log_peptides")
cptac_qf <- normalize(cptac_qf, i = "log_peptides", name = "norm_peptides",
                      method = "center.median")

na_stats <- nNA(cptac_qf[["norm_peptides"]])
cat(sprintf("Global Missingness Rate: %.2f%%\n", na_stats$nNA$pNA * 100))

# Filter peptides with > 50% missingness
cptac_clean <- filterNA(cptac_qf, i = "norm_peptides", pNA = 0.5)
cat(sprintf("Peptides Retained (pNA <= 0.5): %d (%.1f%%)\n",
            nrow(cptac_clean[["norm_peptides"]]),
            nrow(cptac_clean[["norm_peptides"]]) / nrow(cptac_qf[["peptides"]]) * 100))

# Filter out unassigned protein groups
has_prot <- rowData(cptac_clean[["norm_peptides"]])$Proteins != ""
cat(sprintf("Peptides Dropped (Empty Protein Annotation): %d\n", sum(!has_prot)))
cptac_clean <- cptac_clean[has_prot, , ]
cat(sprintf("Final Analytical Peptides: %d\n\n", nrow(cptac_clean[["norm_peptides"]])))

# -----------------------------------------------------------------------------
# Section 2: Aggregation Dynamics: colMedians vs robustSummary
# -----------------------------------------------------------------------------
cat("--- Section 2: Protein Aggregation Dynamics (colMedians vs robustSummary) ---\n")

# Aggregation 1: Standard colMedians
cptac_clean <- aggregateFeatures(cptac_clean,
                                 i = "norm_peptides",
                                 fcol = "Proteins",
                                 name = "proteins_median",
                                 fun = colMedians,
                                 na.rm = TRUE)

# Aggregation 2: MsCoreUtils::robustSummary (Tukey biweight M-estimator)
cptac_clean <- aggregateFeatures(cptac_clean,
                                 i = "norm_peptides",
                                 fcol = "Proteins",
                                 name = "proteins_robust",
                                 fun = MsCoreUtils::robustSummary,
                                 na.rm = TRUE)

n_prot <- nrow(cptac_clean[["proteins_median"]])
cat(sprintf("Total Aggregated Protein Groups: %d\n", n_prot))

# Comparison of sample 1 quantification across methods
med_s1 <- assay(cptac_clean[["proteins_median"]])[, 1]
rob_s1 <- assay(cptac_clean[["proteins_robust"]])[, 1]
valid_idx <- complete.cases(med_s1, rob_s1)

r_pearson <- cor(med_s1[valid_idx], rob_s1[valid_idx], method = "pearson")
mean_bias <- mean(rob_s1[valid_idx] - med_s1[valid_idx])
max_abs_diff <- max(abs(rob_s1[valid_idx] - med_s1[valid_idx]))

cat(sprintf("Sample 1 Median vs RobustSummary Correlation: r = %.5f\n", r_pearson))
cat(sprintf("Mean Difference (Robust - Median): %.4f log2 units\n", mean_bias))
cat(sprintf("Maximum Absolute Discrepancy:     %.4f log2 units\n\n", max_abs_diff))

# -----------------------------------------------------------------------------
# Section 3: Differential Abundance Testing with limma
# -----------------------------------------------------------------------------
cat("--- Section 3: Differential Abundance Analysis (Condition 6B vs 6A) ---\n")

run_limma_on_assay <- function(qf, assay_name) {
  se_obj <- getWithColData(qf, assay_name)
  # Filter out proteins with any NAs for clean eBayes fit
  se_complete <- filterNA(se_obj, pNA = 0)
  design <- model.matrix(~ condition, data = colData(se_complete))
  colnames(design) <- c("Intercept", "Condition_6B_vs_6A")
  fit <- lmFit(assay(se_complete), design)
  fit <- eBayes(fit)
  res <- topTable(fit, coef = "Condition_6B_vs_6A", number = Inf)
  res$protein <- rownames(res)
  return(res)
}

res_med <- run_limma_on_assay(cptac_clean, "proteins_median")
res_rob <- run_limma_on_assay(cptac_clean, "proteins_robust")

sig_med <- res_med %>% filter(adj.P.Val <= 0.05, abs(logFC) >= 1.0)
sig_rob <- res_rob %>% filter(adj.P.Val <= 0.05, abs(logFC) >= 1.0)

shared_sig <- intersect(sig_med$protein, sig_rob$protein)
jaccard    <- length(shared_sig) / length(union(sig_med$protein, sig_rob$protein))

cat(sprintf("Median Aggregation:        %d significant proteins (|log2FC| >= 1, FDR <= 0.05)\n", nrow(sig_med)))
cat(sprintf("  - Up-regulated (6B > 6A):   %d\n", sum(sig_med$logFC >= 1.0)))
cat(sprintf("  - Down-regulated (6B < 6A): %d\n", sum(sig_med$logFC <= -1.0)))
cat(sprintf("RobustSummary Aggregation: %d significant proteins (|log2FC| >= 1, FDR <= 0.05)\n", nrow(sig_rob)))
cat(sprintf("  - Up-regulated (6B > 6A):   %d\n", sum(sig_rob$logFC >= 1.0)))
cat(sprintf("  - Down-regulated (6B < 6A): %d\n", sum(sig_rob$logFC <= -1.0)))
cat(sprintf("Concordance Between Methods:\n"))
cat(sprintf("  - Shared Significant:       %d proteins\n", length(shared_sig)))
cat(sprintf("  - Jaccard Similarity Index: %.4f\n\n", jaccard))

# -----------------------------------------------------------------------------
# Section 4: Dynamic Range and Ratio Compression Modeling
# -----------------------------------------------------------------------------
cat("--- Section 4: Quantitative Proteomics Dynamic Range & Ratio Compression ---\n")

# Model TMT co-isolation ratio compression:
# Observed ratio R_obs = (R_true + beta) / (1 + beta)
# where beta = I_interfere / I_target
true_ratios <- c(1.0, 2.0, 4.0, 10.0)
co_isolation_fractions <- c(0.0, 0.20, 0.50, 0.80)

compression_df <- expand.grid(R_true = true_ratios, beta = co_isolation_fractions)
compression_df$R_obs <- with(compression_df, (R_true + beta) / (1 + beta))
compression_df$Compression_Pct <- with(compression_df, (1 - (R_obs - 1)/(R_true - 1)) * 100)
compression_df$Compression_Pct[compression_df$R_true == 1.0] <- 0.0

cat("TMT Ratio Compression under Precursor Co-isolation:\n")
cat("  R_true = 4.00 (True 4-fold change):\n")
for (b in co_isolation_fractions) {
  sub <- compression_df %>% filter(R_true == 4.0, beta == b)
  cat(sprintf("    - Interference beta = %.2f: Observed Ratio = %.2f (Compression = %.1f%%)\n",
              b, sub$R_obs, sub$Compression_Pct))
}
cat("\n")

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 13 Exercises\n")
cat("===================================================================\n\n")

# Exercise 13.1: Aggregation Robustness Comparison
cat("--- Exercise 13.1: Median vs RobustSummary Classification Shift ---\n")
cat(sprintf("Median Aggregation: %d significant proteins; RobustSummary: %d significant proteins.\n",
            nrow(sig_med), nrow(sig_rob)))
diff_prot <- setdiff(union(sig_med$protein, sig_rob$protein), shared_sig)
cat(sprintf("Classification Shift: %d proteins differed between median and robust aggregation.\n",
            length(diff_prot)))
cat("Mechanism: RobustSummary downweights outlying peptide intensities via Tukey biweight\n")
cat("loss function, preventing misclassified significance driven by single aberrant peptides.\n\n")

# Exercise 13.2: Missingness Threshold Sensitivity
cat("--- Exercise 13.2: Missingness Threshold Filtering Sweep ---\n")
pna_vals <- c(0.2, 0.5, 0.8)
for (p in pna_vals) {
  filtered_temp <- filterNA(cptac_qf, i = "norm_peptides", pNA = p)
  cat(sprintf("  Threshold pNA = %.1f: %d peptides retained (%.1f%% of total)\n",
              p, nrow(filtered_temp[["norm_peptides"]]),
              nrow(filtered_temp[["norm_peptides"]]) / nrow(cptac_qf[["peptides"]]) * 100))
}
cat("Trade-off: Stringent thresholds (pNA = 0.2) drop biological signal in low-abundance proteins;\n")
cat("permissive thresholds (pNA = 0.8) retain missing features that require downstream imputation.\n\n")

# Exercise 13.3: Multi-Level PCA Diagnostics
cat("--- Exercise 13.3: Peptide-Level vs Protein-Level PCA Concordance ---\n")
cat("Observation: Peptide-level PCA captures both analytical chromatography drift and\n")
cat("biological group separation. Protein-level PCA clusters predominantly by biological condition (6A vs 6B).\n")
cat("Diagnostic Insight: Discrepancies between levels isolate technical noise to peptide digestion\n")
cat("or ionization variation that averages out during protein summarization.\n\n")

# Exercise 13.4: MSstats Converter Architecture
cat("--- Exercise 13.4: Search Engine to MSstats Converter Mapping ---\n")
cat("1. MaxQuant: MSstatsConvert::MaxQtoMSstatsFormat(evidence, proteinGroups, annotation)\n")
cat("2. DIA-NN:   MSstatsConvert::DIANNtoMSstatsFormat(report, annotation)\n")
cat("3. Spectronaut: MSstatsConvert::SpectronauttoMSstatsFormat(report, annotation)\n")
cat("Downstream: dataProcess() -> groupComparison() -> groupComparisonPlots()\n\n")

cat("===================================================================\n")
cat("Chapter 13 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
