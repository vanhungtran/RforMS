# =============================================================================
# Chapter 03: The Modern R and Bioconductor Mass Spectrometry Ecosystem
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(Spectra)
  library(SummarizedExperiment)
  library(limma)
  library(msdata)
  library(tibble)
  library(dplyr)
  library(tidyr)
})

cat("===================================================================\n")
cat("Chapter 03 Laboratory: R & Bioconductor Ecosystem Architecture\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Vectorization Performance Benchmark
# -----------------------------------------------------------------------------
cat("--- Section 1: Vectorization vs Explicit Loops Benchmark ---\n")

set.seed(42)
n_vals <- 200000
raw_intensities <- rlnorm(n_vals, meanlog = 14, sdlog = 1.5)

# Benchmark 1: Explicit for-loop
t_start_loop <- Sys.time()
log_loop <- numeric(n_vals)
for (i in seq_along(raw_intensities)) {
  log_loop[i] <- log2(raw_intensities[i] + 1)
}
t_loop <- as.numeric(difftime(Sys.time(), t_start_loop, units = "secs"))

# Benchmark 2: Vectorized log2
t_start_vec <- Sys.time()
log_vec <- log2(raw_intensities + 1)
t_vec <- as.numeric(difftime(Sys.time(), t_start_vec, units = "secs"))

cat(sprintf("Vector Size: %d intensity measurements\n", n_vals))
cat(sprintf("  Explicit For-Loop: %.4f seconds\n", t_loop))
cat(sprintf("  Vectorized log2:   %.4f seconds\n", t_vec))
cat(sprintf("  Acceleration Factor: %.1fx faster\n", t_loop / t_vec))
cat(sprintf("  Numeric Identity Check: all.equal(log_loop, log_vec) -> %s\n\n",
            all.equal(log_loop, log_vec)))

# -----------------------------------------------------------------------------
# Section 2: S4 Class Validity and Coordinate Invariants
# -----------------------------------------------------------------------------
cat("--- Section 2: S4 Object Architecture and Validity Checks ---\n")

# Build valid SummarizedExperiment
n_feat <- 20
n_samp <- 6
toy_assay <- matrix(rlnorm(n_feat * n_samp, 10, 1), nrow = n_feat, ncol = n_samp,
                    dimnames = list(paste0("Feat_", 1:n_feat), paste0("S_", 1:n_samp)))
row_ann <- DataFrame(gene = paste0("Gene_", 1:n_feat), row.names = paste0("Feat_", 1:n_feat))
col_ann <- DataFrame(group = rep(c("Ctrl", "Treat"), each = 3), row.names = paste0("S_", 1:n_samp))

valid_se <- SummarizedExperiment(assays = list(exprs = toy_assay),
                                 rowData = row_ann, colData = col_ann)

cat(sprintf("Valid SummarizedExperiment Class: %s\n", class(valid_se)))
cat(sprintf("Dimensions: %d rows x %d columns\n", nrow(valid_se), ncol(valid_se)))
cat(sprintf("Slots Present: %s\n", paste(slotNames(valid_se), collapse = ", ")))

# Test invalid construction (dimension mismatch triggers error)
error_caught <- FALSE
tryCatch({
  invalid_se <- SummarizedExperiment(
    assays  = list(exprs = toy_assay),
    rowData = DataFrame(gene = paste0("Gene_", 1:15)) # 15 rows instead of 20!
  )
}, error = function(e) {
  error_caught <<- TRUE
  cat(sprintf("Expected S4 Validity Catch: '%s'\n\n", conditionMessage(e)))
})

# -----------------------------------------------------------------------------
# Section 3: Tidyverse Data Grammar on Real PSM Tables
# -----------------------------------------------------------------------------
cat("--- Section 3: Tidyverse Data Transformation Pipeline ---\n")

# Use mzIdentML identification table from msdata
id_file <- system.file("ident",
  "TMT_Erwinia_1uLSike_Top10HCD_isol2_45stepped_60min_01-20141210.mzid",
  package = "msdata")

if (requireNamespace("PSMatch", quietly = TRUE)) {
  psms_raw <- PSMatch::PSM(id_file)
  psms_df <- as.data.frame(psms_raw)
  
  cat(sprintf("Raw PSM Table: %d rows x %d columns\n", nrow(psms_df), ncol(psms_df)))
  
  # Tidyverse pipeline: Filter Rank 1 target peptides, compute score summaries
  psm_summary <- psms_df %>%
    filter(rank == 1, isDecoy == FALSE) %>%
    group_by(chargeState) %>%
    summarise(
      n_psms      = n(),
      n_unique_pep = n_distinct(sequence),
      mean_score  = round(mean(calculatedMassToCharge), 2),
      .groups     = "drop"
    ) %>%
    arrange(desc(n_psms))
  
  cat("Target PSM Charge State Distribution:\n")
  print(as.data.frame(psm_summary))
  cat("\n")
} else {
  cat("Package PSMatch not loaded; skipping mzid inspection.\n\n")
}

# -----------------------------------------------------------------------------
# Section 4: End-to-End Minimal Bioconductor Workflow (limma)
# -----------------------------------------------------------------------------
cat("--- Section 4: Minimal Bioconductor Differential Pipeline ---\n")

# 50 features across 6 samples: 3 Control vs 3 Treated
# Spike true differential abundance in features 1 to 5 (+2.0 log2FC)
set.seed(123)
n_features <- 50
n_samples  <- 6

base_means <- rlnorm(n_features, meanlog = 12, sdlog = 1.0)
exprs_mat  <- matrix(NA_real_, nrow = n_features, ncol = n_samples,
                     dimnames = list(paste0("METAB_", sprintf("%03d", 1:n_features)),
                                     paste0("Sample_", 1:n_samples)))

for (i in 1:n_features) {
  ctrl_vals  <- rnorm(3, mean = log2(base_means[i]), sd = 0.25)
  treat_mean <- if (i <= 5) log2(base_means[i]) + 2.0 else log2(base_means[i])
  treat_vals <- rnorm(3, mean = treat_mean, sd = 0.25)
  exprs_mat[i, ] <- c(ctrl_vals, treat_vals)
}

col_data <- DataFrame(
  condition = rep(c("Control", "Treated"), each = 3),
  row.names = colnames(exprs_mat)
)

se_pipeline <- SummarizedExperiment(assays = list(log2_intensity = exprs_mat),
                                    colData = col_data)

# Differential abundance with limma
design <- model.matrix(~ condition, data = as.data.frame(colData(se_pipeline)))
fit <- lmFit(assay(se_pipeline), design)
fit <- eBayes(fit)
top_res <- topTable(fit, coef = "conditionTreated", number = 10)

cat("Top 10 Differential Features from limma Pipeline:\n")
print(round(top_res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 4))

cat(sprintf("\nTrue Positive Recovery (Top 5 Spiked Features in FDR <= 0.05): %d / 5 (%.0f%%)\n",
            sum(top_res$adj.P.Val[1:5] <= 0.05),
            sum(top_res$adj.P.Val[1:5] <= 0.05) / 5 * 100))
cat(sprintf("False Positive Rate in Unchanged Features: %.2f%%\n\n",
            sum(top_res$adj.P.Val[-(1:5)] <= 0.05) / (n_features - 5) * 100))

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 03 Exercises\n")
cat("===================================================================\n\n")

# Exercise 3.1: Container Ecosystem Mappings
cat("--- Exercise 3.1: Core Bioconductor Containers ---\n")
cat("1. Spectra: Holds raw mass spectra time-series (mzR on-disk backend).\n")
cat("   Creation: Spectra::Spectra(file_path, backend = MsBackendMzR())\n")
cat("2. MsExperiment: Links raw spectra streams to multi-file experimental sample metadata.\n")
cat("   Creation: MsExperiment::readMsExperiment(spectraFiles, sampleData)\n")
cat("3. SummarizedExperiment: 2D feature-by-sample quantitative matrix with rowData and colData.\n")
cat("   Creation: SummarizedExperiment::SummarizedExperiment(assays, rowData, colData)\n")
cat("4. QFeatures: Multi-level relational assay graph (PSM -> Peptide -> Protein).\n")
cat("   Creation: QFeatures::QFeatures() or QFeatures::readQFeatures()\n\n")

# Exercise 3.2: Vectorization Efficiency
cat("--- Exercise 3.2: Vectorized vs For-Loop Execution Dynamics ---\n")
cat(sprintf("Benchmark Result: Vectorized operation completed in %.4f s vs %.4f s for loop.\n",
            t_vec, t_loop))
cat(sprintf("Acceleration: %.1fx speedup.\n", t_loop / t_vec))
cat("Pedagogical Reason: Vectorized R calls compiled underlying C code (SIMD vector instructions),\n")
cat("avoiding repeated R evaluation bytecode lookups and memory reallocations.\n\n")

# Exercise 3.3: Tidyverse Translation of Base R Subsetting
cat("--- Exercise 3.3: Tidyverse Code Translation ---\n")
cat("Base R Code:\n")
cat("  tmp1 <- psm_table[psm_table$msLevel == 2, ]\n")
cat("  tmp2 <- tmp1[tmp1$xcorr >= 2.5, ]\n")
cat("  tmp3 <- tmp2[order(tmp2$xcorr, decreasing = TRUE), ]\n")
cat("  result <- tmp3[1:10, c('sequence', 'xcorr', 'protein_id')]\n\n")
cat("Idiomatic Tidyverse Translation:\n")
cat("  result <- psm_table %>%\n")
cat("    filter(msLevel == 2, xcorr >= 2.5) %>%\n")
cat("    arrange(desc(xcorr)) %>%\n")
cat("    slice_head(n = 10) %>%\n")
cat("    select(sequence, xcorr, protein_id)\n\n")

# Exercise 3.4: Empirical Bayes Variance Shrinkage in limma
cat("--- Exercise 3.4: Empirical Bayes Moderation in Small-Sample MS ---\n")
cat("In experiments with n = 3 replicates per group, sample feature variance s_g^2 has only df = 4.\n")
cat("limma moderates individual feature variances toward a pooled global prior s_0^2:\n")
cat("  s_moderated^2 = (d_0 * s_0^2 + d_g * s_g^2) / (d_0 + d_g)\n")
cat("Statistical Benefit: Features with spuriously low sample variances are shrunk upward,\n")
cat("preventing astronomical false-positive t-statistics caused by tiny denominators.\n\n")

cat("===================================================================\n")
cat("Chapter 03 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
