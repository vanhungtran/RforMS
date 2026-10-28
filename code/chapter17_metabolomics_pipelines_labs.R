# =============================================================================
# Chapter 17: Untargeted Metabolomics Workflows with XCMS
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

Sys.setenv(RENV_CONFIG_NAMESPACES_CHECK = "FALSE")

suppressPackageStartupMessages({
  library(BiocParallel)
  register(SerialParam(), default = TRUE)
  library(MsExperiment)
  library(xcms)
  library(Spectra)
  library(limma)
  library(tibble)
  library(dplyr)
  library(faahKO)
})

cat("===================================================================\n")
cat("Chapter 17 Laboratory: XCMS Metabolomics Preprocessing Pipeline\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Ingestion of Cohort CDF Files into MsExperiment
# -----------------------------------------------------------------------------
cat("--- Section 1: Cohort Ingestion & MsExperiment Architecture ---\n")

cdf_all <- dir(system.file("cdf", package = "faahKO"), full.names = TRUE, recursive = TRUE)

# Select balanced teaching cohort: 2 KO and 2 WT samples
cdf_cohort <- cdf_all[c(1, 2, 7, 8)]
sample_names <- sub(".CDF", "", basename(cdf_cohort), fixed = TRUE)
sample_groups <- c("KO", "KO", "WT", "WT")

sample_df <- data.frame(
  sample_name     = sample_names,
  sample_group    = sample_groups,
  injection_order = 1:4,
  row.names       = basename(cdf_cohort)
)

ms_exp <- readMsExperiment(spectraFiles = cdf_cohort, sampleData = sample_df)

cat(sprintf("Cohort Files Ingested: %d runs\n", length(cdf_cohort)))
for (i in seq_along(cdf_cohort)) {
  cat(sprintf("  [%d] %s (Group: %s)\n", i, sample_names[i], sample_groups[i]))
}
cat(sprintf("Underlying Spectra Object: %d total scans\n\n", length(spectra(ms_exp))))

# -----------------------------------------------------------------------------
# Section 2: CentWave Peak Detection
# -----------------------------------------------------------------------------
cat("--- Section 2: CentWave Chromatographic Peak Detection ---\n")

cwp <- CentWaveParam(
  peakwidth = c(20, 50),
  ppm       = 30,
  snthresh  = 10,
  prefilter = c(3, 100)
)

t0_pks <- Sys.time()
ms_exp <- findChromPeaks(ms_exp, param = cwp, BPPARAM = SerialParam())
t_pks  <- as.numeric(difftime(Sys.time(), t0_pks, units = "secs"))

peaks_mat <- chromPeaks(ms_exp)
n_total_peaks <- nrow(peaks_mat)

cat(sprintf("CentWave Peak Detection Completed in %.2f seconds\n", t_pks))
cat(sprintf("Total Detected Peaks Across Cohort: %d\n", n_total_peaks))

peaks_per_sample <- table(peaks_mat[, "sample"])
for (s in seq_along(peaks_per_sample)) {
  cat(sprintf("  Sample %d (%s): %d peaks\n", s, sample_names[s], peaks_per_sample[s]))
}
cat("\n")

# -----------------------------------------------------------------------------
# Section 3: Retention Time Alignment and Correspondence Grouping
# -----------------------------------------------------------------------------
cat("--- Section 3: Obiwarp RT Alignment & Peak Density Grouping ---\n")

# Step 1: OBI-Warp alignment
owp <- ObiwarpParam(binSize = 1.0)
ms_exp <- adjustRtime(ms_exp, param = owp)

rtime_adj <- rtime(ms_exp, adjusted = TRUE)
rtime_raw <- rtime(ms_exp, adjusted = FALSE)
max_drift <- max(abs(rtime_adj - rtime_raw), na.rm = TRUE)

cat(sprintf("OBI-Warp Alignment Completed. Maximum RT Shift: %.2f seconds\n", max_drift))

# Step 2: Peak Density Correspondence Grouping
pdp <- PeakDensityParam(
  sampleGroups = sampleData(ms_exp)$sample_group,
  minFraction  = 0.5,
  bw           = 30,
  binSize      = 0.025
)

ms_exp <- groupChromPeaks(ms_exp, param = pdp)
feat_defs <- featureDefinitions(ms_exp)
n_features <- nrow(feat_defs)

cat(sprintf("Peak Density Correspondence Grouping: %d Consensus Features Defined\n", n_features))
cat(sprintf("Feature m/z Range: %.4f - %.4f Da | RT Range: %.1f - %.1f s\n\n",
            min(feat_defs$mzmed), max(feat_defs$mzmed),
            min(feat_defs$rtmed), max(feat_defs$rtmed)))

# -----------------------------------------------------------------------------
# Section 4: Gap Filling, QC Filtering, and Differential Abundance
# -----------------------------------------------------------------------------
cat("--- Section 4: Gap Filling & Quality Control Feature Filtering ---\n")

# Matrix before gap filling
feat_mat_raw <- featureValues(ms_exp, value = "into", filled = FALSE)
missing_before <- sum(is.na(feat_mat_raw))
missing_before_pct <- mean(is.na(feat_mat_raw)) * 100

# Perform ChromPeakAreaParam gap filling
fcp <- ChromPeakAreaParam()
ms_exp <- fillChromPeaks(ms_exp, param = fcp, BPPARAM = SerialParam())

feat_mat_filled <- featureValues(ms_exp, value = "into", filled = TRUE)
missing_after <- sum(is.na(feat_mat_filled))
missing_after_pct <- mean(is.na(feat_mat_filled)) * 100

cat(sprintf("Missing Matrix Cells Before Gap Filling: %d / %d (%.2f%%)\n",
            missing_before, length(feat_mat_raw), missing_before_pct))
cat(sprintf("Missing Matrix Cells After Gap Filling:  %d / %d (%.2f%%)\n",
            missing_after, length(feat_mat_filled), missing_after_pct))
cat(sprintf("Gaps Successfully Integrated: %d\n", missing_before - missing_after))

# QC Filtering: Retain features with detection rate >= 50% and CV < 30%
feature_detect <- rowMeans(!is.na(feat_mat_filled))
feature_cv <- apply(feat_mat_filled, 1, function(x) {
  sd(x, na.rm = TRUE) / mean(x, na.rm = TRUE) * 100
})

keep <- feature_detect >= 0.5 & feature_cv < 30.0
n_qc_pass <- sum(keep, na.rm = TRUE)

cat(sprintf("\nQC Filtering (Detection Rate >= 50%% & CV < 30%%):\n"))
cat(sprintf("  Features Passing QC Filter: %d / %d (%.1f%%)\n",
            n_qc_pass, n_features, n_qc_pass / n_features * 100))

# Differential Testing with limma on Normalized Log2 Intensities
mat_qc <- feat_mat_filled[keep, ]
mat_log2 <- log2(mat_qc + 1)
# Median centering
sample_meds <- apply(mat_log2, 2, median, na.rm = TRUE)
mat_norm <- sweep(mat_log2, 2, sample_meds, "-")

group_fac <- factor(sampleData(ms_exp)$sample_group, levels = c("WT", "KO"))
design <- model.matrix(~ group_fac)
colnames(design) <- c("Intercept", "KO_vs_WT")

fit <- lmFit(mat_norm, design)
fit <- eBayes(fit)
top_res <- topTable(fit, coef = "KO_vs_WT", number = 10)

cat("\nTop Differential Metabolite Features (KO vs WT):\n")
print(round(top_res[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 4))

sig_metabs <- top_res %>% filter(adj.P.Val <= 0.05, abs(logFC) >= 1.0)
cat(sprintf("\nSignificant Features (|log2FC| >= 1, FDR <= 0.05): %d\n\n", nrow(sig_metabs)))

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 17 Exercises\n")
cat("===================================================================\n\n")

# Exercise 17.1: CentWave Mass Accuracy Parameterization
cat("--- Exercise 17.1: Orbitrap vs TOF Parameter Tuning ---\n")
cat("For an Orbitrap instrument (resolving power 60,000 - 240,000), set ppm = 5.0 - 10.0.\n")
cat("For a standard TOF instrument (resolving power 10,000 - 30,000), set ppm = 25.0 - 35.0.\n")
cat("Consequence: Using ppm = 30 on an Orbitrap merges near-isobaric metabolite ions,\n")
cat("whereas using ppm = 5 on a TOF splits single chromatographic peaks into fragmented slices.\n\n")

# Exercise 17.2: Alignment vs Correspondence Mechanics
cat("--- Exercise 17.2: Alignment vs Correspondence ---\n")
cat("1. Alignment (adjustRtime): Modifies the retention time coordinates of all spectra\n")
cat("   in a sample to eliminate nonlinear chromatographic column drift.\n")
cat("2. Correspondence (groupChromPeaks): Matches peaks across aligned samples into a 2D\n")
cat("   consensus feature definition based on m/z and retention time density kernels.\n")
cat("Sequential Necessity: Correspondence can theoretically run without alignment, but\n")
cat("chromatographic drift > bw causes identical metabolites to split into duplicate features.\n\n")

# Exercise 17.3: Gap Filling Analytical Value
cat("--- Exercise 17.3: Missing Value Recovery via fillChromPeaks ---\n")
cat(sprintf("Laboratory Result: Gap filling reduced missing cells from %d (%.2f%%) to %d (%.2f%%).\n",
            missing_before, missing_before_pct, missing_after, missing_after_pct))
cat("Analytical Value: Recovers low-abundance signal that missed the peak detection threshold\n")
cat("in individual replicates, preventing artificial truncation in downstream linear models.\n\n")

# Exercise 17.4: QC Pool Reproducibility Filtering
cat("--- Exercise 17.4: Biological vs Technical Variance Filtering ---\n")
cat("Filter Rule: Retain features with QC-pool CV < 20% (or < 30% across cohort).\n")
cat("Theoretical Justification: Features with technical CV > 30% have analytical noise\n")
cat("exceeding the biological effect size, inflating false discoveries in multiple hypothesis testing.\n\n")

cat("===================================================================\n")
cat("Chapter 17 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
