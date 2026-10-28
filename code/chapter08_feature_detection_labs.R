# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 08 Laboratory: Chromatographic Feature Detection & Alignment with xcms
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(xcms)
  library(MSnbase)
  library(BiocParallel)
  library(ggplot2)
})

# Enforce serial execution on Windows to eliminate worker socket overhead
register(SerialParam(), default = TRUE)

message("\n=== Step 1: Loading On-Disk MSnExp Data Container ===")
cdf_dir <- system.file("cdf", package = "faahKO")
raw_files <- list.files(cdf_dir, pattern = "\\.CDF$", full.names = TRUE, recursive = TRUE)
selected_files <- raw_files[c(1, 2, 7, 8)]

sample_info <- data.frame(
  sample_name = sub("\\.CDF$", "", basename(selected_files), ignore.case = TRUE),
  condition   = factor(c("KO", "KO", "WT", "WT")),
  row.names   = basename(selected_files)
)

raw_data <- readMSData(
  selected_files,
  pdata = new("AnnotatedDataFrame", sample_info),
  mode  = "onDisk"
)

cat(sprintf("Loaded on-disk MS container: %d samples, %d total spectra.\n",
            length(selected_files), length(raw_data)))
print(table(raw_data$condition))

message("\n=== Step 2: CentWave Peak Picking with Parameter Tuning ===")
cwp <- CentWaveParam(
  ppm       = 25,
  peakwidth = c(20, 50),
  snthresh  = 10,
  noise     = 500,
  prefilter = c(3, 100)
)

data_peaks <- findChromPeaks(raw_data, param = cwp, BPPARAM = SerialParam())
chrom_peaks <- chromPeaks(data_peaks)

cat(sprintf("Total chromatographic peaks detected across 4 samples: %d peaks.\n",
            nrow(chrom_peaks)))
peaks_per_sample <- table(chrom_peaks[, "sample"])
names(peaks_per_sample) <- sample_info$sample_name[as.numeric(names(peaks_per_sample))]
print(peaks_per_sample)

message("\n=== Step 3: Obiwarp Alignment and Density-Based Feature Correspondence ===")
owp <- ObiwarpParam(
  binSize   = 0.6,
  response  = 1,
  distFun   = "cor_opt",
  gapInit   = 0.3,
  gapExtend = 2.4
)

data_aligned <- adjustRtime(data_peaks, param = owp)

pdp <- PeakDensityParam(
  sampleGroups = pData(data_aligned)$condition,
  minFraction  = 0.50,
  bw           = 5,
  binSize      = 0.025
)

data_grouped <- groupChromPeaks(data_aligned, param = pdp)
feature_def <- featureDefinitions(data_grouped)

cat(sprintf("Number of aligned features across cohort: %d features.\n", nrow(feature_def)))
cat(sprintf("Mean chromatographic feature width: %.1f seconds.\n", 
            mean(feature_def$rtmax - feature_def$rtmin)))

message("\n=== Step 4: Raw Ion Re-integration via fillChromPeaks ===")
mat_before <- featureValues(data_grouped, value = "into")
na_count_before <- sum(is.na(mat_before))

data_filled <- fillChromPeaks(data_grouped, param = FillChromPeaksParam(), BPPARAM = SerialParam())
expr_mat <- featureValues(data_filled, value = "into")
na_count_after <- sum(is.na(expr_mat))

cat(sprintf("Missing values before gap-filling: %d (%.1f%% of matrix).\n",
            na_count_before, na_count_before / length(mat_before) * 100))
cat(sprintf("Missing values after fillChromPeaks: %d (%.1f%% of matrix).\n",
            na_count_after, na_count_after / length(expr_mat) * 100))

cat("\n--- Top 4 Aligned and Gap-Filled Features ---\n")
print(head(round(expr_mat, 0), 4))

message("\n=== Step 5: Production-Grade Feature Quality Filter ===")
filter_high_quality_features <- function(feature_mat, sample_classes, 
                                        max_na_fraction = 0.20, max_qc_cv = 0.30) {
  stopifnot(is.matrix(feature_mat), length(sample_classes) == ncol(feature_mat))
  initial_features <- nrow(feature_mat)
  
  na_frac <- apply(feature_mat, 1, function(x) mean(is.na(x) | x == 0))
  pass_na <- na_frac <= max_na_fraction
  
  qc_indices <- which(sample_classes == "QC")
  if (length(qc_indices) >= 3) {
    qc_mat <- feature_mat[, qc_indices, drop = FALSE]
    qc_mean <- rowMeans(qc_mat, na.rm = TRUE)
    qc_sd   <- apply(qc_mat, 1, sd, na.rm = TRUE)
    qc_cv   <- ifelse(qc_mean > 0, qc_sd / qc_mean, Inf)
    pass_cv <- is.finite(qc_cv) & (qc_cv <= max_qc_cv)
  } else {
    message("Fewer than 3 QC samples detected; skipping technical CV filter.")
    pass_cv <- rep(TRUE, initial_features)
  }
  
  retained_idx <- which(pass_na & pass_cv)
  filtered_mat <- feature_mat[retained_idx, , drop = FALSE]
  
  cat(sprintf("[STATUS: COMPLETE] Filtered %d features -> %d retained (%.1f%% pass rate).\n",
              initial_features, nrow(filtered_mat), nrow(filtered_mat) / initial_features * 100))
  cat(sprintf("Excluded for excessive missingness (> %.0f%%): %d features.\n",
              max_na_fraction * 100, sum(!pass_na)))
  return(filtered_mat)
}

set.seed(42)
sim_mat <- matrix(rlnorm(500 * 8, 10, 0.15), nrow = 500)
sim_mat[sample(length(sim_mat), 500)] <- NA
classes <- c(rep("Bio", 5), rep("QC", 3))

filtered_res <- filter_high_quality_features(sim_mat, classes, max_na_fraction = 0.20, max_qc_cv = 0.30)
message("\n[PASS] Chapter 08 Laboratory successfully executed.")
