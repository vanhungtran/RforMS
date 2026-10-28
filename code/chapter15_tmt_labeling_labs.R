# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 15 Laboratory: Multiplexed Isobaric Labeling (TMT) Workflows
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(QFeatures)
  library(MsCoreUtils)
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

message("\n=== Listing 15.1: Environment Setup and Synthetic TMT Multiplex Construction ===")
set.seed(101)
n_psms <- 1000
channels <- c("X126", "X127", "X128", "X129", "X130", "X131")

base_intensity <- 2^rnorm(n_psms, mean = 16, sd = 1.2)
interference <- base_intensity * 0.35

raw_mat <- matrix(0, nrow = n_psms, ncol = 6)
colnames(raw_mat) <- channels

for (c_idx in 1:3) {
  raw_mat[, c_idx] <- (base_intensity * 1.0) + interference + rnorm(n_psms, sd = 500)
}
for (c_idx in 4:6) {
  raw_mat[, c_idx] <- (base_intensity * 4.0) + interference + rnorm(n_psms, sd = 500)
}

psm_df <- as.data.frame(raw_mat)
psm_df$Protein <- paste0("PROT_", rep(1:200, each = 5))
psm_df$PSM_ID  <- paste0("PSM_", 1:n_psms)
rownames(psm_df) <- psm_df$PSM_ID

cat(sprintf("Simulated %d PSMs mapping to 200 proteins across 6 TMT channels.\n", n_psms))
cat(sprintf("True Fold Change: 4.00 (log2FC = 2.00) with 35%% co-isolation interference.\n"))

message("\n=== Listing 15.2: Building the QFeatures Data Infrastructure and Sample Normalization ===")
coldata <- data.frame(
  Channel   = channels,
  Condition = c(rep("Control", 3), rep("Treatment", 3)),
  row.names = channels
)

se <- SummarizedExperiment(
  assays  = list(raw = as.matrix(psm_df[, channels])),
  rowData = psm_df[, c("Protein", "PSM_ID")],
  colData = coldata
)

qf <- QFeatures(list(psms = se), colData = coldata)

# Step 1: Sample-Loading Normalization (equalize total channel sums)
qf <- normalize(qf, i = "psms", name = "psms_sl", method = "sum")

# Step 2: Channel Median Centering on log2 scale
qf <- logTransform(qf, i = "psms_sl", name = "psms_log")
qf <- normalize(qf, i = "psms_log", name = "psms_norm", method = "center.median")
cat("Completed 2-step sample loading and channel median normalization.\n")

message("\n=== Listing 15.3: Protein Aggregation via Robust M-Estimation ===")
qf <- aggregateFeatures(
  qf,
  i    = "psms_norm",
  name = "proteins",
  fcol = "Protein",
  fun  = MsCoreUtils::robustSummary
)

prot_mat <- assay(qf[["proteins"]])
cat(sprintf("Aggregated 1000 PSMs into %d protein groups.\n", nrow(prot_mat)))

message("\n=== Listing 15.4: Differential Abundance Testing and Ratio Compression Audit ===")
design <- model.matrix(~ Condition, data = as.data.frame(colData(qf)))
fit    <- lmFit(prot_mat, design)
fit    <- eBayes(fit)
res    <- topTable(fit, coef = "ConditionTreatment", number = Inf)

mean_obs_log2fc <- mean(res$logFC)
mean_obs_fc     <- 2^mean_obs_log2fc

cat(sprintf("Theoretical Biological Fold Change: 4.00 (log2FC = 2.00)\n"))
cat(sprintf("Empirically Observed Fold Change:   %.2f (log2FC = %.2f)\n", 
            mean_obs_fc, mean_obs_log2fc))
cat(sprintf("Ratio Compression Magnitude:       %.1f%% signal dampening!\n", 
            (1 - mean_obs_fc / 4.0) * 100))

message("\n=== Listing 15.5: Automated TMT Channel Coefficient of Variation Audit ===")
audit_tmt_channels <- function(tmt_matrix, cv_threshold = 10.0) {
  cv_results <- data.frame(
    Channel = colnames(tmt_matrix),
    Mean_Intensity = colMeans(tmt_matrix, na.rm = TRUE),
    SD_Intensity   = apply(tmt_matrix, 2, sd, na.rm = TRUE)
  ) %>%
    mutate(
      CV_Percent = (SD_Intensity / Mean_Intensity) * 100,
      Status     = ifelse(CV_Percent > cv_threshold, "FLAGGED (High CV)", "PASS")
    )
  return(cv_results)
}

test_data <- matrix(rnorm(600, mean = 20000, sd = 1500), nrow = 100, ncol = 6)
test_data[, 4] <- test_data[, 4] * rnorm(100, mean = 1, sd = 0.35)
colnames(test_data) <- paste0("TMT_", 1:6)

qc_summary <- audit_tmt_channels(test_data)
print(qc_summary[, c("Channel", "CV_Percent", "Status")])
message("\n[PASS] Chapter 15 Laboratory successfully executed.")
