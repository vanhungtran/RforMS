# =============================================================================
# Chapter 07: Raw Spectral Quality Control and Instrument Diagnostics
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(Spectra)
  library(msdata)
  library(tibble)
  library(dplyr)
  library(tidyr)
})

cat("===================================================================\n")
cat("Chapter 07 Laboratory: Mass Spectrometry Quality Control Auditing\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Experimental Design & Metadata Confounding Audit
# -----------------------------------------------------------------------------
cat("--- Section 1: Experimental Design & Confounding Audit ---\n")

# Simulate a 12-sample study design: evaluate Confounded vs Randomized Design
set.seed(42)
confounded_meta <- tibble(
  sample_id        = paste0("Sample_", sprintf("%02d", 1:12)),
  biological_group = rep(c("Control", "Disease"), each = 6),
  batch            = c(rep("Batch_1", 6), rep("Batch_2", 6)), # Completely confounded!
  injection_order  = 1:12
)

cat("Confounded Experimental Design Matrix:\n")
tbl_conf <- table(confounded_meta$batch, confounded_meta$biological_group)
print(tbl_conf)

# Test of independence
is_confounded <- any(rowSums(tbl_conf > 0) == 1)
cat(sprintf("Confounding Flag: %s (Batch effect completely collinear with Disease status)\n",
            if (is_confounded) "CRITICAL_FAIL" else "PASS"))

# Balanced randomized design
randomized_meta <- tibble(
  sample_id        = paste0("Sample_", sprintf("%02d", 1:12)),
  biological_group = rep(c("Control", "Disease"), each = 6),
  batch            = rep(c("Batch_1", "Batch_2"), 6), # Balanced across batches
  injection_order  = sample(1:12)                     # Randomized run sequence
)

cat("\nBalanced & Randomized Experimental Design Matrix:\n")
tbl_rand <- table(randomized_meta$batch, randomized_meta$biological_group)
print(tbl_rand)
cat("Status: PASS (Equal group representation per batch; zero collinearity)\n\n")

# -----------------------------------------------------------------------------
# Section 2: Real Raw File TIC and BPC Extraction with Spectra
# -----------------------------------------------------------------------------
cat("--- Section 2: Raw File TIC & BPC Auditing via Spectra ---\n")

f_raw <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
sps <- Spectra(f_raw, backend = MsBackendMzR())

# Filter to MS1 survey scans
sps_ms1 <- filterMsLevel(sps, 1L)
rts <- rtime(sps_ms1)
tics <- tic(sps_ms1)
bpis <- sps_ms1$basePeakIntensity

cat(sprintf("Auditing File: %s\n", basename(f_raw)))
cat(sprintf("Survey Scans: %d | RT Range: [%.2f, %.2f] s ([%.2f, %.2f] min)\n",
            length(sps_ms1), min(rts), max(rts), min(rts)/60, max(rts)/60))
cat(sprintf("Area Under TIC (AUC): %.2e counts\n", sum(tics * c(diff(rts), median(diff(rts))))))
cat(sprintf("TIC CV across survey scans: %.2f%%\n", sd(tics) / mean(tics) * 100))
cat(sprintf("Median Base Peak Intensity: %.2e counts\n\n", median(bpis)))

# -----------------------------------------------------------------------------
# Section 3: Mass Accuracy Drift and Lock Mass Monitoring
# -----------------------------------------------------------------------------
cat("--- Section 3: Mass Accuracy Drift across Injection Sequence ---\n")

# Reference compound: theoretical m/z = 508.1234 (e.g., standard hexapeptide)
theor_mz <- 508.1234
n_injections <- 10

# Simulate progressive calibration drift across 10 sequential injections
inj_seq <- 1:n_injections
drift_ppm_true <- 0.65 * (inj_seq - 1) - 1.2 # drift from -1.2 ppm to +4.65 ppm
obs_mz <- theor_mz * (1 + drift_ppm_true * 1e-6) + rnorm(n_injections, 0, 0.0001)

mass_acc_df <- tibble(
  Injection_Order = inj_seq,
  Theoretical_MZ  = theor_mz,
  Observed_MZ     = round(obs_mz, 5),
  Error_Da        = round(obs_mz - theor_mz, 5),
  Error_PPM       = round((obs_mz - theor_mz) / theor_mz * 1e6, 2)
)

print(as.data.frame(mass_acc_df))

fit_drift <- lm(Error_PPM ~ Injection_Order, data = mass_acc_df)
cat(sprintf("\nMass Drift Rate: %.3f ppm per injection (p = %.2e)\n",
            coef(fit_drift)[2], summary(fit_drift)$coefficients[2, 4]))
cat(sprintf("Tolerance Limit: |ppm| <= 5.0 ppm -> All injections meet criteria: %s\n\n",
            all(abs(mass_acc_df$Error_PPM) <= 5.0)))

# -----------------------------------------------------------------------------
# Section 4: Blank Contamination and Carryover Audit
# -----------------------------------------------------------------------------
cat("--- Section 4: Solvent Blank Contamination & Carryover Audit ---\n")

# High-concentration calibrator (500 nM) followed by Solvent Blank
analytes <- c("Tryptophan", "Kynurenine", "Phenylalanine", "Tyrosine")
high_std_area <- c(4500000, 2800000, 6200000, 3100000)
blank_post_area <- c(4200, 18500, 3100, 89000) # Tyrosine shows 2.87% carryover!

carryover_df <- tibble(
  Analyte           = analytes,
  High_Std_Area     = high_std_area,
  Post_Blank_Area   = blank_post_area,
  Carryover_Pct     = round(blank_post_area / high_std_area * 100, 3),
  FDA_Criteria_Pass = blank_post_area / high_std_area * 100 <= 1.0 # Standard threshold <= 1%
)

print(as.data.frame(carryover_df))
cat("\nCarryover Audit Findings:\n")
cat("  - Tryptophan, Kynurenine, and Phenylalanine exhibit < 0.7% carryover (PASS).\n")
cat("  - Tyrosine exhibits 2.87% carryover, exceeding the 1.0% bioanalytical threshold (FLAG_CARRYOVER).\n\n")

# -----------------------------------------------------------------------------
# Section 5: Pooled QC Precision and Minimum Reliable Fold-Change
# -----------------------------------------------------------------------------
cat("--- Section 5: Pooled QC Coefficient of Variation (CV) Audit ---\n")

# 100 features measured across 6 pooled QC injections
n_feat <- 100
n_qc   <- 6

# Feature intensities with varying technical noise
qc_mat <- matrix(
  round(rlnorm(n_feat * n_qc, meanlog = 14, sdlog = 1.0) *
          matrix(rnorm(n_feat * n_qc, mean = 1, sd = rep(c(0.05, 0.12, 0.28, 0.45), each = 25)),
                 nrow = n_feat, ncol = n_qc)),
  nrow = n_feat, ncol = n_qc
)

# Compute feature CV%
qc_means <- rowMeans(qc_mat)
qc_sds   <- apply(qc_mat, 1, sd)
qc_cvs   <- round(qc_sds / qc_means * 100, 2)

cv_binned <- tibble(
  Tier_CV = c("<= 15% (High Precision)", "15% - 20% (Acceptable)", "20% - 30% (Metabolomics Limit)", "> 30% (Unreliable Noise)"),
  Feature_Count = c(
    sum(qc_cvs <= 15),
    sum(qc_cvs > 15 & qc_cvs <= 20),
    sum(qc_cvs > 20 & qc_cvs <= 30),
    sum(qc_cvs > 30)
  ),
  Percentage = round(c(
    sum(qc_cvs <= 15),
    sum(qc_cvs > 15 & qc_cvs <= 20),
    sum(qc_cvs > 20 & qc_cvs <= 30),
    sum(qc_cvs > 30)
  ) / n_feat * 100, 1)
)

print(as.data.frame(cv_binned))

# Minimum reliable biological fold-change derived from technical noise (2 * CV)
median_cv <- median(qc_cvs)
min_reliable_fc <- 1 + 2 * (median_cv / 100)
cat(sprintf("\nMedian QC Feature CV: %.2f%%\n", median_cv))
cat(sprintf("Statistical Lower Bound for Biological Fold-Change (1 + 2*CV): %.2fx\n",
            min_reliable_fc))
cat("Pedagogical Implication: Any biological fold-change below this threshold is indistinguishable from technical instrument noise.\n\n")

# -----------------------------------------------------------------------------
# Section 6: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 07 Exercises\n")
cat("===================================================================\n\n")

# Exercise 7.1: Confounding and Degrees of Freedom
cat("--- Exercise 7.1: Confounded Experimental Design Degrees of Freedom ---\n")
# In an experiment with 20 patients (10 Control, 10 Treated) run on 2 batches:
# Batch 1 contains 10 Controls; Batch 2 contains 10 Treated.
# Model: Intensity ~ Group + Batch
cat("Experimental Model: y ~ Condition + Batch\n")
cat("Design Matrix Rank Audit:\n")
cat("  - Vector for Condition (Treated = 1, Control = 0): [0, 0, ..., 1, 1]\n")
cat("  - Vector for Batch (Batch 2 = 1, Batch 1 = 0):     [0, 0, ..., 1, 1]\n")
cat("  - Linear Dependence: Condition == Batch (Pearson correlation r = 1.000)\n")
cat("  - Rank of Design Matrix: rank([X_Condition, X_Batch]) = 1 (Rank Deficient)\n")
cat("Conclusion: The condition effect and batch effect are mathematically inseparable.\n")
cat("Zero degrees of freedom remain to estimate technical batch differences.\n\n")

# Exercise 7.2: Mass Accuracy Drift Calculation
cat("--- Exercise 7.2: Mass Accuracy PPM Formula and Drift Analysis ---\n")
theor_pep <- 1296.6853 # Bovine serum albumin peptide [M+H]+
obs_pep_1 <- 1296.6885 # Injection 1
obs_pep_20 <- 1296.6742 # Injection 20

ppm_1 <- (obs_pep_1 - theor_pep) / theor_pep * 1e6
ppm_20 <- (obs_pep_20 - theor_pep) / theor_pep * 1e6
delta_ppm <- ppm_20 - ppm_1

cat(sprintf("BSA Peptide Theoretical m/z: %.4f\n", theor_pep))
cat(sprintf("  Injection 01 Observed: %.4f -> Error: %+.2f ppm\n", obs_pep_1, ppm_1))
cat(sprintf("  Injection 20 Observed: %.4f -> Error: %+.2f ppm\n", obs_pep_20, ppm_20))
cat(sprintf("  Total Sequence Drift: %+.2f ppm\n", delta_ppm))
cat("Evaluation: Injection 20 exceeds 5 ppm tolerance (|ppm| = 8.56 > 5.0), requiring post-acquisition recalibration.\n\n")

# Exercise 7.3: Carryover Determination and ULOQ Contamination
cat("--- Exercise 7.3: Autosampler Carryover Ratio and LLOQ Comparison ---\n")
uloq_peak_area  <- 8500000
blank_post_uloq <- 34000
lloq_peak_area  <- 25000

carryover_ratio <- blank_post_uloq / uloq_peak_area * 100
carryover_vs_lloq <- blank_post_uloq / lloq_peak_area * 100

cat(sprintf("ULOQ Peak Area:          %.0f counts\n", uloq_peak_area))
cat(sprintf("Post-ULOQ Blank Area:    %.0f counts\n", blank_post_uloq))
cat(sprintf("LLOQ Peak Area:          %.0f counts\n", lloq_peak_area))
cat(sprintf("Carryover relative to ULOQ: %.2f%%\n", carryover_ratio))
cat(sprintf("Carryover relative to LLOQ: %.1f%% (FDA criteria requires <= 20.0%% of LLOQ)\n",
            carryover_vs_lloq))
cat(sprintf("Regulatory Verdict: %s (Carryover equals 136%% of LLOQ signal, causing false positives at low concentrations)\n\n",
            if (carryover_vs_lloq > 20) "FAIL" else "PASS"))

# Exercise 7.4: Pooled QC CV and Minimum Resolvable Effect Size
cat("--- Exercise 7.4: Pooled QC Coefficient of Variation Distribution ---\n")
cat("Given a metabolomics QC pool where a feature has mean intensity 450,000 and SD 99,000:\n")
feat_cv <- 99000 / 450000 * 100
fc_lower <- 1 + 2 * (feat_cv / 100)
fc_upper <- 1 / (1 - 2 * (feat_cv / 100))

cat(sprintf("Feature QC CV: %.1f%%\n", feat_cv))
cat(sprintf("Minimum Detectable Up-Regulation:   +%.1f%% (Fold-Change >= %.2fx)\n",
            (fc_lower - 1) * 100, fc_lower))
cat(sprintf("Minimum Detectable Down-Regulation: -%.1f%% (Fold-Change <= %.2fx)\n",
            (1 - 1/fc_upper) * 100, 1/fc_upper))
cat("Scientific Conclusion: Reporting a 1.25x fold-change for this metabolite is invalid,\n")
cat("as 1.25 falls well within the 95% technical noise interval [0.56x, 1.44x].\n\n")

cat("===================================================================\n")
cat("Chapter 07 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
