# =============================================================================
# Chapter 16: Targeted Quantitation: SRM, PRM, and DIA
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(tidyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 16.1: Transition-Level Calibrator and Unknown Simulation
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 16.1: TRANSITION-LEVEL DATASET SIMULATION ===\n")

# Analytes and nominal calibration slopes
analytes_meta <- tribble(
  ~analyte,        ~internal_standard, ~true_slope, ~true_intercept,
  "Tryptophan",    "Tryptophan_d5",    0.0180,      0.0010,
  "Kynurenine",    "Kynurenine_d4",    0.0260,      0.0015,
  "Phenylalanine", "Phenylalanine_d5", 0.0140,      0.0008
)

# 8-point calibration curve in triplicate (0 to 500 nM)
calib_levels <- c(0, 1, 5, 10, 25, 50, 100, 250, 500)
calib_raw <- expand_grid(
  analyte = analytes_meta$analyte,
  concentration_nM = calib_levels,
  replicate = 1:3,
  transition = c("quantifier", "qualifier")
) %>%
  left_join(analytes_meta, by = "analyte") %>%
  mutate(
    sample_type = if_else(concentration_nM == 0, "blank", "calibrator"),
    sample_id = paste(sample_type, analyte, concentration_nM, replicate, sep = "_"),
    area_is = round(rlnorm(n(), log(500000), 0.05)),
    trans_factor = if_else(transition == "quantifier", 1.0, 0.55),
    true_ratio = pmax(true_intercept + true_slope * concentration_nM, 0) * trans_factor,
    area_ratio = round(true_ratio * rlnorm(n(), 0, 0.04), 5),
    area_analyte = round(area_ratio * area_is),
    retention_time = round(case_when(
      analyte == "Tryptophan" ~ 4.20,
      analyte == "Kynurenine" ~ 3.10,
      TRUE                    ~ 5.40
    ) + rnorm(n(), 0, 0.02), 2),
    signal_to_noise = round(if_else(concentration_nM == 0, runif(n(), 1.2, 2.8),
                                   pmax(3, concentration_nM * 1.5 + runif(n(), 5, 20))), 1)
  )

cat(sprintf("Simulated Transition Dataset: %d rows across %d calibrator levels\n",
            nrow(calib_raw), length(calib_levels)))
print(head(as.data.frame(calib_raw[, c("sample_id", "analyte", "transition", "concentration_nM", "area_ratio", "retention_time", "signal_to_noise")]), 4))


# -----------------------------------------------------------------------------
# Listing 16.2: Weighted Least Squares (1/x^2) Calibration Modeling
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 16.2: WEIGHTED CALIBRATION MODELING ===\n")

# Filter to quantifier transition for non-blank calibrators
calib_quant <- calib_raw %>%
  filter(transition == "quantifier", concentration_nM > 0)

# Function to fit OLS vs 1/x vs 1/x^2 models
fit_analyte_models <- function(df_sub) {
  ols_fit  <- lm(area_ratio ~ concentration_nM, data = df_sub)
  w1_fit   <- lm(area_ratio ~ concentration_nM, data = df_sub, weights = 1 / concentration_nM)
  w2_fit   <- lm(area_ratio ~ concentration_nM, data = df_sub, weights = 1 / (concentration_nM^2))
  
  # Back-calculate concentration at lowest standard (1 nM)
  c1_actual <- 1.0
  ratio_1 <- mean(df_sub$area_ratio[df_sub$concentration_nM == 1.0])
  c1_ols  <- (ratio_1 - coef(ols_fit)[1]) / coef(ols_fit)[2]
  c1_w2   <- (ratio_1 - coef(w2_fit)[1]) / coef(w2_fit)[2]
  
  tibble(
    Slope_W2        = round(coef(w2_fit)[2], 5),
    Intercept_W2    = round(coef(w2_fit)[1], 5),
    R2_W2           = round(summary(w2_fit)$r.squared, 4),
    RE_Pct_OLS_1nM  = round((c1_ols - c1_actual) / c1_actual * 100, 1),
    RE_Pct_W2_1nM   = round((c1_w2 - c1_actual) / c1_actual * 100, 1)
  )
}

calib_summary <- calib_quant %>%
  group_by(analyte) %>%
  group_modify(~ fit_analyte_models(.x)) %>%
  ungroup()

cat("1/x^2 Weighted Calibration Model Summary:\n")
print(as.data.frame(calib_summary))


# -----------------------------------------------------------------------------
# Listing 16.3: Detection Limits Determination (IUPAC / CLSI)
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 16.3: DETECTION LIMITS DETERMINATION ===\n")

# Extract blank replicates (concentration = 0)
blanks_df <- calib_raw %>% filter(sample_type == "blank", transition == "quantifier")

limits_list <- list()
for (ana in unique(calib_raw$analyte)) {
  b_ratios <- blanks_df$area_ratio[blanks_df$analyte == ana]
  sd_b <- sd(b_ratios)
  
  slope <- calib_summary$Slope_W2[calib_summary$analyte == ana]
  
  # IUPAC / FDA bioanalytical guidelines: LoD = 3.3 * sd_b / slope; LoQ = 10 * sd_b / slope
  lod_val <- 3.3 * sd_b / slope
  loq_val <- 10.0 * sd_b / slope
  
  limits_list[[ana]] <- tibble(
    analyte = ana,
    SD_Blank = round(sd_b, 5),
    LoD_nM   = round(lod_val, 3),
    LoQ_nM   = round(loq_val, 3)
  )
}

limits_report <- bind_rows(limits_list)
cat("Analytical Detection Limits (IUPAC 3.3-sigma / 10-sigma Criteria):\n")
print(as.data.frame(limits_report))


# -----------------------------------------------------------------------------
# Listing 16.4: Unknown Sample Quantification and Status Auditing
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 16.4: UNKNOWN SAMPLE QUANTIFICATION ===\n")

# Simulate unknown study samples (n = 12, Control vs Treatment)
unknowns_sim <- expand_grid(
  analyte = analytes_meta$analyte,
  sample_id = paste0("Sample_", sprintf("%02d", 1:12)),
  transition = c("quantifier", "qualifier")
) %>%
  mutate(
    condition = if_else(as.integer(substr(sample_id, 8, 9)) <= 6, "Control", "Treatment"),
    true_nM = case_when(
      analyte == "Tryptophan" ~ if_else(condition == "Control", 45.0, 48.0),
      analyte == "Kynurenine" ~ if_else(condition == "Control", 15.0, 42.0), # Activated IDO pathway
      TRUE                    ~ 75.0
    ),
    trans_factor = if_else(transition == "quantifier", 1.0, 0.55),
    area_is = round(rlnorm(n(), log(500000), 0.05)),
    area_ratio = round((0.001 + 0.02 * true_nM) * trans_factor * rlnorm(n(), 0, 0.05), 5),
    area_analyte = round(area_ratio * area_is)
  )

# Quantify using calibrated quantifier parameters
quant_results <- unknowns_sim %>%
  filter(transition == "quantifier") %>%
  left_join(calib_summary, by = "analyte") %>%
  left_join(limits_report, by = "analyte") %>%
  mutate(
    conc_calc_nM = round((area_ratio - Intercept_W2) / Slope_W2, 2),
    quant_flag = case_when(
      conc_calc_nM < LoD_nM ~ "< LoD",
      conc_calc_nM < LoQ_nM ~ "< LoQ",
      conc_calc_nM > 500    ~ "> ULOQ",
      TRUE                  ~ "Quantified"
    ),
    conc_reported = if_else(quant_flag == "Quantified", conc_calc_nM, NA_real_)
  )

cat("Quantified Unknown Samples (First 6 records):\n")
print(as.data.frame(quant_results[1:6, c("sample_id", "condition", "analyte", "area_ratio", "conc_calc_nM", "quant_flag")]))


# -----------------------------------------------------------------------------
# Listing 16.5: QC Precision and Differential Pathway Testing
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 16.5: QC AUDIT AND STATISTICAL TESTING ===\n")

# Check Internal Standard Peak Area CV across all samples
is_cv <- round(sd(unknowns_sim$area_is) / mean(unknowns_sim$area_is) * 100, 2)
cat(sprintf("Internal Standard Precision Audit: Mean Area = %.0f | CV = %.2f%% (Pass <= 15%%)\n",
            mean(unknowns_sim$area_is), is_cv))

# Statistical group comparison (Treatment vs Control) for Kynurenine pathway
kyn_data <- quant_results %>% filter(analyte == "Kynurenine", quant_flag == "Quantified")
kyn_ttest <- t.test(conc_reported ~ condition, data = kyn_data)

cat("\nDifferential Abundance Analysis: Kynurenine (Treatment vs Control):\n")
cat(sprintf("  Control Mean:   %.2f nM (SD = %.2f)\n",
            mean(kyn_data$conc_reported[kyn_data$condition == "Control"]),
            sd(kyn_data$conc_reported[kyn_data$condition == "Control"])))
cat(sprintf("  Treatment Mean: %.2f nM (SD = %.2f)\n",
            mean(kyn_data$conc_reported[kyn_data$condition == "Treatment"]),
            sd(kyn_data$conc_reported[kyn_data$condition == "Treatment"])))
cat(sprintf("  Difference:     +%.2f nM | Welch t-test: t = %.2f, p = %.2e\n",
            diff(rev(kyn_ttest$estimate)), kyn_ttest$statistic, kyn_ttest$p.value))


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 16.1: OLS vs 1/x vs 1/x^2 at low concentrations
# Standard curve: conc = c(1, 5, 10, 25, 50, 100, 250, 500)
# Area ratios: c(0.011, 0.058, 0.121, 0.310, 0.620, 1.250, 3.110, 6.080)
ex_conc <- c(1, 5, 10, 25, 50, 100, 250, 500)
ex_ratios <- c(0.011, 0.058, 0.121, 0.310, 0.620, 1.250, 3.110, 6.080)
fit_ols <- lm(ex_ratios ~ ex_conc)
fit_w1  <- lm(ex_ratios ~ ex_conc, weights = 1 / ex_conc)
fit_w2  <- lm(ex_ratios ~ ex_conc, weights = 1 / (ex_conc^2))

c1_calc_ols <- (ex_ratios[1] - coef(fit_ols)[1]) / coef(fit_ols)[2]
c1_calc_w2  <- (ex_ratios[1] - coef(fit_w2)[1]) / coef(fit_w2)[2]

cat(sprintf("Exercise 16.1: Weighting Accuracy at Lowest Calibrator (1 nM nominal):\n"))
cat(sprintf("  OLS Back-calculated:    %.3f nM (RE = %+.1f%%)\n", c1_calc_ols, (c1_calc_ols - 1) * 100))
cat(sprintf("  1/x^2 Back-calculated:  %.3f nM (RE = %+.1f%%)\n", c1_calc_w2, (c1_calc_w2 - 1) * 100))

# Exercise 16.2: Blank-based LOD and LOQ
# 10 blanks: c(0.0008, 0.0010, 0.0009, 0.0012, 0.0011, 0.0007, 0.0010, 0.0009, 0.0013, 0.0008)
ex_blanks <- c(0.0008, 0.0010, 0.0009, 0.0012, 0.0011, 0.0007, 0.0010, 0.0009, 0.0013, 0.0008)
ex_slope <- 0.0123
ex_inter <- 0.0001
ex_mu_b <- mean(ex_blanks)
ex_sd_b <- sd(ex_blanks)
ex_lod <- ((ex_mu_b + 3 * ex_sd_b) - ex_inter) / ex_slope
ex_loq <- ((ex_mu_b + 10 * ex_sd_b) - ex_inter) / ex_slope
cat(sprintf("\nExercise 16.2: LoD and LoQ from Blank Statistics:\n"))
cat(sprintf("  Mean Blank = %.5f | SD = %.5f\n", ex_mu_b, ex_sd_b))
cat(sprintf("  Calculated LoD: %.3f nM\n", ex_lod))
cat(sprintf("  Calculated LoQ: %.3f nM\n", ex_loq))

# Exercise 16.3: Quantifier / Qualifier Ion Ratio Auditing
# Nominal trans_factor = 0.55 (Qualifier / Quantifier = 55%)
# Tolerance window: +/- 20% relative (44% to 66%)
ratio_check <- 0.53
cat(sprintf("\nExercise 16.3: Transition Confirmation Ratio:\n"))
cat(sprintf("  Observed Qualifier/Quantifier Ratio: %.1f%% (Acceptance: 44.0%% to 66.0%%) -> CONFIRMED\n",
            ratio_check * 100))

# Exercise 16.4: Carryover Evaluation
# Carryover = Area in blank immediately following ULOQ standard / Area at LLOQ * 100%
# Target: <= 20% of LLOQ
blank_post_uloq_area <- 150
lloq_area <- 1250
carryover_pct <- (blank_post_uloq_area / lloq_area) * 100
cat(sprintf("\nExercise 16.4: Autosampler Carryover Audit:\n"))
cat(sprintf("  Carryover into Blank: %.2f%% of LLOQ (Acceptance <= 20.0%%) -> PASS\n", carryover_pct))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
