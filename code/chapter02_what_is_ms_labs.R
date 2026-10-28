# =============================================================================
# Chapter 02: Mass Spectrometry Principles and Instrument Physics
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(tidyr)
})

cat("===================================================================\n")
cat("Chapter 02 Laboratory: Mass Spectrometry Physics & Signal Processing\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Resolving Charge State from Isotopic Peak Spacing
# -----------------------------------------------------------------------------
cat("--- Section 1: Charge State Deconvolution from Isotope Spacing ---\n")

# Proton mass constant (Da)
m_proton <- 1.007276466879
c13_diff <- 1.003354835 # Mass difference between 13C and 12C

# Simulated experimental isotopic clusters for different charge states
# Peptide 1 (z = 1): Monoisotopic m/z = 450.2500
# Peptide 2 (z = 2): Monoisotopic m/z = 450.2500
# Peptide 3 (z = 3): Monoisotopic m/z = 450.2500
# Intact Protein (z = 20): Monoisotopic m/z = 848.5824 (Myoglobin, ~16,951.5 Da)

test_charges <- c(1L, 2L, 3L, 20L)
mz_monos <- c(450.2500, 450.2500, 450.2500, 848.5824)

deconv_list <- list()
for (i in seq_along(test_charges)) {
  z_true <- test_charges[i]
  mz_0 <- mz_monos[i]
  delta_mz_true <- c13_diff / z_true
  mz_1 <- mz_0 + delta_mz_true
  
  # Measured spacing with small instrument noise
  measured_delta <- round(mz_1 - mz_0, 4)
  inferred_z <- round(c13_diff / measured_delta)
  inferred_mass <- inferred_z * mz_0 - inferred_z * m_proton
  
  deconv_list[[i]] <- tibble(
    True_Z         = z_true,
    Observed_MZ_0  = mz_0,
    Observed_MZ_1  = round(mz_1, 4),
    Delta_MZ       = measured_delta,
    Inferred_Z     = as.integer(inferred_z),
    Neutral_Mass   = round(inferred_mass, 2)
  )
}

deconv_df <- bind_rows(deconv_list)
print(as.data.frame(deconv_df))
cat("\nKey Takeaway: As charge z increases, isotope spacing shrinks proportionally to 1/z.\n")
cat("For z = 1: delta = 1.0034 Da; for z = 2: delta = 0.5017 Da; for z = 3: delta = 0.3345 Da.\n\n")

# -----------------------------------------------------------------------------
# Section 2: Mass Resolving Power and Isobaric Separation
# -----------------------------------------------------------------------------
cat("--- Section 2: Resolving Power and Peak Separation ---\n")

# Formula: R = m / Delta_m_FWHM  =>  Delta_m_FWHM = m / R
analyzers <- tibble(
  Analyzer = c("Unit Quadrupole", "Time-of-Flight (TOF)", "Orbitrap Standard", "Orbitrap Ultra-High"),
  Nominal_R = c(1000, 30000, 60000, 240000)
)

test_mz <- 500.0
analyzers <- analyzers %>%
  mutate(
    FWHM_at_500 = round(test_mz / Nominal_R, 5),
    PPM_Tolerance = round((test_mz / Nominal_R) / test_mz * 1e6 / 2, 1) # ~ half-window
  )

print(as.data.frame(analyzers))

# Isobaric doublet challenge:
# Glutamine: C5H10N2O3 = 146.0691 Da
# Lysine:    C6H14N2O2 = 146.1055 Da
m_gln <- 146.0691
m_lys <- 146.1055
delta_isobaric <- abs(m_lys - m_gln)
req_resolving_power <- round(m_gln / delta_isobaric)

cat(sprintf("\nIsobaric Doublet: Glutamine (%.4f) vs Lysine (%.4f) | Delta = %.4f Da\n",
            m_gln, m_lys, delta_isobaric))
cat(sprintf("Required Resolving Power for Baseline Separation: R >= %d\n", req_resolving_power))
cat(sprintf("Can Quadrupole separate them? %s (R = 1,000 < %d)\n",
            if (1000 >= req_resolving_power) "YES" else "NO", req_resolving_power))
cat(sprintf("Can TOF separate them?        %s (R = 30,000 >= %d)\n\n",
            if (30000 >= req_resolving_power) "YES" else "NO", req_resolving_power))

# -----------------------------------------------------------------------------
# Section 3: Stochastic Precursor Sampling in DDA vs DIA
# -----------------------------------------------------------------------------
cat("--- Section 3: Stochastic Sampling Simulation in DDA vs DIA ---\n")

set.seed(42)
n_precursors <- 300 # Precursors eluting in a 1-minute window
n_replicates <- 5   # 5 technical replicate runs

# True abundance spanning 4 orders of magnitude (log-normal)
true_abundance <- rlnorm(n_precursors, meanlog = 12, sdlog = 2.0)
true_ranks <- rank(-true_abundance)

# Top15 DDA Selection: Top 15 highest intensity ions selected per cycle
# Lower abundance ions suffer stochastic selection driven by ion arrival Poisson noise
dda_selection_matrix <- matrix(FALSE, nrow = n_precursors, ncol = n_replicates)
for (rep_i in 1:n_replicates) {
  # Add 15% technical sampling variance to intensity
  observed_intensity <- true_abundance * rlnorm(n_precursors, 0, 0.15)
  top_15_idx <- order(observed_intensity, decreasing = TRUE)[1:25] # Top 25 selected
  dda_selection_matrix[top_15_idx, rep_i] <- TRUE
}

# Missingness in DDA
dda_detections <- rowSums(dda_selection_matrix)
dda_missing_rate <- mean(dda_detections == 0) * 100
dda_full_reprod  <- mean(dda_detections == n_replicates) * 100

cat(sprintf("Simulated DDA Run (%d precursors, Top25 selection, %d replicates):\n",
            n_precursors, n_replicates))
cat(sprintf("  Precursors Detected in All 5 Replicates: %.1f%%\n", dda_full_reprod))
cat(sprintf("  Precursors Completely Missing (0 of 5):   %.1f%%\n", dda_missing_rate))
cat("DIA Mode Comparison:\n")
cat("  - Precursors Interrogated in All 5 Replicates: 100.0% (Zero stochastic missingness)\n")
cat("  - Trade-off: DIA spectra are chimeric and require deconvolution / spectral libraries.\n\n")

# -----------------------------------------------------------------------------
# Section 4: Ion Suppression Modeling in ESI
# -----------------------------------------------------------------------------
cat("--- Section 4: ESI Ion Suppression Kinetics ---\n")

# Analyte concentration range: 1 to 100 nM
conc <- c(1, 5, 10, 25, 50, 100)
# Clean solvent response: Linear calibration slope = 10,000 counts/nM
resp_clean <- 10000 * conc + rnorm(length(conc), 0, 5000)

# Plasma matrix with co-eluting phospholipid: suppresses ionization by 65%
suppression_factor <- 0.35 # Only 35% of ion beam reaches detector
resp_matrix <- 10000 * suppression_factor * conc + rnorm(length(conc), 0, 3000)

fit_clean  <- lm(resp_clean ~ conc)
fit_matrix <- lm(resp_matrix ~ conc)

cat(sprintf("Clean Solvent Calibration Slope: %.0f counts/nM (R2 = %.4f)\n",
            coef(fit_clean)[2], summary(fit_clean)$r.squared))
cat(sprintf("Matrix Suppression Slope:        %.0f counts/nM (R2 = %.4f)\n",
            coef(fit_matrix)[2], summary(fit_matrix)$r.squared))
cat(sprintf("Apparent Matrix Effect:          %.1f%% of clean signal (%.1f%% Signal Loss)\n\n",
            coef(fit_matrix)[2] / coef(fit_clean)[2] * 100,
            (1 - coef(fit_matrix)[2] / coef(fit_clean)[2]) * 100))

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 02 Exercises\n")
cat("===================================================================\n\n")

# Exercise 2.1: Resolving Charge State and Mass from Observed m/z
cat("--- Exercise 2.1: Resolving Charge State and Neutral Mass ---\n")
# Triply-charged peptide observed with monoisotopic peak at m/z = 624.8348
obs_mz0 <- 624.8348
obs_mz1 <- 625.1693
obs_delta <- obs_mz1 - obs_mz0
calc_z <- round(c13_diff / obs_delta)
calc_mass <- calc_z * obs_mz0 - calc_z * m_proton

cat(sprintf("Observed Monoisotopic Peak m/z: %.4f\n", obs_mz0))
cat(sprintf("Observed First 13C Isotope m/z: %.4f\n", obs_mz1))
cat(sprintf("Measured Peak Spacing Delta(m/z): %.4f Da\n", obs_delta))
cat(sprintf("Inferred Charge State: z = %d\n", calc_z))
cat(sprintf("Deconvoluted Neutral Mass: M = %.4f Da\n\n", calc_mass))

# Exercise 2.2: Isobaric Separation Resolving Power
cat("--- Exercise 2.2: Resolving Power for Isobaric Metabolites ---\n")
# Leucine metabolite vs Asparagine derivative:
# Nominal mass = 132 Da; exact masses: 132.0773 vs 132.1025 Da
m1 <- 132.0773
m2 <- 132.1025
dm <- m2 - m1
r_needed <- round(m1 / dm)

cat(sprintf("Metabolite 1: %.4f Da | Metabolite 2: %.4f Da\n", m1, m2))
cat(sprintf("Mass Difference Delta(m): %.4f Da (%.1f mDa)\n", dm, dm * 1000))
cat(sprintf("Minimum Resolving Power Required: R >= %d\n", r_needed))
cat("Instrumentation Verdict: Quadrupole (R=1000) fails; TOF or Orbitrap (R >= 10,000) resolves the doublet.\n\n")

# Exercise 2.3: DDA Stochastic Sampling and Missing Values
cat("--- Exercise 2.3: Stochastic Sampling and Missing Data Mechanism ---\n")
cat("In DDA, precursors are selected based on instantaneous MS1 intensity:\n")
cat("  - Abundant peptides: Consistently in the TopN candidates across all runs (Observed in 100% replicates).\n")
cat("  - Low-abundance peptides: Hover near the dynamic exclusion / selection threshold.\n")
cat("    Poisson counting fluctuations determine whether they are fragmented in Run A or dropped in Run B.\n")
cat("Classification: This missingness is Missing Not At Random (MNAR) / Censored at low abundances,\n")
cat("accompanied by stochastic Missing at Random (MAR) near the selection boundary.\n\n")

# Exercise 2.4: Internal Standards for Ion Suppression Correction
cat("--- Exercise 2.4: Stable Isotope Labeled (SIL) Internal Standards ---\n")
cat("Ion suppression attenuates the ESI droplet evaporation rate equally for an analyte and its SIL analog:\n")
cat("  - Analyte signal with 60% suppression:  I_analyte = 0.40 * k * C_analyte\n")
cat("  - SIL standard with 60% suppression:   I_SIL     = 0.40 * k * C_SIL\n")
cat("  - Measured Area Ratio:                 Ratio     = I_analyte / I_SIL = C_analyte / C_SIL\n")
cat("Mathematical Conclusion: The suppression coefficient (0.40) cancels algebraically,\n")
cat("preserving accurate absolute quantification despite severe matrix effects.\n\n")

cat("===================================================================\n")
cat("Chapter 02 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
