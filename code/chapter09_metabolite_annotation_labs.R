# =============================================================================
# Chapter 09: Metabolite Annotation and Adduct Deconvolution
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(MetaboCoreUtils)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 9.1: Exact Mass and Adduct Inversion Arithmetic
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 9.1: EXACT MASS AND ADDUCT ARITHMETIC ===\n")

# Neutral formulas for key biological metabolites
metabolites <- c(
  Glucose     = "C6H12O6",
  Dopamine    = "C8H11NO2",
  Caffeine    = "C8H10N4O2",
  Tryptophan  = "C11H12N2O2"
)

# Calculate exact neutral monoisotopic masses
exact_masses <- MetaboCoreUtils::calculateMass(metabolites)

# Positive-mode adducts to evaluate
target_adducts_pos <- c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M+H-H2O]+")

# Compute theoretical m/z for Dopamine across adducts
dopamine_mz <- MetaboCoreUtils::mass2mz(exact_masses["Dopamine"], adduct = target_adducts_pos)

adduct_summary <- data.frame(
  Compound       = rep("Dopamine (C8H11NO2)", length(target_adducts_pos)),
  Neutral_Mass   = rep(unname(round(exact_masses["Dopamine"], 4)), length(target_adducts_pos)),
  Adduct         = target_adducts_pos,
  Theoretical_mz = round(as.numeric(dopamine_mz[1, ]), 4),
  stringsAsFactors = FALSE
)

# Demonstrate neutral mass inversion from observed [M+Na]+ ion
obs_mz <- 176.0682   # Observed Dopamine [M+Na]+
inv_mass <- MetaboCoreUtils::mz2mass(obs_mz, adduct = "[M+Na]+")
mass_err_ppm <- ((inv_mass - exact_masses["Dopamine"]) / exact_masses["Dopamine"]) * 1e6

cat("Dopamine Electrospray Adduct Profile (ESI+):\n")
print(adduct_summary)
cat(sprintf("\nAdduct Inversion Audit: Observed [M+Na]+ = %.4f Da\n", obs_mz))
cat(sprintf("  Inferred Neutral Mass: %.4f Da | True Mass: %.4f Da | Delta: %+.4f Da (%+.2f ppm)\n",
            inv_mass, exact_masses["Dopamine"], inv_mass - exact_masses["Dopamine"], mass_err_ppm))


# -----------------------------------------------------------------------------
# Listing 9.2: Peak Shape Correlation and Multiplet Deconvolution
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 9.2: PEAK SHAPE CORRELATION DECONVOLUTION ===\n")

# Simulate chromatographic time series across elution window (t = 2.0 to 2.8 min)
time_scans <- seq(2.0, 2.8, by = 0.02)
rt_apex <- 2.40
peak_sigma <- 0.06

# Simulate co-eluting features from a single metabolite:
# Feature 1: Monoisotopic [M+H]+ (m/z = 181.0707)
# Feature 2: In-source water loss [M+H-H2O]+ (m/z = 163.0601)
# Feature 3: Sodium adduct [M+Na]+ (m/z = 203.0526)
# Feature 4: Unrelated co-eluting contaminant peak (rt_apex = 2.44)
int_mh   <- 120000 * exp(-0.5 * ((time_scans - rt_apex) / peak_sigma)^2) + rnorm(length(time_scans), 0, 800)
int_loss <-  24000 * exp(-0.5 * ((time_scans - rt_apex) / peak_sigma)^2) + rnorm(length(time_scans), 0, 300)
int_mna  <-  48000 * exp(-0.5 * ((time_scans - (rt_apex + 0.005)) / peak_sigma)^2) + rnorm(length(time_scans), 0, 500)
int_cont <-  35000 * exp(-0.5 * ((time_scans - 2.46) / peak_sigma)^2) + rnorm(length(time_scans), 0, 400)

eic_matrix <- cbind(
  FT_181_MH   = pmax(0, int_mh),
  FT_163_Loss = pmax(0, int_loss),
  FT_203_MNa  = pmax(0, int_mna),
  FT_unrel    = pmax(0, int_cont)
)

cor_mat <- round(cor(eic_matrix), 3)
cat("Chromatographic Peak Shape Correlation Matrix (Pearson r):\n")
print(cor_mat)

# Deconvolution grouping threshold: r >= 0.85
cluster_members <- colnames(eic_matrix)[cor_mat["FT_181_MH", ] >= 0.85]
cat(sprintf("\nDeconvoluted Pseudo-Spectrum Group 1 (r >= 0.85 with [M+H]+):\n"))
print(cluster_members)


# -----------------------------------------------------------------------------
# Listing 9.3: Tolerance-Aware Database Query and PPM Error Control
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 9.3: TOLERANCE-AWARE DATABASE MATCHING ===\n")

# Curated reference database snippet
ref_database <- data.frame(
  Compound_Name = c("Glucose", "Fructose", "myo-Inositol", "Dopamine", "Caffeine", "Citric acid"),
  Formula       = c("C6H12O6", "C6H12O6", "C6H12O6", "C8H11NO2", "C8H10N4O2", "C6H8O7"),
  Exact_Mass    = c(180.06339, 180.06339, 180.06339, 153.07898, 194.08038, 192.02700),
  stringsAsFactors = FALSE
)

# Experimental observed features (neutral mass deconvoluted)
query_features <- data.frame(
  Feature_ID     = c("FT_001", "FT_002", "FT_003"),
  Observed_Mass  = c(180.0639, 153.0784, 312.1245),
  Retention_Time = c(2.45, 1.82, 8.40),
  stringsAsFactors = FALSE
)

# Search with tolerance tau = 5.0 ppm
ppm_tolerance <- 5.0

matches_list <- list()
for (i in seq_len(nrow(query_features))) {
  obs_m <- query_features$Observed_Mass[i]
  delta_limit <- obs_m * (ppm_tolerance / 1e6)
  candidates <- ref_database %>%
    filter(abs(Exact_Mass - obs_m) <= delta_limit)
  
  if (nrow(candidates) > 0) {
    candidates$Feature_ID     <- query_features$Feature_ID[i]
    candidates$Observed_Mass  <- obs_m
    candidates$Retention_Time <- query_features$Retention_Time[i]
    candidates$PPM_Error      <- round(abs(candidates$Exact_Mass - obs_m) / candidates$Exact_Mass * 1e6, 2)
    matches_list[[length(matches_list) + 1]] <- candidates
  }
}

matched_df <- do.call(rbind, matches_list) %>%
  select(Feature_ID, Observed_Mass, Retention_Time, Compound_Name, Formula, PPM_Error)

cat(sprintf("Exact Mass Database Search Results (Tolerance = +/- %.1f ppm):\n", ppm_tolerance))
print(matched_df)


# -----------------------------------------------------------------------------
# Listing 9.4: High-Resolution Isotopic Fine Structure Modeling
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 9.4: ISOTOPIC FINE STRUCTURE MODELING ===\n")

# Contrast Palmitic Acid (C16H32O2) vs isobaric peptide fragment with 1 Sulfur:
# Palmitic acid: C16 H32 O2 (Exact mass = 256.24023 Da)
# Isobaric model: C12 H20 N2 O2 S (Exact mass = 256.12455 Da - nominal match at m/z 257)

# Theoretical M+2 fine structure shifts:
# Carbon-13 double substitution (13C2): Delta m = 2 * 1.003355 = 2.006710 Da
# Sulfur-34 substitution (34S):        Delta m = 1.995796 Da
fine_mass_diff <- abs(2.006710 - 1.995796)
resolving_power_req <- 257 / fine_mass_diff

cat("High-Resolution Isotopic Fine Structure Mass Differences:\n")
cat(sprintf("  13C2 Isotope Mass Shift:  +%.6f Da\n", 2.006710))
cat(sprintf("  34S  Isotope Mass Shift:  +%.6f Da\n", 1.995796))
cat(sprintf("  Mass Separation (Delta m): %.6f Da (%.2f mDa)\n", fine_mass_diff, fine_mass_diff * 1000))
cat(sprintf("  Resolving Power Required to baseline-resolve 34S from 13C2 at m/z 257: R >= %.0f\n",
            resolving_power_req))


# -----------------------------------------------------------------------------
# Listing 9.5: Structured Annotation Assembly and MSI Level Classification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 9.5: ANNOTATION ASSEMBLY AND MSI CLASSIFICATION ===\n")

# Aggregate multiple candidate annotations into unambiguous report
annotation_report <- query_features %>%
  left_join(matched_df, by = c("Feature_ID", "Observed_Mass", "Retention_Time")) %>%
  group_by(Feature_ID, Observed_Mass, Retention_Time) %>%
  summarise(
    Matched_Count  = sum(!is.na(Compound_Name)),
    Candidates     = if (sum(!is.na(Compound_Name)) > 0) paste(Compound_Name, collapse = " | ") else "Unmatched",
    Formulas       = if (sum(!is.na(Formula)) > 0) paste(unique(Formula), collapse = " | ") else "None",
    Min_PPM_Error  = if (sum(!is.na(PPM_Error)) > 0) min(PPM_Error) else NA_real_,
    .groups        = "drop"
  ) %>%
  mutate(
    MSI_Confidence = case_when(
      Matched_Count == 0 ~ "Level 4 (Unidentified Feature)",
      Matched_Count == 1 ~ "Level 3 (Putatively Characterized Class)",
      Matched_Count > 1  ~ "Level 3 (Ambiguous Isomeric Set)"
    )
  )

cat("Final Structured Metabolite Annotation Report:\n")
print(as.data.frame(annotation_report))


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 9.1: Dopamine adduct arithmetic
dopamine_mass <- MetaboCoreUtils::calculateMass("C8H11NO2")
dop_mh_pos   <- MetaboCoreUtils::mass2mz(dopamine_mass, adduct = "[M+H]+")
dop_mna_pos  <- MetaboCoreUtils::mass2mz(dopamine_mass, adduct = "[M+Na]+")
dop_loss_pos <- MetaboCoreUtils::mass2mz(dopamine_mass, adduct = "[M+H-H2O]+")
dop_mh_neg   <- MetaboCoreUtils::mass2mz(dopamine_mass, adduct = "[M-H]-")

cat(sprintf("Exercise 9.1: Dopamine Neutral Mass = %.4f Da\n", dopamine_mass))
cat(sprintf("  [M+H]+:     %.4f Da\n", dop_mh_pos))
cat(sprintf("  [M+Na]+:    %.4f Da\n", dop_mna_pos))
cat(sprintf("  [M+H-H2O]+: %.4f Da\n", dop_loss_pos))
cat(sprintf("  [M-H]-:     %.4f Da\n", dop_mh_neg))

# Exercise 9.2: PPM Error
obs_sucrose <- 341.1082
theo_sucrose <- 341.1078
sucrose_err_da <- obs_sucrose - theo_sucrose
sucrose_err_ppm <- (sucrose_err_da / theo_sucrose) * 1e6
cat(sprintf("\nExercise 9.2: Sucrose Error = %+.4f Da (%.2f ppm) -> Within 5 ppm: %s\n",
            sucrose_err_da, sucrose_err_ppm, ifelse(abs(sucrose_err_ppm) <= 5.0, "YES", "NO")))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
