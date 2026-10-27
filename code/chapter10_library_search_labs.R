# =============================================================================
# Chapter 10: Tandem Spectral Library Searching
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(Spectra)
  library(MsCoreUtils)
  library(MetaboAnnotation)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 10.1: Constructing Query and Library Spectra Objects
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 10.1: CONSTRUCTING SPECTRA OBJECTS ===\n")

# Query spectrum (observed experimental MS/MS for unknown analyte)
# Precursor: m/z = 205.0972 [M+H]+ (Tryptophan candidate)
query_df <- DataFrame(
  msLevel          = 2L,
  precursorMz      = 205.0972,
  precursorCharge  = 1L,
  rtime            = 312.4, # seconds (5.20 min)
  name             = "Query_Feature_205"
)
query_df$mz <- list(c(118.0651, 130.0651, 144.0808, 146.0600, 159.0917, 188.0706, 205.0972))
query_df$intensity <- list(c(12400, 35600, 89200, 15300, 120000, 48200, 18500))

query_sps <- Spectra(query_df)

# Reference library spectra: Curated standards
ref_df <- DataFrame(
  msLevel         = c(2L, 2L, 2L),
  precursorMz     = c(205.0972, 205.0972, 182.0812),
  precursorCharge = c(1L, 1L, 1L),
  rtime           = c(310.2, 425.0, 180.5), # Tryptophan authentic standard rtime = 310.2 s (5.17 min)
  name            = c("L-Tryptophan", "D-Tryptophan_Isomer", "Tyrosine")
)
ref_df$mz <- list(
  # L-Tryptophan reference MS/MS
  c(118.0653, 130.0651, 144.0808, 146.0601, 159.0918, 188.0707, 205.0970),
  # Near-isomer with shifted fragment intensities
  c(118.0653, 130.0651, 144.0808, 170.0600, 188.0707, 205.0970),
  # Tyrosine reference MS/MS (off-target precursor)
  c(119.0491, 136.0757, 147.0441, 165.0546, 182.0812)
)
ref_df$intensity <- list(
  c(15000, 42000, 95000, 18000, 135000, 52000, 22000),
  c(45000, 12000, 15000, 85000, 21000, 35000),
  c(25000, 98000, 45000, 15000, 82000)
)
ref_sps <- Spectra(ref_df)

cat(sprintf("Constructed Query: %d spectrum | Reference Library: %d spectra\n",
            length(query_sps), length(ref_sps)))
cat(sprintf("Query Precursor: %.4f Da | Library Precursors: %s\n",
            precursorMz(query_sps)[1],
            paste(round(precursorMz(ref_sps), 4), collapse = ", ")))


# -----------------------------------------------------------------------------
# Listing 10.2: Normalized Dot-Product (Cosine Similarity) Evaluation
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 10.2: COSINE SIMILARITY EVALUATION ===\n")

# Calculate pairwise cosine similarity between query and all reference spectra
# Parameters: fragment tolerance = 0.02 Da, square-root intensity weighting (b = 0.5)
sim_scores <- compareSpectra(
  query_sps, ref_sps,
  tolerance = 0.02,
  ppm = 10,
  FUN = MsCoreUtils::ndotproduct
)

scores_df <- data.frame(
  Library_Compound = ref_df$name,
  Precursor_mz     = ref_df$precursorMz,
  RT_Ref_sec       = ref_df$rtime,
  Cosine_Score     = round(as.numeric(sim_scores), 4),
  stringsAsFactors = FALSE
)

cat("Pairwise Spectral Cosine Similarities against Query:\n")
print(scores_df)


# -----------------------------------------------------------------------------
# Listing 10.3: High-Throughput Matching with MetaboAnnotation
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 10.3: MATCHING WITH METABOANNOTATION ===\n")

# Define matching parameter: 15 ppm precursor tolerance, 0.02 Da fragment tolerance, cosine >= 0.70
param <- CompareSpectraParam(
  ppm              = 15,
  tolerance        = 0.02,
  requirePrecursor = TRUE,
  score            = "ndotproduct",
  THRESHFUN        = function(x) which(x >= 0.70)
)

match_res <- matchSpectra(query_sps, ref_sps, param = param)
matched_df <- matchedData(match_res)

cat(sprintf("MetaboAnnotation match results (Threshold >= 0.70):\n"))
cat(sprintf("  Query spectrum matched to: %s\n", matched_df$target_name))
cat(sprintf("  Cosine score: %.4f | Precursor ppm error: %.2f ppm\n",
            matched_df$score,
            abs(matched_df$query_precursorMz - matched_df$target_precursorMz) / matched_df$target_precursorMz * 1e6))


# -----------------------------------------------------------------------------
# Listing 10.4: Decoy Spectral Library and Empirical Null Distribution
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 10.4: EMPIRICAL NULL AND P-VALUE CALCULATION ===\n")

# Generate 1,000 decoy spectra by mass-shuffling reference peaks
n_decoys <- 1000
decoy_scores <- numeric(n_decoys)
true_peaks <- peaksData(ref_sps)[[1]]  # L-Tryptophan peaks

for (k in seq_len(n_decoys)) {
  # Shuffle fragment m/z values uniformly across mass range [100, 210]
  shuffled_mz <- sort(runif(nrow(true_peaks), min = 100, max = 210))
  decoy_df <- DataFrame(msLevel = 2L, precursorMz = 205.0972)
  decoy_df$mz <- list(shuffled_mz)
  decoy_df$intensity <- list(true_peaks[, "intensity"])
  decoy_sp <- Spectra(decoy_df)
  decoy_scores[k] <- compareSpectra(query_sps, decoy_sp, tolerance = 0.02, FUN = MsCoreUtils::ndotproduct)[1]
}

# Empirical match p-value for the true hit (score = 0.9996)
obs_score <- scores_df$Cosine_Score[1]
p_val_empirical <- sum(decoy_scores >= obs_score) / n_decoys
score_99th <- quantile(decoy_scores, 0.99)

cat(sprintf("Decoy Null Model (N = %d shuffled spectra):\n", n_decoys))
cat(sprintf("  Mean Null Cosine Score:    %.4f (SD = %.4f)\n", mean(decoy_scores), sd(decoy_scores)))
cat(sprintf("  99th Percentile Threshold: %.4f\n", score_99th))
cat(sprintf("  Observed Hit Cosine Score: %.4f\n", obs_score))
cat(sprintf("  Empirical False Match p-value: p < %.4f\n", max(1 / n_decoys, p_val_empirical)))


# -----------------------------------------------------------------------------
# Listing 10.5: Multi-Criteria Filtering and MSI Level Assignment
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 10.5: MULTI-CRITERIA MSI CLASSIFICATION ===\n")

# Combine Cosine score, Precursor delta, and Retention Time deviation
query_rt <- query_df$rtime[1]

curated_report <- scores_df %>%
  mutate(
    Precursor_PPM = round(abs(Precursor_mz - query_df$precursorMz[1]) / Precursor_mz * 1e6, 2),
    Delta_RT_sec  = round(abs(RT_Ref_sec - query_rt), 1),
    MSI_Level     = case_when(
      Cosine_Score >= 0.85 & Precursor_PPM <= 10 & Delta_RT_sec <= 10 ~ "Level 1 (Confirmed via Standard)",
      Cosine_Score >= 0.80 & Precursor_PPM <= 15                      ~ "Level 2 (Putatively Annotated Compound)",
      Cosine_Score >= 0.60 & Precursor_PPM <= 25                      ~ "Level 3 (Putative Compound Class)",
      TRUE                                                            ~ "Level 4 (Unassigned Feature)"
    )
  )

cat("Curated Spectral Match and MSI Classification Report:\n")
print(curated_report)


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 10.1: Normalized Dot Product Manual Derivation
# Query: m/z 100 (int 100), 150 (int 400), 200 (int 900)
# Target: m/z 100 (int 144), 150 (int 256), 200 (int 1024)
# Square-root intensity weights: u = c(10, 20, 30), v = c(12, 16, 32)
u <- c(10, 20, 30)
v <- c(12, 16, 32)
dot_prod <- sum(u * v)
norm_u <- sqrt(sum(u^2))
norm_v <- sqrt(sum(v^2))
cos_val <- dot_prod / (norm_u * norm_v)

cat(sprintf("Exercise 10.1: Manual Cosine Similarity:\n"))
cat(sprintf("  Numerator (dot product):   %d\n", dot_prod))
cat(sprintf("  Denominator: sqrt(%d) * sqrt(%d) = %.4f * %.4f = %.4f\n",
            sum(u^2), sum(v^2), norm_u, norm_v, norm_u * norm_v))
cat(sprintf("  Calculated Cosine Score:   %.4f\n", cos_val))

# Exercise 10.2: Hypergeometric match probability
# N = 1000 bins, m = 8 library peaks, n = 6 query peaks, k >= 4 matches
p_hyper <- sum(dhyper(4:6, m = 8, n = 1000 - 8, k = 6))
cat(sprintf("\nExercise 10.2: Hypergeometric Probability (k >= 4 of 6 matches in 1000 bins):\n"))
cat(sprintf("  P(K >= 4 | H0) = %.4e (p < 0.0001)\n", p_hyper))

# Exercise 10.3: Reverse Cosine with unassigned contaminant peaks
# Query has 3 shared peaks plus 2 intense contaminant peaks
# Forward dot product dilutes denominator, reverse dot product considers only library peaks
w_lib <- c(12, 16, 32) # library peaks
w_qry_shared <- c(10, 20, 30)
w_qry_all <- c(10, 20, 30, 40, 50) # includes 2 contaminant peaks (int 1600, 2500 -> sqrt 40, 50)
cos_forward <- sum(w_qry_shared * w_lib) / (sqrt(sum(w_qry_all^2)) * sqrt(sum(w_lib^2)))
cos_reverse <- sum(w_qry_shared * w_lib) / (sqrt(sum(w_qry_shared^2)) * sqrt(sum(w_lib^2)))
cat(sprintf("\nExercise 10.3: Forward vs Reverse Cosine under Co-eluting Contamination:\n"))
cat(sprintf("  Forward Cosine: %.4f (severely penalized by contaminants)\n", cos_forward))
cat(sprintf("  Reverse Cosine: %.4f (immune to unassigned query peaks)\n", cos_reverse))

# Exercise 10.4: Precursor Filtering Speedup
n_total_lib <- 500000
ppm_win <- 20
frac_window <- (2 * ppm_win * 1e-6) # window fraction of mass spectrum assuming uniform distribution
n_filtered <- n_total_lib * frac_window * 5 # typical candidate density factor
cat(sprintf("\nExercise 10.4: Candidate Space Reduction:\n"))
cat(sprintf("  Full Library Comparisons: %d\n", n_total_lib))
cat(sprintf("  Window-Filtered Candidates (20 ppm window): ~%d (speedup > %dx)\n",
            as.integer(n_filtered), as.integer(n_total_lib / n_filtered)))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
