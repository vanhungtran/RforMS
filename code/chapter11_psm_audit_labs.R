# =============================================================================
# Chapter 11: Peptide-Spectrum Match Auditing and Target-Decoy FDR
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(PSMatch)
  library(msdata)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 11.1: Loading mzIdentML Data and Inspecting PSM Metadata
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 11.1: LOADING MZIDENTML DATA ===\n")

# Locate bundled Erwinia proteome Top10-HCD mzIdentML file
id_file <- msdata::ident(full.names = TRUE)
cat(sprintf("Source mzIdentML file: %s\n", basename(id_file)))

# Read PSM dataset
psm_raw <- readPSMs(id_file)
cat(sprintf("Initial PSM dimensions: %d rows x %d columns\n", nrow(psm_raw), ncol(psm_raw)))

# Target vs Decoy counts in raw search results
target_decoy_counts <- table(psm_raw$isDecoy)
cat(sprintf("Target PSMs: %d (%.1f%%) | Decoy PSMs: %d (%.1f%%)\n",
            target_decoy_counts["FALSE"],
            target_decoy_counts["FALSE"] / nrow(psm_raw) * 100,
            target_decoy_counts["TRUE"],
            target_decoy_counts["TRUE"] / nrow(psm_raw) * 100))


# -----------------------------------------------------------------------------
# Listing 11.2: Empirical Target-Decoy FDR Evaluation
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 11.2: TARGET-DECOY FDR EVALUATION ===\n")

# Filter to rank 1 matches before FDR evaluation
psm_r1 <- filterPsmRank(psm_raw)

# Calculate empirical FDR across score thresholds
fdr_at_score <- function(thresh, psm_df) {
  above <- psm_df$MS.GF.RawScore >= thresh
  n_d <- sum(above & psm_df$isDecoy)
  n_t <- sum(above & !psm_df$isDecoy)
  fdr <- if (n_t == 0) 0 else n_d / n_t
  c(Score_Threshold = thresh, N_Target = n_t, N_Decoy = n_d, FDR = fdr)
}

score_grid <- seq(0, 100, by = 10)
fdr_table <- as.data.frame(t(sapply(score_grid, fdr_at_score, psm_df = psm_r1)))
fdr_table$FDR_Pct <- round(fdr_table$FDR * 100, 2)

cat("Empirical Target-Decoy FDR across MS-GF+ RawScore thresholds:\n")
print(fdr_table[, c("Score_Threshold", "N_Target", "N_Decoy", "FDR_Pct")])

# Exact score threshold achieving FDR < 1.0%
fine_scores <- seq(70, 100, by = 1)
fine_fdr <- as.data.frame(t(sapply(fine_scores, fdr_at_score, psm_df = psm_r1)))
threshold_1pct <- min(fine_fdr$Score_Threshold[fine_fdr$FDR <= 0.01])
cat(sprintf("\nStrict 1.0%% FDR Score Threshold: MS-GF+ RawScore >= %.0f (Retained Targets: %d, Decoys: %d)\n",
            threshold_1pct,
            fine_fdr$N_Target[fine_fdr$Score_Threshold == threshold_1pct],
            fine_fdr$N_Decoy[fine_fdr$Score_Threshold == threshold_1pct]))


# -----------------------------------------------------------------------------
# Listing 11.3: Stepwise Filtering with PSMatch
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 11.3: STEPWISE FILTERING WITH PSMATCH ===\n")

n_initial <- nrow(psm_raw)

# Step 1: Filter to rank 1
psm_step1 <- filterPsmRank(psm_raw)
n_rank1 <- nrow(psm_step1)
cat(sprintf("Step 1 (Rank 1 filter):     Retained %d PSMs (removed %d suboptimal matches)\n",
            n_rank1, n_initial - n_rank1))

# Step 2: Remove decoy hits
psm_step2 <- filterPsmDecoy(psm_step1)
n_target <- nrow(psm_step2)
cat(sprintf("Step 2 (Decoy removal):     Retained %d target PSMs (removed %d decoy hits)\n",
            n_target, n_rank1 - n_target))

# Step 3: Remove shared (degenerate) peptides
psm_step3 <- filterPsmShared(psm_step2)
n_unique <- nrow(psm_step3)
cat(sprintf("Step 3 (Shared removal):    Retained %d unique-protein PSMs (removed %d shared)\n",
            n_unique, n_target - n_unique))


# -----------------------------------------------------------------------------
# Listing 11.4: Systematic Failure Mode Diagnostics
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 11.4: SYSTEMATIC FAILURE MODE DIAGNOSTICS ===\n")

# A. Missed Cleavage Detection (Trypsin rule: K/R not followed by P)
count_missed_cleavages <- function(seq) {
  sites <- gregexpr("[KR](?!P)", seq, perl = TRUE)[[1]]
  n <- if (sites[1] == -1) 0L else length(sites)
  if (substr(seq, nchar(seq), nchar(seq)) %in% c("K", "R")) n <- n - 1L
  max(n, 0L)
}

psm_step2$missed_cleavages <- sapply(psm_step2$sequence, count_missed_cleavages)
mc_table <- table(psm_step2$missed_cleavages)
cat("Missed Cleavages Distribution (Target PSMs):\n")
print(mc_table)

# B. Precursor Charge-State Distribution
charge_table <- table(psm_step2$chargeState)
cat("\nPrecursor Charge-State Distribution (Target PSMs):\n")
print(charge_table)

# C. Common Contaminants Audit (Trypsin, Keratin, BSA)
contam_indices <- grep("trypsin|keratin|albumin", psm_step2$DatabaseDescription, ignore.case = TRUE)
cat(sprintf("\nContaminant PSMs Detected: %d hits (%.2f%% of targets)\n",
            length(contam_indices), length(contam_indices) / nrow(psm_step2) * 100))
print(table(psm_step2$DatabaseDescription[contam_indices]))


# -----------------------------------------------------------------------------
# Listing 11.5: Assembling Clean Identification Report
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 11.5: ASSEMBLING IDENTIFICATION REPORT ===\n")

# Compute missed cleavages for confident unique PSMs
psm_step3$missed_cleavages <- sapply(psm_step3$sequence, count_missed_cleavages)

# Combine confident unique PSMs into report data frame
ident_report <- as.data.frame(psm_step3) %>%
  select(
    Spectrum = spectrumID,
    Sequence = sequence,
    Charge = chargeState,
    Protein = DatabaseAccess,
    Score = MS.GF.RawScore,
    QValue = MS.GF.QValue,
    Missed_Cleavages = missed_cleavages
  ) %>%
  arrange(desc(Score))

cat(sprintf("Clean Identification Report Summary:\n"))
cat(sprintf("  Confident PSMs:    %d\n", nrow(ident_report)))
cat(sprintf("  Unique Peptides:   %d\n", length(unique(ident_report$Sequence))))
cat(sprintf("  Distinct Proteins: %d\n", length(unique(ident_report$Protein))))
cat("\nTop 5 Confident Identifications by MS-GF+ Score:\n")
print(head(ident_report, 5))


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 11.1: SpecEValue vs RawScore FDR thresholding
# Lower SpecEValue is more significant
fdr_at_evalue <- function(eval_thresh, psm_df) {
  below <- psm_df$MS.GF.SpecEValue <= eval_thresh
  n_d <- sum(below & psm_df$isDecoy)
  n_t <- sum(below & !psm_df$isDecoy)
  fdr <- if (n_t == 0) 0 else n_d / n_t
  c(EValue = eval_thresh, N_Target = n_t, N_Decoy = n_d, FDR = fdr)
}
eval_candidates <- 10^seq(-15, -5, by = 1)
eval_res <- as.data.frame(t(sapply(eval_candidates, fdr_at_evalue, psm_df = psm_r1)))
eval_1pct <- max(eval_res$EValue[eval_res$FDR <= 0.01])
cat(sprintf("Exercise 11.1: SpecEValue 1.0%% FDR cutoff: SpecEValue <= %.2e (Targets: %d, Decoys: %d)\n",
            eval_1pct,
            eval_res$N_Target[eval_res$EValue == eval_1pct],
            eval_res$N_Decoy[eval_res$EValue == eval_1pct]))

# Exercise 11.2: Association between missed cleavages and precursor charge
mean_charge_by_mc <- aggregate(chargeState ~ missed_cleavages, data = psm_step2, FUN = mean)
cat("\nExercise 11.2: Mean Precursor Charge by Missed Cleavages:\n")
print(round(mean_charge_by_mc, 3))

# Exercise 11.3: Order of Operations Sensitivity
# Shared -> Decoy -> Rank
psm_alt1 <- filterPsmShared(psm_raw)
psm_alt2 <- filterPsmDecoy(psm_alt1)
psm_alt3 <- filterPsmRank(psm_alt2)
cat(sprintf("\nExercise 11.3: Filter Order Sensitivity:\n"))
cat(sprintf("  Standard Order (Rank -> Decoy -> Shared): %d PSMs\n", nrow(psm_step3)))
cat(sprintf("  Alternative Order (Shared -> Decoy -> Rank): %d PSMs (Delta = %d)\n",
            nrow(psm_alt3), nrow(psm_step3) - nrow(psm_alt3)))

# Exercise 11.4: Connected Components Topology
cc <- ConnectedComponents(psm_step2)
adj <- adjacencyMatrix(cc)
cat(sprintf("\nExercise 11.4: Protein-Peptide Graph Topology:\n"))
cat(sprintf("  Total Peptides in Network: %d\n", nrow(adj)))
cat(sprintf("  Total Proteins in Network: %d\n", ncol(adj)))
cat(sprintf("  Total Connected Components: %d\n", length(cc)))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
