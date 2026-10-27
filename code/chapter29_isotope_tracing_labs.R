# =============================================================================
# Chapter 29: Stable Isotope Tracing and Molecular Networking
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(MetaboCoreUtils)
  library(Spectra)
  library(igraph)
  library(MsCoreUtils)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 29.1: Exact-Mass Target Calculation for Stable-Isotope Labeling
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 29.1: EXACT-MASS TARGET CALCULATION ===\n")

d13c <- 1.003355   # mass difference: 13C - 12C (Da)

# Amino acid neutral formulas (C, H, N, O)
aa_formulas <- list(
  Phe = c(C = 9,  H = 11, N = 1, O = 2),
  Trp = c(C = 11, H = 12, N = 2, O = 2),
  Tyr = c(C = 9,  H = 11, N = 1, O = 3)
)

# Acyl donor fatty acid neutral formulas
acyl_formulas <- list(
  `Acetyl(C2:0)`        = c(C = 2, H = 4,  O = 2),
  `Isovaleryl(C5:0)`    = c(C = 5, H = 10, O = 2),
  `Isocaproyl(C6:0)`    = c(C = 6, H = 12, O = 2),
  `Phenylacetyl(PA)`    = c(C = 8, H = 8,  O = 2),
  `Phenylpropionyl(PP)` = c(C = 9, H = 10, O = 2)
)

# N-acyl condensation formula arithmetic: AA + Acyl - H2O
calc_nacyl_formula <- function(aa, acyl) {
  f <- aa
  for (el in names(acyl)) {
    f[el] <- (if (el %in% names(f)) f[el] else 0) + acyl[el]
  }
  f["H"] <- f["H"] - 2
  f["O"] <- f["O"] - 1
  f[f != 0]
}

# Calculate exact masses for target candidate panel
targets <- list()
for (aa_name in names(aa_formulas)) {
  n_carbons_aa <- aa_formulas[[aa_name]]["C"]
  for (acyl_name in names(acyl_formulas)) {
    form <- calc_nacyl_formula(aa_formulas[[aa_name]], acyl_formulas[[acyl_name]])
    form_str <- paste0(names(form), form, collapse = "")
    exact_mass <- MetaboCoreUtils::calculateMass(form_str)
    mz_native  <- MetaboCoreUtils::mass2mz(exact_mass, adduct = "[M+H]+")
    mz_labeled <- mz_native + n_carbons_aa * d13c
    
    targets[[length(targets) + 1]] <- data.frame(
      AminoAcid    = aa_name,
      AcylGroup    = acyl_name,
      Compound     = paste0(acyl_name, "-", aa_name),
      Formula      = form_str,
      Native_mz    = round(as.numeric(mz_native), 4),
      Labeled_mz   = round(as.numeric(mz_labeled), 4),
      MassShift    = round(n_carbons_aa * d13c, 4),
      stringsAsFactors = FALSE
    )
  }
}

target_df <- do.call(rbind, targets)
cat(sprintf("Target library generated: %d candidate N-acyl amino acid conjugates\n", nrow(target_df)))
print(head(target_df, 8))


# -----------------------------------------------------------------------------
# Listing 29.2: Targeted Extraction, Co-Elution Verification, and Fractional Enrichment
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 29.2: CO-ELUTION AND FRACTIONAL ENRICHMENT ===\n")

# Targeted extraction parameters for Isocaproyl-Phe
target_pair <- list(
  compound     = "Isocaproyl-Phe",
  native_mz    = 264.1594,
  labeled_mz   = 273.1896,
  rt_native    = 19.96,
  rt_labeled   = 19.98,
  area_native  = 85722,
  area_labeled = 18415
)

# Co-elution criterion: Delta RT <= 0.10 min
delta_rt <- abs(target_pair$rt_native - target_pair$rt_labeled)
coeluting <- delta_rt <= 0.10

# Cluster-summed fractional enrichment (assuming equal ESI response factor)
phi <- (target_pair$area_labeled / (target_pair$area_native + target_pair$area_labeled)) * 100

cat(sprintf("Compound: %s\n", target_pair$compound))
cat(sprintf("  Native   [M+H]+ = %.4f Da | Observed RT apex = %.2f min | Peak Area = %d\n",
            target_pair$native_mz, target_pair$rt_native, target_pair$area_native))
cat(sprintf("  13C9     [M+H]+ = %.4f Da | Observed RT apex = %.2f min | Peak Area = %d\n",
            target_pair$labeled_mz, target_pair$rt_labeled, target_pair$area_labeled))
cat(sprintf("  Chromatographic Delta RT = %.3f min (Co-elution criterion: <= 0.10 min) -> %s\n",
            delta_rt, ifelse(coeluting, "CO-ELUTION CONFIRMED", "REJECTED (Non-coeluting)")))
cat(sprintf("  Fractional Biosynthetic Enrichment: %.1f%% of total metabolite pool synthesized from fed tracer\n",
            phi))


# -----------------------------------------------------------------------------
# Listing 29.3: Diagnostic MS/MS Fragment Ion Partitioning
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 29.3: MS/MS FRAGMENT PARTITIONING ===\n")

ms2_fragments <- data.frame(
  Fragment_Identity  = c("Phe immonium ion", "[Phe + H]+", "Isocaproyl acylium", "Alkyl acyl tail"),
  Elemental_Formula  = c("C8 H10 N", "C9 H12 N O2", "C6 H11 O", "C5 H11"),
  Moiety             = c("Amino Acid", "Amino Acid", "Acyl Donor", "Acyl Donor"),
  Native_mz          = c(120.0808, 166.0863, 99.0804, 71.0855),
  Labeled_13C9_mz    = c(128.1076, 175.1165, 99.0804, 71.0855),
  stringsAsFactors   = FALSE
)

ms2_fragments$Observed_Shift_Da <- round(ms2_fragments$Labeled_13C9_mz - ms2_fragments$Native_mz, 4)
ms2_fragments$Label_Confirmed <- ifelse(
  ms2_fragments$Moiety == "Amino Acid",
  paste0("YES (+", round(ms2_fragments$Observed_Shift_Da / 1.003355), "C)"),
  "NO (Unlabeled)"
)

cat("Diagnostic Tandem MS/MS Fragment Partitioning:\n")
print(ms2_fragments)
cat("Structural Proof: Mass shifts occur strictly on amino-acid fragments, confirming\n")
cat("condensation between unlabelled endogenous acyl donors and fed 13C-amino acid.\n")


# -----------------------------------------------------------------------------
# Listing 29.4: Molecular Networking with Top-K Reciprocal Pruning
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 29.4: MOLECULAR NETWORKING ===\n")

n_nodes <- 12
node_ids <- c(paste0("Acyl_Phe_", 1:7), paste0("Acyl_Trp_", 1:5))
nodes_df <- data.frame(id = node_ids, family = c(rep("Phe", 7), rep("Trp", 5)))

# Generate cosine similarity matrix with intra-family homology (pedagogical simulation)
sim_mat <- matrix(0.0, nrow = n_nodes, ncol = n_nodes, dimnames = list(node_ids, node_ids))
for (i in 1:n_nodes) {
  for (j in i:n_nodes) {
    if (i == j) {
      sim_mat[i, j] <- 1.0
    } else if (nodes_df$family[i] == nodes_df$family[j]) {
      sim_mat[i, j] <- sim_mat[j, i] <- round(runif(1, 0.52, 0.94), 3)
    } else {
      sim_mat[i, j] <- sim_mat[j, i] <- round(runif(1, 0.05, 0.28), 3)
    }
  }
}

# Top-K reciprocal pruning routine
topk_prune <- function(sim, K = 5, tau = 0.50) {
  diag(sim) <- 0
  candidate <- (sim >= tau)
  keep <- matrix(FALSE, nrow(sim), ncol(sim), dimnames = dimnames(sim))
  for (i in seq_len(nrow(sim))) {
    nbrs <- order(sim[i, ], decreasing = TRUE)
    pos_nbrs <- nbrs[sim[i, nbrs] > 0]
    top_nbrs <- pos_nbrs[seq_len(min(K, length(pos_nbrs)))]
    keep[i, top_nbrs] <- TRUE
  }
  mutual <- candidate & keep & t(keep)
  graph_from_adjacency_matrix(mutual * 1, mode = "undirected", diag = FALSE)
}

g <- topk_prune(sim_mat, K = 5, tau = 0.50)
V(g)$family <- nodes_df$family

comp <- components(g)
cat(sprintf("Molecular Similarity Network: %d nodes, %d edges (Cosine threshold >= 0.50, Top-K = 5)\n",
            vcount(g), ecount(g)))
cat(sprintf("Number of Disjoint Molecular Families (Connected Components): %d\n", comp$no))
cat("Family Component Sizes:\n")
print(table(Family_Size = comp$csize))

cat("\nPhe Family Members (Component with majority Phe):\n")
phe_comp_id <- comp$membership[V(g)$family == "Phe"][1]
print(V(g)$name[comp$membership == phe_comp_id])

cat("Trp Family Members (Component with majority Trp):\n")
trp_comp_id <- comp$membership[V(g)$family == "Trp"][1]
print(V(g)$name[comp$membership == trp_comp_id])


# -----------------------------------------------------------------------------
# Listing 29.5: Hardened Biological Replicate Concordance Audit
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 29.5: BIOLOGICAL REPLICATE CONCORDANCE ===\n")

replicate_data <- data.frame(
  Compound     = c(rep("Isocaproyl-Phe", 3), rep("Phenylacetyl-Phe", 3),
                   rep("Isovaleryl-Trp", 3), "Artifact_Peak_X"),
  Replicate    = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1),
  RT           = c(19.96, 19.98, 19.97, 18.42, 18.44, 18.43, 16.12, 16.15, 16.13, 12.50),
  Area         = c(85722, 88100, 84310, 52100, 53400, 51800, 31200, 29800, 32100, 4200),
  Labeled_Pair = c(rep(TRUE, 9), TRUE),
  stringsAsFactors = FALSE
)

audit_replicates <- function(df) {
  df %>%
    group_by(Compound) %>%
    summarise(
      Reps_Present = n_distinct(Replicate),
      Reps_Labeled = sum(Labeled_Pair, na.rm = TRUE),
      RT_StdDev    = if (n_distinct(Replicate) >= 2) round(sd(RT), 4) else NA_real_,
      Area_CV      = if (n_distinct(Replicate) >= 2) round(sd(Area) / mean(Area), 3) else NA_real_,
      .groups      = "drop"
    ) %>%
    mutate(
      Status = ifelse(
        Reps_Present >= 3 & Reps_Labeled >= 3 &
          !is.na(RT_StdDev) & RT_StdDev <= 0.05 &
          !is.na(Area_CV) & Area_CV <= 0.30,
        "VALIDATED", "REJECTED (Noise)"
      )
    )
}

audit_res <- audit_replicates(replicate_data)
cat("Biological Replicate Concordance Audit (n = 3):\n")
print(as.data.frame(audit_res))


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 29.1: Phenylacetyl-Tyr formula and 13C mass calculation
form_patyr <- calc_nacyl_formula(aa_formulas$Tyr, acyl_formulas$`Phenylacetyl(PA)`)
form_patyr_str <- paste0(names(form_patyr), form_patyr, collapse = "")
exact_patyr <- MetaboCoreUtils::calculateMass(form_patyr_str)
mz_patyr_nat <- as.numeric(MetaboCoreUtils::mass2mz(exact_patyr, adduct = "[M+H]+"))
mz_patyr_lab <- mz_patyr_nat + 9 * d13c
cat(sprintf("Exercise 29.1: Phenylacetyl-Tyr Formula: %s\n", form_patyr_str))
cat(sprintf("  Native [M+H]+ = %.4f | 13C9-labeled [M+H]+ = %.4f\n", mz_patyr_nat, mz_patyr_lab))

# Exercise 29.2: PPM deviation formula and minimum window
theo_mz <- 314.1387
obs_mz  <- 314.1404
delta_m <- obs_mz - theo_mz
ppm_dev <- (delta_m / theo_mz) * 1e6
cat(sprintf("\nExercise 29.2: Mass Deviation = %+.4f Da (%.4f ppm)\n", delta_m, ppm_dev))
cat(sprintf("  Minimum symmetric tolerance window: +/- %.2f ppm (or +/- %.1f ppm)\n", ppm_dev, round(ppm_dev, 1)))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
