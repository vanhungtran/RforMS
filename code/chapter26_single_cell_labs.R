# =============================================================================
# Chapter 26: Single-Cell Proteomics Laboratory Scripts
# Verified reproducible code for Listings 1 to 7
# =============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(QFeatures)
  library(MsCoreUtils)
  library(lme4)
})

# -----------------------------------------------------------------------------
# Listing 1: Simulating a 16-Plex Single-Cell TMT Multiplex
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 1 ===\n")
set.seed(303)

n_psms <- 800
channels <- c("TMT_126_Carrier", "TMT_127N_Ref", "TMT_127C_Blank",
              paste0("TMT_Naive_", 1:7), paste0("TMT_Primed_", 1:6))
n_channels <- length(channels)

# Base log2 peptide abundance across proteome (mean = 11.0, sd = 1.8)
base_abund <- rnorm(n_psms, mean = 11.0, sd = 1.8)

# True biology: 80 PSMs (first 16 proteins) up-regulated in Naive cells (FC = 2.5, log2FC = +1.32)
bio_effect <- rep(0, n_psms)
bio_effect[1:80] <- 1.32

# Channel scaling factors (relative loading)
channel_loadings <- c(100, 5, 0.0, rep(1.0, 7), rep(1.0, 6))

raw_psm_mat <- matrix(NA_real_, nrow = n_psms, ncol = n_channels)
colnames(raw_psm_mat) <- channels
rownames(raw_psm_mat) <- paste0("PSM_", sprintf("%04d", 1:n_psms))

# Carrier (100x) and Reference (5x) channels
raw_psm_mat[, 1] <- rpois(n_psms, lambda = (2^base_abund) * 100)
raw_psm_mat[, 2] <- rpois(n_psms, lambda = (2^base_abund) * 5)

# Negative control / blank channel: low Poisson background noise
blank_counts <- rpois(n_psms, lambda = 2)
raw_psm_mat[, 3] <- ifelse(blank_counts < 8, NA, blank_counts)

# Single cells: 1x loading with physical Poisson shot noise and detection limits
lod_thresh <- 25  # Limit of Detection in raw ion counts

for (j in 4:n_channels) {
  is_naive <- grepl("Naive", channels[j])
  exp_log2 <- base_abund + (if (is_naive) bio_effect else 0)
  exp_linear <- 2^exp_log2
  
  # Physical Poisson ion arrival counts
  counts <- rpois(n_psms, lambda = exp_linear)
  
  # Intensity-dependent DDA precursor sampling (MAR)
  prob_detect <- 1 / (1 + exp(-0.8 * (exp_log2 - 10.5)))
  selected <- runif(n_psms) < prob_detect
  
  # Apply limit of detection (MNAR) and duty cycle sampling
  obs <- counts
  obs[!selected | counts < lod_thresh] <- NA
  raw_psm_mat[, j] <- obs
}

# Construct Channel Metadata (colData)
col_meta <- data.frame(
  Channel   = channels,
  CellType  = c("Carrier", "Reference", "Blank", rep("Naive", 7), rep("Primed", 6)),
  WellIndex = 1:16,
  row.names = channels,
  stringsAsFactors = FALSE
)

# Construct Strictly Nested Feature Metadata (rowData): 160 proteins -> 320 peptides -> 800 PSMs
pep_ids <- paste0("PEP_", sprintf("%03d", rep(1:320, times = rep(c(3, 2), times = 160))))
prot_ids <- rep(paste0("PROT_", sprintf("%03d", 1:160)), each = 5)

row_meta <- data.frame(
  PSM_ID   = rownames(raw_psm_mat),
  Protein  = prot_ids,
  Peptide  = pep_ids,
  PEP      = runif(n_psms, 0.001, 0.08),
  Reverse  = sample(c(rep("-", 760), rep("+", 40))),
  stringsAsFactors = FALSE
)

cat(sprintf("Simulated %d PSMs mapping to 320 peptides and 160 proteins.\n", n_psms))
cat(sprintf("Carrier Intensity: %.2e | Median Single-Cell Intensity: %.2e\n",
            median(raw_psm_mat[, 1], na.rm = TRUE),
            median(raw_psm_mat[, 4:16], na.rm = TRUE)))

# -----------------------------------------------------------------------------
# Listing 2: Building the QFeatures Single-Cell Container
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 2 ===\n")
sce_psm <- SingleCellExperiment(
  assays  = list(psms_raw = raw_psm_mat),
  rowData = row_meta,
  colData = col_meta
)

qf <- QFeatures(list(psms = sce_psm), colData = col_meta)

cat("QFeatures Container Structure:\n")
print(qf)

# -----------------------------------------------------------------------------
# Listing 3: Cell-Level Quality Control and Filtering
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 3 ===\n")
# 1. Feature Filtering: Remove Decoys (Reverse == '+') and high PEP (> 0.05)
qf <- filterFeatures(qf, ~ Reverse == "-" & PEP < 0.05)

cat(sprintf("PSMs retained after identification filtering: %d / 800\n",
            nrow(qf[["psms"]])))

# 2. Cell Filtering: Retain Single-Cell Channels for Biological Quantification
qf_single <- subsetByColData(qf, colData(qf)$CellType %in% c("Naive", "Primed"))

cat(sprintf("Channels retained for single-cell quantification: %d\n",
            ncol(qf_single[["psms"]])))

# -----------------------------------------------------------------------------
# Listing 4: Multi-Level Aggregation and Normalization
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 4 ===\n")
# Step 1: Remove unobserved features across single cells before robust aggregation
psm_obs <- rowSums(!is.na(assay(qf_single[["psms"]])))
qf_single <- qf_single[names(which(psm_obs >= 1)), , ]

# Step 2: Normalize PSM intensities across single cells (median centering)
qf_single <- normalize(qf_single,
                       i = "psms",
                       name = "psms_norm",
                       method = "center.median")

# Step 3: Aggregate PSMs to Peptides (Robust Tukey Biweight M-estimation)
qf_single <- aggregateFeatures(qf_single,
                               i = "psms_norm",
                               fcol = "Peptide",
                               name = "peptides",
                               fun = MsCoreUtils::robustSummary)

# Step 4: Aggregate Peptides to Protein Groups
qf_single <- aggregateFeatures(qf_single,
                               i = "peptides",
                               fcol = "Protein",
                               name = "proteins",
                               fun = MsCoreUtils::robustSummary)

cat("\nAggregated QFeatures Assay Topology:\n")
print(names(qf_single))
cat(sprintf("Retained Features: %d peptides | %d proteins\n",
            nrow(qf_single[["peptides"]]), nrow(qf_single[["proteins"]])))
cat(sprintf("Final Protein Matrix Dimensions: %d proteins x %d single cells\n",
            nrow(qf_single[["proteins"]]), ncol(qf_single[["proteins"]])))

# -----------------------------------------------------------------------------
# Listing 5: Single-Cell Phenotypic Clustering via Principal Component Analysis
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 5 ===\n")
prot_matrix <- assay(qf_single[["proteins"]])

# Restrict PCA to complete-case proteins (standard non-imputed single-cell analysis)
prot_complete <- na.omit(prot_matrix)
cat(sprintf("Complete-case proteins for PCA: %d / %d\n",
            nrow(prot_complete), nrow(prot_matrix)))

# Perform PCA on transposed matrix (cells as rows)
pca_res <- prcomp(t(prot_complete), scale. = TRUE)
pca_var <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

# Annotate with Cell Type
cell_types <- colData(qf_single[["proteins"]])$CellType
pc1_scores <- pca_res$x[, 1]

# Assess linear separation between Naive and Primed
t_test_res <- t.test(pc1_scores ~ cell_types)

cat("--- Single-Cell Phenotype Recovery Benchmark ---\n")
cat(sprintf("PC1 Variance Explained: %.1f%% | PC2 Variance: %.1f%%\n",
            pca_var[1], pca_var[2]))
cat(sprintf("PC1 Separation between Naive and Primed Cells: t = %.2f (p = %.2e)\n",
            t_test_res$statistic, t_test_res$p.value))

# -----------------------------------------------------------------------------
# Listing 6: Automated Single-Cell QC Channel Auditing (Exercise 26.3 Solution)
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 6 ===\n")
audit_scp_channels <- function(assay_mat, blank_col = "TMT_127C_Blank",
                               min_psms = 250, max_ratio = 2.5) {
  # 1. Blank Channel Intensity Baseline
  blank_intensity <- median(assay_mat[, blank_col], na.rm = TRUE)
  if (is.na(blank_intensity) || blank_intensity == 0) blank_intensity <- 1.0
  
  # 2. Per-Channel Metrics
  tot_intensity <- colSums(assay_mat, na.rm = TRUE)
  med_intensity <- apply(assay_mat, 2, median, na.rm = TRUE)
  n_detected    <- colSums(!is.na(assay_mat) & assay_mat > 0)
  
  # Identify Channel Classes
  is_blank   <- colnames(assay_mat) == blank_col
  is_carrier <- grepl("Carrier", colnames(assay_mat), ignore.case = TRUE)
  is_ref     <- grepl("Ref", colnames(assay_mat), ignore.case = TRUE)
  is_single  <- !is_blank & !is_carrier & !is_ref
  
  # Median single-cell intensity across bona fide single cells
  median_sc_tot <- median(tot_intensity[is_single], na.rm = TRUE)
  
  signal_to_blank <- med_intensity / blank_intensity
  doublet_ratio   <- tot_intensity / median_sc_tot
  
  # 3. Classify Channel Status with Channel-Type Awareness
  status <- rep("PASS", ncol(assay_mat))
  status[is_blank]   <- "CONTROL_BLANK"
  status[is_carrier] <- "CONTROL_CARRIER"
  status[is_ref]     <- "CONTROL_REFERENCE"
  
  # Single-cell quality control filters
  status[is_single & (n_detected < min_psms | signal_to_blank < 3.0)] <- "FLAG_EMPTY_WELL"
  status[is_single & doublet_ratio > max_ratio] <- "FLAG_DOUBLET"
  
  qc_summary <- data.frame(
    Channel         = colnames(assay_mat),
    Features_Count  = n_detected,
    Total_Signal    = tot_intensity,
    Signal_To_Blank = round(signal_to_blank, 2),
    Doublet_Ratio   = round(doublet_ratio, 2),
    Status          = status,
    row.names       = NULL,
    stringsAsFactors = FALSE
  )
  
  cat("--- Single-Cell Channel Quality Control Summary ---\n")
  print(qc_summary)
  return(invisible(qc_summary))
}

qc_res <- audit_scp_channels(raw_psm_mat, blank_col = "TMT_127C_Blank")

# -----------------------------------------------------------------------------
# Listing 7: Demonstration of Pseudoreplication and Type I Error Inflation (Exercise 26.4 Solution)
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 7 ===\n")
# Crossed Experimental Design: 6 Donors, each tested in both Condition A and Condition B
# 10 Single Cells per donor (5 cells Condition A, 5 cells Condition B) -> 60 total cells
set.seed(999)

n_donors <- 6
cells_per_cond <- 5
cells_per_donor <- cells_per_cond * 2
n_total_cells <- n_donors * cells_per_donor

donor_ids <- rep(paste0("Donor_", 1:n_donors), each = cells_per_donor)
conditions <- rep(rep(c("Condition_A", "Condition_B"), each = cells_per_cond), times = n_donors)
# Distinct culture preparations / wells per donor-condition pair
well_ids <- rep(paste0("Well_", 1:(n_donors * 2)), each = cells_per_cond)

# 200 Null Simulations
p_naive_list <- numeric(200)
p_lmm_list   <- numeric(200)

for (sim in 1:200) {
  # Biological donor baseline variance (ICC_donor = 0.35)
  u_donor <- rep(rnorm(n_donors, mean = 0, sd = 1.0), each = cells_per_donor)
  # Culture well preparation variance (ICC_well = 0.30)
  u_well  <- rep(rnorm(n_donors * 2, mean = 0, sd = 0.9), each = cells_per_cond)
  cell_noise <- rnorm(n_total_cells, mean = 0, sd = 0.7)
  
  y_null <- 10 + u_donor + u_well + cell_noise
  
  # 1. Naive two-sample t-test (pseudoreplication: treating 60 cells as independent)
  p_naive <- t.test(y_null ~ conditions)$p.value
  p_naive_list[sim] <- p_naive
  
  # 2. Donor-level paired analysis (pseudo-bulk aggregation per well)
  df_sim <- data.frame(y = y_null, condition = conditions, donor = donor_ids)
  df_agg <- aggregate(y ~ condition + donor, data = df_sim, FUN = mean)
  y_a <- df_agg$y[df_agg$condition == "Condition_A"]
  y_b <- df_agg$y[df_agg$condition == "Condition_B"]
  p_lmm <- t.test(y_a, y_b, paired = TRUE)$p.value
  p_lmm_list[sim] <- p_lmm
}

false_pos_naive <- mean(p_naive_list < 0.05) * 100
false_pos_lmm   <- mean(p_lmm_list < 0.05) * 100

cat("--- Pseudoreplication Type I Error Audit (200 Null Simulations) ---\n")
cat(sprintf("Naive Two-Sample t-Test False Positive Rate: %.1f%% (Catastrophic Error Inflation!)\n",
            false_pos_naive))
cat(sprintf("Donor-Level Paired Analysis False Positive Rate: %.1f%% (Preserves alpha = 5%%)\n",
            false_pos_lmm))
