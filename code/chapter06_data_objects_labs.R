# =============================================================================
# Chapter 06: Modern Mass Spectrometry Data Containers
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(Spectra)
  library(MsExperiment)
  library(SummarizedExperiment)
  library(QFeatures)
  library(MsCoreUtils)
  library(tibble)
  library(dplyr)
  library(tidyr)
})

cat("===================================================================\n")
cat("Chapter 06 Laboratory: Modern Mass Spectrometry Data Containers\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Container 1 - Spectra Scans and SpectraVariables
# -----------------------------------------------------------------------------
cat("--- Section 1: Spectra Container Mechanics ---\n")

# Create a synthetic Spectra container with 6 scans (MS1 and MS2)
mz_list <- list(
  c(200.1, 400.2, 600.3, 800.4),
  c(150.1, 250.2, 350.3),
  c(180.1, 290.2),
  c(200.1, 400.25, 600.35, 800.45),
  c(145.1, 245.2, 345.3),
  c(175.1, 285.2)
)
intensity_list <- list(
  c(10000, 50000, 80000, 20000),
  c(500, 2500, 1200),
  c(800, 1600),
  c(12000, 55000, 85000, 22000),
  c(600, 2800, 1400),
  c(900, 1750)
)

sp_df <- DataFrame(
  msLevel = c(1L, 2L, 2L, 1L, 2L, 2L),
  rtime   = c(10.5, 11.2, 11.8, 14.5, 15.1, 15.7),
  precursorMz = c(NA, 400.2, 600.3, NA, 400.25, 600.35),
  precursorCharge = c(NA, 2L, 3L, NA, 2L, 3L)
)
sp_df$mz <- mz_list
sp_df$intensity <- intensity_list

sps <- Spectra(sp_df)
cat(sprintf("Constructed Spectra: %d scans (MS1: %d, MS2: %d)\n",
            length(sps), sum(msLevel(sps) == 1), sum(msLevel(sps) == 2)))
cat(sprintf("Precursor m/z values in MS2 scans: %s\n\n",
            paste(round(precursorMz(filterMsLevel(sps, 2L)), 2), collapse = ", ")))

# -----------------------------------------------------------------------------
# Section 2: Container 2 - MsExperiment Experimental Design Binding
# -----------------------------------------------------------------------------
cat("--- Section 2: MsExperiment Design Binding ---\n")

# Bind Spectra to sample metadata table using readMsExperiment
tmt_files <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
sample_tbl <- data.frame(
  sample_id       = "Run_01",
  condition       = "Control",
  batch           = "Batch_A",
  injection_order = 1L,
  stringsAsFactors = FALSE
)

mse <- readMsExperiment(spectraFiles = tmt_files, sampleData = sample_tbl)

cat(sprintf("MsExperiment constructed: %d samples, %d total spectra\n",
            nrow(sampleData(mse)), length(spectra(mse))))
cat(sprintf("Sample conditions: %s\n",
            paste(sampleData(mse)$condition, collapse = ", ")))
cat(sprintf("Spectra MS levels: MS1 (%d), MS2 (%d), MS3 (%d)\n\n",
            sum(msLevel(spectra(mse)) == 1),
            sum(msLevel(spectra(mse)) == 2),
            sum(msLevel(spectra(mse)) == 3)))

# -----------------------------------------------------------------------------
# Section 3: Container 3 - SummarizedExperiment for Feature Matrices
# -----------------------------------------------------------------------------
cat("--- Section 3: SummarizedExperiment Synchronized Subsetting ---\n")

set.seed(42)
n_features <- 8
n_samples  <- 6

assay_counts <- matrix(
  round(rlnorm(n_features * n_samples, meanlog = 12, sdlog = 1.2)),
  nrow = n_features,
  ncol = n_samples,
  dimnames = list(
    paste0("PEP_", sprintf("%03d", 1:n_features)),
    paste0("Sample_", sprintf("%02d", 1:n_samples))
  )
)

row_meta <- DataFrame(
  peptide_id  = paste0("PEP_", sprintf("%03d", 1:n_features)),
  protein_id  = c("PROT_A", "PROT_A", "PROT_A", "PROT_B", "PROT_B", "PROT_C", "PROT_C", "PROT_D"),
  gene_symbol = c("ACTB", "ACTB", "ACTB", "GAPDH", "GAPDH", "ALB", "ALB", "INS"),
  length_aa   = c(12, 15, 9, 14, 18, 11, 16, 8),
  row.names   = paste0("PEP_", sprintf("%03d", 1:n_features))
)

col_meta <- DataFrame(
  sample_id   = paste0("Sample_", sprintf("%02d", 1:n_samples)),
  condition   = rep(c("Control", "Treated"), each = 3),
  batch       = rep(c("B1", "B2", "B1", "B2", "B1", "B2")),
  row.names   = paste0("Sample_", sprintf("%02d", 1:n_samples))
)

se <- SummarizedExperiment(
  assays  = list(raw_intensity = assay_counts),
  rowData = row_meta,
  colData = col_meta
)

cat(sprintf("SummarizedExperiment: %d features x %d samples\n", nrow(se), ncol(se)))
cat(sprintf("Assay Name: %s | Dimensions: [%d, %d]\n",
            assayNames(se)[1], nrow(assay(se)), ncol(assay(se))))

# Coordinate subsetting verification: Filter to PROT_A peptides in Treated samples
se_sub <- se[rowData(se)$protein_id == "PROT_A", se$condition == "Treated"]
cat(sprintf("Subset se[protein_id == 'PROT_A', condition == 'Treated']:\n"))
cat(sprintf("  Features Retained: %d (%s)\n", nrow(se_sub), paste(rownames(se_sub), collapse = ", ")))
cat(sprintf("  Samples Retained:  %d (%s)\n", ncol(se_sub), paste(colnames(se_sub), collapse = ", ")))
cat(sprintf("  RowData Matches Assay Rows: %s\n", identical(rownames(rowData(se_sub)), rownames(assay(se_sub)))))
cat(sprintf("  ColData Matches Assay Cols: %s\n\n", identical(rownames(colData(se_sub)), colnames(assay(se_sub)))))

# -----------------------------------------------------------------------------
# Section 4: Container 4 - QFeatures Hierarchical Aggregation Pipeline
# -----------------------------------------------------------------------------
cat("--- Section 4: QFeatures Hierarchical Aggregation ---\n")

# Build QFeatures starting from PSM-level SummarizedExperiment
data(feat1, package = "QFeatures")
cat(sprintf("Built-in QFeatures 'feat1': Assays = %s\n", paste(names(feat1), collapse = ", ")))
cat(sprintf("PSM Assay Dimensions: %d rows x %d samples\n",
            nrow(feat1[["psms"]]), ncol(feat1[["psms"]])))

# Aggregate PSMs to Peptides using robustSummary (Tukey median polish)
feat1 <- aggregateFeatures(
  feat1,
  i    = "psms",
  fcol = "Sequence",
  name = "peptides",
  fun  = MsCoreUtils::robustSummary
)

# Aggregate Peptides to Proteins using robustSummary
feat1 <- aggregateFeatures(
  feat1,
  i    = "peptides",
  fcol = "Protein",
  name = "proteins",
  fun  = MsCoreUtils::robustSummary
)

cat(sprintf("Aggregated QFeatures: %d Assays (%s)\n",
            length(feat1), paste(names(feat1), collapse = " -> ")))
cat(sprintf("  Level 1 (psms):     %d rows\n", nrow(feat1[["psms"]])))
cat(sprintf("  Level 2 (peptides): %d rows\n", nrow(feat1[["peptides"]])))
cat(sprintf("  Level 3 (proteins): %d rows\n", nrow(feat1[["proteins"]])))

# Inspect assayLink parentage
al_prot <- assayLink(feat1, "proteins")
cat(sprintf("AssayLink 'proteins' -> hits from parent assay '%s':\n", al_prot@from))
hits_df <- as.data.frame(al_prot@hits)
print(hits_df[, c("names_from", "names_to")])
cat("\n")

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 06 Exercises\n")
cat("===================================================================\n\n")

# Exercise 6.1: Subsetting Synchrony and Dimensions
cat("--- Exercise 6.1: Subsetting Invariant Verification ---\n")
# Question: In the 8x6 SummarizedExperiment, a researcher subsets:
# se_filtered <- se[rowData(se)$length_aa >= 12, se$batch == "B1"]
se_ex1 <- se[rowData(se)$length_aa >= 12, se$batch == "B1"]
cat(sprintf("Original Matrix: %d features x %d samples\n", nrow(se), ncol(se)))
cat(sprintf("Filtered Matrix: %d features x %d samples\n", nrow(se_ex1), ncol(se_ex1)))
cat(sprintf("Retained Features (length_aa >= 12): %s\n", paste(rownames(se_ex1), collapse = ", ")))
cat(sprintf("Retained Samples (batch == 'B1'):     %s\n", paste(colnames(se_ex1), collapse = ", ")))
cat(sprintf("Integrity Check: nrow(assay) == nrow(rowData) [%s] | ncol(assay) == ncol(colData) [%s]\n\n",
            nrow(assay(se_ex1)) == nrow(rowData(se_ex1)),
            ncol(assay(se_ex1)) == ncol(colData(se_ex1))))

# Exercise 6.2: Aggregation Summary Statistics: Median Polish vs Mean vs Sum
cat("--- Exercise 6.2: Aggregation Estimator Sensitivity under Outliers ---\n")
# Protein with 3 peptides measured across 2 samples:
# Sample 1: [100, 110, 1000] (peptide 3 is an unremoved contaminant / miscleaved outlier)
# Sample 2: [105, 115, 120]
toy_mat <- matrix(
  c(100, 110, 1000, 105, 115, 120),
  nrow = 3, ncol = 2,
  dimnames = list(c("Pep_1", "Pep_2", "Pep_3"), c("S1", "S2"))
)

agg_mean   <- colMeans(toy_mat)
agg_sum    <- colSums(toy_mat)
agg_median <- apply(toy_mat, 2, median)
agg_robust <- MsCoreUtils::robustSummary(toy_mat)

cat("Protein Abundance Estimators for Toy Peptide Matrix:\n")
cat(sprintf("  Sample 1 [100, 110, 1000]: Sum = %.1f | Mean = %.1f | Median = %.1f | RobustSummary = %.2f\n",
            agg_sum["S1"], agg_mean["S1"], agg_median["S1"], agg_robust["S1"]))
cat(sprintf("  Sample 2 [105, 115, 120]:  Sum = %.1f | Mean = %.1f | Median = %.1f | RobustSummary = %.2f\n",
            agg_sum["S2"], agg_mean["S2"], agg_median["S2"], agg_robust["S2"]))
cat(sprintf("Apparent Fold Change (S1 / S2):\n"))
cat(sprintf("  Sum:           %.2fx (Outlier inflates fold change 3.6-fold)\n", agg_sum["S1"] / agg_sum["S2"]))
cat(sprintf("  Mean:          %.2fx (Linear distortion)\n", agg_mean["S1"] / agg_mean["S2"]))
cat(sprintf("  Median:        %.2fx (Resistant)\n", agg_median["S1"] / agg_median["S2"]))
cat(sprintf("  RobustSummary: %.2fx (Tukey biweight damping)\n\n", agg_robust["S1"] / agg_robust["S2"]))

# Exercise 6.3: Memory Footprint: QFeatures Tree vs Flat Denormalized Data Frame
cat("--- Exercise 6.3: Memory Scalability of QFeatures vs Denormalized Tables ---\n")
n_psms <- 50000
n_peps <- 15000
n_prots <- 3000
n_cols  <- 10

# In a flat table, 50,000 rows duplicate peptide and protein metadata across all 10 sample columns
flat_df <- data.frame(
  psm_id = paste0("PSM_", 1:n_psms),
  sequence = sample(paste0("PEP_", 1:n_peps), n_psms, replace = TRUE),
  protein  = sample(paste0("PROT_", 1:n_prots), n_psms, replace = TRUE),
  matrix(runif(n_psms * n_cols), nrow = n_psms, ncol = n_cols)
)
size_flat_mb <- as.numeric(object.size(flat_df)) / (1024^2)

cat(sprintf("Flat Denormalized PSM Table (50,000 rows x 13 columns): %.2f MB\n", size_flat_mb))
cat("QFeatures Normalization Advantage:\n")
cat("  - Normalizes metadata into relational assays connected via lightweight integer indices.\n")
cat("  - Eliminates string redundancy for peptide sequences and protein accessions.\n")
cat("  - Preserves exact parentage links with zero redundant floating-point array duplication.\n\n")

# Exercise 6.4: Cross-Object Interoperability: Extracting SE from QFeatures for Limma
cat("--- Exercise 6.4: Seamless Extraction of SummarizedExperiment for Modeling ---\n")
prot_se <- feat1[["proteins"]]
cat(sprintf("Extracted Standalone Assay: %s\n", class(prot_se)))
cat(sprintf("Dimensions: %d proteins x %d samples\n", nrow(prot_se), ncol(prot_se)))
cat(sprintf("Extracting Assay Matrix: is.matrix(assay(prot_se)) -> %s\n",
            is.matrix(assay(prot_se))))
cat(sprintf("Sample Metadata Preserved: colnames = [%s]\n",
            paste(colnames(colData(prot_se)), collapse = ", ")))
cat("Pedagogical Implication: Any package accepting SummarizedExperiment (limma, DEqMS, proDA)\n")
cat("interoperates directly with feat1[['proteins']] without flattening to CSV or disk serialization.\n\n")

cat("===================================================================\n")
cat("Chapter 06 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
