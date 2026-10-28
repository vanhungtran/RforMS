# =============================================================================
# Chapter 31: Capstone Multi-Omics Clinical Case Studies
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

Sys.setenv(RENV_CONFIG_NAMESPACES_CHECK = "FALSE")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(limma)
  library(DEP)
  library(tibble)
  library(dplyr)
})

cat("===================================================================\n")
cat("Chapter 31 Laboratory: Multi-Omics Capstone Case Studies\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Case Study 1: Label-Free Quantitative Proteomics (DEP::UbiLength)
# -----------------------------------------------------------------------------
cat("--- Case Study 1: Label-Free Proteomics (UbiLength HeLa Perturbation) ---\n")

data("UbiLength", package = "DEP")
cat(sprintf("Raw MaxQuant ProteinGroups: %d rows x %d columns\n",
            nrow(UbiLength), ncol(UbiLength)))

# Step 1: Filter reverse hits and potential contaminants
clean_data <- UbiLength %>%
  filter(Reverse != "+", Potential.contaminant != "+")

# Ensure unique gene symbols
clean_data <- DEP::make_unique(clean_data, names = "Gene.names", ids = "Protein.IDs")
cat(sprintf("Filtered Analytical Proteins: %d (Dropped: %d)\n",
            nrow(clean_data), nrow(UbiLength) - nrow(clean_data)))

# Step 2: Extract LFQ intensity columns and construct SummarizedExperiment
lfq_cols <- grep("^LFQ.intensity", colnames(clean_data))
exp_design <- DEP::UbiLength_ExpDesign

se_prot <- DEP::make_se(clean_data, lfq_cols, exp_design)
cat(sprintf("Constructed SummarizedExperiment: %d proteins x %d samples\n",
            nrow(se_prot), ncol(se_prot)))

# Step 3: Variance stabilizing normalization (VSN) & missing-value handling
se_vsn <- DEP::normalize_vsn(se_prot)

# Filter proteins with zero missing values across conditions, then impute
se_filt <- DEP::filter_missval(se_vsn, thr = 0)
cat(sprintf("Proteins Retained Without Structural Missingness: %d\n", nrow(se_filt)))

# MinProb imputation for left-censored missing data (MNAR)
set.seed(42)
se_imp <- DEP::impute(se_filt, fun = "MinProb", q = 0.01)

# Step 4: Differential abundance testing via limma
diff_prot <- DEP::test_diff(se_imp, type = "control", control = "Ctrl")
dep_prot  <- DEP::add_rejections(diff_prot, alpha = 0.05, lfc = 1.0)
res_prot  <- DEP::get_results(dep_prot)

# Focus on primary contrast: Ubi6_vs_Ctrl
sig_proteins <- res_prot %>%
  filter(significant) %>%
  arrange(desc(abs(Ubi6_vs_Ctrl_ratio)))

cat(sprintf("Differentially Abundant Proteins (Any Contrast, |log2FC| >= 1.0, FDR < 0.05): %d\n",
            nrow(sig_proteins)))
cat(sprintf("  - Ubi6 vs Ctrl Significant (|log2FC| >= 1.0): %d\n",
            sum(sig_proteins$Ubi6_vs_Ctrl_significant)))
cat(sprintf("    * Upregulated:   %d\n", sum(sig_proteins$Ubi6_vs_Ctrl_ratio >= 1.0 & sig_proteins$Ubi6_vs_Ctrl_significant)))
cat(sprintf("    * Downregulated: %d\n", sum(sig_proteins$Ubi6_vs_Ctrl_ratio <= -1.0 & sig_proteins$Ubi6_vs_Ctrl_significant)))

cat("\nTop 5 Differential Proteins (HeLa Ubi6 vs Ctrl):\n")
top_5_prot <- sig_proteins %>%
  filter(Ubi6_vs_Ctrl_significant) %>%
  slice_head(n = 5) %>%
  select(name, ID, Ubi6_vs_Ctrl_ratio, Ubi6_vs_Ctrl_p.adj)
print(as.data.frame(top_5_prot))
cat("\n")

# -----------------------------------------------------------------------------
# Case Study 2: Untargeted Metabolomics Pipeline & Adduct Annotation
# -----------------------------------------------------------------------------
cat("--- Case Study 2: Untargeted Metabolomics (faahKO Murine Cohort) ---\n")

# Simulated feature matrix derived from faahKO OBI-Warp aligned consensus features
set.seed(123)
n_metabs <- 827
n_samples <- 4

base_intensity <- rlnorm(n_metabs, meanlog = 16, sdlog = 1.5)
metab_mat <- matrix(NA_real_, nrow = n_metabs, ncol = n_samples,
                    dimnames = list(paste0("FT_", sprintf("%04d", 1:n_metabs)),
                                    c("KO_15", "KO_16", "WT_15", "WT_16")))

for (i in 1:n_metabs) {
  ko_mean <- if (i <= 10) base_intensity[i] * 4.0 else base_intensity[i] # Spiked true KO accumulation
  ko_vals <- rlnorm(2, meanlog = log(ko_mean), sdlog = 0.20)
  wt_vals <- rlnorm(2, meanlog = log(base_intensity[i]), sdlog = 0.20)
  metab_mat[i, ] <- c(ko_vals, wt_vals)
}

# Construct SummarizedExperiment
col_data <- DataFrame(
  sample_name = colnames(metab_mat),
  genotype    = rep(c("KO", "WT"), each = 2),
  row.names   = colnames(metab_mat)
)

row_data <- DataFrame(
  mzmed = seq(200.1, 599.9, length.out = n_metabs),
  rtmed = seq(2520, 4480, length.out = n_metabs),
  row.names = rownames(metab_mat)
)

se_metab <- SummarizedExperiment(
  assays  = list(intensity = metab_mat),
  rowData = row_data,
  colData = col_data
)

cat(sprintf("Metabolomics SummarizedExperiment: %d features x %d samples\n",
            nrow(se_metab), ncol(se_metab)))

# Normalization: Log2 + Median Centering
mat_log2 <- log2(assay(se_metab) + 1)
sample_meds <- apply(mat_log2, 2, median, na.rm = TRUE)
mat_norm <- sweep(mat_log2, 2, sample_meds - mean(sample_meds), "-")

# Differential Analysis with limma
design_metab <- model.matrix(~ 0 + factor(colData(se_metab)$genotype, levels = c("WT", "KO")))
colnames(design_metab) <- c("WT", "KO")
fit_m <- lmFit(mat_norm, design_metab)
contrast_m <- makeContrasts(KO - WT, levels = design_metab)
fit_m <- contrasts.fit(fit_m, contrast_m)
fit_m <- eBayes(fit_m)

top_metab <- topTable(fit_m, number = 10, adjust.method = "BH")
cat("Top Differential Metabolite Features (KO vs WT):\n")
print(round(top_metab[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 4))

sig_metab_count <- sum(top_metab$adj.P.Val <= 0.05 & abs(top_metab$logFC) >= 1.0)
cat(sprintf("\nSignificant Features (|log2FC| >= 1.0, FDR <= 0.05): %d / 10\n", sig_metab_count))

# Candidate Adduct Calculation for FAAH Substrate (Oleamide, C18H35NO, 281.2719 Da)
m_proton <- 1.007276
m_sodium <- 22.989218
mono_oleamide <- 281.271864

adduct_pos_h  <- mono_oleamide + m_proton
adduct_pos_na <- mono_oleamide + m_sodium

cat("\nFAAH Substrate (Oleamide, C18H35NO) Theoretical Adduct Masses:\n")
cat(sprintf("  Neutral Monoisotopic Mass: %.4f Da\n", mono_oleamide))
cat(sprintf("  [M+H]+  Adduct m/z:        %.4f Th\n", adduct_pos_h))
cat(sprintf("  [M+Na]+ Adduct m/z:        %.4f Th\n\n", adduct_pos_na))

# -----------------------------------------------------------------------------
# Section 3: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 31 Exercises\n")
cat("===================================================================\n\n")

# Exercise 31.1: Multi-Omics Shared Container Architecture
cat("--- Exercise 31.1: Shared Computational Grammar Across Omics ---\n")
cat("Shared Components: SummarizedExperiment 2D coordinate invariant, limma empirical\n")
cat("Bayes variance shrinkage, Benjamini-Hochberg FDR control, and Quarto parameterized reporting.\n")
cat("Divergent Components: Raw tier data structures (QFeatures PSM->peptide->protein DAG\n")
cat("vs xcms chromatographic peak density correspondence and OBI-Warp alignment).\n\n")

# Exercise 31.2: Blocking Factor Adjustment
cat("--- Exercise 31.2: Blocking Factor Incorporation in limma ---\n")
cat("Model: lmFit(assay, design = model.matrix(~ replicate + condition))\n")
cat("Analytical Benefit: Adjusting for replicate batch effects absorbs paired technical variance,\n")
cat("increasing residual degrees of freedom precision and enhancing statistical power.\n\n")

# Exercise 31.3: Oleamide Adduct Matching
cat("--- Exercise 31.3: Candidate Adduct Identification ---\n")
cat(sprintf("Matching candidate Oleamide [M+H]+ (%.4f) against detected feature m/z values\n", adduct_pos_h))
cat("within +/- 15 ppm confirms putative MSI Level 2 identity based on exact mass and retention time.\n\n")

cat("===================================================================\n")
cat("Chapter 31 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
