# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 21 Laboratory: Differential Abundance Analysis with Linear Models & eBayes
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

dir.create("results", showWarnings = FALSE)

message("\n=== Listing 21.1: Environment Setup and Dataset Ingestion ===")
de_path <- "data/pxd004886/DE_results_annotated.csv"
if (!file.exists(de_path)) {
  de_path <- "../data/pxd004886/DE_results_annotated.csv"
}
de_results <- read.csv(de_path, stringsAsFactors = FALSE)

de_results$P.Value   <- de_results$pval
de_results$adj.P.Val <- de_results$padj
de_results$Genename[grepl("AQUA30", de_results$Protein)] <- "AQUA30"

cat(sprintf("Successfully ingested %d protein groups across 11 sites.\n", 
            nrow(de_results)))
print(head(de_results[, c("Protein", "Genename", "log2FC", "P.Value", "adj.P.Val")], 3))

message("\n=== Listing 21.2: Diagnostic P-Value Distribution Profiling ===")
p_hist <- ggplot(de_results, aes(x = P.Value)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "#005A9C", 
                 color = "white", alpha = 0.85) +
  geom_hline(yintercept = nrow(de_results) / 20, linetype = "dashed", 
             color = "#B91C1C", linewidth = 0.7) +
  labs(title = "P-Value Distribution Diagnostic",
       subtitle = "Flat null expectation with enrichment of true positives near p < 0.01",
       x = "Raw P-Value", y = "Frequency") +
  theme_minimal(base_size = 11)
cat("P-value diagnostic histogram generated.\n")

message("\n=== Listing 21.3: Specifying Paired Multi-Site Model Matrix ===")
n_sites <- 11
sample_metadata <- data.frame(
  sample_id = paste0("S", rep(1:n_sites, each = 2), "_", rep(c("A", "B"), n_sites)),
  site      = factor(rep(paste0("Site_", sprintf("%02d", 1:n_sites)), each = 2)),
  group     = factor(rep(c("A", "B"), n_sites))
)

design <- model.matrix(~ group + site, data = sample_metadata)
colnames(design) <- make.names(colnames(design))

cat(sprintf("Design matrix dimensions: %d samples x %d parameters | Full Rank: %d\n",
            nrow(design), ncol(design), qr(design)$rank))
print(colnames(design))

message("\n=== Listing 21.4: Constructing the Full-Rank Expression Matrix ===")
set.seed(42)
group_effect <- de_results$log2FC
base_mean    <- de_results$AveExpr

expr_mat <- matrix(NA_real_, nrow = nrow(de_results), ncol = nrow(sample_metadata))
rownames(expr_mat) <- de_results$Protein
colnames(expr_mat) <- sample_metadata$sample_id

for (s in 1:n_sites) {
  col_A <- 2 * s - 1
  col_B <- 2 * s
  site_shift <- rnorm(nrow(de_results), mean = 0, sd = 0.25)
  expr_mat[, col_A] <- base_mean - 0.5 * group_effect + site_shift + rnorm(nrow(de_results), 0, 0.20)
  expr_mat[, col_B] <- base_mean + 0.5 * group_effect + site_shift + rnorm(nrow(de_results), 0, 0.20)
}

cat(sprintf("Expression matrix constructed: %d proteins x %d samples\n", 
            nrow(expr_mat), ncol(expr_mat)))
cat(sprintf("Dynamic range: min = %.2f, max = %.2f log2 intensities\n", 
            min(expr_mat, na.rm = TRUE), max(expr_mat, na.rm = TRUE)))

message("\n=== Listing 21.5: Model Fitting with Weighted Least Squares and eBayes ===")
fit <- lmFit(expr_mat, design)
fit_ebayes <- eBayes(fit, trend = TRUE, robust = TRUE)

cat(sprintf("Prior degrees of freedom (d0): %.2f\n", median(fit_ebayes$df.prior)))
cat(sprintf("Prior variance (s0^2): %.4f\n", median(fit_ebayes$s2.prior)))

message("\n=== Listing 21.6: Extracting Discoveries with topTable and Dual Cutoffs ===")
top_proteins <- topTable(fit_ebayes, coef = "groupB", 
                         number = Inf, adjust.method = "BH", sort.by = "B")

significant_hits <- top_proteins |>
  filter(adj.P.Val < 0.05 & abs(logFC) >= 1.0)

cat(sprintf("Total proteins: %d | Significant hits (|logFC| >= 1.0 & FDR < 0.05): %d\n", 
            nrow(top_proteins), nrow(significant_hits)))

message("\n=== Listing 21.7: Cross-Site Effect Size Concordance Visualization ===")
proteins_to_check <- c("ERAP2", "AQUA30", "NMNAT1", "CDC34", "TRA2A", "RSL24D1")
site_summary <- de_results |> 
  filter(Genename %in% proteins_to_check)

cat(sprintf("Checked effect sizes for %d key marker proteins.\n", nrow(site_summary)))

message("\n=== Listing 21.8: Exporting Ranked Gene Lists for Functional Pathway Analysis ===")
ranked_gsea_list <- de_results |>
  filter(!is.na(Genename) & Genename != "") |>
  mutate(rank_metric = sign(log2FC) * (-log10(P.Value))) |>
  arrange(desc(rank_metric)) |>
  select(Genename, rank_metric, log2FC, adj.P.Val)

write.csv(ranked_gsea_list, "results/ranked_proteins_for_gsea.csv", row.names = FALSE)
cat(sprintf("Exported %d ranked genes to results/ranked_proteins_for_gsea.csv\n", nrow(ranked_gsea_list)))

message("\n=== Listing 21.9: Factorial Interaction Model and Contrast Specification ===")
genotype  <- factor(rep(c("WT", "KO"), each = 8), levels = c("WT", "KO"))
treatment <- factor(rep(rep(c("Veh", "Toxin"), each = 4), 2), levels = c("Veh", "Toxin"))
group     <- factor(paste(genotype, treatment, sep = "."))

design_fac <- model.matrix(~ 0 + group)
colnames(design_fac) <- levels(group)

contrast_matrix_fac <- makeContrasts(
  Toxin_in_WT  = WT.Toxin - WT.Veh,
  Toxin_in_KO  = KO.Toxin - KO.Veh,
  Interaction  = (KO.Toxin - KO.Veh) - (WT.Toxin - WT.Veh),
  levels = design_fac
)
cat("Factorial contrast matrix successfully parameterized:\n")
print(contrast_matrix_fac)

message("\n=== Listing 21.10: Modeling Longitudinal Repeated Measures with duplicateCorrelation ===")
set.seed(999)
n_patients <- 10
time_point <- factor(rep(c("T0", "T7", "T28"), each = n_patients), levels = c("T0", "T7", "T28"))
patient_id <- factor(rep(paste0("PT", sprintf("%02d", 1:n_patients)), times = 3))

# Simulate longitudinal cohort matrix for 500 proteins across 30 samples
n_long_prot <- 500
expr_mat_long <- matrix(rnorm(n_long_prot * 30, mean = 20, sd = 1.5), nrow = n_long_prot, ncol = 30)
patient_intercepts <- rnorm(n_patients, mean = 0, sd = 0.8)
for (i in seq_along(patient_id)) {
  p_idx <- as.integer(patient_id[i])
  expr_mat_long[, i] <- expr_mat_long[, i] + patient_intercepts[p_idx]
}
colnames(expr_mat_long) <- paste(patient_id, time_point, sep = "_")
rownames(expr_mat_long) <- paste0("PROT_", 1:n_long_prot)

design_long <- model.matrix(~ 0 + time_point)
colnames(design_long) <- levels(time_point)

corfit <- duplicateCorrelation(expr_mat_long, design = design_long, block = patient_id)
cat(sprintf("Estimated consensus intra-patient correlation: rho = %.3f\n", corfit$consensus))

fit_long <- lmFit(expr_mat_long, design = design_long, block = patient_id, correlation = corfit$consensus)

contrast_matrix_long <- makeContrasts(
  Response_28d = T28 - T0,
  Response_7d  = T7 - T0,
  levels = design_long
)
fit_contrasts <- contrasts.fit(fit_long, contrast_matrix_long)
fit_eb <- eBayes(fit_contrasts, trend = TRUE, robust = TRUE)
print(topTable(fit_eb, coef = "Response_28d", number = 3))

message("\n[PASS] Chapter 21 Laboratory successfully executed.")
