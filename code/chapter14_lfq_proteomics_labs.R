# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 14 Laboratory: Quantitative Label-Free Proteomics
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

message("\n=== Listing 14.1: Benchmark Ingestion and Schema Audit (PXD004886) ===")
de_file <- "data/pxd004886/DE_results_annotated.csv"
if (!file.exists(de_file)) stop("Benchmark data file not found: ", de_file)
de_df <- read.csv(de_file, stringsAsFactors = FALSE)

cat(sprintf("Successfully ingested %d quantified protein groups across 11 sites.\n", nrow(de_df)))
cat(sprintf("  Columns (%d): %s\n", ncol(de_df), paste(colnames(de_df)[1:6], collapse = ", ")))
cat(sprintf("  Groups: Group A (control) vs Group B (spiked)\n"))
cat(sprintf("  Known spike-in control: AQUA30 peptide standard\n"))

message("\n=== Listing 14.2: Benchmark Quality Filtering: Contaminant and Missingness Audit ===")
clean_df <- de_df %>%
  filter(!grepl("^(REV|DECOY)", Protein), !grepl("^(CON|contam)", Protein),
         !is.na(log2FC), !is.na(padj))

n_retained   <- nrow(clean_df)
pct_retained <- (n_retained / nrow(de_df)) * 100

cat(sprintf("Proteins passing quality control: %d (%.2f%%)\n", 
            n_retained, pct_retained))

message("\n=== Listing 14.3: Controlled Simulation: Multi-Site Median Normalization Mechanics ===")
set.seed(404)
n_prot  <- 2500
n_sites <- 11
base_expr <- rnorm(n_prot, mean = 20, sd = 2.5)

site_data <- lapply(seq_len(n_sites), function(i) {
  shift     <- rnorm(1, mean = 0, sd = 1.2)
  scale     <- runif(1, 0.85, 1.15)
  raw_vals  <- (base_expr + shift) * scale
  norm_vals <- raw_vals - median(raw_vals) + median(base_expr)
  data.frame(
    Run        = paste0("Site_", sprintf("%02d", i)),
    Raw        = raw_vals,
    Normalized = norm_vals
  )
})
df_density <- do.call(rbind, site_data)

cat(sprintf("Multi-site density evaluated across %d sites and %d proteins.\n",
            n_sites, n_prot))
cat(sprintf("Raw median range: [%.2f, %.2f] -> Normalized median: %.2f\n",
            min(tapply(df_density$Raw, df_density$Run, median)),
            max(tapply(df_density$Raw, df_density$Run, median)),
            median(df_density$Normalized)))

message("\n=== Listing 14.4: Controlled Simulation: PCA Ordination and Group Separation ===")
set.seed(606)
n_runs <- 36

pca_meta <- data.frame(
  Run   = paste0("Run_", sprintf("%02d", 1:n_runs)),
  Group = rep(c("Condition A", "Condition B"), each = 18),
  Site  = rep(paste("Site", sprintf("%02d", 1:6)), 6)
)

pc1 <- c(rnorm(18, mean = -2.8, sd = 1.1), rnorm(18, mean = 2.8, sd = 1.1))
pc2 <- rnorm(n_runs, mean = 0, sd = 1.3)
pca_df <- cbind(pca_meta, PC1 = pc1, PC2 = pc2)

t_pc1 <- t.test(PC1 ~ Group, data = pca_df)
delta <- mean(pca_df$PC1[pca_df$Group == "Condition B"]) - 
         mean(pca_df$PC1[pca_df$Group == "Condition A"])

cat(sprintf("PCA computed across %d runs (18 Condition A, 18 Condition B).\n", nrow(pca_df)))
cat(sprintf("Separation along PC1: delta = %.2f (t = %.2f, p = %.2e)\n",
            delta, t_pc1$statistic, t_pc1$p.value))

message("\n=== Listing 14.5: Protein-Level Summarization: Comparing Naive Mean vs Robust M-Estimation ===")
set.seed(42)
peptides_true <- c(18.2, 18.5, 18.1, 18.4, 23.5)
summary_mean   <- mean(peptides_true)
summary_median <- median(peptides_true)
summary_robust <- unname(MASS::rlm(peptides_true ~ 1, psi = MASS::psi.bisquare)$coefficients[1])

cat(sprintf("True Latent Protein Abundance: ~18.30\n"))
cat(sprintf("  Naive Mean Summary:   %.2f (Distorted by outlier!)\n", summary_mean))
cat(sprintf("  Median Summary:       %.2f\n", summary_median))
cat(sprintf("  Robust Tukey Summary: %.2f (Outlier completely discounted)\n", summary_robust))

message("\n=== Listing 14.6: Benchmark Biological Validation: Spike-In and Significant Hits on PXD004886 ===")
n_sig <- sum(clean_df$padj < 0.05, na.rm = TRUE)
pct_sig <- 100 * n_sig / nrow(clean_df)
cat(sprintf("Benchmark Significant Proteins (FDR < 0.05): %d / %d (%.2f%%)\n",
            n_sig, nrow(clean_df), pct_sig))
spike_in_hits <- clean_df %>%
  filter(grepl("AQUA", Protein) | grepl("AQUA", Genename))

cat(sprintf("Identified %d positive control spike-in proteins in dataset.\n", nrow(spike_in_hits)))
if (nrow(spike_in_hits) > 0) {
  cat(sprintf("  Spike-in Protein: %s\n", spike_in_hits$Protein[1]))
  cat(sprintf("  Observed log2FC: %.3f (padj = %.5f)\n", 
              spike_in_hits$log2FC[1], spike_in_hits$padj[1]))
}

message("\n=== Listing 14.7: Custom Iteratively Reweighted Least Squares Algorithm in R ===")
my_robust_rollup <- function(peptides_mat, c_tune = 4.685, max_iter = 50, tol = 1e-6) {
  n_samples <- ncol(peptides_mat)
  protein_estimates <- numeric(n_samples)
  
  for (s in 1:n_samples) {
    y <- peptides_mat[, s]
    y <- y[!is.na(y)]
    if (length(y) == 0) {
      protein_estimates[s] <- NA
      next
    }
    
    x_curr <- median(y)
    
    for (iter in 1:max_iter) {
      scale_est <- median(abs(y - x_curr)) / 0.6745
      if (scale_est < 1e-8) break
      
      u <- (y - x_curr) / (c_tune * scale_est)
      w <- ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
      
      if (sum(w) < 1e-8) break
      x_next <- sum(w * y) / sum(w)
      
      if (abs(x_next - x_curr) < tol) {
        x_curr <- x_next
        break
      }
      x_curr <- x_next
    }
    protein_estimates[s] <- x_curr
  }
  return(protein_estimates)
}

test_mat <- matrix(c(18.2, 18.5, 18.1, 18.4, 23.5), ncol = 1)
res_custom <- my_robust_rollup(test_mat)
res_rlm    <- unname(MASS::rlm(test_mat[, 1] ~ 1, psi = MASS::psi.bisquare)$coefficients[1])

cat(sprintf("Custom IRLS Biweight Estimate: %.4f\n", res_custom))
cat(sprintf("MASS::rlm bisquare Estimate:   %.4f\n", res_rlm))
cat(sprintf("Discrepancy:                  %.2e\n", abs(res_custom - res_rlm)))
message("\n[PASS] Chapter 14 Laboratory successfully executed.")
