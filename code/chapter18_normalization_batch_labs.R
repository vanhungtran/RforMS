# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 18 Laboratory: Feature Normalization and Batch Correction
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(sva)
  library(dplyr)
  library(ggplot2)
})

message("\n=== Listing 18.1: Simulating Multi-Batch Proteomic Matrix with Confounding Risk ===")
set.seed(42)
n_features <- 500
n_samples  <- 30

condition <- factor(rep(c("Control", "Treatment"), each = 15))
batch     <- factor(rep(rep(c("B1", "B2", "B3"), each = 5), 2))

true_baseline <- matrix(rnorm(n_features * n_samples, mean = 20, sd = 1.5), 
                        nrow = n_features, ncol = n_samples)
de_indices <- 1:50
for (j in 16:30) {
  true_baseline[de_indices, j] <- true_baseline[de_indices, j] + 1.8
}

batch_loc <- c(B1 = 0.0, B2 = 2.2, B3 = -1.8)
batch_scl <- c(B1 = 1.0, B2 = 1.6, B3 = 0.7)

raw_mat <- true_baseline
for (s_idx in 1:n_samples) {
  b <- as.character(batch[s_idx])
  raw_mat[, s_idx] <- (raw_mat[, s_idx] + batch_loc[b]) * batch_scl[b]
}
colnames(raw_mat) <- paste(condition, batch, 1:n_samples, sep = "_")
rownames(raw_mat) <- paste0("PROT_", 1:n_features)

cat(sprintf("Simulated %d proteins across %d samples in 3 distinct batches.\n", 
            n_features, n_samples))

message("\n=== Listing 18.2: Empirical Bayes Batch Adjustment with ComBat ===")
mod_bio <- model.matrix(~ condition)

combat_mat <- ComBat(
  dat = raw_mat,
  batch = batch,
  mod = mod_bio,
  par.prior = TRUE,
  prior.plots = FALSE
)

pca_raw <- prcomp(t(raw_mat), scale. = TRUE)
pca_adj <- prcomp(t(combat_mat), scale. = TRUE)

r2_raw_batch <- summary(lm(pca_raw$x[, 1] ~ batch))$r.squared
r2_raw_cond  <- summary(lm(pca_raw$x[, 1] ~ condition))$r.squared

r2_adj_batch <- summary(lm(pca_adj$x[, 1] ~ batch))$r.squared
r2_adj_cond  <- summary(lm(pca_adj$x[, 1] ~ condition))$r.squared

cat(sprintf("PCA Diagnostics (PC1 Variance Attribution):\n"))
cat(sprintf("BEFORE ComBat: Batch R^2 = %.3f | Condition R^2 = %.3f\n", 
            r2_raw_batch, r2_raw_cond))
cat(sprintf("AFTER  ComBat: Batch R^2 = %.3f | Condition R^2 = %.3f\n", 
            r2_adj_batch, r2_adj_cond))

message("\n=== Listing 18.3: Confounding Audit via the Design Matrix Condition Number ===")
audit_batch_confounding <- function(metadata, condition_col, batch_col) {
  cond  <- factor(metadata[[condition_col]])
  batch <- factor(metadata[[batch_col]])
  
  cross_tab <- table(Condition = cond, Batch = batch)
  cat("--- Experimental Cross-Tabulation Matrix ---\n")
  print(cross_tab)
  
  X_joint <- model.matrix(~ cond + batch)
  svd_decomp <- svd(X_joint)
  singular_vals <- svd_decomp$d
  kappa <- max(singular_vals) / min(singular_vals)
  
  cat(sprintf("\nDesign Matrix Condition Number: %.2f\n", kappa))
  if (kappa > 30) {
    cat("[CRITICAL WARNING] Severe batch-condition confounding detected (kappa > 30)!\n")
  } else {
    cat("[STATUS: PASS] Design is well-balanced across batches (kappa <= 30).\n")
  }
  return(invisible(list(cross_tab = cross_tab, kappa = kappa)))
}

demo_meta <- data.frame(
  Condition = rep(c("Control", "Treatment"), each = 15),
  Batch     = rep(rep(c("B1", "B2", "B3"), each = 5), 2)
)

audit_batch_confounding(demo_meta, "Condition", "Batch")

message("\n=== Listing 18.4: Dilution Series Benchmarking of Normalization Linearity ===")
set.seed(123)
dilutions <- c(1, 2, 5, 10)
n_reps <- 3
dilution_vector <- rep(dilutions, each = n_reps)
n_cols <- length(dilution_vector)

base_conc <- runif(100, 1e4, 1e7)
synth_mat <- matrix(0, nrow = 100, ncol = n_cols)

for (col_idx in 1:n_cols) {
  dil <- dilution_vector[col_idx]
  mult_err <- rnorm(100, mean = 1, sd = 0.12)
  add_err  <- rnorm(100, mean = 0, sd = 5000)
  synth_mat[, col_idx] <- pmax(10, base_conc * dil * mult_err + add_err)
}

norm_log2 <- log2(synth_mat)
norm_med <- sweep(norm_log2, 2, apply(norm_log2, 2, median) - median(norm_log2))

cor_raw <- apply(norm_log2, 1, function(x) cor(x, log2(dilution_vector)))
cor_med <- apply(norm_med, 1, function(x) cor(x, log2(dilution_vector)))

cat("--- Dilution Series Linearity Preservation (Mean Pearson r) ---\n")
cat(sprintf("  Raw Log2:         r = %.4f (SD = %.4f)\n", mean(cor_raw), sd(cor_raw)))
cat(sprintf("  Median Centered:  r = %.4f (SD = %.4f)\n", mean(cor_med), sd(cor_med)))
message("\n[PASS] Chapter 18 Laboratory successfully executed.")
