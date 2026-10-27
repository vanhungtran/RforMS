# =============================================================================
# Chapter 20: Experimental Design, Replication, and Power in Mass Spectrometry
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(pwr)
})

# Set seed for reproducibility across all simulations
set.seed(42)

# -----------------------------------------------------------------------------
# Listing 20.1: Randomized Complete Block Design (RCBD) and Injection Sequence
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 20.1: RCBD AND INJECTION SEQUENCE ===\n")

# Define biological cohort: 24 Control and 24 Disease subjects
n_ctrl <- 24
n_dise <- 24
n_total <- n_ctrl + n_dise

samples <- data.frame(
  subject_id = paste0("SUBJ_", sprintf("%03d", 1:n_total)),
  condition  = rep(c("Control", "Disease"), c(n_ctrl, n_dise)),
  stringsAsFactors = FALSE
)

# Allocate into 2 balanced digestion plates (RCBD: 12 Control, 12 Disease per plate)
ctrl_idx <- which(samples$condition == "Control")
dise_idx <- which(samples$condition == "Disease")

plate1_idx <- c(sample(ctrl_idx, 12), sample(dise_idx, 12))
plate2_idx <- setdiff(1:n_total, plate1_idx)

samples$plate <- NA_character_
samples$plate[plate1_idx] <- "Plate_1"
samples$plate[plate2_idx] <- "Plate_2"

# Randomize injection order within each plate block
seq_list <- list()
qc_idx <- 1

# 1. Five pre-run column conditioning injections
for (k in 1:5) {
  seq_list[[length(seq_list) + 1]] <- data.frame(
    run_order   = length(seq_list) + 1,
    sample_id   = paste0("QC_COND_", k),
    sample_type = "QC_Conditioning",
    condition   = "QC",
    plate       = "None",
    stringsAsFactors = FALSE
  )
}

# 2. Sequential block processing: Plate 1 then Plate 2
for (plt in c("Plate_1", "Plate_2")) {
  plt_samples <- samples[samples$plate == plt, ]
  plt_samples <- plt_samples[sample(nrow(plt_samples)), ]  # randomize within plate
  
  for (i in seq_len(nrow(plt_samples))) {
    seq_list[[length(seq_list) + 1]] <- data.frame(
      run_order   = length(seq_list) + 1,
      sample_id   = plt_samples$subject_id[i],
      sample_type = "Biological",
      condition   = plt_samples$condition[i],
      plate       = plt,
      stringsAsFactors = FALSE
    )
    
    # Interleave 1 pooled QC every 8 biological injections
    if (i %% 8 == 0) {
      seq_list[[length(seq_list) + 1]] <- data.frame(
        run_order   = length(seq_list) + 1,
        sample_id   = paste0("QC_POOL_", sprintf("%02d", qc_idx)),
        sample_type = "QC_Pool",
        condition   = "QC",
        plate       = plt,
        stringsAsFactors = FALSE
      )
      qc_idx <- qc_idx + 1
    }
  }
}

injection_schedule <- do.call(rbind, seq_list)
cat(sprintf("Total injections scheduled: %d (Biological: %d, QC Pool: %d, Conditioning: %d)\n",
            nrow(injection_schedule),
            sum(injection_schedule$sample_type == "Biological"),
            sum(injection_schedule$sample_type == "QC_Pool"),
            sum(injection_schedule$sample_type == "QC_Conditioning")))
print(head(injection_schedule, 10))
print(table(injection_schedule$condition, injection_schedule$plate))


# -----------------------------------------------------------------------------
# Listing 20.2: Pre-Acquisition Confounder Audit and Balance Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 20.2: CONFOUNDER AUDIT ===\n")

# Extract biological samples only for statistical confounding audit
bio_sched <- injection_schedule[injection_schedule$sample_type == "Biological", ]

# 1. Contingency table check (Plate vs Condition)
cont_tbl <- table(bio_sched$plate, bio_sched$condition)
cat("Contingency Table (Plate x Condition):\n")
print(cont_tbl)

# Chi-squared test of independence between plate and condition
chisq_res <- chisq.test(cont_tbl)
cat(sprintf("Independence Chi-Square: X2 = %.3f, df = %d, p-value = %.4f\n",
            chisq_res$statistic, chisq_res$parameter, chisq_res$p.value))

# 2. Pearson correlation between biological condition and run order
cond_binary <- as.numeric(bio_sched$condition == "Disease")
cor_res <- cor.test(cond_binary, bio_sched$run_order)
cat(sprintf("Run-Order Correlation: r = %.4f, 95%% CI: [%.4f, %.4f], p-value = %.4f\n",
            cor_res$estimate, cor_res$conf.int[1], cor_res$conf.int[2], cor_res$p.value))

# 3. Two-way ANOVA: Run order explained by Condition and Plate
aov_fit <- aov(run_order ~ condition + plate, data = bio_sched)
cat("ANOVA of Run Order on Experimental Factors:\n")
print(summary(aov_fit))

# Automated assertion gates: pass/fail criteria
audit_passed <- (abs(cor_res$estimate) < 0.15) && (summary(aov_fit)[[1]]["condition", "Pr(>F)"] > 0.05)
cat(sprintf("Confounder Audit Status: %s\n", if (audit_passed) "PASSED (No Confounding)" else "FAILED (Confounded)"))


# -----------------------------------------------------------------------------
# Listing 20.3: Analytical Sample Size and FDR Multiplicity Adjustments
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 20.3: ANALYTICAL POWER CALCULATIONS ===\n")

# Parameters
delta <- 1.0       # log2 fold change = 1.0 (2-fold difference)
sigma <- 0.70      # log2 standard deviation
d_cohen <- delta / sigma
target_power <- 0.80

# 1. Single-hypothesis testing at nominal alpha = 0.05
res_single <- pwr.t.test(d = d_cohen, sig.level = 0.05, power = target_power, type = "two.sample")
n_single <- ceiling(res_single$n)
cat(sprintf("Single-test (alpha = 0.05): Cohen's d = %.3f -> Required n = %d per group (exact: %.2f)\n",
            d_cohen, n_single, res_single$n))

# 2. Multiple comparisons adjustment for M = 5,000 features
m_features <- 5000
alpha_bonf <- 0.05 / m_features
res_bonf <- pwr.t.test(d = d_cohen, sig.level = alpha_bonf, power = target_power, type = "two.sample")
n_bonf <- ceiling(res_bonf$n)
cat(sprintf("Multiplicity (M = %d, alpha* = %.1e): Required n = %d per group (exact: %.2f)\n",
            m_features, alpha_bonf, n_bonf, res_bonf$n))

# 3. Missing data inflation adjustment (p_miss = 0.25)
p_miss <- 0.25
n_inflated <- ceiling(n_bonf / (1 - p_miss))
cat(sprintf("Missingness-inflated (p_miss = %.2f): Enrolled n* = %d per group\n",
            p_miss, n_inflated))

# 4. Multi-parameter power table across effect sizes
deltas <- c(0.5, 0.75, 1.0, 1.25, 1.5, 2.0)
power_table <- data.frame(
  Delta_log2FC = deltas,
  Fold_Change  = round(2^deltas, 2),
  Cohens_d     = round(deltas / sigma, 3),
  n_Nominal_05 = sapply(deltas, function(dl) ceiling(pwr.t.test(d = dl/sigma, sig.level = 0.05, power = 0.80)$n)),
  n_FDR_5000   = sapply(deltas, function(dl) ceiling(pwr.t.test(d = dl/sigma, sig.level = 0.05/5000, power = 0.80)$n)),
  n_With_25Pct_Miss = sapply(deltas, function(dl) ceiling(ceiling(pwr.t.test(d = dl/sigma, sig.level = 0.05/5000, power = 0.80)$n) / (1 - 0.25)))
)
print(power_table)


# -----------------------------------------------------------------------------
# Listing 20.4: Monte Carlo Power Simulation with Empirical Bayes Moderation
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 20.4: MONTE CARLO POWER SIMULATION ===\n")

simulate_eb_power <- function(n_per_group, n_features = 2000, n_de = 100, 
                              log2fc = 1.0, sigma_mean = 0.60, prior_df = 4, 
                              n_sim = 40) {
  power_limma <- numeric(n_sim)
  power_ttest <- numeric(n_sim)
  
  for (s in seq_len(n_sim)) {
    group <- factor(rep(c("Control", "Treated"), each = n_per_group))
    design <- model.matrix(~ group)
    
    # Feature-specific variances drawn from scaled inverse chi-square
    feature_sd <- sqrt(sigma_mean^2 * prior_df / rchisq(n_features, df = prior_df))
    
    # Simulate null log2 intensity matrix
    mat <- matrix(rnorm(n_features * 2 * n_per_group, mean = 20.0, 
                        sd = rep(feature_sd, 2 * n_per_group)), 
                  nrow = n_features, byrow = FALSE)
    
    # Inject true differences into the first n_de features
    mat[1:n_de, group == "Treated"] <- mat[1:n_de, group == "Treated"] + log2fc
    
    # 1. limma empirical Bayes moderated fit
    fit_eb <- eBayes(lmFit(mat, design))
    p_limma <- fit_eb$p.value[, 2]
    q_limma <- p.adjust(p_limma, method = "BH")
    power_limma[s] <- mean(q_limma[1:n_de] < 0.05)
    
    # 2. Ordinary unmoderated two-sample t-tests
    # Ordinary t-statistic: coef / (stdev.unscaled * sigma)
    t_ordinary <- fit_eb$coefficients[, 2] / (fit_eb$stdev.unscaled[, 2] * fit_eb$sigma)
    df_ordinary <- fit_eb$df.residual[1]
    p_ordinary <- 2 * pt(-abs(t_ordinary), df = df_ordinary)
    q_ordinary <- p.adjust(p_ordinary, method = "BH")
    power_ttest[s] <- mean(q_ordinary[1:n_de] < 0.05)
  }
  
  c(limma_power = mean(power_limma), ttest_power = mean(power_ttest))
}

sample_sizes <- c(3, 5, 8, 12)
cat("Simulating empirical power across sample sizes (n = 3, 5, 8, 12)...\n")
sim_matrix <- sapply(sample_sizes, function(n) {
  simulate_eb_power(n_per_group = n, n_features = 1500, n_de = 75, log2fc = 1.0, 
                    sigma_mean = 0.60, prior_df = 4, n_sim = 30)
})

sim_df <- data.frame(
  Replicates_Per_Group = sample_sizes,
  limma_eBayes_Power   = round(sim_matrix["limma_power", ], 3),
  Ordinary_tTest_Power = round(sim_matrix["ttest_power", ], 3),
  Power_Gain_Pct       = round((sim_matrix["limma_power", ] - sim_matrix["ttest_power", ]) * 100, 1)
)
print(sim_df)


# -----------------------------------------------------------------------------
# Listing 20.5: Duplicate Correlation Modeling of Technical Replicates
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 20.5: DUPLICATE CORRELATION MODELING ===\n")

# Simulate 6 biological subjects per group with 3 technical replicates each (N = 36)
n_bio <- 6
n_tech <- 3
n_groups <- 2
total_samples <- n_bio * n_groups * n_tech  # 36

bio_subject <- factor(rep(paste0("Subject_", 1:(n_bio * n_groups)), each = n_tech))
bio_group   <- factor(rep(c("Control", "Treated"), each = n_bio * n_tech))

# Simulating 500 features
n_feat <- 500
sigma_bio  <- 0.50
sigma_tech <- 0.20
true_icc   <- sigma_bio^2 / (sigma_bio^2 + sigma_tech^2)
cat(sprintf("Theoretical Biological SD = %.2f, Technical SD = %.2f -> True ICC = %.4f\n",
            sigma_bio, sigma_tech, true_icc))

# Simulate intensity data with subject random effects
subj_rand_eff <- matrix(rnorm(n_feat * (n_bio * n_groups), mean = 0, sd = sigma_bio), 
                        nrow = n_feat, ncol = n_bio * n_groups)
tech_noise     <- matrix(rnorm(n_feat * total_samples, mean = 0, sd = sigma_tech), 
                        nrow = n_feat, ncol = total_samples)

intensity_mat <- matrix(20.0, nrow = n_feat, ncol = total_samples)
for (s in seq_len(n_bio * n_groups)) {
  col_indices <- which(bio_subject == paste0("Subject_", s))
  intensity_mat[, col_indices] <- intensity_mat[, col_indices] + subj_rand_eff[, s]
}
intensity_mat <- intensity_mat + tech_noise

# Inject treatment effect into first 50 features
intensity_mat[1:50, bio_group == "Treated"] <- intensity_mat[1:50, bio_group == "Treated"] + 0.80

# Design matrix for treatment group
design <- model.matrix(~ bio_group)

# 1. Estimate consensus technical correlation across features
dup_cor <- duplicateCorrelation(intensity_mat, design, block = bio_subject)
cat(sprintf("limma Estimated Consensus Intraclass Correlation: rho = %.4f (atanh = %.4f)\n",
            dup_cor$consensus.correlation, atanh(dup_cor$consensus.correlation)))

# 2. Fit linear model with blocking on biological subject
fit_blocked <- lmFit(intensity_mat, design, block = bio_subject, correlation = dup_cor$consensus)
fit_blocked <- eBayes(fit_blocked)

# 3. Fit naive linear model treating technical replicates as independent (PSEUDOREPLICATION)
fit_naive <- lmFit(intensity_mat, design)
fit_naive <- eBayes(fit_naive)

# Compare degrees of freedom and standard error estimates
cat(sprintf("Blocked Model:   Residual df = %.1f, Mean SE(coef) = %.4f\n",
            fit_blocked$df.residual[1], mean(fit_blocked$stdev.unscaled[, 2] * fit_blocked$sigma)))
cat(sprintf("Naive Model:     Residual df = %.1f, Mean SE(coef) = %.4f\n",
            fit_naive$df.residual[1], mean(fit_naive$stdev.unscaled[, 2] * fit_naive$sigma)))

# Number of features called significant at 5% FDR
q_blocked <- p.adjust(fit_blocked$p.value[, 2], method = "BH")
q_naive   <- p.adjust(fit_naive$p.value[, 2], method = "BH")

cat(sprintf("Blocked Model: %d features significant at 5%% FDR (True Positives: %d / 50)\n",
            sum(q_blocked < 0.05), sum(q_blocked[1:50] < 0.05)))
cat(sprintf("Naive Model:   %d features significant at 5%% FDR (True Positives: %d / 50, False Positives: %d)\n",
            sum(q_naive < 0.05), sum(q_naive[1:50] < 0.05), sum(q_naive[51:500] < 0.05)))


# -----------------------------------------------------------------------------
# Listing 20.6: Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISE SOLUTIONS VERIFICATION ===\n")

# Exercise 20.2: Analytical Sample Size and FDR Inflation Function
calc_ms_sample_size <- function(delta, sigma, m_features = 8000, power = 0.80, fdr = 0.05, p_miss = 0.20) {
  d_eff <- delta / sigma
  alpha_nominal <- 0.05
  alpha_bonf    <- fdr / m_features
  
  n_single <- ceiling(pwr.t.test(d = d_eff, sig.level = alpha_nominal, power = power)$n)
  n_mcomp  <- ceiling(pwr.t.test(d = d_eff, sig.level = alpha_bonf, power = power)$n)
  n_final  <- ceiling(n_mcomp / (1 - p_miss))
  
  data.frame(
    Delta = delta,
    Sigma = sigma,
    Cohens_d = round(d_eff, 3),
    n_Single_Hypothesis = n_single,
    n_Bonferroni_FDR = n_mcomp,
    n_Enrolled_With_Missingness = n_final
  )
}

# Exercise 20.2 Part 1 & 2: delta = log2(1.5) = 0.585, sigma = 0.65, M = 8000
ex2_res <- calc_ms_sample_size(delta = log2(1.5), sigma = 0.65, m_features = 8000, power = 0.80, p_miss = 0.25)
cat("Exercise 20.2 Output:\n")
print(ex2_res)

# Exercise 20.3: ICC and Variance Allocation
# Cost: c_bio = $150, c_tech = $50. Total budget C = $3,000 per group.
# sigma_bio = 0.80, sigma_tech = 0.30.
# Optimal m* = (sigma_tech / sigma_bio) * sqrt(c_bio / c_tech)
sigma_b <- 0.80
sigma_t <- 0.30
c_b <- 150
c_t <- 50
m_opt <- (sigma_t / sigma_b) * sqrt(c_b / c_t)
cat(sprintf("\nExercise 20.3 Output:\nOptimal technical replicates m* = %.3f -> round to m = %d\n",
            m_opt, round(m_opt)))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
