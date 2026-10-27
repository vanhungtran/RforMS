# =============================================================================
# Chapter 24: Biomarker Modeling Without Data Leakage for MS-Based Omics
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# =============================================================================

suppressPackageStartupMessages({
  library(randomForest)
  library(survival)
  library(pROC)
  library(limma)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 24.1: Simulating High-Dimensional Clinical MS Proteomics Cohort
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 24.1: SIMULATING CLINICAL MS COHORT ===\n")

n_samples  <- 120
n_features <- 500
n_true_de  <- 15

# 1. Simulate null log2 intensity matrix
feature_means <- rnorm(n_features, mean = 22.0, sd = 1.5)
feature_sds   <- runif(n_features, min = 0.5, max = 0.9)

intensity_mat <- matrix(NA_real_, nrow = n_samples, ncol = n_features)
colnames(intensity_mat) <- paste0("Prot_", sprintf("%03d", 1:n_features))
rownames(intensity_mat) <- paste0("Sample_", sprintf("%03d", 1:n_samples))

for (j in seq_len(n_features)) {
  intensity_mat[, j] <- rnorm(n_samples, mean = feature_means[j], sd = feature_sds[j])
}

# 2. Simulate clinical metadata
classes <- factor(rep(c("Control", "Disease"), each = n_samples / 2))
ages    <- round(rnorm(n_samples, mean = 60, sd = 9))
batches <- factor(rep(paste0("Batch_", 1:4), length.out = n_samples))

# Inject subtle biological signal into first 15 proteins for Disease group
# log2 fold-changes range from 0.35 to 0.75
true_log2fc <- seq(0.35, 0.75, length.out = n_true_de)
for (j in seq_len(n_true_de)) {
  intensity_mat[classes == "Disease", j] <- intensity_mat[classes == "Disease", j] + true_log2fc[j]
}

# Simulate realistic survival time linked to the underlying disease status
base_hazard <- ifelse(classes == "Disease", 0.04, 0.015)
surv_time   <- round(rexp(n_samples, rate = base_hazard), 1)
surv_event  <- rbinom(n_samples, size = 1, prob = ifelse(classes == "Disease", 0.70, 0.35))

clinical_df <- data.frame(
  sample_id = rownames(intensity_mat),
  class     = classes,
  age       = ages,
  batch     = batches,
  time      = surv_time,
  status    = surv_event,
  stringsAsFactors = FALSE
)

cat(sprintf("Cohort dimensions: %d samples x %d features (%d true biomarker signals)\n",
            nrow(intensity_mat), ncol(intensity_mat), n_true_de))
print(head(clinical_df, 6))


# -----------------------------------------------------------------------------
# Listing 24.2: Leakage-Free Partitioning and Preprocessing
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 24.2: PARTITIONING AND PREPROCESSING ===\n")

# Split 75% training (90 samples) and 25% held-out test (30 samples)
train_idx <- sample(seq_len(n_samples), size = round(0.75 * n_samples))
test_idx  <- setdiff(seq_len(n_samples), train_idx)

train_raw  <- intensity_mat[train_idx, ]
test_raw   <- intensity_mat[test_idx,  ]
train_meta <- clinical_df[train_idx, ]
test_meta  <- clinical_df[test_idx,  ]

cat(sprintf("Training partition: %d samples (Control: %d, Disease: %d)\n",
            nrow(train_raw), sum(train_meta$class == "Control"), sum(train_meta$class == "Disease")))
cat(sprintf("Test partition:     %d samples (Control: %d, Disease: %d)\n",
            nrow(test_raw), sum(test_meta$class == "Control"), sum(test_meta$class == "Disease")))

# Feature scaling: parameters learned STRICTLY on training partition
train_centers <- colMeans(train_raw)
train_scales  <- apply(train_raw, 2, sd)
train_scales[train_scales == 0] <- 1

train_scaled <- scale(train_raw, center = train_centers, scale = train_scales)
test_scaled  <- scale(test_raw,  center = train_centers, scale = train_scales)

cat("Preprocessing verified: Test partition standardized using training mean and SD.\n")


# -----------------------------------------------------------------------------
# Listing 24.3: Nested CV vs Leaky Global Selection Benchmark
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 24.3: NESTED CV VS LEAKY BENCHMARK ===\n")

# A. Honest Nested Cross-Validation (Feature selection inside folds)
k_folds  <- 5
n_train  <- nrow(train_scaled)
fold_ids <- sample(rep(seq_len(k_folds), length.out = n_train))
nested_aucs <- numeric(k_folds)

for (k in seq_len(k_folds)) {
  in_fold  <- (fold_ids != k)
  val_fold <- (fold_ids == k)
  
  fold_train_x <- train_scaled[in_fold, ]
  fold_train_y <- train_meta$class[in_fold]
  fold_val_x   <- train_scaled[val_fold, ]
  fold_val_y   <- train_meta$class[val_fold]
  
  # Selection inside fold using two-sample t-tests
  p_vals <- apply(fold_train_x, 2, function(x) {
    t.test(x ~ fold_train_y)$p.value
  })
  selected_feats <- names(sort(p_vals))[1:15]
  
  # Fit Random Forest on selected features
  rf_fold <- randomForest(x = fold_train_x[, selected_feats], y = fold_train_y, ntree = 200)
  pred_probs <- predict(rf_fold, fold_val_x[, selected_feats], type = "prob")[, "Disease"]
  
  roc_fold <- roc(fold_val_y, pred_probs, quiet = TRUE)
  nested_aucs[k] <- as.numeric(auc(roc_fold))
}

cat(sprintf("Honest Nested CV AUC across 5 folds: Mean = %.3f, SD = %.3f\n",
            mean(nested_aucs), sd(nested_aucs)))

# B. Leaky Global Selection on Pure Noise (Catastrophic False Positive Demonstration)
# Generate pure Gaussian noise with NO true signal (n = 90, p = 1,000)
null_noise <- matrix(rnorm(n_train * 1000), nrow = n_train, ncol = 1000)
colnames(null_noise) <- paste0("Noise_", sprintf("%04d", 1:1000))

# Leaky step: select top 15 noise features on the entire dataset
p_leaky <- apply(null_noise, 2, function(x) t.test(x ~ train_meta$class)$p.value)
leaky_top_feats <- names(sort(p_leaky))[1:15]

# Cross-validate those pre-selected features
leaky_aucs <- numeric(k_folds)
for (k in seq_len(k_folds)) {
  in_fold  <- (fold_ids != k)
  val_fold <- (fold_ids == k)
  rf_leaky <- randomForest(x = null_noise[in_fold, leaky_top_feats], y = train_meta$class[in_fold], ntree = 200)
  pred_leaky <- predict(rf_leaky, null_noise[val_fold, leaky_top_feats], type = "prob")[, "Disease"]
  leaky_aucs[k] <- as.numeric(auc(roc(train_meta$class[val_fold], pred_leaky, quiet = TRUE)))
}
cat(sprintf("Catastrophic Leaky Selection on Pure Noise: Apparent AUC = %.3f (True AUC = 0.500)\n",
            mean(leaky_aucs)))


# -----------------------------------------------------------------------------
# Listing 24.4: Final Model Training, Test Evaluation, and ROC Analysis
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 24.4: FINAL MODEL EVALUATION ON HELD-OUT TEST SET ===\n")

# Feature selection on the COMPLETE training set
p_train_full <- apply(train_scaled, 2, function(x) t.test(x ~ train_meta$class)$p.value)
final_biomarkers <- names(sort(p_train_full))[1:15]

cat("Top 10 Biomarkers Selected on Training Set:\n")
print(head(data.frame(
  Feature = final_biomarkers,
  pValue  = round(p_train_full[final_biomarkers], 6),
  row.names = NULL
), 10))

# Train final classifier
final_rf <- randomForest(
  x = train_scaled[, final_biomarkers],
  y = train_meta$class,
  ntree = 500,
  importance = TRUE
)

# Predict on held-out test set
test_probs <- predict(final_rf, test_scaled[, final_biomarkers], type = "prob")[, "Disease"]
test_preds <- predict(final_rf, test_scaled[, final_biomarkers], type = "response")

# Confusion matrix and clinical accuracy metrics
conf_mat <- table(Predicted = test_preds, Actual = test_meta$class)
cat("Test Set Confusion Matrix:\n")
print(conf_mat)

tp <- conf_mat["Disease", "Disease"]
tn <- conf_mat["Control", "Control"]
fp <- conf_mat["Disease", "Control"]
fn <- conf_mat["Control", "Disease"]

sensitivity <- tp / (tp + fn)
specificity <- tn / (tn + fp)
accuracy    <- (tp + tn) / sum(conf_mat)
brier_score <- mean((test_probs - as.numeric(test_meta$class == "Disease"))^2)

cat(sprintf("Clinical Accuracy: %.1f%% | Sensitivity: %.1f%% | Specificity: %.1f%%\n",
            accuracy * 100, sensitivity * 100, specificity * 100))
cat(sprintf("Probability Brier Score: %.4f (lower is better, 0.25 is random chance)\n",
            brier_score))

# ROC and DeLong 95% Confidence Interval
roc_test <- roc(test_meta$class, test_probs, quiet = TRUE)
auc_ci   <- ci.auc(roc_test, method = "delong")
cat(sprintf("Held-Out Test AUC: %.3f (DeLong 95%% CI: [%.3f, %.3f])\n",
            as.numeric(auc(roc_test)), auc_ci[1], auc_ci[3]))


# -----------------------------------------------------------------------------
# Listing 24.5: Clinical Survival Prognostication and Multivariable Cox Model
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 24.5: SURVIVAL AND PROGNOSTIC STRATIFICATION ===\n")

# Augment test metadata with model score
test_meta$model_score <- test_probs
test_meta$risk_group  <- factor(
  ifelse(test_probs >= median(test_probs), "High_Risk", "Low_Risk"),
  levels = c("Low_Risk", "High_Risk")
)

# Kaplan-Meier survival curves
surv_obj <- Surv(time = test_meta$time, event = test_meta$status)
km_fit   <- survfit(surv_obj ~ risk_group, data = test_meta)

# Log-rank test
logrank_res <- survdiff(surv_obj ~ risk_group, data = test_meta)
logrank_p   <- 1 - pchisq(logrank_res$chisq, df = 1)
cat(sprintf("Log-Rank Test for Survival Stratification: ChiSq = %.3f, p-value = %.4f\n",
            logrank_res$chisq, logrank_p))

# Multivariable Cox Proportional Hazards Model adjusting for age and batch
cox_fit <- coxph(surv_obj ~ model_score + age + batch, data = test_meta)
cox_sum <- summary(cox_fit)
cat("Multivariable Cox Regression Results (Test Set):\n")
print(round(cox_sum$coefficients[, c("coef", "exp(coef)", "se(coef)", "Pr(>|z|)")], 4))
cat(sprintf("Model Score Hazard Ratio: HR = %.2f (95%% CI: [%.2f, %.2f], p = %.4f)\n",
            cox_sum$conf.int["model_score", "exp(coef)"],
            cox_sum$conf.int["model_score", "lower .95"],
            cox_sum$conf.int["model_score", "upper .95"],
            cox_sum$coefficients["model_score", "Pr(>|z|)"]))


# -----------------------------------------------------------------------------
# Listing 24.6: Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 24.1: Quantifying optimistic bias from global t-tests on training set
p_all <- apply(intensity_mat, 2, function(x) t.test(x ~ clinical_df$class)$p.value)
leaky_feats_all <- names(sort(p_all))[1:15]

# Leaky model fit
rf_leaky_full <- randomForest(x = intensity_mat[train_idx, leaky_feats_all], 
                              y = clinical_df$class[train_idx], ntree = 300)
pred_leaky_test <- predict(rf_leaky_full, intensity_mat[test_idx, leaky_feats_all], type = "prob")[, "Disease"]
leaky_auc_test  <- as.numeric(auc(roc(clinical_df$class[test_idx], pred_leaky_test, quiet = TRUE)))

cat(sprintf("Exercise 24.1: Honest Test AUC = %.3f vs Leaky Selection Test AUC = %.3f (Delta AUC = %+.3f)\n",
            as.numeric(auc(roc_test)), leaky_auc_test, leaky_auc_test - as.numeric(auc(roc_test))))

# Exercise 24.2: limma empirical Bayes selection inside CV
design_cv <- model.matrix(~ train_meta$class)
fit_cv <- eBayes(lmFit(t(train_scaled), design_cv))
limma_selected <- rownames(topTable(fit_cv, coef = 2, number = 15, sort.by = "P"))
overlap_count  <- length(intersect(final_biomarkers, limma_selected))
cat(sprintf("Exercise 24.2: Overlap between Student's t and limma eBayes top 15: %d / 15 features (%.1f%%)\n",
            overlap_count, overlap_count / 15 * 100))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
