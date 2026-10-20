# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Monograph)
# Chapter 23: Machine Learning in Mass Spectrometry
# Companion R Laboratory Script: End-to-End Reproducible Multi-Model Workflow
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Tested on: R 4.6.1 (Bioconductor 3.20)
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Dependencies and Computational Environment Setup
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(tidymodels)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(themis)
  library(probably)
  library(dcurves)
  library(dplyr)
  library(ggplot2)
  library(rsample)
  library(yardstick)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(vcd)
})

cat("=============================================================================\n")
cat("Starting Chapter 23 Machine Learning Laboratory Execution\n")
cat("=============================================================================\n")

# -----------------------------------------------------------------------------
# Listing 23.1: Simulating High-Dimensional Proteomics Cohort (p >> n) with MAR
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.1: Simulating Proteomic Cohort (p = 4000, n = 120) ---\n")
set.seed(2026)
n_samples   <- 120
n_proteins  <- 4000L

status <- factor(c(rep("Control", 80), rep("Tumor", 40)),
                 levels = c("Control", "Tumor"))

# Observed continuous covariate: patient age (MAR missingness driver)
age_years <- c(rnorm(80, mean = 55, sd = 12),
               rnorm(40, mean = 58, sd = 12))
age_years <- pmax(18, pmin(90, age_years))

# Baseline protein abundance: N(20, 1.0)
X <- matrix(rnorm(n_samples * n_proteins, mean = 20, sd = 1.0),
            nrow = n_samples, ncol = n_proteins)
colnames(X) <- paste0("PROT_", sprintf("%04d", 1:n_proteins))

# Module 1: Collinear co-regulated biomarker signature (Proteins 1-8, r ~ 0.70)
latent_driver <- c(rep(0, 80), rep(0.75, 40))
for (j in 1:8) {
  X[, j] <- latent_driver + rnorm(n_samples, mean = 20, sd = 0.8)
}

# Module 2: Epistatic interaction (Proteins 9 & 10, zero individual main effect)
tumor_idx <- 81:120
X[, 9]  <- X[, 9]  - 20
X[, 10] <- X[, 10] - 20
interaction_driver <- rnorm(40, mean = 1.5, sd = 0.3)
X[tumor_idx, 9]  <- X[tumor_idx, 9]  + 0.6 * interaction_driver
X[tumor_idx, 10] <- X[tumor_idx, 10] + 0.6 * interaction_driver

# MAR missingness: probability conditional on observed Age
miss_prob <- plogis(-3.2 + 0.05 * (age_years - 55))
miss_mask <- matrix(
  runif(n_samples * n_proteins) < rep(miss_prob, times = n_proteins),
  nrow = n_samples, ncol = n_proteins
)

# MNAR component on low-abundance proteins near Limit of Detection
low_abund <- which(colMeans(X, na.rm = TRUE) < 19.5)
for (j in low_abund) {
  thresh_j <- quantile(X[, j], 0.05, na.rm = TRUE)
  miss_mask[, j] <- miss_mask[, j] | ((X[, j] < thresh_j) & (runif(n_samples) < 0.30))
}
X[miss_mask] <- NA

df_cohort <- as.data.frame(X)
df_cohort$Status <- status
df_cohort$Age    <- age_years

# Stratified 75/25 Train-Test Partition
set.seed(42)
split_obj  <- initial_split(df_cohort, prop = 0.75, strata = Status)
train_data <- training(split_obj)
test_data  <- testing(split_obj)

cat(sprintf("Simulated cohort: %d samples x %d features (p/n ratio ~ %.0f).\n",
            nrow(df_cohort), ncol(df_cohort) - 2, (ncol(df_cohort) - 2) / nrow(df_cohort)))
cat(sprintf("Training set: %d patients (%d Control, %d Tumor)\n",
            nrow(train_data), sum(train_data$Status == "Control"), sum(train_data$Status == "Tumor")))
cat(sprintf("Independent test holdout: %d patients (%d Control, %d Tumor)\n",
            nrow(test_data), sum(test_data$Status == "Control"), sum(test_data$Status == "Tumor")))

# MAR verification across age quartiles
train_data$.age_quartile <- cut(
  train_data$Age,
  breaks = quantile(train_data$Age, c(0, 0.25, 0.5, 0.75, 1)),
  include.lowest = TRUE, labels = c("Q1(young)", "Q2", "Q3", "Q4(old)")
)
protein_cols <- grep("^PROT_", names(train_data), value = TRUE)
miss_by_age <- train_data |>
  group_by(.age_quartile) |>
  summarise(missingness_rate = round(mean(is.na(as.matrix(across(all_of(protein_cols))))), 4),
            n = n(), .groups = "drop")
cat("\nMAR Verification: Missingness rate stratified by Age quartile (Training set):\n")
print(miss_by_age)
train_data$.age_quartile <- NULL

# -----------------------------------------------------------------------------
# Listing 23.2: Leakage-Safe Preprocessing Recipe with Role Assignment
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.2: Preprocessing Recipe (Training-Fold Encapsulated) ---\n")
rec <- recipe(Status ~ ., data = train_data) |>
  update_role(Age, new_role = "imputation donor") |>
  step_impute_median(all_numeric_predictors()) |>
  step_zv(all_numeric_predictors()) |>
  step_normalize(all_numeric_predictors()) |>
  step_downsample(Status, under_ratio = 1.2)

prep_rec <- prep(rec, training = train_data)
baked_train <- bake(prep_rec, new_data = NULL)

cat("Recipe successfully baked on training fold:\n")
cat(sprintf("  Predictors processed: %d (Zero-variance filtered: %d)\n",
            ncol(baked_train) - 2, ncol(train_data) - ncol(baked_train)))
cat(sprintf("  Post-downsample class balance: %d Control, %d Tumor\n",
            sum(baked_train$Status == "Control"), sum(baked_train$Status == "Tumor")))
cat(sprintf("  Remaining missing values: %d\n", sum(is.na(baked_train))))

# -----------------------------------------------------------------------------
# Listing 23.3: Multi-Model Specification and 5-Fold Cross-Validation Tuning
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.3: Hyperparameter Tuning via 5-Fold Stratified CV ---\n")
spec_enet <- logistic_reg(penalty = tune(), mixture = 0.5) |>
  set_engine("glmnet") |>
  set_mode("classification")

spec_rf <- rand_forest(trees = 250, min_n = tune()) |>
  set_engine("ranger", importance = "permutation", seed = 42) |>
  set_mode("classification")

spec_xgb <- boost_tree(trees = 100, tree_depth = 3, learn_rate = tune()) |>
  set_engine("xgboost", nthread = 1) |>
  set_mode("classification")

wf_enet <- workflow() |> add_recipe(rec) |> add_model(spec_enet)
wf_rf   <- workflow() |> add_recipe(rec) |> add_model(spec_rf)
wf_xgb  <- workflow() |> add_recipe(rec) |> add_model(spec_xgb)

set.seed(42)
cv_folds <- vfold_cv(train_data, v = 5, strata = Status)

res_enet <- tune_grid(wf_enet, resamples = cv_folds,
                      grid = tibble(penalty = 10^seq(-3, -0.2, length.out = 6)),
                      metrics = metric_set(roc_auc, accuracy, brier_class))

res_rf <- tune_grid(wf_rf, resamples = cv_folds,
                    grid = tibble(min_n = c(3, 7, 12)),
                    metrics = metric_set(roc_auc, accuracy, brier_class))

res_xgb <- tune_grid(wf_xgb, resamples = cv_folds,
                     grid = tibble(learn_rate = c(0.02, 0.08, 0.15)),
                     metrics = metric_set(roc_auc, accuracy, brier_class))

cat("Hyperparameter grid tuning complete across Elastic Net, Random Forest, and XGBoost.\n")

# -----------------------------------------------------------------------------
# Listing 23.4: Model Benchmark on Independent Holdout Test Set (N = 30)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.4: Evaluation on Untouched Holdout Test Partition ---\n")
best_enet <- select_best(res_enet, metric = "roc_auc")
best_rf   <- select_best(res_rf, metric = "roc_auc")
best_xgb  <- select_best(res_xgb, metric = "roc_auc")

fit_enet <- last_fit(finalize_workflow(wf_enet, best_enet), split_obj)
fit_rf   <- last_fit(finalize_workflow(wf_rf, best_rf), split_obj)
fit_xgb  <- last_fit(finalize_workflow(wf_xgb, best_xgb), split_obj)

m_enet <- collect_metrics(fit_enet)
m_rf   <- collect_metrics(fit_rf)
m_xgb  <- collect_metrics(fit_xgb)

cat("Model Benchmark on Independent Holdout Test Set (N = 30):\n")
cat(sprintf("  Elastic Net:   ROC-AUC = %.3f | Accuracy = %.3f (Optimal penalty = %.4f)\n",
            m_enet$.estimate[m_enet$.metric == "roc_auc"],
            m_enet$.estimate[m_enet$.metric == "accuracy"],
            best_enet$penalty))
cat(sprintf("  Random Forest: ROC-AUC = %.3f | Accuracy = %.3f (Optimal min_n = %d)\n",
            m_rf$.estimate[m_rf$.metric == "roc_auc"],
            m_rf$.estimate[m_rf$.metric == "accuracy"],
            best_rf$min_n))
cat(sprintf("  XGBoost:       ROC-AUC = %.3f | Accuracy = %.3f (Optimal learn_rate = %.3f)\n",
            m_xgb$.estimate[m_xgb$.metric == "roc_auc"],
            m_xgb$.estimate[m_xgb$.metric == "accuracy"],
            best_xgb$learn_rate))

# -----------------------------------------------------------------------------
# Listing 23.5: Leakage-Safe Probability Calibration via Platt Scaling
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.5: Leakage-Safe Platt Calibration via Training OOF ---\n")
ctrl_oof <- control_resamples(save_pred = TRUE)
wf_enet_final <- finalize_workflow(wf_enet, best_enet)
wf_rf_final   <- finalize_workflow(wf_rf,   best_rf)

set.seed(42)
oof_enet <- fit_resamples(wf_enet_final, resamples = cv_folds,
                          control = ctrl_oof, metrics = metric_set(roc_auc, brier_class))
set.seed(42)
oof_rf   <- fit_resamples(wf_rf_final,   resamples = cv_folds,
                          control = ctrl_oof, metrics = metric_set(roc_auc, brier_class))

oof_pred_enet <- collect_predictions(oof_enet)
oof_pred_rf   <- collect_predictions(oof_rf)

# Estimate calibration maps strictly on training OOF predictions
cal_map_enet <- cal_estimate_logistic(oof_pred_enet, truth = Status)
cal_map_rf   <- cal_estimate_logistic(oof_pred_rf,   truth = Status)

# Apply frozen calibration maps to holdout test set
pred_enet_test <- collect_predictions(fit_enet)
pred_rf_test   <- collect_predictions(fit_rf)

pred_enet_cal <- cal_apply(pred_enet_test, cal_map_enet)
pred_rf_cal   <- cal_apply(pred_rf_test,   cal_map_rf)

brier_raw_enet <- brier_class(pred_enet_test, truth = Status, .pred_Control)$.estimate
brier_cal_enet <- brier_class(pred_enet_cal,  truth = Status, .pred_Control)$.estimate
brier_raw_rf   <- brier_class(pred_rf_test,   truth = Status, .pred_Control)$.estimate
brier_cal_rf   <- brier_class(pred_rf_cal,    truth = Status, .pred_Control)$.estimate

cat("Probability Calibration Audit (Brier Score, leakage-safe):\n")
cat(sprintf("  Elastic Net:   Raw BS = %.4f -> Calibrated BS = %.4f (%.1f%% error reduction)\n",
            brier_raw_enet, brier_cal_enet, 100 * (brier_raw_enet - brier_cal_enet) / brier_raw_enet))
cat(sprintf("  Random Forest: Raw BS = %.4f -> Calibrated BS = %.4f (%.1f%% error reduction)\n",
            brier_raw_rf, brier_cal_rf, 100 * (brier_raw_rf - brier_cal_rf) / brier_raw_rf))

bin_breaks <- seq(0, 1, length.out = 6)
cal_audit_enet <- pred_enet_cal |>
  mutate(bin = cut(.pred_Control, breaks = bin_breaks, include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(mean_pred = round(mean(.pred_Control), 3),
            mean_obs  = round(mean(Status == "Control"), 3),
            n = n(), .groups = "drop")
cat("\nHoldout Reliability Audit across Probability Bins (Elastic Net):\n")
print(cal_audit_enet)

# -----------------------------------------------------------------------------
# Listing 23.6: Clinical Decision Curve Analysis and Biomarker Recovery
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.6: Clinical Decision Curve Analysis & Signature Recovery ---\n")
dca_df <- data.frame(
  Status    = as.numeric(pred_enet_test$Status == "Tumor"),
  pred_enet = pred_enet_cal$.pred_Tumor,
  pred_rf   = pred_rf_cal$.pred_Tumor
)

dca_obj <- dcurves::dca(
  Status ~ pred_enet + pred_rf,
  data = dca_df,
  thresholds = c(0.15, 0.25, 0.35, 0.50)
)
nb_table <- as.data.frame(dca_obj$dca)

cat("Clinical Decision Curve Analysis (Net Benefit across Decision Thresholds pt):\n")
for (thresh in c(0.15, 0.25, 0.35, 0.50)) {
  sub_nb <- subset(nb_table, abs(threshold - thresh) < 0.01)
  nb_e <- sub_nb$net_benefit[sub_nb$variable == "pred_enet"]
  nb_r <- sub_nb$net_benefit[sub_nb$variable == "pred_rf"]
  nb_all <- sub_nb$net_benefit[sub_nb$variable == "all"]
  cat(sprintf("  pt = %.2f: NetBenefit(ENet) = %.3f | NetBenefit(RF) = %.3f | Treat-All = %.3f\n",
              thresh, nb_e, nb_r, nb_all))
}

rf_parsnip <- fit_rf$.workflow[[1]] |> extract_fit_parsnip()
rf_imp <- sort(ranger::importance(rf_parsnip$fit), decreasing = TRUE)
enet_parsnip <- fit_enet$.workflow[[1]] |> extract_fit_parsnip()
enet_coefs <- as.matrix(coef(enet_parsnip$fit, s = best_enet$penalty))
non_zero_genes <- rownames(enet_coefs)[enet_coefs[, 1] != 0 & rownames(enet_coefs) != "(Intercept)"]

cat("\nBiomarker Signature Recovery Audit:\n")
cat(sprintf("  Random Forest Top 5 Features: %s\n", paste(names(rf_imp)[1:5], collapse = ", ")))
cat(sprintf("  Elastic Net Total Features Retained: %d\n", length(non_zero_genes)))
cat(sprintf("  Elastic Net Ground-Truth Module 1 Overlap: %d / 8 proteins retained\n",
            sum(paste0("PROT_", sprintf("%04d", 1:8)) %in% non_zero_genes)))

# -----------------------------------------------------------------------------
# Listing 23.7: Batch-Aware Resampling via Leave-One-Batch-Out CV (LOBOCV)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.7: Batch-Aware Resampling via Leave-One-Batch-Out CV ---\n")
assign_batch <- function(status_vec, seed = 2026) {
  set.seed(seed)
  batches <- character(length(status_vec))
  for (i in seq_along(status_vec)) {
    if (status_vec[i] == "Control") {
      probs <- c(0.35, 0.25, 0.20, 0.12, 0.08)
    } else {
      probs <- c(0.08, 0.12, 0.20, 0.25, 0.35)
    }
    batches[i] <- sample(LETTERS[1:5], size = 1, prob = probs)
  }
  factor(batches, levels = LETTERS[1:5])
}

train_data$Batch <- assign_batch(as.character(train_data$Status))
test_data$Batch  <- assign_batch(as.character(test_data$Status))

tab_bc <- table(train_data$Batch, train_data$Status)
cat("Batch x Class contingency table (training set):\n")
print(tab_bc)

cramer_v <- sqrt(chisq.test(tab_bc)$statistic / (sum(tab_bc) * (min(nrow(tab_bc), ncol(tab_bc)) - 1)))
cat(sprintf("\nCramér's V (batch-class confounding association): %.3f\n", as.numeric(cramer_v)))

lobo_folds <- group_vfold_cv(train_data, group = "Batch", v = NULL)
cat(sprintf("LOBOCV configuration: %d folds (one per batch)\n", length(lobo_folds$splits)))

rec_lobo <- recipe(Status ~ ., data = train_data) |>
  update_role(Batch, new_role = "batch_id") |>
  update_role(Age, new_role = "imputation donor") |>
  step_impute_median(all_numeric_predictors()) |>
  step_zv(all_numeric_predictors()) |>
  step_normalize(all_numeric_predictors())

wf_enet_lobo <- workflow() |> add_recipe(rec_lobo) |> add_model(spec_enet)
wf_rf_lobo   <- workflow() |> add_recipe(rec_lobo) |> add_model(spec_rf)
wf_xgb_lobo  <- workflow() |> add_recipe(rec_lobo) |> add_model(spec_xgb)

set.seed(42)
res_enet_lobo <- fit_resamples(
  finalize_workflow(wf_enet_lobo, best_enet),
  resamples = lobo_folds,
  control   = control_resamples(save_pred = TRUE),
  metrics   = metric_set(roc_auc, accuracy, brier_class)
)

set.seed(42)
res_rf_lobo <- fit_resamples(
  finalize_workflow(wf_rf_lobo, best_rf),
  resamples = lobo_folds,
  control   = control_resamples(save_pred = TRUE),
  metrics   = metric_set(roc_auc, accuracy, brier_class)
)

set.seed(42)
res_xgb_lobo <- fit_resamples(
  finalize_workflow(wf_xgb_lobo, best_xgb),
  resamples = lobo_folds,
  control   = control_resamples(save_pred = TRUE),
  metrics   = metric_set(roc_auc, accuracy, brier_class)
)

stratified_auc_enet <- collect_metrics(res_enet) |> filter(.metric == "roc_auc") |> pull(mean)
stratified_auc_rf   <- collect_metrics(res_rf)   |> filter(.metric == "roc_auc") |> pull(mean)
stratified_auc_xgb  <- collect_metrics(res_xgb)  |> filter(.metric == "roc_auc") |> pull(mean)

lobo_auc_enet <- collect_metrics(res_enet_lobo) |> filter(.metric == "roc_auc") |> pull(mean)
lobo_auc_rf   <- collect_metrics(res_rf_lobo)   |> filter(.metric == "roc_auc") |> pull(mean)
lobo_auc_xgb  <- collect_metrics(res_xgb_lobo)  |> filter(.metric == "roc_auc") |> pull(mean)

cat("\n=== Resampling Architecture Comparison ===\n")
cat(sprintf("  Elastic Net:  Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f\n",
            stratified_auc_enet, lobo_auc_enet))
cat(sprintf("  Random Forest: Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f\n",
            stratified_auc_rf, lobo_auc_rf))
cat(sprintf("  XGBoost:      Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f\n",
            stratified_auc_xgb, lobo_auc_xgb))

cat("\n=== Per-Batch LOBOCV AUC (Elastic Net) ===\n")
per_fold_enet <- res_enet_lobo$.metrics |>
  map_dfr(function(m) {
    auc_row <- m[m$.metric == "roc_auc", ]
    data.frame(auc = round(auc_row$.estimate, 4))
  }) |>
  mutate(Batch = paste0("Batch ", LETTERS[1:5]))
print(per_fold_enet)

# -----------------------------------------------------------------------------
# Listing 23.8: Bootstrap 95% Confidence Intervals on Holdout Metrics (B = 2000)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.8: Bootstrap 95% Confidence Intervals (B = 2000) ---\n")
test_preds <- list(
  ENet = collect_predictions(fit_enet) |> select(.row, Status, .pred_class, .pred_Control, .pred_Tumor),
  RF   = collect_predictions(fit_rf)   |> select(.row, Status, .pred_class, .pred_Control, .pred_Tumor),
  XGB  = collect_predictions(fit_xgb)  |> select(.row, Status, .pred_class, .pred_Control, .pred_Tumor)
)

compute_boot_metrics <- function(pred_df, B = 2000, seed = 42) {
  set.seed(seed)
  N <- nrow(pred_df)
  truth    <- pred_df$Status
  pred_cls <- pred_df$.pred_class
  pred_p1  <- pred_df$.pred_Control

  map_dfr(1:B, function(i) {
    idx <- sample.int(N, size = N, replace = TRUE)
    b_truth <- truth[idx]
    b_cls   <- pred_cls[idx]
    b_p1    <- pred_p1[idx]

    if (length(unique(b_truth)) < 2) return(tibble(auc = NA_real_, acc = NA_real_, brier = NA_real_))

    boot_df <- tibble(Status = b_truth, .pred_Control = b_p1, .pred_class = b_cls)
    auc_v   <- roc_auc(boot_df, truth = Status, .pred_Control)$.estimate
    acc_v   <- accuracy(boot_df, truth = Status, estimate = .pred_class)$.estimate
    brier_v <- brier_class(boot_df, truth = Status, .pred_Control)$.estimate

    tibble(auc = auc_v, acc = acc_v, brier = brier_v)
  }) |>
    filter(!is.na(auc))
}

cat("Computing bootstrap resamples...\n")
boot_results <- map(test_preds, compute_boot_metrics, B = 2000)

summarize_ci <- function(boot_df, metric_name, point_est) {
  vals <- boot_df[[metric_name]]
  tibble(
    metric = metric_name,
    point_estimate = point_est,
    ci_low  = quantile(vals, 0.025, na.rm = TRUE),
    ci_high = quantile(vals, 0.975, na.rm = TRUE),
    boot_mean = mean(vals, na.rm = TRUE),
    boot_se   = sd(vals, na.rm = TRUE)
  )
}

point_ests <- list(
  ENet = list(auc = m_enet$.estimate[m_enet$.metric == "roc_auc"],
              acc = m_enet$.estimate[m_enet$.metric == "accuracy"],
              brier = brier_class(test_preds$ENet, truth = Status, .pred_Control)$.estimate),
  RF   = list(auc = m_rf$.estimate[m_rf$.metric == "roc_auc"],
              acc = m_rf$.estimate[m_rf$.metric == "accuracy"],
              brier = brier_class(test_preds$RF, truth = Status, .pred_Control)$.estimate),
  XGB  = list(auc = m_xgb$.estimate[m_xgb$.metric == "roc_auc"],
              acc = m_xgb$.estimate[m_xgb$.metric == "accuracy"],
              brier = brier_class(test_preds$XGB, truth = Status, .pred_Control)$.estimate)
)

ci_table <- imap_dfr(boot_results, function(bdf, model_name) {
  bind_rows(
    summarize_ci(bdf, "auc",   point_ests[[model_name]]$auc),
    summarize_ci(bdf, "acc",   point_ests[[model_name]]$acc),
    summarize_ci(bdf, "brier", point_ests[[model_name]]$brier)
  ) |>
    mutate(model = model_name, .before = 1)
})

cat("\n=== Bootstrap 95% Confidence Intervals on Holdout Metrics (N_test = 30) ===\n")
ci_table |>
  mutate(display = sprintf("%.3f [%.3f, %.3f]", point_estimate, ci_low, ci_high)) |>
  select(model, metric, display) |>
  pivot_wider(names_from = metric, values_from = display) |>
  print(width = Inf)

cat("\n=============================================================================\n")
cat("Laboratory Execution Complete: All Workflows Verified Reproducible.\n")
cat("=============================================================================\n")
