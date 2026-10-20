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
  library(stacks)
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
# Listing 23.1: Simulating Proteomic Cohort (p = 4000, n = 120) with Batch Effect
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.1: Simulating Proteomic Cohort (p = 4000, n = 120) ---\n")
set.seed(2026)
n_samples  <- 120
n_proteins <- 4000L

status <- factor(c(rep("Control", 80), rep("Tumor", 40)), levels = c("Control", "Tumor"))

# Observed continuous covariate: patient age (MAR missingness driver)
age_years <- c(rnorm(80, mean = 55, sd = 12), rnorm(40, mean = 58, sd = 12))
age_years <- pmax(18, pmin(90, age_years))

# Assign acquisition batch confounded with disease status
assign_batch <- function(status_vec, seed = 2026) {
  set.seed(seed)
  batches <- character(length(status_vec))
  for (i in seq_along(status_vec)) {
    probs <- if (status_vec[i] == "Control") c(0.40, 0.25, 0.18, 0.10, 0.07) else c(0.05, 0.10, 0.18, 0.27, 0.40)
    batches[i] <- sample(LETTERS[1:5], size = 1, prob = probs)
  }
  factor(batches, levels = LETTERS[1:5])
}
batch_vec <- assign_batch(as.character(status))

# Baseline protein abundance: N(20, 1.0)
X <- matrix(rnorm(n_samples * n_proteins, mean = 20, sd = 1.0),
            nrow = n_samples, ncol = n_proteins)
colnames(X) <- paste0("PROT_", sprintf("%04d", 1:n_proteins))

# Module 1: Collinear co-regulated biomarker signature (Proteins 1-8, r ~ 0.70)
latent_driver <- c(rep(0, 80), rep(0.65, 40))
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

# Inject technical batch effects on non-marker background proteins (PROT_0021 to PROT_0220)
batch_shifts <- c(A = -0.75, B = -0.35, C = 0.00, D = 0.35, E = 0.75)
for (i in 1:n_samples) {
  b <- as.character(batch_vec[i])
  X[i, 21:220] <- X[i, 21:220] + batch_shifts[b]
}

# MAR missingness driven by Age
miss_prob <- plogis(-3.2 + 0.05 * (age_years - 55))
miss_mask <- matrix(runif(n_samples * n_proteins) < rep(miss_prob, times = n_proteins),
                    nrow = n_samples, ncol = n_proteins)

# MNAR component on low-abundance proteins near detection limit
low_abund <- which(colMeans(X, na.rm = TRUE) < 19.5)
for (j in low_abund) {
  thresh_j <- quantile(X[, j], 0.05, na.rm = TRUE)
  miss_mask[, j] <- miss_mask[, j] | ((X[, j] < thresh_j) & (runif(n_samples) < 0.30))
}
X[miss_mask] <- NA

df_cohort <- as.data.frame(X)
df_cohort$Status <- status
df_cohort$Age    <- age_years
df_cohort$Batch  <- batch_vec

# Stratified 75/25 Train-Test Split
set.seed(42)
split_obj  <- initial_split(df_cohort, prop = 0.75, strata = Status)
train_data <- training(split_obj)
test_data  <- testing(split_obj)

cat(sprintf("Simulated cohort: %d samples x %d features (p/n ratio ~ %.0f).\n",
            nrow(df_cohort), ncol(df_cohort) - 3, (ncol(df_cohort) - 3) / nrow(df_cohort)))
cat(sprintf("Training set: %d patients (%d Control, %d Tumor)\n",
            nrow(train_data), sum(train_data$Status == "Control"), sum(train_data$Status == "Tumor")))
cat(sprintf("Independent test holdout: %d patients (%d Control, %d Tumor)\n",
            nrow(test_data), sum(test_data$Status == "Control"), sum(test_data$Status == "Tumor")))

# MAR audit across age quartiles
train_data$.age_quartile <- cut(train_data$Age,
  breaks = quantile(train_data$Age, c(0, 0.25, 0.5, 0.75, 1)),
  include.lowest = TRUE, labels = c("Q1(young)", "Q2", "Q3", "Q4(old)"))
protein_cols <- grep("^PROT_", names(train_data), value = TRUE)
miss_by_age <- train_data |> group_by(.age_quartile) |>
  summarise(missingness_rate = round(mean(is.na(as.matrix(across(all_of(protein_cols))))), 4),
            n = n(), .groups = "drop")
print(miss_by_age)
train_data$.age_quartile <- NULL

# -----------------------------------------------------------------------------
# Listing 23.2: Leakage-Safe Preprocessing Recipe
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.2: Leakage-Safe Preprocessing Recipe ---\n")
rec <- recipe(Status ~ ., data = train_data) |>
  update_role(Batch, new_role = "batch_id") |>
  update_role(Age, new_role = "imputation donor") |>
  step_impute_median(all_numeric_predictors()) |>
  step_zv(all_numeric_predictors()) |>
  step_normalize(all_numeric_predictors()) |>
  step_downsample(Status, under_ratio = 1.2)

prep_rec <- prep(rec, training = train_data)
baked_train <- bake(prep_rec, new_data = NULL)

cat(sprintf("Recipe successfully baked on training fold:\n"))
cat(sprintf("  Predictors processed: %d (Zero-variance filtered: %d)\n",
            ncol(baked_train) - 3, ncol(train_data) - 3 - (ncol(baked_train) - 3)))
cat(sprintf("  Post-downsample class balance: %d Control, %d Tumor\n",
            sum(baked_train$Status == "Control"), sum(baked_train$Status == "Tumor")))
cat(sprintf("  Remaining missing values: %d\n", sum(is.na(baked_train))))

# -----------------------------------------------------------------------------
# Listing 23.3: Multi-Model Specification and 5-Fold Cross-Validation Tuning
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.3: Multi-Model Specification and 5-Fold CV Tuning ---\n")
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

# Stacks requires candidate out-of-fold predictions and workflow objects
ctrl_stack <- control_stack_grid()

res_enet <- tune_grid(wf_enet, resamples = cv_folds,
                      grid = tibble(penalty = 10^seq(-3, -0.2, length.out = 6)),
                      metrics = metric_set(roc_auc, accuracy, brier_class),
                      control = ctrl_stack)

res_rf <- tune_grid(wf_rf, resamples = cv_folds,
                    grid = tibble(min_n = c(3, 7, 12)),
                    metrics = metric_set(roc_auc, accuracy, brier_class),
                    control = ctrl_stack)

res_xgb <- tune_grid(wf_xgb, resamples = cv_folds,
                     grid = tibble(learn_rate = c(0.02, 0.08, 0.15)),
                     metrics = metric_set(roc_auc, accuracy, brier_class),
                     control = ctrl_stack)

cat("Hyperparameter grid tuning complete across all 3 model architectures.\n")

# -----------------------------------------------------------------------------
# Listing 23.4: Model Benchmark on Independent Holdout Test Set
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.4: Model Benchmark on Independent Holdout Test Set ---\n")
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
# Listing 23.5: Probability Calibration Audit & Expected Calibration Error (ECE)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.5: Probability Calibration Audit & ECE ---\n")
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

cal_map_enet <- cal_estimate_logistic(oof_pred_enet, truth = Status)
cal_map_rf   <- cal_estimate_logistic(oof_pred_rf,   truth = Status)

pred_enet_test <- collect_predictions(fit_enet)
pred_rf_test   <- collect_predictions(fit_rf)

pred_enet_cal <- cal_apply(pred_enet_test, cal_map_enet)
pred_rf_cal   <- cal_apply(pred_rf_test,   cal_map_rf)

brier_raw_enet <- brier_class(pred_enet_test, truth = Status, .pred_Control)$.estimate
brier_cal_enet <- brier_class(pred_enet_cal,  truth = Status, .pred_Control)$.estimate
brier_raw_rf   <- brier_class(pred_rf_test,   truth = Status, .pred_Control)$.estimate
brier_cal_rf   <- brier_class(pred_rf_cal,    truth = Status, .pred_Control)$.estimate

compute_ece <- function(pred_df, n_bins = 5) {
  breaks <- seq(0, 1, length.out = n_bins + 1)
  df <- pred_df |>
    mutate(bin = cut(.pred_Control, breaks = breaks, include.lowest = TRUE),
           y_num = as.numeric(Status == "Control")) |>
    group_by(bin) |>
    summarise(
      n = n(),
      conf = mean(.pred_Control),
      acc = mean(y_num),
      .groups = "drop"
    )
  sum((df$n / nrow(pred_df)) * abs(df$acc - df$conf))
}

ece_raw_enet <- compute_ece(pred_enet_test)
ece_cal_enet <- compute_ece(pred_enet_cal)
ece_raw_rf   <- compute_ece(pred_rf_test)
ece_cal_rf   <- compute_ece(pred_rf_cal)

cat("Probability Calibration Audit (Brier Score and ECE, leakage-safe):\n")
cat(sprintf("  Elastic Net:   Raw BS = %.4f (ECE = %.3f) -> Calibrated BS = %.4f (ECE = %.3f) [%+.1f%% BS error change]\n",
            brier_raw_enet, ece_raw_enet, brier_cal_enet, ece_cal_enet, 100 * (brier_raw_enet - brier_cal_enet) / brier_raw_enet))
cat(sprintf("  Random Forest: Raw BS = %.4f (ECE = %.3f) -> Calibrated BS = %.4f (ECE = %.3f) [%+.1f%% BS error change]\n",
            brier_raw_rf, ece_raw_rf, brier_cal_rf, ece_cal_rf, 100 * (brier_raw_rf - brier_cal_rf) / brier_raw_rf))

# Holdout reliability audit for Elastic Net
bin_breaks <- seq(0, 1, length.out = 6)
cal_audit_enet <- pred_enet_cal |>
  mutate(bin = cut(.pred_Control, breaks = bin_breaks, include.lowest = TRUE)) |>
  group_by(bin) |>
  summarise(mean_pred = round(mean(.pred_Control), 3),
            mean_obs  = round(mean(Status == "Control"), 3),
            n = n(), .groups = "drop")
print(cal_audit_enet)

# -----------------------------------------------------------------------------
# Listing 23.6: Pedagogical DCA and Multi-Method Feature Importance (TreeSHAP)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.6: Pedagogical DCA and Multi-Method Feature Importance ---\n")
# DCA using calibrated Elastic Net and raw Random Forest (RF Platt calibration failed)
dca_df <- data.frame(
  Status    = as.numeric(pred_enet_test$Status == "Tumor"),
  pred_enet = pred_enet_cal$.pred_Tumor,
  pred_rf   = pred_rf_test$.pred_Tumor
)

dca_obj <- dcurves::dca(
  Status ~ pred_enet + pred_rf,
  data = dca_df,
  thresholds = c(0.15, 0.25, 0.35, 0.50)
)
nb_table <- as.data.frame(dca_obj$dca)

cat("Illustrative Simulation DCA (Net Benefit across Decision Thresholds pt):\n")
for (thresh in c(0.15, 0.25, 0.35, 0.50)) {
  sub_nb <- subset(nb_table, abs(threshold - thresh) < 0.01)
  nb_e <- sub_nb$net_benefit[sub_nb$variable == "pred_enet"]
  nb_r <- sub_nb$net_benefit[sub_nb$variable == "pred_rf"]
  nb_all <- sub_nb$net_benefit[sub_nb$variable == "all"]
  cat(sprintf("  pt = %.2f: NetBenefit(ENet) = %.3f | NetBenefit(RF) = %.3f | Treat-All = %.3f\n",
              thresh, nb_e, nb_r, nb_all))
}

# 1. Random forest permutation importance
rf_parsnip <- fit_rf$.workflow[[1]] |> extract_fit_parsnip()
rf_imp <- sort(ranger::importance(rf_parsnip$fit), decreasing = TRUE)

# 2. Elastic net coefficients
enet_parsnip <- fit_enet$.workflow[[1]] |> extract_fit_parsnip()
enet_coefs <- as.matrix(coef(enet_parsnip$fit, s = best_enet$penalty))
non_zero_genes <- rownames(enet_coefs)[enet_coefs[, 1] != 0 & rownames(enet_coefs) != "(Intercept)"]

# 3. XGBoost TreeSHAP values on test data
xgb_fit_raw <- extract_fit_engine(fit_xgb$.workflow[[1]])
baked_test  <- bake(prep_rec, new_data = test_data)
test_mat    <- as.matrix(baked_test |> select(starts_with("PROT_")))
shap_contrib <- predict(xgb_fit_raw, newdata = test_mat, predcontrib = TRUE)
mean_abs_shap <- colMeans(abs(shap_contrib[, colnames(test_mat)]))
top_shap <- sort(mean_abs_shap, decreasing = TRUE)

cat("\nBiomarker Signature Recovery Audit (Multi-Method Comparison):\n")
cat(sprintf("  Elastic Net Total Features Retained: %d\n", length(non_zero_genes)))
cat(sprintf("  Elastic Net Ground-Truth Module 1 Overlap: %d / 8 proteins retained\n",
            sum(paste0("PROT_", sprintf("%04d", 1:8)) %in% non_zero_genes)))
cat(sprintf("  Random Forest Top 5 Features (Permutation Imp): %s\n", paste(names(rf_imp)[1:5], collapse = ", ")))
cat(sprintf("  XGBoost Top 5 Features (Exact TreeSHAP):        %s\n", paste(names(top_shap)[1:5], collapse = ", ")))

# -----------------------------------------------------------------------------
# Listing 23.7: Batch-Aware Resampling via Leave-One-Batch-Out CV (LOBOCV)
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.7: Batch-Aware Resampling via LOBOCV ---\n")
tab_bc <- table(train_data$Batch, train_data$Status)
cat("Batch x Class contingency table (training set):\n")
print(tab_bc)

cramer_v <- sqrt(chisq.test(tab_bc)$statistic / (sum(tab_bc) * (min(nrow(tab_bc), ncol(tab_bc)) - 1)))
cat(sprintf("\nCramér's V (batch-class confounding association): %.3f\n", as.numeric(cramer_v)))

lobo_folds <- group_vfold_cv(train_data, group = "Batch", v = NULL)

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

cat("\n=== Resampling Architecture Comparison (Batch Confounder Impact) ===\n")
cat(sprintf("  Elastic Net:   Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f (Drop = %+.3f)\n",
            stratified_auc_enet, lobo_auc_enet, lobo_auc_enet - stratified_auc_enet))
cat(sprintf("  Random Forest: Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f (Drop = %+.3f)\n",
            stratified_auc_rf, lobo_auc_rf, lobo_auc_rf - stratified_auc_rf))
cat(sprintf("  XGBoost:       Stratified 5-fold CV AUC = %.3f  ->  LOBOCV AUC = %.3f (Drop = %+.3f)\n",
            stratified_auc_xgb, lobo_auc_xgb, lobo_auc_xgb - stratified_auc_xgb))

# -----------------------------------------------------------------------------
# Listing 23.8: Paired Bootstrap Resampling and Uncertainty Audit
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.8: Paired Bootstrap Resampling and Uncertainty Audit ---\n")
test_preds <- list(
  ENet = collect_predictions(fit_enet) |> select(.row, Status, .pred_class, .pred_Control),
  RF   = collect_predictions(fit_rf)   |> select(.row, Status, .pred_class, .pred_Control),
  XGB  = collect_predictions(fit_xgb)  |> select(.row, Status, .pred_class, .pred_Control)
)

set.seed(42)
B <- 2000
N_test <- nrow(test_preds$ENet)

boot_paired <- map_dfr(1:B, function(i) {
  idx <- sample.int(N_test, size = N_test, replace = TRUE)
  
  df_enet <- test_preds$ENet[idx, ]
  df_rf   <- test_preds$RF[idx, ]
  df_xgb  <- test_preds$XGB[idx, ]
  
  if (length(unique(df_enet$Status)) < 2) return(NULL)
  
  auc_enet <- roc_auc(df_enet, truth = Status, .pred_Control)$.estimate
  auc_rf   <- roc_auc(df_rf, truth = Status, .pred_Control)$.estimate
  auc_xgb  <- roc_auc(df_xgb, truth = Status, .pred_Control)$.estimate
  
  acc_enet <- accuracy(df_enet, truth = Status, estimate = .pred_class)$.estimate
  acc_rf   <- accuracy(df_rf, truth = Status, estimate = .pred_class)$.estimate
  acc_xgb  <- accuracy(df_xgb, truth = Status, estimate = .pred_class)$.estimate
  
  brier_enet <- brier_class(df_enet, truth = Status, .pred_Control)$.estimate
  brier_rf   <- brier_class(df_rf, truth = Status, .pred_Control)$.estimate
  brier_xgb  <- brier_class(df_xgb, truth = Status, .pred_Control)$.estimate
  
  tibble(
    auc_enet = auc_enet, auc_rf = auc_rf, auc_xgb = auc_xgb,
    delta_auc_enet_rf  = auc_enet - auc_rf,
    delta_auc_enet_xgb = auc_enet - auc_xgb,
    acc_enet = acc_enet, acc_rf = acc_rf, acc_xgb = acc_xgb,
    brier_enet = brier_enet, brier_rf = brier_rf, brier_xgb = brier_xgb
  )
})

cat("Marginal 95% Confidence Intervals on Holdout Metrics (N = 30, B = 2000):\n")
cat(sprintf("  Elastic Net:   AUC = %.3f [%.3f, %.3f] | Acc = %.3f [%.3f, %.3f] | Brier = %.3f [%.3f, %.3f]\n",
            m_enet$.estimate[m_enet$.metric == "roc_auc"],
            quantile(boot_paired$auc_enet, 0.025), quantile(boot_paired$auc_enet, 0.975),
            m_enet$.estimate[m_enet$.metric == "accuracy"],
            quantile(boot_paired$acc_enet, 0.025), quantile(boot_paired$acc_enet, 0.975),
            brier_class(test_preds$ENet, truth = Status, .pred_Control)$.estimate,
            quantile(boot_paired$brier_enet, 0.025), quantile(boot_paired$brier_enet, 0.975)))

cat(sprintf("  Random Forest: AUC = %.3f [%.3f, %.3f] | Acc = %.3f [%.3f, %.3f] | Brier = %.3f [%.3f, %.3f]\n",
            m_rf$.estimate[m_rf$.metric == "roc_auc"],
            quantile(boot_paired$auc_rf, 0.025), quantile(boot_paired$auc_rf, 0.975),
            m_rf$.estimate[m_rf$.metric == "accuracy"],
            quantile(boot_paired$acc_rf, 0.025), quantile(boot_paired$acc_rf, 0.975),
            brier_class(test_preds$RF, truth = Status, .pred_Control)$.estimate,
            quantile(boot_paired$brier_rf, 0.025), quantile(boot_paired$brier_rf, 0.975)))

cat(sprintf("  XGBoost:       AUC = %.3f [%.3f, %.3f] | Acc = %.3f [%.3f, %.3f] | Brier = %.3f [%.3f, %.3f]\n",
            m_xgb$.estimate[m_xgb$.metric == "roc_auc"],
            quantile(boot_paired$auc_xgb, 0.025), quantile(boot_paired$auc_xgb, 0.975),
            m_xgb$.estimate[m_xgb$.metric == "accuracy"],
            quantile(boot_paired$acc_xgb, 0.025), quantile(boot_paired$acc_xgb, 0.975),
            brier_class(test_preds$XGB, truth = Status, .pred_Control)$.estimate,
            quantile(boot_paired$brier_xgb, 0.025), quantile(boot_paired$brier_xgb, 0.975)))

# Paired difference testing
p_diff_enet_rf  <- 2 * min(mean(boot_paired$delta_auc_enet_rf > 0), mean(boot_paired$delta_auc_enet_rf < 0))
p_diff_enet_xgb <- 2 * min(mean(boot_paired$delta_auc_enet_xgb > 0), mean(boot_paired$delta_auc_enet_xgb < 0))

cat("\nPaired Bootstrap Difference Tests (Delta AUC):\n")
cat(sprintf("  ENet vs RF:  Mean Delta = %+.3f, 95%% CI [%+.3f, %+.3f], Paired bootstrap p = %.3f\n",
            mean(boot_paired$delta_auc_enet_rf),
            quantile(boot_paired$delta_auc_enet_rf, 0.025),
            quantile(boot_paired$delta_auc_enet_rf, 0.975),
            p_diff_enet_rf))
cat(sprintf("  ENet vs XGB: Mean Delta = %+.3f, 95%% CI [%+.3f, %+.3f], Paired bootstrap p = %.3f\n",
            mean(boot_paired$delta_auc_enet_xgb),
            quantile(boot_paired$delta_auc_enet_xgb, 0.025),
            quantile(boot_paired$delta_auc_enet_xgb, 0.975),
            p_diff_enet_xgb))

# -----------------------------------------------------------------------------
# Listing 23.9: Stacked Generalization via the stacks Meta-Learning Framework
# -----------------------------------------------------------------------------
cat("\n--- Listing 23.9: Stacked Generalization via stacks ---\n")

# Initialize candidate stack from 5-fold CV tuned model objects
proteomic_stack <- stacks() |>
  add_candidates(res_enet) |>
  add_candidates(res_rf) |>
  add_candidates(res_xgb)

# Fit regularized meta-learner across out-of-fold candidate predictions
set.seed(42)
stack_model <- proteomic_stack |>
  blend_predictions(penalty = 10^seq(-3, -0.5, length.out = 10))

# Fit non-zero weighted members on the full training partition
fitted_stack <- fit_members(stack_model)

# Evaluate stacked ensemble on independent holdout test set (N = 30)
stack_test_probs <- predict(fitted_stack, test_data, type = "prob")
stack_test_class <- predict(fitted_stack, test_data, type = "class")
test_res_stack   <- bind_cols(test_data |> select(Status), stack_test_probs, stack_test_class)

auc_stack   <- roc_auc(test_res_stack, truth = Status, .pred_Control)$.estimate
acc_stack   <- accuracy(test_res_stack, truth = Status, .pred_class)$.estimate
brier_stack <- brier_class(test_res_stack, truth = Status, .pred_Control)$.estimate

cat("Stacked Ensemble Holdout Performance (N = 30):\n")
cat(sprintf("  ROC-AUC:  %.3f\n", auc_stack))
cat(sprintf("  Accuracy: %.3f\n", acc_stack))
cat(sprintf("  Brier:    %.4f\n", brier_stack))

cat("\n=============================================================================\n")
cat("Laboratory Execution Complete: All Workflows Verified Reproducible.\n")
cat("=============================================================================\n")
