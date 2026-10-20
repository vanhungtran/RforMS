# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Monograph Series)
# Chapter 19: Publication-Quality Figures via figures4papers House Style
#
# Style Guide: C:/Users/lucas/OneDrive/figures4papers-main/scientific-figure-making-r
# Journal Preset: "nature" (183 mm double-column, 7 pt base, 5 pt text floor)
# Author: Lucas VHH TRAN (CK-CARE Davos)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(scales)
})

# Source figures4papers theme & helpers
theme_pub_path <- "C:/Users/lucas/OneDrive/figures4papers-main/scientific-figure-making-r/R/theme_pub.R"
source(theme_pub_path)

# Source companion laboratory script to obtain exact data objects
source("code/chapter19_missing_data_labs.R")

J <- "nature"
W <- pub_width(J, "double") # 183 mm

cat("\n=============================================================================\n")
cat("Generating Chapter 19 Figures in figures4papers House Style (Nature Preset)\n")
cat("=============================================================================\n\n")

# -----------------------------------------------------------------------------
# Figure 19.1: Empirical Missingness Diagnosis and Completeness Audit
# -----------------------------------------------------------------------------
cat("--- Generating Figure 19.1: Missingness Diagnosis ---\n")

mean_obs_vec <- rowMeans(obs_filtered, na.rm = TRUE)
miss_rate_vec <- rowSums(is.na(obs_filtered)) / n

df_diag <- tibble(
  Mean_Abundance = mean_obs_vec,
  Missing_Rate   = miss_rate_vec,
  Is_Low         = mean_obs_vec < 20.8
)

grid_x <- seq(min(mean_obs_vec), max(mean_obs_vec), length.out = 200)
pred_y <- plogis(coef_log[1, 1] + coef_log[2, 1] * grid_x)
df_pred <- tibble(Mean_Abundance = grid_x, Missing_Rate = pred_y)

p1a <- ggplot() +
  geom_point(data = df_diag, aes(x = Mean_Abundance, y = Missing_Rate, color = Is_Low),
             alpha = 0.55, size = 1.3) +
  geom_line(data = df_pred, aes(x = Mean_Abundance, y = Missing_Rate),
            color = PALETTE[["red_strong"]], linewidth = 0.8) +
  geom_vline(xintercept = 20.2, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  annotate("text", x = 20.4, y = 0.80, label = "LOD = 20.2",
           hjust = 0, size = 6.0 / .pt, color = "gray20") +
  annotate("text", x = 22.5, y = 0.55,
           label = "Slope = -0.358 (p = 1.8e-07)",
           size = 6.0 / .pt, color = PALETTE[["red_strong"]], hjust = 0) +
  scale_color_manual(values = c("FALSE" = PALETTE[["blue_main"]], "TRUE" = PALETTE[["red_strong"]])) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.85)) +
  labs(x = expression("Mean Observed Abundance ("*log[2]*" Intensity)"),
       y = "Missingness Fraction (%)") +
  theme_pub(J) +
  theme(legend.position = "none")

sample_completeness <- colMeans(!is.na(obs_filtered)) * 100
df_samples <- tibble(
  Sample = factor(colnames(obs_filtered), levels = colnames(obs_filtered)),
  Group  = group,
  Completeness = sample_completeness
)

p1b <- ggplot(df_samples, aes(x = Sample, y = Completeness, fill = Group)) +
  geom_col(width = 0.7) +
  geom_hline(yintercept = 70, linetype = "dashed", color = PALETTE[["red_strong"]], linewidth = 0.5) +
  annotate("text", x = 1.5, y = 85, label = "70% filter threshold",
           size = 6.0 / .pt, color = PALETTE[["red_strong"]], hjust = 0) +
  scale_fill_manual(values = c("Control" = PALETTE[["blue_main"]], "Disease" = PALETTE[["blue_secondary"]])) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 25, 50, 70, 100)) +
  labs(x = "Mass Spectrometry Run",
       y = "Quantified Features (%)") +
  theme_pub(J) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.5),
    legend.position = "top",
    legend.title = element_blank()
  )

fig1 <- (p1a + p1b) + pub_tags(J)

save_pub(fig1, "latex/figures/fig_19_1_missingness_diagnosis",
         width = W, height = 75, units = "mm", journal = J)
save_pub(fig1, "docs/latex/figures/fig_19_1_missingness_diagnosis",
         width = W, height = 75, units = "mm", journal = J)

# -----------------------------------------------------------------------------
# Figure 19.2: Post-Imputation Empirical Density Distributions
# -----------------------------------------------------------------------------
cat("\n--- Generating Figure 19.2: Imputation Distributions ---\n")

df_dens <- bind_rows(
  tibble(Intensity = as.vector(truth_filtered), Method = "Ground Truth (Oracle)"),
  tibble(Intensity = as.vector(obs_filtered[!is.na(obs_filtered)]), Method = "Observed (Truncated)"),
  tibble(Intensity = as.vector(mat_mindet), Method = "MinDet (Deterministic)"),
  tibble(Intensity = as.vector(mat_perseus), Method = "Perseus (Downshifted Gaussian)"),
  tibble(Intensity = as.vector(mat_qrilc), Method = "QRILC (Quantile Regression)"),
  tibble(Intensity = as.vector(mat_knn), Method = "k-NN (Donor Matching)")
)

method_levels <- c(
  "Ground Truth (Oracle)", "Observed (Truncated)", "MinDet (Deterministic)",
  "Perseus (Downshifted Gaussian)", "QRILC (Quantile Regression)", "k-NN (Donor Matching)"
)
df_dens$Method <- factor(df_dens$Method, levels = method_levels)

method_colors <- c(
  "Ground Truth (Oracle)"          = PALETTE[["blue_main"]],
  "Observed (Truncated)"           = "gray40",
  "MinDet (Deterministic)"         = PALETTE[["red_strong"]],
  "Perseus (Downshifted Gaussian)" = PALETTE[["blue_secondary"]],
  "QRILC (Quantile Regression)"    = PALETTE[["teal"]],
  "k-NN (Donor Matching)"          = PALETTE[["violet"]]
)

method_lines <- c(
  "Ground Truth (Oracle)"          = "solid",
  "Observed (Truncated)"           = "dashed",
  "MinDet (Deterministic)"         = "solid",
  "Perseus (Downshifted Gaussian)" = "solid",
  "QRILC (Quantile Regression)"    = "solid",
  "k-NN (Donor Matching)"          = "solid"
)

fig2 <- ggplot(df_dens, aes(x = Intensity, color = Method, linetype = Method)) +
  geom_density(linewidth = 0.7, trim = FALSE) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = method_lines) +
  scale_x_continuous(limits = c(15, 28), breaks = seq(16, 28, 2)) +
  annotate("segment", x = 17.52, xend = 17.52, y = 0.28, yend = 0.18,
           arrow = arrow(length = unit(1.5, "mm")), color = PALETTE[["red_strong"]], linewidth = 0.5) +
  annotate("text", x = 17.52, y = 0.30, label = "MinDet point-mass spike (17.52)",
           size = 5.5 / .pt, color = PALETTE[["red_strong"]], hjust = 0.5) +
  labs(x = expression("Protein Abundance ("*log[2]*" Intensity)"),
       y = "Empirical Density",
       color = "", linetype = "") +
  theme_pub(J) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.5),
    legend.key.width = unit(6, "mm")
  )

save_pub(fig2, "latex/figures/fig_19_2_imputation_distributions",
         width = W, height = 80, units = "mm", journal = J)
save_pub(fig2, "docs/latex/figures/fig_19_2_imputation_distributions",
         width = W, height = 80, units = "mm", journal = J)

# -----------------------------------------------------------------------------
# Figure 19.3: Differential Sensitivity and Variance Inflation Audit
# -----------------------------------------------------------------------------
cat("\n--- Generating Figure 19.3: Differential Sensitivity & Variance Inflation ---\n")

make_volcano_df <- function(logfc, pval, label_method, prot_names = rownames(obs_filtered)) {
  tibble(
    Protein = prot_names,
    logFC   = logfc,
    p_val   = pval,
    negLogP = -log10(pval),
    padj    = p.adjust(pval, method = "BH"),
    Status  = case_when(
      prot_names %in% true_up_names   ~ "Upregulated (+2.0)",
      prot_names %in% true_down_names ~ "Downregulated (-2.5)",
      TRUE                            ~ "True Null"
    ),
    Method  = label_method
  )
}

df_volcano <- bind_rows(
  make_volcano_df(res_oracle$logFC,  res_oracle$pval,  "Oracle (Ground Truth)"),
  make_volcano_df(res_mindet$logFC,  res_mindet$pval,  "MinDet (Deterministic)"),
  make_volcano_df(res_qrilc$logFC,   res_qrilc$pval,   "QRILC (Stochastic Left-Censored)"),
  make_volcano_df(res_proda$diff,    res_proda$pval,   "proDA (Dropout Likelihood)", prot_names = res_proda$name)
)

df_volcano$Method <- factor(df_volcano$Method, levels = c(
  "Oracle (Ground Truth)", "MinDet (Deterministic)",
  "QRILC (Stochastic Left-Censored)", "proDA (Dropout Likelihood)"
))

p3a <- ggplot(df_volcano, aes(x = logFC, y = negLogP, color = Status)) +
  geom_point(alpha = 0.7, size = 1.0) +
  facet_wrap(~ Method, ncol = 2) +
  scale_color_manual(values = c(
    "Upregulated (+2.0)"   = PALETTE[["red_strong"]],
    "Downregulated (-2.5)" = PALETTE[["blue_main"]],
    "True Null"            = PALETTE[["neutral"]]
  )) +
  scale_x_continuous(limits = c(-5.5, 5.5), breaks = seq(-4, 4, 2)) +
  scale_y_continuous(limits = c(0, 8.5), breaks = seq(0, 8, 2)) +
  labs(x = expression(log[2]*" Fold Change (Disease vs. Control)"),
       y = expression(-log[10]*"("*italic(p)*"-value)"),
       color = "") +
  theme_pub(J) +
  theme(
    strip.text = element_text(size = 6, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 5.5)
  )

df_rubin_audit <- tibble(
  Missing_Pct = (rowSums(is.na(obs_filtered)) / n) * 100,
  Var_Inflation = (T_var - U_bar) / U_bar * 100,
  Status = case_when(
    rownames(obs_filtered) %in% true_up_names   ~ "Upregulated",
    rownames(obs_filtered) %in% true_down_names ~ "Downregulated",
    TRUE                                        ~ "True Null"
  )
)

p3b <- ggplot(df_rubin_audit, aes(x = Missing_Pct, y = Var_Inflation)) +
  geom_point(aes(color = Status), alpha = 0.65, size = 1.3) +
  geom_smooth(method = "loess", color = PALETTE[["blue_main"]], fill = "gray85",
              linewidth = 0.7, se = TRUE) +
  scale_color_manual(values = c(
    "Upregulated"   = PALETTE[["red_strong"]],
    "Downregulated" = PALETTE[["blue_secondary"]],
    "True Null"     = PALETTE[["neutral"]]
  )) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 15)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10)) +
  labs(x = "Feature Missingness Rate (%)",
       y = "Rubin Variance Inflation Factor (%)",
       color = "") +
  theme_pub(J) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.5)
  )

fig3 <- wrap_plots(p3a, p3b, ncol = 2, widths = c(1.3, 1)) + pub_tags(J)

save_pub(fig3, "latex/figures/fig_19_3_differential_sensitivity",
         width = W, height = 95, units = "mm", journal = J)
save_pub(fig3, "docs/latex/figures/fig_19_3_differential_sensitivity",
         width = W, height = 95, units = "mm", journal = J)

# -----------------------------------------------------------------------------
# Figure 19.4: Templ Missingness Matrix Map, Margin Plot & Precision/ADCS Tradeoff
# -----------------------------------------------------------------------------
cat("\n--- Generating Figure 19.4: Templ Diagnostics ---\n")

# Panel A: VIM Sorted Matrix Plot
prot_order <- order(rowMeans(obs_filtered, na.rm = TRUE), decreasing = TRUE)
sorted_mat <- obs_filtered[prot_order, ]

mat_df <- as.data.frame(sorted_mat) %>%
  mutate(Protein_Rank = row_number()) %>%
  pivot_longer(cols = -Protein_Rank, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Is_Missing = is.na(Intensity))

p4a <- ggplot(mat_df, aes(x = Sample, y = Protein_Rank, fill = Is_Missing)) +
  geom_raster() +
  scale_fill_manual(values = c("FALSE" = "gray85", "TRUE" = PALETTE[["red_strong"]]),
                    labels = c("FALSE" = "Observed", "TRUE" = "Missing (Dropout)"),
                    name = "") +
  scale_y_reverse(breaks = c(1, 100, 200, 300, 386)) +
  labs(x = "Mass Spectrometry Sample",
       y = "Protein Rank (Abundance)") +
  theme_pub(J) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.0),
    legend.position = "top",
    legend.text = element_text(size = 5.5)
  )

# Panel B: VIM Margin Scatterplot
p_df <- tibble(
  Prot_1 = obs_filtered[1, ],
  Prot_2 = obs_filtered[2, ],
  Status = ifelse(is.na(obs_filtered[2, ]), "Missing in Prot 2", "Observed in Prot 2")
)

p4b <- ggplot(p_df, aes(x = Prot_1, y = ifelse(is.na(Prot_2), 18.0, Prot_2),
                        color = Status, shape = Status)) +
  geom_point(size = 2.0, alpha = 0.85) +
  geom_hline(yintercept = 18.3, linetype = "dashed", color = PALETTE[["red_strong"]], linewidth = 0.4) +
  scale_color_manual(values = c("Observed in Prot 2" = PALETTE[["blue_main"]],
                                "Missing in Prot 2"  = PALETTE[["red_strong"]])) +
  scale_shape_manual(values = c("Observed in Prot 2" = 16, "Missing in Prot 2" = 17)) +
  scale_y_continuous(limits = c(17.5, 27), breaks = seq(18, 26, 2)) +
  labs(x = expression("Protein 1 ("*log[2]*" Intensity)"),
       y = expression("Protein 2 ("*log[2]*" Intensity)"),
       color = "", shape = "") +
  theme_pub(J) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 5.5)
  )

# Panel C: Templ NRMSE vs ADCS Tradeoff
calc_templ_metrics_simple <- function(imputed_mat, name) {
  y_true <- truth_filtered[missing_mask]
  y_imp  <- imputed_mat[missing_mask]
  rmse <- sqrt(mean((y_imp - y_true)^2))
  nrmse <- rmse / sd(y_true)
  sigma_true <- cov(t(truth_filtered))
  sigma_imp  <- cov(t(imputed_mat))
  adcs <- (1 / p_filt) * norm(sigma_imp - sigma_true, type = "F")
  tibble(Method = name, NRMSE = nrmse, ADCS = adcs)
}

templ_df <- bind_rows(
  calc_templ_metrics_simple(mat_mindet, "MinDet"),
  calc_templ_metrics_simple(mat_perseus, "Perseus"),
  calc_templ_metrics_simple(mat_qrilc, "QRILC"),
  calc_templ_metrics_simple(mat_knn, "k-NN"),
  calc_templ_metrics_simple(mat_hybrid, "Hybrid")
)

p4c <- ggplot(templ_df, aes(x = NRMSE, y = ADCS, label = Method, color = Method)) +
  geom_point(size = 2.5) +
  geom_text_repel(fontface = "bold", size = 5.5 / .pt, box.padding = 0.3) +
  scale_color_manual(values = c(
    "MinDet"  = PALETTE[["red_strong"]],
    "Perseus" = PALETTE[["blue_main"]],
    "QRILC"   = PALETTE[["teal"]],
    "k-NN"    = PALETTE[["violet"]],
    "Hybrid"  = PALETTE[["blue_secondary"]]
  )) +
  scale_x_continuous(limits = c(1.15, 1.70), breaks = seq(1.2, 1.6, 0.1)) +
  scale_y_continuous(limits = c(0.28, 0.46), breaks = seq(0.30, 0.45, 0.05)) +
  labs(x = "Normalized Root Mean Squared Error (NRMSE)",
       y = "Covariance Difference (ADCS)") +
  theme_pub(J) +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(l = 8, r = 4))
  )

fig4 <- wrap_plots(p4a, p4b, p4c, design = "AB\nCC") + pub_tags(J)

save_pub(fig4, "latex/figures/fig_19_4_templ_diagnostics",
         width = W, height = 110, units = "mm", journal = J)
save_pub(fig4, "docs/latex/figures/fig_19_4_templ_diagnostics",
         width = W, height = 110, units = "mm", journal = J)

cat("\n=============================================================================\n")
cat("All Four Chapter 19 Figures Successfully Generated in figures4papers Style!\n")
cat("=============================================================================\n")
