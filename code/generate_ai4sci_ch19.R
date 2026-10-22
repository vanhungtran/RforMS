# =============================================================================
# generate_ai4sci_ch19.R
# Publication Figure Generation for Chapter 19: Missing Data in Mass Spectrometry
# Using the AI4SciFig Standards (theme_scifigure)
#
# Specifications:
# - Full 4-spine closed box border (0.5 pt solid black)
# - Outward-pointing tick marks (2.5 pt length)
# - Subtle dashed reference grid (#d9d9d9, 0.3 pt)
# - Publication typography: 6.5 - 8.0 pt base sizes
# - Nature Biomedical color palette (#0F4D92, #D95F02, #1B9E77, #7570B3, #E7298A)
# - Exact page margin fit: 165 mm width (matches LaTeX \textwidth at margin = 1 in)
# - Zero legend overlap: Legends placed outside plotting regions (bottom/right)
# - Zero overflow: Headroom on axes, ample panel heights, no clipping
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(patchwork)
  library(ggrepel)
  library(scales)
})

# Source AI4SciFig theme engine
source("R/theme_scifig.R")
source("R/scifig_plots.R")

# Load companion lab data
cat("Loading Chapter 19 Laboratory Environment...\n")
source("code/chapter19_missing_data_labs.R")

# Target directories
latex_fig_dir <- "latex/figures"
docs_fig_dir  <- "docs/latex/figures"
scifig_dir    <- "figures/scifig"

for (d in c(latex_fig_dir, docs_fig_dir, scifig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Master save function: 165 mm exact textwidth fit
save_ai4sci <- function(plot_obj, base_name, width_mm = 165, height_mm = 72, dpi = 300) {
  targets <- c(
    file.path(latex_fig_dir, base_name),
    file.path(docs_fig_dir, base_name),
    file.path(scifig_dir, base_name)
  )
  for (tgt in targets) {
    for (fmt in c("png", "pdf")) {
      ggsave(
        filename = paste0(tgt, ".", fmt),
        plot = plot_obj,
        width = width_mm,
        height = height_mm,
        units = "mm",
        dpi = dpi
      )
    }
  }
  cat(sprintf("  [SAVED] %s (PNG & PDF @ %.1f x %.1f mm)\n", base_name, width_mm, height_mm))
}

# =============================================================================
# Figure 19.1: Empirical Missingness Diagnosis and Group Completeness
# =============================================================================
cat("\n--- Generating Figure 19.1: Missingness Diagnosis (AI4Sci Style) ---\n")

p_filt <- nrow(obs_filtered)
n_samples <- ncol(obs_filtered)

mean_abund <- rowMeans(obs_filtered, na.rm = TRUE)
frac_miss  <- rowMeans(is.na(obs_filtered))
diag_df <- tibble(
  Protein = rownames(obs_filtered),
  Mean_Abundance = mean_abund,
  Missing_Fraction = frac_miss
)

b0 <- coef(log_fit)[1]
b1 <- coef(log_fit)[2]
curve_df <- tibble(
  x = seq(18.0, 27.5, length.out = 200),
  y = 1 / (1 + exp(-(b0 + b1 * x)))
)

# Panel A: Logistic diagnosis curve with clear corner label
p1a <- ggplot() +
  geom_point(data = diag_df, aes(x = Mean_Abundance, y = Missing_Fraction),
             colour = "#0F4D92", alpha = 0.45, size = 1.2) +
  geom_line(data = curve_df, aes(x = x, y = y),
            colour = "#D95F02", linewidth = 0.8) +
  geom_vline(xintercept = 20.2, linetype = "dashed", colour = "#444444", linewidth = 0.4) +
  annotate("label", x = 20.4, y = 0.64,
           label = "LOD = 20.2",
           size = 2.1, fill = "white", alpha = 0.9, hjust = 0, label.padding = unit(0.15, "lines")) +
  annotate("label", x = 27.2, y = 0.55,
           label = "Slope = -0.358 (p = 1.8e-07)\nMNAR Left-Censoring",
           size = 2.1, fill = "white", alpha = 0.9, hjust = 1, label.padding = unit(0.15, "lines")) +
  scale_x_continuous(limits = c(18.0, 27.5), breaks = seq(18, 26, 2)) +
  scale_y_continuous(limits = c(0, 0.70), breaks = seq(0, 0.60, 0.15),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "(a) Logistic Missingness Diagnosis",
    subtitle = "Dropout probability increases sharply below limit of detection",
    x = expression(Mean~Observed~Log[2]~Intensity),
    y = "Missingness Fraction"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    legend.position = "none",
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel B: Cohort Completeness Profile with bottom legend and ample headroom
sample_completeness <- colMeans(!is.na(obs_filtered)) * 100
sample_groups <- ifelse(grepl("^Ctrl", colnames(obs_filtered)), "Control", "Disease")

sample_df <- tibble(
  Sample = factor(colnames(obs_filtered), levels = colnames(obs_filtered)),
  Completeness = sample_completeness,
  Group = factor(sample_groups, levels = c("Control", "Disease"))
)

p1b <- ggplot(sample_df, aes(x = Sample, y = Completeness, fill = Group)) +
  geom_col(width = 0.68, alpha = 0.9, colour = "black", linewidth = 0.25) +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "#222222", linewidth = 0.5) +
  annotate("label", x = 20, y = 74,
           label = "70% Retention Threshold",
           size = 2.0, fill = "white", alpha = 0.9, hjust = 1, fontface = "bold", label.padding = unit(0.15, "lines")) +
  scale_fill_manual(values = c("Control" = "#0F4D92", "Disease" = "#D95F02")) +
  scale_y_continuous(limits = c(0, 105), breaks = seq(0, 100, 25),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "(b) Sample Completeness by Cohort",
    subtitle = "70% rule retains 96.5% of proteins with true conditional signatures",
    x = "Mass Spectrometry Run",
    y = "Quantified Features (%)",
    fill = "Cohort"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.2),
    legend.position = "bottom",
    legend.text = element_text(size = 6.0),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

fig1 <- (p1a + p1b) + plot_layout(widths = c(1, 1.15))
save_ai4sci(fig1, "fig_19_1_missingness_diagnosis", width_mm = 165, height_mm = 72)

# =============================================================================
# Figure 19.2: Post-Imputation Empirical Density Distributions
# =============================================================================
cat("\n--- Generating Figure 19.2: Imputation Distributions (AI4Sci Style) ---\n")

density_df <- bind_rows(
  tibble(Value = as.vector(truth_filtered), Strategy = "Ground Truth (Oracle)"),
  tibble(Value = as.vector(na.omit(as.vector(obs_filtered))), Strategy = "Observed (with Dropouts)"),
  tibble(Value = as.vector(mat_mindet), Strategy = "MinDet (Deterministic)"),
  tibble(Value = as.vector(mat_perseus), Strategy = "Perseus (Stochastic)"),
  tibble(Value = as.vector(mat_qrilc), Strategy = "QRILC (Quantile Reg.)"),
  tibble(Value = as.vector(mat_knn), Strategy = "k-NN (Donor Matching)")
) %>%
  mutate(Strategy = factor(Strategy, levels = c(
    "Ground Truth (Oracle)",
    "Observed (with Dropouts)",
    "MinDet (Deterministic)",
    "Perseus (Stochastic)",
    "QRILC (Quantile Reg.)",
    "k-NN (Donor Matching)"
  )))

strat_colors <- c(
  "Ground Truth (Oracle)"    = "#111111",
  "Observed (with Dropouts)" = "#777777",
  "MinDet (Deterministic)"   = "#D95F02",
  "Perseus (Stochastic)"     = "#0F4D92",
  "QRILC (Quantile Reg.)"    = "#1B9E77",
  "k-NN (Donor Matching)"    = "#7570B3"
)

strat_linetypes <- c(
  "Ground Truth (Oracle)"    = "dashed",
  "Observed (with Dropouts)" = "solid",
  "MinDet (Deterministic)"   = "solid",
  "Perseus (Stochastic)"     = "solid",
  "QRILC (Quantile Reg.)"    = "solid",
  "k-NN (Donor Matching)"    = "solid"
)

strat_linewidths <- c(
  "Ground Truth (Oracle)"    = 0.75,
  "Observed (with Dropouts)" = 0.65,
  "MinDet (Deterministic)"   = 0.75,
  "Perseus (Stochastic)"     = 0.75,
  "QRILC (Quantile Reg.)"    = 0.75,
  "k-NN (Donor Matching)"    = 0.75
)

# Position legend on RIGHT so it NEVER covers density peaks or spike annotation
p2 <- ggplot(density_df, aes(x = Value, colour = Strategy, linetype = Strategy, linewidth = Strategy)) +
  geom_density(adjust = 1.1, trim = FALSE) +
  scale_colour_manual(values = strat_colors) +
  scale_linetype_manual(values = strat_linetypes) +
  scale_linewidth_manual(values = strat_linewidths) +
  annotate("segment", x = 17.52, xend = 17.52, y = 0.38, yend = 0.22,
           arrow = arrow(length = unit(1.8, "mm"), type = "closed"),
           colour = "#D95F02", linewidth = 0.65) +
  annotate("label", x = 17.52, y = 0.42,
           label = "MinDet Point-Mass Spike\n(Variance Deflation)",
           size = 2.1, fill = "white", alpha = 0.9, colour = "#D95F02",
           fontface = "bold", label.padding = unit(0.15, "lines")) +
  scale_x_continuous(limits = c(14.5, 29.0), breaks = seq(16, 28, 2)) +
  scale_y_continuous(limits = c(0, 0.50), breaks = seq(0, 0.50, 0.10),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Empirical Intensity Density Distributions Post-Imputation",
    subtitle = "Stochastic models restore continuous left tails; MinDet introduces an artificial point-mass spike",
    x = expression(Protein~Log[2]~Intensity),
    y = "Empirical Probability Density",
    colour = "Strategy",
    linetype = "Strategy",
    linewidth = "Strategy"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 5.8),
    legend.key.width = unit(5, "mm"),
    plot.margin = margin(t = 4, r = 6, b = 4, l = 4, unit = "pt")
  )

save_ai4sci(p2, "fig_19_2_imputation_distributions", width_mm = 165, height_mm = 74)

# =============================================================================
# Figure 19.3: Differential Sensitivity and Variance Inflation
# =============================================================================
cat("\n--- Generating Figure 19.3: Differential Sensitivity & Variance Inflation (AI4Sci Style) ---\n")

get_de_df <- function(mat) {
  res <- run_limma(mat)
  tibble(
    Protein = rownames(mat),
    logFC = as.numeric(res$logFC),
    P_Value = as.numeric(res$pval),
    adj_P_Val = p.adjust(res$pval, method = "BH")
  )
}

de_truth  <- get_de_df(truth_filtered)
de_mindet <- get_de_df(mat_mindet)
de_qrilc  <- get_de_df(mat_qrilc)
de_proda  <- as_tibble(res_proda) %>%
  transmute(
    Protein = name,
    logFC = diff,
    P_Value = pval,
    adj_P_Val = adj_pval
  )

annotate_volcano <- function(df) {
  df %>%
    mutate(
      Status = case_when(
        Protein %in% paste0("PROT_", sprintf("%04d", 1:20))   ~ "True Upregulated (+2.0)",
        Protein %in% paste0("PROT_", sprintf("%04d", 31:60))  ~ "True Downregulated (-2.5)",
        TRUE                                                   ~ "True Null"
      ),
      Significant = adj_P_Val < 0.05 & abs(logFC) > 1.0,
      Status_Plot = factor(Status, levels = c("True Upregulated (+2.0)", "True Downregulated (-2.5)", "True Null"))
    )
}

plot_v_truth  <- annotate_volcano(de_truth)
plot_v_mindet <- annotate_volcano(de_mindet)
plot_v_qrilc  <- annotate_volcano(de_qrilc)
plot_v_proda  <- annotate_volcano(de_proda)

# 4 Volcano subplots arranged side-by-side in 1 row (clean, no overflow)
make_clean_volcano <- function(df, panel_tag, title_text, tp, fp, fdr, show_y = FALSE) {
  sub_text <- sprintf("TP = %d | FP = %d | FDR = %.1f%%", tp, fp, fdr)
  ggplot(df, aes(x = logFC, y = -log10(P_Value), colour = Status_Plot)) +
    geom_point(alpha = 0.65, size = 0.9, show.legend = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#555555", linewidth = 0.3) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", colour = "#555555", linewidth = 0.3) +
    scale_colour_manual(values = c(
      "True Upregulated (+2.0)"   = "#0F4D92",
      "True Downregulated (-2.5)" = "#D95F02",
      "True Null"                 = "#A5A5A5"
    )) +
    scale_x_continuous(limits = c(-5.5, 5.5), breaks = seq(-4, 4, 2)) +
    scale_y_continuous(limits = c(0, 8.5), breaks = seq(0, 8, 2)) +
    labs(
      title = paste0(panel_tag, " ", title_text),
      subtitle = sub_text,
      x = expression(Log[2]~FC),
      y = if (show_y) expression(-Log[10]~italic(P)) else ""
    ) +
    theme_scifigure(base_size = 6.2, grid = "both") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 6.8, face = "bold"),
      plot.subtitle = element_text(size = 5.5, colour = "#333333"),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
    )
}

v1 <- make_clean_volcano(plot_v_truth,  "(a)", "Oracle Truth",  48, 0, 0.0, show_y = TRUE)
v2 <- make_clean_volcano(plot_v_mindet, "(b)", "MinDet",        25, 2, 7.4, show_y = FALSE)
v3 <- make_clean_volcano(plot_v_qrilc,  "(c)", "QRILC",         29, 3, 9.4, show_y = FALSE)
v4 <- make_clean_volcano(plot_v_proda,  "(d)", "proDA",         20, 0, 0.0, show_y = FALSE)

volcano_row <- (v1 | v2 | v3 | v4)

# Panel B: Rubin Variance Inflation Factor
mi_plot_df <- rubin_summary %>%
  mutate(
    Missing_Pct = rowMeans(is.na(obs_filtered)) * 100,
    Status = case_when(
      Protein %in% paste0("PROT_", sprintf("%04d", 1:20))   ~ "True Upregulated (+2.0)",
      Protein %in% paste0("PROT_", sprintf("%04d", 31:60))  ~ "True Downregulated (-2.5)",
      TRUE                                                   ~ "True Null"
    ),
    Status_Plot = factor(Status, levels = c("True Upregulated (+2.0)", "True Downregulated (-2.5)", "True Null"))
  )

p3b <- ggplot(mi_plot_df, aes(x = Missing_Pct, y = Variance_Inflation, colour = Status_Plot)) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_smooth(method = "loess", colour = "#111111", linewidth = 0.7, se = TRUE, fill = "gray85") +
  scale_colour_manual(values = c(
    "True Upregulated (+2.0)"   = "#0F4D92",
    "True Downregulated (-2.5)" = "#D95F02",
    "True Null"                 = "#A5A5A5"
  )) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 15)) +
  scale_y_continuous(limits = c(0, 48), breaks = seq(0, 40, 10)) +
  labs(
    title = "(e) Rubin Variance Inflation Factor",
    subtitle = "Between-imputation variance expands confidence bounds as missingness increases",
    x = "Feature Missingness Rate (%)",
    y = "Rubin Variance Inflation (%)",
    colour = "Truth Status"
  ) +
  theme_scifigure(base_size = 6.8, grid = "both") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.8),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 4, r = 4, b = 2, l = 4, unit = "pt")
  )

# Stack row of 4 volcanoes above the Rubin plot; clean bottom legend from p3b
fig3 <- (volcano_row / p3b) +
  plot_layout(heights = c(1.1, 1.0))

save_ai4sci(fig3, "fig_19_3_differential_sensitivity", width_mm = 165, height_mm = 96)

# =============================================================================
# Figure 19.4: Multivariate Templ Diagnostics (Templ 2023)
# =============================================================================
cat("\n--- Generating Figure 19.4: Templ Diagnostics (AI4Sci Style) ---\n")

# Panel A: VIM Abundance-Sorted Matrix Map with bottom legend
prot_order <- order(rowMeans(obs_filtered, na.rm = TRUE), decreasing = TRUE)
sorted_mat <- obs_filtered[prot_order, ]

mat_df <- as.data.frame(sorted_mat) %>%
  mutate(Protein_Rank = row_number()) %>%
  pivot_longer(cols = -Protein_Rank, names_to = "Sample", values_to = "Intensity") %>%
  mutate(Is_Missing = is.na(Intensity))

p4a <- ggplot(mat_df, aes(x = Sample, y = Protein_Rank, fill = Is_Missing)) +
  geom_raster() +
  scale_fill_manual(values = c("FALSE" = "#E8E8E8", "TRUE" = "#D95F02"),
                    labels = c("FALSE" = "Observed", "TRUE" = "Missing (Dropout)"),
                    name = "") +
  scale_y_reverse(breaks = c(1, 100, 200, 300, 386)) +
  labs(
    title = "(a) VIM Missingness Matrix Map",
    subtitle = "Dropouts cluster heavily at lower abundance ranks",
    x = "Mass Spectrometry Run (Sample)",
    y = "Protein Rank (Abundance Order)"
  ) +
  theme_scifigure(base_size = 6.8, grid = "none") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4.8),
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel B: Bivariate Margin Scatterplot with bottom legend
p_df <- tibble(
  Prot_1 = obs_filtered[1, ],
  Prot_2 = obs_filtered[2, ],
  Status = ifelse(is.na(obs_filtered[2, ]), "Missing in Prot 2", "Observed in Prot 2")
)

p4b <- ggplot(p_df, aes(x = Prot_1, y = ifelse(is.na(Prot_2), 18.0, Prot_2),
                        colour = Status, shape = Status)) +
  geom_point(size = 2.0, alpha = 0.85) +
  geom_hline(yintercept = 18.3, linetype = "dashed", colour = "#D95F02", linewidth = 0.45) +
  scale_colour_manual(values = c("Observed in Prot 2" = "#0F4D92",
                                "Missing in Prot 2"  = "#D95F02")) +
  scale_shape_manual(values = c("Observed in Prot 2" = 16, "Missing in Prot 2" = 17)) +
  scale_y_continuous(limits = c(17.5, 27.0), breaks = seq(18, 26, 2)) +
  labs(
    title = "(b) Bivariate Margin Scatterplot",
    subtitle = "Marginal dots below line expose conditional dropouts",
    x = expression(Protein~1~(Log[2]~Intensity)),
    y = expression(Protein~2~(Log[2]~Intensity)),
    colour = "",
    shape = ""
  ) +
  theme_scifigure(base_size = 6.8, grid = "both") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel C: Templ NRMSE vs ADCS Tradeoff
calc_templ_metrics <- function(imputed_mat, name) {
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
  calc_templ_metrics(mat_mindet, "MinDet"),
  calc_templ_metrics(mat_perseus, "Perseus"),
  calc_templ_metrics(mat_qrilc, "QRILC"),
  calc_templ_metrics(mat_knn, "k-NN"),
  calc_templ_metrics(mat_hybrid, "Hybrid")
)

p4c <- ggplot(templ_df, aes(x = NRMSE, y = ADCS, label = Method, colour = Method)) +
  geom_point(size = 2.8) +
  geom_text_repel(fontface = "bold", size = 2.6, box.padding = 0.35, point.padding = 0.25) +
  annotate("label", x = 1.20, y = 0.44,
           label = "Optimal Frontier: Low NRMSE & Low ADCS\n(Perseus balances individual error and covariance preservation)",
           size = 2.0, fill = "#F0F7FF", colour = "#0F4D92", fontface = "italic",
           hjust = 0, label.padding = unit(0.2, "lines")) +
  scale_colour_manual(values = c(
    "MinDet"  = "#D95F02",
    "Perseus" = "#0F4D92",
    "QRILC"   = "#1B9E77",
    "k-NN"    = "#7570B3",
    "Hybrid"  = "#E7298A"
  )) +
  scale_x_continuous(limits = c(1.15, 1.70), breaks = seq(1.2, 1.6, 0.1)) +
  scale_y_continuous(limits = c(0.28, 0.46), breaks = seq(0.30, 0.45, 0.05)) +
  labs(
    title = "(c) Multi-Criterion Quality Tradeoff (Templ 2023)",
    subtitle = "Perseus achieves lowest NRMSE (1.218) and lowest covariance distortion (0.313)",
    x = "Normalized Root Mean Squared Error (NRMSE)",
    y = "Covariance Difference (ADCS)"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    legend.position = "none",
    axis.title.y = element_text(margin = margin(l = 6, r = 4)),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1.05, 1.0))
save_ai4sci(fig4, "fig_19_4_templ_diagnostics", width_mm = 165, height_mm = 112)

# =============================================================================
# Figure 19.5: Parameter Sensitivity Auditing and Cross-Method Concordance
# =============================================================================
cat("\n--- Generating Figure 19.5: Sensitivity & Concordance (AI4Sci Style) ---\n")

# Panel A: Perseus Hyperparameter Tradeoff
p5a_df <- perseus_sens_df %>%
  mutate(Omega_Label = paste0("omega == ", omega))

p5a <- ggplot(p5a_df, aes(x = NRMSE, y = Bias_Down, colour = factor(omega), shape = factor(omega))) +
  geom_line(aes(group = omega), linewidth = 0.5, linetype = "dashed", alpha = 0.7) +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = sprintf("delta==%.1f", delta)), parse = TRUE,
                  size = 2.1, box.padding = 0.25, point.padding = 0.2, show.legend = FALSE) +
  scale_colour_manual(values = c("0.2" = "#1B9E77", "0.3" = "#0F4D92", "0.4" = "#D95F02"),
                      name = expression(Width~shrinkage~(omega))) +
  scale_shape_manual(values = c("0.2" = 15, "0.3" = 16, "0.4" = 17),
                     name = expression(Width~shrinkage~(omega))) +
  labs(
    title = "(a) Perseus Hyperparameter Tradeoff",
    subtitle = "Higher downshift delta reduces fold-change bias but increases NRMSE error",
    x = "Normalized Root Mean Squared Error (NRMSE)",
    y = expression(Downward~Fold-Change~Bias~(Log[2]))
  ) +
  theme_scifigure(base_size = 6.8, grid = "both") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    legend.title = element_text(size = 5.8, face = "bold"),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel B: Hybrid Segmentation Cutoff Transition
p5b <- ggplot(hybrid_sens_df, aes(x = Cutoff_tau)) +
  geom_line(aes(y = Bias_Down, colour = "Downward Bias (log2 FC)"), linewidth = 0.6) +
  geom_point(aes(y = Bias_Down, colour = "Downward Bias (log2 FC)"), size = 2.0) +
  geom_line(aes(y = FDR_pct / 10, colour = "False Discovery Rate (/10%)"), linewidth = 0.6, linetype = "dashed") +
  geom_point(aes(y = FDR_pct / 10, colour = "False Discovery Rate (/10%)"), size = 2.0, shape = 17) +
  annotate("rect", xmin = 21.8, xmax = 22.2, ymin = 0, ymax = 1.8,
           alpha = 0.12, fill = "#1B9E77") +
  annotate("text", x = 22.0, y = 1.75, label = "Optimal Cutoff: tau = 22.0\n(0.0% FDR, Bias = +0.906)",
           size = 2.0, fontface = "italic", colour = "#1B9E77", hjust = 0.5) +
  scale_y_continuous(
    name = expression(Downward~Fold-Change~Bias~(Log[2])),
    sec.axis = sec_axis(~ . * 10, name = "Empirical False Discovery Rate (%)")
  ) +
  scale_x_continuous(breaks = c(21.2, 21.6, 22.0, 22.4, 22.8)) +
  scale_colour_manual(
    values = c("Downward Bias (log2 FC)" = "#0F4D92", "False Discovery Rate (/10%)" = "#D95F02"),
    name = ""
  ) +
  labs(
    title = "(b) Hybrid Cutoff Optimization",
    subtitle = "tau = 22.0 balances MAR donor matching and MNAR tail imputation",
    x = expression(Abundance~Segmentation~Cutoff~(tau))
  ) +
  theme_scifigure(base_size = 6.8, grid = "both") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel C: Cross-Method Biomarker Concordance Heatmap
jaccard_df <- as.data.frame(as.table(jaccard_mat)) %>%
  rename(Method_A = Var1, Method_B = Var2, Jaccard = Freq) %>%
  mutate(
    Method_A = factor(Method_A, levels = rev(colnames(jaccard_mat))),
    Method_B = factor(Method_B, levels = colnames(jaccard_mat))
  )

p5c <- ggplot(jaccard_df, aes(x = Method_B, y = Method_A, fill = Jaccard)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", Jaccard),
                colour = ifelse(Jaccard > 0.6, "white", "black")),
            size = 2.0, fontface = "bold", show.legend = FALSE) +
  scale_fill_gradient(low = "#F7FBFF", high = "#084594", limits = c(0, 1),
                      name = "Jaccard Index") +
  scale_colour_manual(values = c("white" = "white", "black" = "#202020")) +
  labs(
    title = "(c) Cross-Method Biomarker Concordance",
    subtitle = "Single imputation methods cluster together (J > 0.75), sharing false positives",
    x = "",
    y = ""
  ) +
  theme_scifigure(base_size = 6.8, grid = "none") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5),
    axis.text.y = element_text(size = 5.5),
    legend.position = "bottom",
    legend.key.width = unit(15, "pt"),
    legend.key.height = unit(5, "pt"),
    legend.text = element_text(size = 5.2),
    legend.title = element_text(size = 5.5, face = "bold"),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Panel D: Rubin MI Information Loss & Variance Inflation
p5d_df <- fmi_binned %>%
  mutate(Miss_Bin = factor(Miss_Bin, levels = c("0%", "1-10%", "11-20%", "21-30%", ">30%")))

p5d <- ggplot(p5d_df, aes(x = Miss_Bin)) +
  geom_col(aes(y = Mean_RVI_Pct, fill = "Variance Inflation (%)"), width = 0.55, alpha = 0.85) +
  geom_errorbar(aes(ymin = Mean_RVI_Pct, ymax = Max_RVI_Pct), width = 0.2, colour = "#0F4D92", linewidth = 0.5) +
  geom_line(aes(y = Mean_FMI_Pct * 2, group = 1, colour = "Fraction Missing Info (x2 %)"), linewidth = 0.6) +
  geom_point(aes(y = Mean_FMI_Pct * 2, colour = "Fraction Missing Info (x2 %)"), size = 2.0, shape = 17) +
  geom_text(aes(y = Max_RVI_Pct + 1.2, label = sprintf("nu[adj]==%.1f", Mean_df_adj)),
            parse = TRUE, size = 2.0, colour = "#202020") +
  scale_y_continuous(
    name = "Mean & Max Relative Variance Increase (RVI %)",
    limits = c(0, 24),
    breaks = seq(0, 20, 5),
    sec.axis = sec_axis(~ . / 2, name = "Mean Fraction Missing Info (FMI %)")
  ) +
  scale_fill_manual(values = c("Relative Variance Increase (%)" = "#4292C6"), name = "") +
  scale_colour_manual(values = c("Fraction Missing Info (x2 %)" = "#D95F02"), name = "") +
  labs(
    title = "(d) Rubin MI Information Loss & Relative Variance Increase",
    subtitle = "RVI expands up to 20.8% and df shrinks from 18.0 to 14.6 as missingness increases",
    x = "Protein Missingness Rate Strata"
  ) +
  theme_scifigure(base_size = 6.8, grid = "both") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 5.5),
    legend.margin = margin(t = -2, b = 0),
    plot.margin = margin(t = 3, r = 4, b = 3, l = 3, unit = "pt")
  )

# Assemble Figure 19.5 (2x2 layout, 165 mm width x 125 mm height)
fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_layout(heights = c(1.0, 1.15))

save_ai4sci(fig5, "fig_19_5_sensitivity_concordance", width_mm = 165, height_mm = 125)

cat("\n=============================================================================\n")
cat("All Five Chapter 19 Figures Successfully Generated in AI4Sci Style (165 mm Margin Fit)!\n")
cat("=============================================================================\n")
