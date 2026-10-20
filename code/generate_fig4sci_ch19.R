# =============================================================================
# generate_fig4sci_ch19.R
# Publication Figure Generation for Chapter 19: Missing Data in Mass Spectrometry
# Using the AI4SciFig / Fig4sci Scientific Figure Standards (theme_scifigure)
#
# Standards:
# - Full 4-spine closed box border (0.5 pt solid black)
# - Outward-pointing tick marks (2.5 pt length)
# - Subtle dashed reference grid (#d9d9d9, 0.3 pt)
# - Publication typography: 7.0 - 8.5 pt font sizes
# - Nature Biomedical colorblind-safe palette (#0F4D92, #D95F02, #1B9E77, etc.)
# - True physical dimensions: Double column (180 mm)
# - Dual export: 300/600 DPI PNG + vector PDF
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

# Source AI4SciFig theme and plot engine
source("R/theme_scifig.R")
source("R/scifig_plots.R")

# Ensure simulation and lab data are generated
cat("Loading Chapter 19 Laboratory Environment...\n")
source("code/chapter19_missing_data_labs.R")

# Output directory targets
latex_fig_dir <- "latex/figures"
docs_fig_dir  <- "docs/latex/figures"
scifig_dir    <- "figures/scifig"

for (d in c(latex_fig_dir, docs_fig_dir, scifig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Helper to save to all target locations
save_fig4sci_all <- function(plot_obj, base_name, width_mm = 180, height_mm = 75, dpi = 300) {
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
  cat(sprintf("  [SAVED] %s (PNG & PDF @ %d x %d mm)\n", base_name, round(width_mm), round(height_mm)))
}

# =============================================================================
# Figure 19.1: Empirical Missingness Diagnosis and Group Completeness
# =============================================================================
cat("\n--- Generating Figure 19.1: Missingness Diagnosis (Fig4sci Style) ---\n")

# Panel A: Logistic missingness diagnosis curve
p_filt <- nrow(obs_filtered)
n_samples <- ncol(obs_filtered)

mean_abund <- rowMeans(obs_filtered, na.rm = TRUE)
frac_miss  <- rowMeans(is.na(obs_filtered))
diag_df <- tibble(
  Protein = rownames(obs_filtered),
  Mean_Abundance = mean_abund,
  Missing_Fraction = frac_miss
)

# Logistic model prediction curve
b0 <- coef(log_fit)[1]
b1 <- coef(log_fit)[2]
curve_df <- tibble(
  x = seq(18.0, 27.5, length.out = 200),
  y = 1 / (1 + exp(-(b0 + b1 * x)))
)

p1a <- ggplot() +
  geom_point(data = diag_df, aes(x = Mean_Abundance, y = Missing_Fraction),
             colour = "#0F4D92", alpha = 0.45, size = 1.3) +
  geom_line(data = curve_df, aes(x = x, y = y),
            colour = "#D95F02", linewidth = 0.85) +
  geom_vline(xintercept = 20.2, linetype = "dashed", colour = "#444444", linewidth = 0.4) +
  annotate("label", x = 20.4, y = 0.62,
           label = "Detection Limit\n(LOD = 20.2)",
           size = 2.1, fill = "white", alpha = 0.85, hjust = 0, label.padding = unit(0.15, "lines")) +
  annotate("label", x = 24.2, y = 0.52,
           label = "Logit Slope = -0.358\n(p = 1.8e-07)\nMNAR Left-Censoring",
           size = 2.1, fill = "white", alpha = 0.85, label.padding = unit(0.15, "lines")) +
  scale_x_continuous(limits = c(18.0, 27.5), breaks = seq(18, 26, 2)) +
  scale_y_continuous(limits = c(0, 0.70), breaks = seq(0, 0.60, 0.15),
                     labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "(a) Logistic Missingness Diagnosis",
    subtitle = "Abundance-dependent dropout confirms left-censored MNAR",
    x = expression(Mean~Observed~Log[2]~Intensity),
    y = "Missingness Fraction Across Cohort"
  ) +
  theme_scifigure(base_size = 7.5, grid = "both")

# Panel B: Cohort Completeness Profile
sample_completeness <- colMeans(!is.na(obs_filtered)) * 100
sample_groups <- ifelse(grepl("^Ctrl", colnames(obs_filtered)), "Control", "Disease")

sample_df <- tibble(
  Sample = factor(colnames(obs_filtered), levels = colnames(obs_filtered)),
  Completeness = sample_completeness,
  Group = factor(sample_groups, levels = c("Control", "Disease"))
)

p1b <- ggplot(sample_df, aes(x = Sample, y = Completeness, fill = Group)) +
  geom_col(width = 0.68, alpha = 0.9, colour = "black", linewidth = 0.3) +
  geom_hline(yintercept = 70, linetype = "dashed", colour = "#222222", linewidth = 0.5) +
  annotate("label", x = 10.5, y = 73,
           label = "70% Group-Wise Retention Threshold",
           size = 2.2, fill = "white", alpha = 0.9, fontface = "bold", label.padding = unit(0.2, "lines")) +
  scale_fill_manual(values = c("Control" = "#0F4D92", "Disease" = "#D95F02")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                     expand = expansion(mult = c(0, 0.04))) +
  labs(
    title = "(b) Sample Completeness by Cohort",
    subtitle = "Group-wise rule retains 96.5% of proteins with conditional signatures",
    x = "Mass Spectrometry Run (Sample)",
    y = "Observed Features (%)",
    fill = "Cohort"
  ) +
  theme_scifigure(base_size = 7.5, grid = "both") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6.0),
    legend.position = "top",
    legend.margin = margin(b = -2)
  )

fig1 <- (p1a + p1b) + plot_layout(widths = c(1, 1.15))
save_fig4sci_all(fig1, "fig_19_1_missingness_diagnosis", width_mm = 180, height_mm = 74)

# =============================================================================
# Figure 19.2: Post-Imputation Empirical Density Distributions
# =============================================================================
cat("\n--- Generating Figure 19.2: Imputation Distributions (Fig4sci Style) ---\n")

density_df <- bind_rows(
  tibble(Value = as.vector(truth_filtered), Strategy = "Ground Truth (Oracle)"),
  tibble(Value = as.vector(na.omit(as.vector(obs_filtered))), Strategy = "Observed (with Dropouts)"),
  tibble(Value = as.vector(mat_mindet), Strategy = "MinDet (Deterministic)"),
  tibble(Value = as.vector(mat_perseus), Strategy = "Perseus (Stochastic)"),
  tibble(Value = as.vector(mat_qrilc), Strategy = "QRILC (Quantile Reg.)"),
  tibble(Value = as.vector(mat_knn), Strategy = "k-NN (Donor Matching)"),
  tibble(Value = as.vector(mat_hybrid), Strategy = "Hybrid MAR/MNAR")
) %>%
  mutate(Strategy = factor(Strategy, levels = c(
    "Ground Truth (Oracle)",
    "Observed (with Dropouts)",
    "MinDet (Deterministic)",
    "Perseus (Stochastic)",
    "QRILC (Quantile Reg.)",
    "k-NN (Donor Matching)",
    "Hybrid MAR/MNAR"
  )))

strat_colors <- c(
  "Ground Truth (Oracle)"    = "#111111",
  "Observed (with Dropouts)" = "#777777",
  "MinDet (Deterministic)"   = "#D95F02",
  "Perseus (Stochastic)"     = "#0F4D92",
  "QRILC (Quantile Reg.)"    = "#1B9E77",
  "k-NN (Donor Matching)"    = "#7570B3",
  "Hybrid MAR/MNAR"          = "#E7298A"
)

strat_linetypes <- c(
  "Ground Truth (Oracle)"    = "dashed",
  "Observed (with Dropouts)" = "solid",
  "MinDet (Deterministic)"   = "solid",
  "Perseus (Stochastic)"     = "solid",
  "QRILC (Quantile Reg.)"    = "solid",
  "k-NN (Donor Matching)"    = "solid",
  "Hybrid MAR/MNAR"          = "solid"
)

strat_linewidths <- c(
  "Ground Truth (Oracle)"    = 0.80,
  "Observed (with Dropouts)" = 0.65,
  "MinDet (Deterministic)"   = 0.75,
  "Perseus (Stochastic)"     = 0.75,
  "QRILC (Quantile Reg.)"    = 0.75,
  "k-NN (Donor Matching)"    = 0.75,
  "Hybrid MAR/MNAR"          = 0.75
)

p2 <- ggplot(density_df, aes(x = Value, colour = Strategy, linetype = Strategy, linewidth = Strategy)) +
  geom_density(adjust = 1.1) +
  scale_colour_manual(values = strat_colors) +
  scale_linetype_manual(values = strat_linetypes) +
  scale_linewidth_manual(values = strat_linewidths) +
  # Arrow pointing down to the MinDet spike at x = 17.52
  annotate("segment", x = 17.52, xend = 17.52, y = 0.23, yend = 0.16,
           arrow = arrow(length = unit(2.2, "mm"), type = "closed"),
           colour = "#D95F02", linewidth = 0.7) +
  annotate("label", x = 17.52, y = 0.26,
           label = "MinDet Point-Mass Spike\n(Variance Deflation)",
           size = 2.2, fill = "white", alpha = 0.9, colour = "#D95F02",
           fontface = "bold", label.padding = unit(0.2, "lines")) +
  scale_x_continuous(limits = c(15.5, 28.0), breaks = seq(16, 28, 2)) +
  scale_y_continuous(limits = c(0, 0.30), breaks = seq(0, 0.30, 0.05)) +
  labs(
    title = "Empirical Intensity Density Distributions Post-Imputation",
    subtitle = "Stochastic models (Perseus, QRILC) restore continuous left tails; MinDet introduces an artificial spike",
    x = expression(Protein~Log[2]~Intensity),
    y = "Empirical Probability Density",
    colour = "Strategy",
    linetype = "Strategy",
    linewidth = "Strategy"
  ) +
  theme_scifigure(base_size = 7.5, grid = "both") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6.2),
    legend.margin = margin(b = -2)
  )

save_fig4sci_all(p2, "fig_19_2_imputation_distributions", width_mm = 180, height_mm = 72)

# =============================================================================
# Figure 19.3: Differential Sensitivity and Variance Inflation
# =============================================================================
cat("\n--- Generating Figure 19.3: Differential Sensitivity & Variance Inflation (Fig4sci Style) ---\n")

# Helper to extract differential abundance using the chapter's moderated eBayes model
get_de_df <- function(mat) {
  res <- run_limma(mat)
  tibble(
    Protein = rownames(mat),
    logFC = as.numeric(res$logFC),
    P_Value = as.numeric(res$pval),
    adj_P_Val = p.adjust(res$pval, method = "BH")
  )
}

# DE dataframes
de_truth  <- get_de_df(truth_filtered)
de_mindet <- get_de_df(mat_mindet)
de_qrilc  <- get_de_df(mat_qrilc)

# proDA results
de_proda <- as_tibble(res_proda) %>%
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
        Protein %in% paste0("PROT_", sprintf("%04d", 1:20))   ~ "True Upregulated",
        Protein %in% paste0("PROT_", sprintf("%04d", 31:60))  ~ "True Downregulated",
        TRUE                                                   ~ "True Null"
      ),
      Significant = adj_P_Val < 0.05 & abs(logFC) > 1.0,
      Status_Plot = factor(Status, levels = c("True Upregulated", "True Downregulated", "True Null"))
    )
}

plot_v_truth  <- annotate_volcano(de_truth)
plot_v_mindet <- annotate_volcano(de_mindet)
plot_v_qrilc  <- annotate_volcano(de_qrilc)
plot_v_proda  <- annotate_volcano(de_proda)

# Volcano subplot helper
make_mini_volcano <- function(df, panel_tag, title_text, tp, fp, fdr) {
  sub_text <- sprintf("TP = %d | FP = %d | FDR = %.1f%%", tp, fp, fdr)
  ggplot(df, aes(x = logFC, y = -log10(P_Value), colour = Status_Plot)) +
    geom_point(alpha = 0.65, size = 1.0) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "#555555", linewidth = 0.35) +
    geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed", colour = "#555555", linewidth = 0.35) +
    scale_colour_manual(values = c(
      "True Upregulated"   = "#0F4D92",
      "True Downregulated" = "#D95F02",
      "True Null"          = "#A5A5A5"
    )) +
    scale_x_continuous(limits = c(-4.8, 4.8), breaks = seq(-4, 4, 2)) +
    scale_y_continuous(limits = c(0, 7.5), breaks = seq(0, 6, 2)) +
    labs(
      title = paste0(panel_tag, " ", title_text),
      subtitle = sub_text,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~italic(P)-value)
    ) +
    theme_scifigure(base_size = 6.5, grid = "both") +
    theme(
      legend.position = "none",
      plot.title = element_text(size = 7.0, face = "bold"),
      plot.subtitle = element_text(size = 6.0, colour = "#333333")
    )
}

v1 <- make_mini_volcano(plot_v_truth,  "(a)", "Truth (Oracle)",  48, 0, 0.0)
v2 <- make_mini_volcano(plot_v_mindet, "(b)", "MinDet",          25, 2, 7.4)
v3 <- make_mini_volcano(plot_v_qrilc,  "(c)", "QRILC",           29, 3, 9.4)
v4 <- make_mini_volcano(plot_v_proda,  "(d)", "proDA (Dropout)", 20, 0, 0.0)

volcano_grid <- (v1 | v2) / (v3 | v4)

# Panel B: Rubin Variance Inflation Factor
mi_plot_df <- rubin_summary %>%
  mutate(
    Missing_Pct = rowMeans(is.na(obs_filtered)) * 100,
    Status = case_when(
      Protein %in% paste0("PROT_", sprintf("%04d", 1:20))   ~ "True Upregulated",
      Protein %in% paste0("PROT_", sprintf("%04d", 31:60))  ~ "True Downregulated",
      TRUE                                                   ~ "True Null"
    ),
    Status_Plot = factor(Status, levels = c("True Upregulated", "True Downregulated", "True Null"))
  )

p3b <- ggplot(mi_plot_df, aes(x = Missing_Pct, y = Variance_Inflation, colour = Status_Plot)) +
  geom_point(alpha = 0.65, size = 1.3) +
  geom_smooth(method = "loess", colour = "#111111", linewidth = 0.75, se = FALSE) +
  scale_colour_manual(values = c(
    "True Upregulated"   = "#0F4D92",
    "True Downregulated" = "#D95F02",
    "True Null"          = "#A5A5A5"
  )) +
  scale_x_continuous(limits = c(0, 65), breaks = seq(0, 60, 15)) +
  scale_y_continuous(limits = c(0, 45), breaks = seq(0, 40, 10)) +
  labs(
    title = "(e) Rubin Variance Inflation",
    subtitle = "Between-imputation variance widens confidence bounds",
    x = "Feature Missingness Rate (%)",
    y = "Rubin Variance Inflation Factor (%)",
    colour = "Truth Status"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 6.0),
    legend.margin = margin(b = -2)
  )

fig3 <- (volcano_grid | p3b) + plot_layout(widths = c(1.35, 1.0))
save_fig4sci_all(fig3, "fig_19_3_differential_sensitivity", width_mm = 180, height_mm = 95)

# =============================================================================
# Figure 19.4: Multivariate Templ Diagnostics (Templ 2023)
# =============================================================================
cat("\n--- Generating Figure 19.4: Templ Diagnostics (Fig4sci Style) ---\n")

# Panel A: VIM Abundance-Sorted Matrix Map
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
                    name = "Status") +
  scale_y_reverse(breaks = c(1, 100, 200, 300, 386)) +
  labs(
    title = "(a) VIM Missingness Matrix Map",
    subtitle = "Dropouts cluster heavily at lower abundance ranks",
    x = "Mass Spectrometry Run (Sample)",
    y = "Protein Rank (Abundance Order)"
  ) +
  theme_scifigure(base_size = 7.0, grid = "none") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5.2),
    legend.position = "top",
    legend.margin = margin(b = -2)
  )

# Panel B: Bivariate Margin Scatterplot
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
    subtitle = "Marginal distribution below dashed line exposes conditional dropouts",
    x = expression(Protein~1~(Log[2]~Intensity)),
    y = expression(Protein~2~(Log[2]~Intensity)),
    colour = "Observation",
    shape = "Observation"
  ) +
  theme_scifigure(base_size = 7.0, grid = "both") +
  theme(
    legend.position = "top",
    legend.margin = margin(b = -2)
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
  geom_point(size = 3.0) +
  geom_text_repel(fontface = "bold", size = 2.7, box.padding = 0.35, point.padding = 0.3) +
  annotate("label", x = 1.25, y = 0.30,
           label = "Optimal Frontier\n(Low NRMSE & Low ADCS)",
           size = 2.2, fill = "#F0F7FF", colour = "#0F4D92", fontface = "italic",
           label.padding = unit(0.2, "lines")) +
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
    subtitle = "Perseus balances element-wise precision (NRMSE) and covariance preservation (ADCS)",
    x = "Normalized Root Mean Squared Error (NRMSE)",
    y = "Covariance Structure Difference (ADCS)"
  ) +
  theme_scifigure(base_size = 7.5, grid = "both") +
  theme(
    legend.position = "none"
  )

fig4 <- (p4a | p4b) / p4c + plot_layout(heights = c(1.1, 1.0))
save_fig4sci_all(fig4, "fig_19_4_templ_diagnostics", width_mm = 180, height_mm = 118)

cat("\n=============================================================================\n")
cat("All Four Chapter 19 Figures Successfully Generated in Fig4sci Style!\n")
cat("=============================================================================\n")
