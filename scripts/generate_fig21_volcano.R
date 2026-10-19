# =============================================================================
# scripts/generate_fig21_volcano.R
# Publication-Grade Volcano Plot Generator (AI4SciFig / Nature Standards)
# For Figure 21.1: 11-Site DIA Reproducibility Study (PXD004886)
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(scales)
})

# Ingest differential abundance dataset
de_file <- "data/pxd004886/DE_results_annotated.csv"
if (!file.exists(de_file)) {
  stop("Input file not found: ", de_file)
}
de_df <- read.csv(de_file, stringsAsFactors = FALSE)

# Clean annotations and classify status
de_df <- de_df %>%
  mutate(
    # Clean gene name
    Gene_Clean = case_when(
      grepl("AQUA30", Genename) | grepl("AQUA30", Protein) ~ "AQUA30 (Spike-in)",
      is.na(Genename) | Genename == "" ~ Protein,
      TRUE ~ Genename
    ),
    neg_log10_padj = -log10(padj),
    # Categorical classification
    Status = case_when(
      grepl("AQUA30", Gene_Clean) ~ "Spike-In Standard (Group B)",
      padj < 0.05 & log2FC > 1.0  ~ "Up in Group B (Endogenous)",
      padj < 0.05 & log2FC < -1.0 ~ "Up in Group A (Endogenous)",
      TRUE                        ~ "Not Significant"
    ),
    Status = factor(Status, levels = c(
      "Up in Group A (Endogenous)",
      "Up in Group B (Endogenous)",
      "Spike-In Standard (Group B)",
      "Not Significant"
    ))
  )

# Separate subsets for layered rendering
df_ns   <- de_df %>% filter(Status == "Not Significant")
df_endo <- de_df %>% filter(Status %in% c("Up in Group A (Endogenous)", "Up in Group B (Endogenous)"))
df_aqua <- de_df %>% filter(Status == "Spike-In Standard (Group B)")
df_sig  <- de_df %>% filter(Status != "Not Significant")

# Construct individual label positions and texts
df_sig <- df_sig %>%
  mutate(
    Label_Text = sprintf("%s\nlog2FC = %+.2f, FDR = %.3f", Gene_Clean, log2FC, padj),
    # Explicit custom nudge targets to ensure zero collisions
    nudge_x = case_when(
      Gene_Clean == "TRA2A"             ~ -0.90,
      Gene_Clean == "ERAP2"             ~ -0.55,
      Gene_Clean == "CDC34"             ~ -0.95,
      Gene_Clean == "RSL24D1"           ~ -0.95,
      Gene_Clean == "NMNAT1"            ~  0.85,
      grepl("AQUA30", Gene_Clean)       ~  0.85,
      TRUE                              ~  0.00
    ),
    nudge_y = case_when(
      Gene_Clean == "TRA2A"             ~ -0.22,
      Gene_Clean == "ERAP2"             ~  0.35,
      Gene_Clean == "CDC34"             ~  0.25,
      Gene_Clean == "RSL24D1"           ~ -0.42,
      Gene_Clean == "NMNAT1"            ~  0.18,
      grepl("AQUA30", Gene_Clean)       ~ -0.10,
      TRUE                              ~  0.00
    )
  )

# Define color palette (Nature Publishing Group aesthetic)
status_colors <- c(
  "Up in Group A (Endogenous)"   = "#1D4ED8",  # Royal Cobalt Blue
  "Up in Group B (Endogenous)"   = "#DC2626",  # Vibrant Ruby Crimson
  "Spike-In Standard (Group B)"  = "#D97706",  # Luminous Amber/Gold
  "Not Significant"              = "#94A3B8"   # Slate Grey
)

# Plot limits
x_lim <- c(-3.4, 3.4)
y_max <- 3.4

# Build ggplot2 canvas
p_volcano <- ggplot(de_df, aes(x = log2FC, y = neg_log10_padj)) +
  # 1. Subtle shaded significance quadrant backgrounds
  annotate("rect", xmin = -Inf, xmax = -1.0, ymin = -log10(0.05), ymax = Inf,
           fill = "#EFF6FF", alpha = 0.55) +
  annotate("rect", xmin = 1.0, xmax = Inf, ymin = -log10(0.05), ymax = Inf,
           fill = "#FFF1F2", alpha = 0.55) +

  # 2. Subtle 2D density contour lines for the bulk null distribution
  stat_density_2d(
    data = df_ns,
    aes(x = log2FC, y = neg_log10_padj),
    colour = "#CBD5E1",
    linewidth = 0.3,
    alpha = 0.70,
    bins = 7
  ) +

  # 3. Non-significant background proteins (plume points)
  geom_point(
    data = df_ns,
    aes(x = log2FC, y = neg_log10_padj, colour = Status),
    alpha = 0.30,
    size = 1.1,
    stroke = 0
  ) +

  # 4. Endogenous significant proteins
  geom_point(
    data = df_endo,
    aes(x = log2FC, y = neg_log10_padj, colour = Status),
    size = 3.3,
    alpha = 0.95
  ) +

  # 5. Positive spike-in standard (highlighted diamond with dark gold border)
  geom_point(
    data = df_aqua,
    aes(x = log2FC, y = neg_log10_padj, colour = Status),
    shape = 23,
    fill = "#F59E0B",
    size = 4.6,
    stroke = 1.3
  ) +

  # 6. Threshold boundary lines
  geom_vline(xintercept = c(-1.0, 1.0), linetype = "dashed",
             colour = "#64748B", linewidth = 0.45) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             colour = "#64748B", linewidth = 0.45) +

  # 7. Line threshold indicators placed in clean whitespace
  annotate("label", x = 2.45, y = -log10(0.05), label = "FDR = 0.05 threshold",
           size = 2.4, fill = "white", colour = "#475569", fontface = "italic",
           label.padding = unit(0.15, "lines"), linewidth = 0.25) +
  annotate("label", x = -1.0, y = 0.18, label = "2-fold down",
           size = 2.2, fill = "white", colour = "#475569",
           label.padding = unit(0.12, "lines"), linewidth = 0.2) +
  annotate("label", x = 1.0, y = 0.18, label = "2-fold up",
           size = 2.2, fill = "white", colour = "#475569",
           label.padding = unit(0.12, "lines"), linewidth = 0.2) +

  # 8. Quadrant header callout badges
  annotate("label", x = -2.35, y = y_max * 0.94, label = "Group A Enriched (n = 4)",
           fill = "#DBEAFE", colour = "#1E40AF", fontface = "bold", size = 2.7,
           label.padding = unit(0.24, "lines"), linewidth = 0.3) +
  annotate("label", x = 2.35, y = y_max * 0.94, label = "Group B Enriched (n = 2)",
           fill = "#FEE2E2", colour = "#991B1B", fontface = "bold", size = 2.7,
           label.padding = unit(0.24, "lines"), linewidth = 0.3) +

  # 9. Repelled labels for endogenous hits
  geom_label_repel(
    data = df_sig %>% filter(!grepl("AQUA30", Gene_Clean)),
    aes(label = Label_Text, colour = Status),
    nudge_x = (df_sig %>% filter(!grepl("AQUA30", Gene_Clean)))$nudge_x,
    nudge_y = (df_sig %>% filter(!grepl("AQUA30", Gene_Clean)))$nudge_y,
    size = 2.5,
    fontface = "bold",
    fill = alpha("white", 0.94),
    label.padding = unit(0.20, "lines"),
    box.padding = 0.2,
    point.padding = 0.3,
    min.segment.length = 0,
    segment.colour = "#64748B",
    linewidth = 0.3,
    segment.size = 0.35,
    show.legend = FALSE,
    max.overlaps = 50
  ) +

  # 10. Dedicated callout label for AQUA30 Spike-In Standard
  geom_label_repel(
    data = df_sig %>% filter(grepl("AQUA30", Gene_Clean)),
    aes(label = Label_Text),
    colour = "#78350F",
    fill = "#FEF3C7",
    size = 2.7,
    fontface = "bold",
    nudge_x = (df_sig %>% filter(grepl("AQUA30", Gene_Clean)))$nudge_x,
    nudge_y = (df_sig %>% filter(grepl("AQUA30", Gene_Clean)))$nudge_y,
    label.padding = unit(0.24, "lines"),
    box.padding = 0.2,
    point.padding = 0.4,
    min.segment.length = 0,
    segment.colour = "#B45309",
    linewidth = 0.4,
    segment.size = 0.45,
    show.legend = FALSE
  ) +

  # Aesthetic scales & unified legend
  scale_colour_manual(
    name = NULL,
    values = status_colors,
    labels = c(
      "Up in Group A (n = 4)",
      "Up in Group B (n = 1)",
      "AQUA30 Spike-in (n = 1)",
      "Not Significant (n = 4,511)"
    )
  ) +
  guides(
    colour = guide_legend(
      override.aes = list(
        shape = c(16, 16, 23, 16),
        fill  = c(NA, NA, "#F59E0B", NA),
        colour = c("#1D4ED8", "#DC2626", "#78350F", "#94A3B8"),
        size  = c(3.2, 3.2, 3.8, 1.8),
        alpha = c(1, 1, 1, 0.6)
      ),
      nrow = 1,
      byrow = TRUE
    )
  ) +
  scale_x_continuous(limits = x_lim, breaks = seq(-3, 3, by = 1)) +
  scale_y_continuous(limits = c(0, y_max), breaks = seq(0, 3.5, by = 0.5), expand = c(0, 0)) +

  # Titles and axis labeling
  labs(
    title    = "Differential Protein Abundance — PXD004886 (11-Site DIA)",
    subtitle = "4,517 proteins tested across 11 laboratories | limma paired empirical Bayes | AQUA30 spike-in validated",
    x        = expression(bold(log[2]~"Fold Change (Group B / Group A)")),
    y        = expression(bold(-log[10]~(adjusted~italic(P)~value)))
  ) +

  # Closed 4-spine Nature / Science full-box theme
  theme_classic(base_size = 8.5, base_family = "sans") +
  theme(
    # Closed box border
    panel.border       = element_rect(colour = "#0F172A", fill = NA, linewidth = 0.6),
    axis.line          = element_blank(),
    # Outward ticks
    axis.ticks         = element_line(colour = "#0F172A", linewidth = 0.45),
    axis.ticks.length  = unit(2.5, "pt"),
    # Subtle dashed grid
    panel.grid.major   = element_line(colour = "#E2E8F0", linewidth = 0.35, linetype = "dashed"),
    panel.grid.minor   = element_blank(),
    # Typography
    plot.title         = element_text(size = 9.5, face = "bold", colour = "#0F172A", margin = margin(b = 2)),
    plot.subtitle      = element_text(size = 8.0, colour = "#475569", margin = margin(b = 6)),
    axis.title.x       = element_text(size = 8.5, face = "bold", colour = "#0F172A", margin = margin(t = 4)),
    axis.title.y       = element_text(size = 8.5, face = "bold", colour = "#0F172A", margin = margin(r = 4)),
    axis.text          = element_text(size = 8.0, colour = "#1E293B"),
    # Legend formatting
    legend.position    = "bottom",
    legend.box         = "horizontal",
    legend.margin      = margin(t = 2, b = 2),
    legend.text        = element_text(size = 7.5, colour = "#1E293B"),
    legend.key         = element_blank(),
    legend.background  = element_rect(fill = alpha("white", 0.8), colour = NA),
    plot.margin        = margin(6, 12, 6, 10)
  )

# Export at physical publication dimensions
out_dirs <- c(
  "figures/scifig",
  "data/figures",
  "docs/figures/scifig",
  "docs/data/figures"
)

for (d in out_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# 160 mm width x 100 mm height (ideal for single-page monograph width)
export_targets <- list(
  list(path = "figures/scifig/fig11_volcano_differential", w = 160, h = 100),
  list(path = "data/figures/fig11_volcano_differential", w = 160, h = 100),
  list(path = "data/figures/05_volcano", w = 160, h = 100),
  list(path = "docs/figures/scifig/fig11_volcano_differential", w = 160, h = 100),
  list(path = "docs/data/figures/fig11_volcano_differential", w = 160, h = 100)
)

cat("Exporting AI4SciFig volcano plot across all targets...\n")
for (target in export_targets) {
  stem <- target$path
  w    <- target$w
  h    <- target$h

  # Vector PDF
  ggsave(paste0(stem, ".pdf"), p_volcano, width = w, height = h, units = "mm")
  # 300 DPI PNG
  ggsave(paste0(stem, ".png"), p_volcano, width = w, height = h, units = "mm", dpi = 300)
  cat(sprintf("  ✓ Saved: %s (.pdf + .png)\n", stem))
}

cat("Volcano plot generation completed successfully.\n")
