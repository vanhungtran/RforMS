#' Generate publication-ready Volcano Plot with AI4SciFig styling
#'
#' @param df Data frame containing differential abundance results
#' @param log2fc_col Name of log2 fold-change column (default: "logFC")
#' @param padj_col Name of adjusted p-value column (default: "adj.P.Val")
#' @param label_col Optional column name for feature/gene labels
#' @param n_top_labels Number of top hits to label (default: 8)
#' @param log2fc_thresh Numeric fold-change cutoff (default: 1.0)
#' @param padj_thresh Numeric FDR cutoff (default: 0.05)
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_volcano <- function(df, log2fc_col = "logFC", padj_col = "adj.P.Val",
                           label_col = NULL, n_top_labels = 8,
                           log2fc_thresh = 1.0, padj_thresh = 0.05,
                           title = "Differential Abundance") {
  plot_df <- make_volcano_data(df, log2fc_col, padj_col, log2fc_thresh, padj_thresh)

  cols <- c(
    "Upregulated"   = "#0F4D92", # Deep Nature Blue
    "Downregulated" = "#D95F02", # Deep Vermilion
    "Not Significant" = "#B0B0B0" # Muted Gray
  )

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data[[log2fc_col]], y = .data[["neg_log10_padj"]])) +
    ggplot2::geom_point(ggplot2::aes(colour = .data[["significance"]]), alpha = 0.7, size = 1.2) +
    ggplot2::scale_colour_manual(values = cols) +
    ggplot2::geom_vline(xintercept = c(-log2fc_thresh, log2fc_thresh), linetype = "dashed", colour = "#666666", linewidth = 0.3) +
    ggplot2::geom_hline(yintercept = -log10(padj_thresh), linetype = "dashed", colour = "#666666", linewidth = 0.3) +
    ggplot2::labs(
      title = title,
      x = expression(Log[2]~Fold~Change),
      y = expression(-Log[10]~Adjusted~italic(P)-value),
      colour = "Regulation"
    ) +
    theme_scifigure(base_size = 7, grid = "both")

  if (!is.null(label_col) && label_col %in% names(plot_df) && n_top_labels > 0) {
    top_hits <- plot_df[plot_df$significance != "Not Significant", ]
    top_hits <- top_hits[order(top_hits[[padj_col]]), ]
    top_hits <- utils::head(top_hits, n_top_labels)

    if (nrow(top_hits) > 0) {
      p <- p + ggplot2::geom_text(
        data = top_hits,
        ggplot2::aes(label = .data[[label_col]]),
        size = 2.2, fontface = "bold", vjust = -0.5, check_overlap = TRUE
      )
    }
  }

  p
}

#' Generate publication-ready PCA score plot with AI4SciFig styling
#'
#' @param matrix Numeric matrix (rows = features, columns = samples)
#' @param metadata Data frame containing sample metadata
#' @param color_by Column name in metadata to color points by
#' @param shape_by Optional column name in metadata for point shapes
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_pca <- function(matrix, metadata, color_by, shape_by = NULL, title = "Principal Component Analysis") {
  complete_mat <- stats::na.omit(matrix)
  pca_res <- stats::prcomp(t(complete_mat), scale. = TRUE)

  var_explained <- (pca_res$sdev^2) / sum(pca_res$sdev^2) * 100
  pc1_lab <- sprintf("PC1 (%.1f%%)", var_explained[1])
  pc2_lab <- sprintf("PC2 (%.1f%%)", var_explained[2])

  df_plot <- as.data.frame(pca_res$x[, 1:2])
  df_plot[[color_by]] <- metadata[[color_by]]

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[["PC1"]], y = .data[["PC2"]]))

  if (!is.null(shape_by) && shape_by %in% names(metadata)) {
    df_plot[[shape_by]] <- metadata[[shape_by]]
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data[[color_by]], shape = .data[[shape_by]]), size = 2)
  } else {
    p <- p + ggplot2::geom_point(ggplot2::aes(colour = .data[[color_by]]), size = 2)
  }

  p +
    scale_color_scifigure(palette = "nature") +
    ggplot2::labs(title = title, x = pc1_lab, y = pc2_lab) +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate missingness diagnostic plot (MCAR vs MNAR)
#'
#' @param matrix Numeric matrix of log2-transformed intensities
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_missingness <- function(matrix, title = "Missingness vs Abundance") {
  diag <- diagnose_missingness(matrix)
  df_plot <- data.frame(
    mean_intensity = diag$mean_abundance,
    fraction_missing = diag$fraction_missing
  )
  df_plot <- df_plot[!is.na(df_plot$mean_intensity), ]

  ggplot2::ggplot(df_plot, ggplot2::aes(x = .data[["mean_intensity"]], y = .data[["fraction_missing"]])) +
    ggplot2::geom_point(alpha = 0.5, size = 1, colour = "#0F4D92") +
    ggplot2::geom_smooth(method = "loess", colour = "#D95F02", linewidth = 0.7, se = FALSE) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Mechanism:", diag$mechanism),
      x = expression(Mean~Log[2]~Intensity),
      y = "Fraction Missing"
    ) +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate publication-ready Biomarker Calibration Plot
#'
#' @param pred_prob Numeric predicted probabilities [0, 1]
#' @param actual True binary outcomes (0/1)
#' @param n_bins Number of probability bins (default: 8)
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_calibration <- function(pred_prob, actual, n_bins = 8, title = "Biomarker Model Calibration") {
  cal_stats <- calc_calibration_metrics(pred_prob, actual)

  bins <- cut(pred_prob, breaks = seq(0, 1, length.out = n_bins + 1), include.lowest = TRUE)
  bin_df <- data.frame(pred = pred_prob, actual = as.numeric(actual), bin = bins)

  agg <- stats::aggregate(cbind(pred, actual) ~ bin, data = bin_df, FUN = mean)

  sub_lab <- sprintf("Slope = %.2f | Intercept = %.2f | Brier = %.3f",
                     cal_stats$calibration_slope, cal_stats$calibration_intercept, cal_stats$brier_score)

  ggplot2::ggplot(agg, ggplot2::aes(x = .data[["pred"]], y = .data[["actual"]])) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#888888", linewidth = 0.4) +
    ggplot2::geom_line(colour = "#0F4D92", linewidth = 0.6) +
    ggplot2::geom_point(colour = "#0F4D92", size = 2) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = title,
      subtitle = sub_lab,
      x = "Predicted Probability",
      y = "Observed Proportion"
    ) +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate publication-ready Density Distribution Plot
#'
#' @param matrix Numeric matrix of log2 intensities (features in rows, samples in cols)
#' @param groups Optional vector of sample groups (length = ncol(matrix))
#' @param title Plot title
#' @param xlab X-axis label
#' @return A ggplot object
#' @export
scifig_density <- function(matrix, groups = NULL, title = "Intensity Density Distribution",
                           xlab = expression(Log[2]~Intensity)) {
  sample_names <- colnames(matrix)
  if (is.null(sample_names)) sample_names <- paste0("Sample_", seq_len(ncol(matrix)))

  df_list <- lapply(seq_len(ncol(matrix)), function(j) {
    vals <- stats::na.omit(matrix[, j])
    data.frame(
      Sample = sample_names[j],
      Group = if (!is.null(groups)) as.character(groups[j]) else "Sample",
      Intensity = vals,
      stringsAsFactors = FALSE
    )
  })
  df_long <- do.call(rbind, df_list)

  p <- ggplot2::ggplot(df_long, ggplot2::aes(x = .data[["Intensity"]], group = .data[["Sample"]]))

  if (!is.null(groups)) {
    p <- p + ggplot2::geom_density(ggplot2::aes(colour = .data[["Group"]]), alpha = 0.4, linewidth = 0.35) +
      scale_color_scifigure(palette = "nature") +
      ggplot2::labs(colour = "Group")
  } else {
    p <- p + ggplot2::geom_density(colour = "#0F4D92", alpha = 0.4, linewidth = 0.35)
  }

  p +
    ggplot2::labs(title = title, x = xlab, y = "Density") +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate publication-ready Chromatographic Peak Plot (EIC)
#'
#' @param rtime Numeric vector of retention times
#' @param intensity Numeric vector of chromatographic intensities
#' @param rt_min Optional numeric retention time start of integrated peak
#' @param rt_max Optional numeric retention time end of integrated peak
#' @param apex_rt Optional numeric retention time of peak apex
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_eic <- function(rtime, intensity, rt_min = NULL, rt_max = NULL, apex_rt = NULL,
                       title = "Extracted Ion Chromatogram") {
  df <- data.frame(RT = rtime, Intensity = intensity)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["RT"]], y = .data[["Intensity"]])) +
    ggplot2::geom_area(fill = "#0F4D92", alpha = 0.15) +
    ggplot2::geom_line(colour = "#0F4D92", linewidth = 0.5)

  if (!is.null(rt_min) && !is.null(rt_max)) {
    peak_sub <- df[df$RT >= rt_min & df$RT <= rt_max, ]
    max_peak_y <- if (nrow(peak_sub) > 0) max(peak_sub$Intensity) else max(df$Intensity)

    p <- p + ggplot2::annotate(
      "rect", xmin = rt_min, xmax = rt_max, ymin = 0, ymax = max_peak_y * 1.05,
      fill = NA, colour = "#D95F02", linetype = "dashed", linewidth = 0.4
    )
  }

  if (!is.null(apex_rt)) {
    apex_int <- df$Intensity[which.min(abs(df$RT - apex_rt))]
    p <- p + ggplot2::annotate("point", x = apex_rt, y = apex_int, colour = "#D95F02", size = 2)
  }

  p +
    ggplot2::labs(title = title, x = "Retention Time (s)", y = "Intensity (counts)") +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate Target-Decoy Score Distribution Plot
#'
#' @param scores Numeric vector of search scores (e.g., SEQUEST XCorr, Mascot, or Andromeda score)
#' @param is_decoy Logical vector indicating decoy PSMs (TRUE = decoy, FALSE = target)
#' @param score_thresh Optional numeric score threshold corresponding to 1% FDR
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_target_decoy <- function(scores, is_decoy, score_thresh = NULL,
                                title = "Target-Decoy PSM Score Distribution") {
  df <- data.frame(
    Score = scores,
    Type = ifelse(is_decoy, "Decoy", "Target")
  )

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Score"]], fill = .data[["Type"]], colour = .data[["Type"]])) +
    ggplot2::geom_density(alpha = 0.35, linewidth = 0.4) +
    ggplot2::scale_fill_manual(values = c("Target" = "#0F4D92", "Decoy" = "#D95F02")) +
    ggplot2::scale_colour_manual(values = c("Target" = "#0F4D92", "Decoy" = "#D95F02"))

  if (!is.null(score_thresh)) {
    p <- p +
      ggplot2::geom_vline(xintercept = score_thresh, linetype = "dashed", colour = "#333333", linewidth = 0.4) +
      ggplot2::annotate("text", x = score_thresh, y = 0, label = " 1% FDR Cutoff",
                        hjust = 0, vjust = -1, size = 2.2, fontface = "bold")
  }

  p +
    ggplot2::labs(title = title, x = "Identification Score", y = "Density", fill = "Database", colour = "Database") +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate Targeted Quantification Calibration Dilution Curve
#'
#' @param concentration Numeric vector of nominal analyte concentrations
#' @param signal Numeric vector of peak area ratios (analyte / internal standard)
#' @param lod Optional numeric limit of detection
#' @param loq Optional numeric limit of quantitation
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_calibration_curve <- function(concentration, signal, lod = NULL, loq = NULL,
                                     title = "PRM / MRM Calibration Curve") {
  df <- data.frame(Concentration = concentration, Signal = signal)
  fit <- stats::lm(Signal ~ Concentration, data = df)
  r2 <- summary(fit)$r.squared
  slope <- coef(fit)[2]

  sub_lab <- sprintf("Linear fit: slope = %.3e | R\u00b2 = %.4f", slope, r2)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[["Concentration"]], y = .data[["Signal"]])) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, colour = "#0F4D92", linewidth = 0.5, se = FALSE) +
    ggplot2::geom_point(colour = "#0F4D92", size = 1.8)

  if (!is.null(lod)) {
    p <- p + ggplot2::geom_vline(xintercept = lod, linetype = "dotted", colour = "#D95F02", linewidth = 0.4) +
      ggplot2::annotate("text", x = lod, y = max(signal) * 0.9, label = "LOD", hjust = -0.1, size = 2, colour = "#D95F02")
  }
  if (!is.null(loq)) {
    p <- p + ggplot2::geom_vline(xintercept = loq, linetype = "dashed", colour = "#D95F02", linewidth = 0.4) +
      ggplot2::annotate("text", x = loq, y = max(signal) * 0.8, label = "LOQ", hjust = -0.1, size = 2, colour = "#D95F02")
  }

  p +
    ggplot2::labs(title = title, subtitle = sub_lab, x = "Nominal Concentration (fmol/\u00b5L)", y = "Peak Area Ratio (PAR)") +
    theme_scifigure(base_size = 7, grid = "both")
}

#' Generate Multi-Cohort ROC Curve Plot
#'
#' @param roc_list Named list of data frames, each with columns 'fpr' and 'tpr'
#' @param auc_values Named numeric vector of AUC values corresponding to each cohort
#' @param title Plot title
#' @return A ggplot object
#' @export
scifig_roc <- function(roc_list, auc_values = NULL, title = "Multi-Cohort Validation ROC") {
  df_list <- lapply(names(roc_list), function(nm) {
    df <- roc_list[[nm]]
    df$Cohort <- nm
    if (!is.null(auc_values) && nm %in% names(auc_values)) {
      df$CohortLabel <- sprintf("%s (AUC = %.3f)", nm, auc_values[[nm]])
    } else {
      df$CohortLabel <- nm
    }
    df
  })
  df_all <- do.call(rbind, df_list)

  ggplot2::ggplot(df_all, ggplot2::aes(x = .data[["fpr"]], y = .data[["tpr"]], colour = .data[["CohortLabel"]])) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", colour = "#888888", linewidth = 0.35) +
    ggplot2::geom_line(linewidth = 0.6) +
    scale_color_scifigure(palette = "nature") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = title,
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      colour = "Cohort"
    ) +
    theme_scifigure(base_size = 7, grid = "both")
}

