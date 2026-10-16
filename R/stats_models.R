#' Format differential abundance results for volcano plotting
#'
#' Classifies features into 'Upregulated', 'Downregulated', or 'Not Significant'
#' based on log2 fold change and adjusted p-value cutoffs.
#'
#' @param df Data frame containing differential abundance results
#' @param log2fc_col Name of the log2 fold-change column
#' @param padj_col Name of the adjusted p-value column
#' @param log2fc_thresh Numeric threshold for log2 fold change (default: 1.0)
#' @param padj_thresh Numeric threshold for adjusted p-value (default: 0.05)
#' @return Data frame with added 'significance' and '-log10(padj)' columns
#' @export
make_volcano_data <- function(df, log2fc_col = "logFC", padj_col = "adj.P.Val",
                              log2fc_thresh = 1.0, padj_thresh = 0.05) {
  res <- df
  lfc <- res[[log2fc_col]]
  padj <- res[[padj_col]]

  res$neg_log10_padj <- -log10(pmax(padj, 1e-300))

  res$significance <- "Not Significant"
  res$significance[lfc >= log2fc_thresh & padj <= padj_thresh] <- "Upregulated"
  res$significance[lfc <= -log2fc_thresh & padj <= padj_thresh] <- "Downregulated"
  res$significance <- factor(res$significance, levels = c("Upregulated", "Downregulated", "Not Significant"))

  res
}

#' Calculate calibration slope and intercept for biomarker models
#'
#' Evaluates whether predicted probabilities agree with actual observed event rates.
#'
#' @param pred_prob Numeric vector of predicted probabilities [0, 1]
#' @param actual Binary numeric or logical vector of true outcomes (0/1)
#' @return A list containing calibration intercept, calibration slope, and Brier score
#' @export
calc_calibration_metrics <- function(pred_prob, actual) {
  actual <- as.numeric(actual)
  # Avoid infinite logit
  p_clamped <- pmin(pmax(pred_prob, 1e-6), 1 - 1e-6)
  logit_p <- stats::qlogis(p_clamped)

  # Logistic regression on logit(p): logit(P(Y=1)) = a + b * logit(p)
  cal_fit <- stats::glm(actual ~ logit_p, family = stats::binomial())
  coefs <- stats::coef(cal_fit)

  # Brier score
  brier <- mean((pred_prob - actual)^2)

  list(
    calibration_intercept = unname(coefs[1]), # Ideal = 0
    calibration_slope = unname(coefs[2]),     # Ideal = 1
    brier_score = brier
  )
}
