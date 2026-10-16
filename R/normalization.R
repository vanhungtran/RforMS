#' Apply median centering normalization across samples
#'
#' Subtracts sample medians and adds global median on log2-transformed data.
#'
#' @param matrix Numeric matrix of log2-transformed intensities (rows = features, cols = samples)
#' @return Normalized numeric matrix
#' @export
normalize_median <- function(matrix) {
  sample_medians <- apply(matrix, 2, stats::median, na.rm = TRUE)
  global_median <- stats::median(sample_medians, na.rm = TRUE)
  sweep(matrix, 2, sample_medians - global_median, "-")
}

#' Apply sample-loading normalization (TMT / isobaric labeling)
#'
#' Scales channel totals to match the mean channel intensity across the plex.
#'
#' @param matrix Numeric matrix of raw reporter ion intensities
#' @return Sample-loading normalized matrix
#' @export
normalize_sample_loading <- function(matrix) {
  channel_sums <- colSums(matrix, na.rm = TRUE)
  target_sum <- mean(channel_sums, na.rm = TRUE)
  factors <- target_sum / channel_sums
  sweep(matrix, 2, factors, "*")
}

#' Diagnose missingness mechanisms (MCAR vs MNAR)
#'
#' Evaluates whether missingness correlates with mean abundance (indicative of
#' left-censored MNAR dropout).
#'
#' @param matrix Numeric matrix of log2-transformed intensities
#' @return A list containing mean abundance per feature, missingness fraction, and correlation
#' @export
diagnose_missingness <- function(matrix) {
  means <- apply(matrix, 1, base::mean, na.rm = TRUE)
  frac_missing <- apply(matrix, 1, function(x) base::mean(is.na(x)))

  # Correlation between abundance and missingness
  valid_idx <- !is.na(means) & frac_missing > 0 & frac_missing < 1
  cor_val <- if (sum(valid_idx) > 5) {
    stats::cor(means[valid_idx], frac_missing[valid_idx], method = "spearman")
  } else {
    NA_real_
  }

  list(
    mean_abundance = means,
    fraction_missing = frac_missing,
    spearman_rho = cor_val,
    mechanism = if (!is.na(cor_val) && cor_val < -0.3) "predominantly MNAR (left-censored)" else "predominantly MAR/MCAR"
  )
}
