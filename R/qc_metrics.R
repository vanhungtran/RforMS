#' Compute Coefficient of Variation (CV) across QC pool samples
#'
#' @param matrix Numeric feature matrix (rows = features, columns = samples)
#' @param qc_cols Character vector or integer indices corresponding to QC pool columns
#' @param na.rm Logical indicating whether to remove NA values (default: TRUE)
#' @return A numeric vector of CV percentages per feature
#' @export
calc_qc_cv <- function(matrix, qc_cols, na.rm = TRUE) {
  if (is.character(qc_cols)) {
    qc_mat <- matrix[, qc_cols, drop = FALSE]
  } else {
    qc_mat <- matrix[, qc_cols, drop = FALSE]
  }

  sds <- apply(qc_mat, 1, stats::sd, na.rm = na.rm)
  means <- apply(qc_mat, 1, base::mean, na.rm = na.rm)

  cv <- (sds / abs(means)) * 100
  cv[!is.finite(cv)] <- NA_real_
  cv
}

#' Audit metadata for batch and condition confounding
#'
#' Evaluates whether the experimental factor of interest is confounded with
#' technical batches or correlated with run order.
#'
#' @param metadata Data frame containing sample annotations
#' @param condition Character name of the condition column
#' @param batch Character name of the batch column
#' @param run_order Optional character name of the run order column
#' @return A list containing the contingency table, chi-square test, and run-order correlation
#' @export
audit_confounding <- function(metadata, condition, batch, run_order = NULL) {
  tab <- table(metadata[[condition]], metadata[[batch]])
  chi <- stats::chisq.test(tab)

  cor_val <- NA_real_
  if (!is.null(run_order) && run_order %in% names(metadata)) {
    cor_val <- stats::cor(
      as.integer(factor(metadata[[condition]])),
      as.numeric(metadata[[run_order]]),
      use = "pairwise.complete.obs"
    )
  }

  list(
    contingency_table = tab,
    chisq_pvalue = chi$p.value,
    run_order_correlation = cor_val,
    confounded = chi$p.value < 0.05 || (!is.na(cor_val) && abs(cor_val) > 0.3)
  )
}
