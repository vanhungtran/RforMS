#' Clean protein accession string
#'
#' Strips standard database prefixes (sp|, tr|, Cont_, etc.) and isolates the primary accession.
#'
#' @param acc Character accession string
#' @return Cleaned accession string
#' @export
clean_accession <- function(acc) {
  if (is.na(acc) || !nzchar(acc)) return(NA_character_)
  # If format is sp|P12345|ENTRY or tr|Q9Y263|ENTRY, extract the accession
  if (grepl("^[a-z]+\\|[^|]+\\|", acc)) {
    acc <- sub("^[a-z]+\\|([^|]+)\\|.*$", "\\1", acc)
  }
  # Strip "N/" prefix common in OpenSWATH/DIA output
  acc <- sub("^\\d+/", "", acc)
  # Keep only the first accession if semicolon-separated
  acc <- trimws(strsplit(acc, ";")[[1]][1])
  # Remove "Cont_" prefix
  acc <- sub("^Cont_", "", acc)
  acc
}

#' Filter PSM table by target-decoy FDR threshold
#'
#' @param psm_df Data frame containing PSM matches
#' @param score_col Name of the search score column (e.g., "hyperscore", "pep")
#' @param is_decoy_col Logical or character name indicating decoy hits
#' @param fdr_target Target false discovery rate threshold (default: 0.01)
#' @param score_decreasing Logical indicating if higher score means higher confidence
#' @return Filtered data frame of confident target PSMs
#' @export
filter_psm_fdr <- function(psm_df, score_col = "score", is_decoy_col = "is_decoy",
                           fdr_target = 0.01, score_decreasing = TRUE) {
  df <- psm_df[order(psm_df[[score_col]], decreasing = score_decreasing), ]
  is_decoy <- as.logical(df[[is_decoy_col]])

  n_decoy <- cumsum(is_decoy)
  n_target <- cumsum(!is_decoy)

  # Target-decoy FDR formula: FDR = n_decoy / n_target
  q_values <- n_decoy / pmax(n_target, 1)

  # Make monotonic
  q_values <- rev(cummin(rev(q_values)))

  keep_idx <- which(!is_decoy & q_values <= fdr_target)
  df[keep_idx, , drop = FALSE]
}
