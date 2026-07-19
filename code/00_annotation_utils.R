# =============================================================================
# annotation_utils.R — Protein annotation for OpenSWATH/DIA data
# Adapted from CardioCare MS pipeline for the R4MS book chapter.
#
# Uses org.Hs.eg.db for local annotation (fast, no network needed).
# Falls back to accession if no gene symbol found.
# =============================================================================

# Ensure the user R library is visible (org.Hs.eg.db lives there)
.libPaths(c("/cluster/home/ltran/R/library", .libPaths()))

# ── Strip OpenSWATH prefix (e.g. "1/P12345", "sp|AQUA30|AQUA30") ─────────────
clean_accession <- function(acc) {
  if (is.na(acc) || !nzchar(acc)) return(NA_character_)
  # Remove "sp|...|" prefix from e.g. "sp|AQUA30|AQUA30"
  acc <- sub("^[a-z]+\\|[^|]+\\|", "", acc)
  # Strip "N/" prefix common in OpenSWATH output
  acc <- sub("^\\d+/", "", acc)
  # Keep only the first accession if semicolon-separated
  acc <- trimws(strsplit(acc, ";")[[1]][1])
  # Remove "Cont_" prefix
  acc <- sub("^Cont_", "", acc)
  acc
}

# ── org.Hs.eg.db annotation ──────────────────────────────────────────────────
annotate_proteins <- function(protein_ids, verbose = TRUE) {
  raw_ids <- protein_ids

  # Clean accessions
  clean_ids <- vapply(raw_ids, clean_accession, character(1))
  unique_clean <- unique(clean_ids[!is.na(clean_ids) & nzchar(clean_ids)])

  if (verbose) cat(sprintf("  Annotating %d unique accessions ...\n", length(unique_clean)))

  # Try org.Hs.eg.db
  has_orgdb <- requireNamespace("AnnotationDbi", quietly = TRUE) &&
               requireNamespace("org.Hs.eg.db", quietly = TRUE)

  if (has_orgdb) {
    if (verbose) cat("  Using org.Hs.eg.db (local) ...\n")
    res <- tryCatch(
      AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = unique_clean,
        columns = c("SYMBOL", "GENENAME"),
        keytype = "UNIPROT"
      ),
      error = function(e) { cat("  org.Hs.eg.db error:", e$message, "\n"); NULL }
    )

    if (!is.null(res) && nrow(res) > 0) {
      # Deduplicate: keep first match per accession
      res <- res[!duplicated(res$UNIPROT), ]
      ann_map <- data.frame(
        Accession = res$UNIPROT,
        Genename = res$SYMBOL,
        Protein.names = res$GENENAME,
        stringsAsFactors = FALSE
      )
    } else {
      ann_map <- NULL
    }
  } else {
    if (verbose) cat("  org.Hs.eg.db not available.\n")
    ann_map <- NULL
  }

  # Build result table
  result <- data.frame(
    Protein = raw_ids,
    Accession = clean_ids,
    stringsAsFactors = FALSE
  )

  if (!is.null(ann_map)) {
    result <- result %>%
      dplyr::left_join(ann_map, by = "Accession")
  } else {
    result$Genename <- NA_character_
    result$Protein.names <- NA_character_
  }

  # Fallback: use accession as gene name where annotation is missing
  result <- result %>%
    dplyr::mutate(
      Genename = dplyr::coalesce(Genename, Accession),
      Protein.names = dplyr::coalesce(Protein.names, Accession)
    )

  n_annotated <- sum(result$Genename != result$Accession, na.rm = TRUE)
  if (verbose) cat(sprintf("  Annotated: %d / %d proteins\n", n_annotated, nrow(result)))

  result
}
