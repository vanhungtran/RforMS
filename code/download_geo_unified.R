#' Download GEO and Proteomics Data — Unified Function
#'
#' Single entry point for downloading transcriptomic (GEO) and proteomic (PRIDE/MassIVE)
#' data. Covers all assay types found across the Cardio-CARE server.
#'
#' @param accessions  Character vector of GEO (GSE/GSM) or PRIDE (PXD) accessions.
#' @param out_dir     Base output directory. Default: "data/geo_downloads".
#' @param data_type   Auto-detected if "auto". One of: "auto", "microarray", "bulk_rnaseq",
#'                    "scrna_seq", "proteomics", "metadata_only".
#' @param download    If FALSE, list files only (dry-run). Default: TRUE.
#' @param skip_raw    Skip raw instrument files. Default: TRUE.
#' @param resume      Resume partial downloads. Default: TRUE.
#' @param max_retries Max retries per file. Default: 3.
#' @param timeout_sec Timeout per download request. Default: 600.
#'
#' @return Invisible tibble with columns: accession, data_type, file, size_mb,
#'         status (downloaded/skipped/failed), path.
#'
#' @examples
#' # List what would be downloaded
#' download_geo_unified(c("GSE121212", "GSE193309", "PXD000547"), download = FALSE)
#'
#' # Download single-cell data with supplement files
#' download_geo_unified("GSE184509", data_type = "scrna_seq")
#'
#' # Download proteomics data
#' download_geo_unified("PXD004886", data_type = "proteomics")
#'
#' # Audit metadata only
#' download_geo_unified(c("GSE32924", "GSE130588"), data_type = "metadata_only")
#'
#' # Batch download multiple types (auto-detected)
#' download_geo_unified(c("GSE121212", "GSE184509", "PXD000547", "GSE32924"))
#'
#' @export
download_geo_unified <- function(
    accessions,
    out_dir       = "data/geo_downloads",
    data_type     = c("auto", "microarray", "bulk_rnaseq", "scrna_seq",
                       "proteomics", "metabolomics", "metadata_only"),
    download      = TRUE,
    skip_raw      = TRUE,
    resume        = TRUE,
    max_retries   = 3,
    timeout_sec   = 600
) {
  data_type <- match.arg(data_type)

  # ── Dependencies ──────────────────────────────────────────────────────────
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    stop("Package 'GEOquery' is required. Install: BiocManager::install('GEOquery')")
  }
  pkg_check <- function(pkg) requireNamespace(pkg, quietly = TRUE)

  # ── Constants ─────────────────────────────────────────────────────────────
  RAW_EXTENSIONS <- c(
    "raw", "wiff", "wiff.scan", "d", "mgf", "mzML", "mzXML",
    "fastq", "fastq.gz", "fq", "fq.gz", "bam", "sam", "cram",
    "sra", "srf", "bcl"
  )

  PRIDE_API  <- "https://www.ebi.ac.uk/pride/ws/archive/v1"
  PRIDE_FTP  <- "ftp://ftp.pride.ebi.ac.uk/pride/data/archive"

  # ── Helpers ───────────────────────────────────────────────────────────────

  fmt_bytes <- function(x) {
    if (length(x) == 0 || is.na(x) || is.null(x) || x == 0) return("0 B")
    units <- c("B", "KB", "MB", "GB", "TB")
    exp   <- floor(log(x, 1024))
    exp   <- min(exp, length(units) - 1)
    sprintf("%.1f %s", x / (1024^exp), units[exp + 1])
  }

  is_raw <- function(fname) {
    ext <- tolower(tools::file_ext(fname))
    # Handle compound extensions (.fastq.gz, .wiff.scan)
    base_lower <- tolower(fname)
    any(vapply(RAW_EXTENSIONS, function(pat) endsWith(base_lower, pat), logical(1)))
  }

  classify_accession <- function(acc) {
    if (data_type != "auto") return(data_type)
    # PRIDE / MassIVE proteomics accessions
    if (grepl("^PXD[0-9]{5,7}$", acc, ignore.case = TRUE)) return("proteomics")
    if (grepl("^MSV[0-9]{5,}$", acc, ignore.case = TRUE))  return("proteomics")
    if (grepl("^RPXD[0-9]{5,}$", acc, ignore.case = TRUE)) return("proteomics")
    # MetaboLights accessions
    if (grepl("^MTBLS[0-9]+$", acc, ignore.case = TRUE)) return("metabolomics")
    if (grepl("^ST[0-9]+$", acc, ignore.case = TRUE)) return("metabolomics")
    # GEO accessions
    if (grepl("^GS[EM][0-9]{4,}$", acc, ignore.case = TRUE)) {
      # Further classify by probing GEO metadata
      return(classify_geo_type(acc))
    }
    warning("Unknown accession format: ", acc, ". Treating as GEO bulk RNA-seq.")
    "bulk_rnaseq"
  }

  # ── Safe GEOquery wrapper with timeout ─────────────────────────────────
  safe_getGEO <- function(acc, ...) {
    setTimeLimit(cpu = 120, elapsed = 120, transient = TRUE)
    on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))
    tryCatch(
      GEOquery::getGEO(acc, ...),
      error = function(e) { cat(sprintf("  GEOquery error for %s: %s\n", acc, e$message)); NULL }
    )
  }

  classify_geo_type <- function(acc) {
    tryCatch({
      gset <- safe_getGEO(acc, GSEMatrix = FALSE, getGPL = FALSE)
      if (is.null(gset)) return("bulk_rnaseq")  # Network/parse failure → safe default

      meta <- tryCatch(GEOquery::Meta(gset), error = function(e) list())
      title       <- tolower(paste(meta$title       %||% "", collapse = " "))
      summary     <- tolower(paste(meta$summary     %||% "", collapse = " "))

      # ── scRNA: supplement files + text evidence combined ──
      supp <- tryCatch({
        GEOquery::getGEOSuppFiles(acc, makeDirectory = FALSE, fetch_files = FALSE)
      }, error = function(e) NULL)

      # Tier 1: supplement has h5/h5ad/loom → definitively scRNA
      has_h5_payload <- !is.null(supp) && nrow(supp) > 0 &&
        any(grepl("\\.(h5ad|h5|loom)(\\.gz)?$", supp$fname, ignore.case = TRUE))
      # Tier 2: supplement has RAW.tar (common scRNA packaging)
      has_raw_tar <- !is.null(supp) && nrow(supp) > 0 &&
        any(grepl("_RAW\\.tar(\\.gz)?$", supp$fname, ignore.case = TRUE))
      # Tier 3: no supplement files at all (older datasets)
      no_supp <- is.null(supp) || nrow(supp) == 0

      sc_markers <- c(
        "single[. ]cell", "single-cell", "scrna", "sc-rna", "scrna-seq",
        "scrnaseq", "10x genomics", "droplet.*seq", "smart-seq",
        "celseq", "dropseq", "indrop"
      )
      text_says_sc <- any(sapply(sc_markers, function(m) grepl(m, paste(title, summary))))

      if (has_h5_payload) return("scrna_seq")                                 # definitive
      if (text_says_sc && (has_raw_tar || no_supp)) return("scrna_seq")       # text + data packaging or no supp

      # ── Platform evidence (hard — takes priority over text) ──
      microarray_gpl <- c(
        # Affymetrix
        "GPL96","GPL97","GPL570","GPL571","GPL1261","GPL1638","GPL3921",
        "GPL6244","GPL6480","GPL6848","GPL6947","GPL8300",
        # Illumina BeadChips
        "GPL2700","GPL6102","GPL6104","GPL6883","GPL6884","GPL10558",
        "GPL13393","GPL16288","GPL17021",
        # Agilent
        "GPL4133","GPL4134","GPL6480","GPL6848","GPL7312","GPL10332",
        "GPL13497","GPL13667","GPL14811","GPL14951","GPL15207","GPL15651",
        "GPL16686","GPL16965","GPL17077",
        # Other known microarray
        "GPL10999","GPL17586","GPL17692","GPL23159",
        "GPL4910","GPL4911","GPL4912","GPL4913","GPL4914","GPL4915","GPL4916",
        "GPL9052","GPL9460"
      )
      rnaseq_gpl <- c(
        "GPL16791","GPL11154","GPL18573","GPL24676","GPL21290","GPL20301",
        "GPL17303","GPL17344","GPL20148","GPL23227","GPL23934","GPL25759",
        "GPL34284","GPL19057","GPL21697"
      )

      platform_type <- tryCatch({
        gpls <- tryCatch(GEOquery::GPLList(gset), error = function(e) list())
        if (length(gpls) == 0) return("unknown")
        types <- vapply(names(gpls), function(id) {
          # Known GPL ID
          if (id %in% microarray_gpl) return("microarray")
          if (id %in% rnaseq_gpl) return("rnaseq")
          # GPL title keywords
          gpl_title <- tryCatch(
            tolower(GEOquery::Meta(gpls[[id]])$title %||% ""), error = function(e) "")
          if (grepl("(affymetrix|illumina.*bead|agilent.*array|genechip|microarray|beadchip|bead ?chip)",
                    gpl_title)) return("microarray")
          if (grepl("(illumina.*seq|nextseq|hiseq|novaseq|rna.seq|sequencing by)",
                    gpl_title)) return("rnaseq")
          # Probe ID patterns in platform table
          tbl <- tryCatch(GEOquery::Table(gpls[[id]]), error = function(e) NULL)
          if (!is.null(tbl)) {
            id_col <- intersect(c("ID", "probe_id", "ProbeName"), colnames(tbl))[1]
            if (!is.na(id_col)) {
              sample_ids <- head(as.character(tbl[[id_col]]), 50)
              if (any(grepl("(_at|_s_at|_x_at|_a_at|_g_at)$", sample_ids))) return("microarray")
              if (any(grepl("^ILMN_", sample_ids))) return("microarray")
              if (any(grepl("^ENSG", sample_ids))) return("rnaseq")
            }
          }
          "unknown"
        }, character(1))
        if (any(types == "rnaseq")) "rnaseq"
        else if (any(types == "microarray")) "microarray"
        else "unknown"
      }, error = function(e) "unknown")

      if (platform_type == "microarray") return("microarray")
      if (platform_type == "rnaseq")     return("bulk_rnaseq")

      # ── Text heuristics (soft — only when platform is unknown) ──
      micro_markers <- c(
        "affymetrix", "illumina beadchip", "microarray", "agilent",
        "expression array", "genechip", "oligonucleotide array"
      )
      if (any(sapply(micro_markers, function(m) grepl(m, paste(title, summary))))) {
        return("microarray")
      }

      rnaseq_markers <- c(
        "rna-seq", "rnaseq", "rna sequencing",
        "rna-sequencing", "next generation sequencing", "high throughput sequencing"
      )
      if (any(sapply(rnaseq_markers, function(m) grepl(m, paste(title, summary))))) {
        return("bulk_rnaseq")
      }

      # Truly ambiguous — default to bulk_rnaseq
      "bulk_rnaseq"
    }, error = function(e) "bulk_rnaseq")
  }

  # ── GEO Download Strategies ───────────────────────────────────────────────

  download_geo_metadata_only <- function(acc, dest) {
    cat(sprintf("\n  [METADATA] %s — fetching SOFT...\n", acc))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    for (attempt in seq_len(max_retries)) {
      result <- tryCatch({
        gset <- GEOquery::getGEO(acc, GSEMatrix = FALSE, getGPL = FALSE)
        rds_path <- file.path(dest, paste0(acc, "_soft.rds"))
        saveRDS(gset, rds_path)

        meta <- GEOquery::Meta(gset)
        meta_tbl <- tibble::tibble(
          accession = acc,
          title     = meta$title %||% NA_character_,
          summary   = substr(meta$summary %||% "", 1, 500),
          n_samples = length(GEOquery::GSMList(gset)),
          platforms = paste(names(GEOquery::GPLList(gset)), collapse = "; "),
          type_guess = classify_geo_type(acc)
        )

        readr::write_csv(meta_tbl, file.path(dest, paste0(acc, "_meta.csv")))
        list(status = "downloaded", file = rds_path, size_mb = file.info(rds_path)$size / 1e6)
      }, error = function(e) {
        if (attempt < max_retries) {
          cat(sprintf("    Retry %d/%d: %s\n", attempt, max_retries, e$message))
          Sys.sleep(2^attempt)
        }
        NULL
      })

      if (!is.null(result)) return(result)
    }
    list(status = "failed", file = NA_character_, size_mb = NA_real_)
  }

  download_geo_supplement <- function(acc, dest) {
    cat(sprintf("\n  [SUPPLEMENT] %s — fetching supplementary files...\n", acc))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    results <- list()

    # Step 1: Get supplement file listing
    supp <- tryCatch({
      GEOquery::getGEOSuppFiles(acc, baseDir = dest, makeDirectory = FALSE,
                                 fetch_files = FALSE)
    }, error = function(e) NULL)

    if (is.null(supp) || nrow(supp) == 0) {
      cat("    No supplement files found.\n")
    } else {
      for (i in seq_len(nrow(supp))) {
        fname   <- supp$fname[i]
        furl    <- supp$url[i]
        fsize_raw <- supp$size[i]
        fsize   <- if (length(fsize_raw) > 0) as.numeric(fsize_raw) else NA_real_
        if (length(fsize) == 0) fsize <- NA_real_
        dest_f  <- file.path(dest, fname)

        if (skip_raw && is_raw(fname)) {
          cat(sprintf("    SKIP (raw): %s (%s)\n", fname, fmt_bytes(fsize)))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6, status = "skipped_raw", path = NA_character_
          )
          next
        }

        if (file.exists(dest_f) && resume) {
          local_sz <- file.info(dest_f)$size
          if (abs(local_sz - fsize) < 1000) {
            cat(sprintf("    SKIP (exists): %s\n", fname))
            results[[length(results) + 1]] <- list(
              file = fname, size_mb = fsize / 1e6, status = "skipped_exists",
              path = dest_f
            )
            next
          }
        }

        if (download) {
          cat(sprintf("    DOWNLOAD: %s (%s)\n", fname, fmt_bytes(fsize)))
          ok <- download_with_retry(furl, dest_f, fsize, max_retries, timeout_sec)
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6,
            status = if (ok) "downloaded" else "failed", path = dest_f
          )
        } else {
          cat(sprintf("    WOULD: %s (%s)\n", fname, fmt_bytes(fsize)))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6, status = "dry_run", path = NA_character_
          )
        }
      }
    }

    # Step 2: Try getting series matrix for metadata
    cat("    Fetching series matrix for metadata...\n")
    tryCatch({
      gse <- GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)
      gse <- Filter(function(x) methods::is(x, "ExpressionSet"),
                    if (is.list(gse)) gse else list(gse))
      if (length(gse) > 0) {
        pdata <- Biobase::pData(gse[[1]])
        readr::write_tsv(pdata, file.path(dest, paste0(acc, "_pdata.tsv")))
        results[[length(results) + 1]] <- list(
          file = paste0(acc, "_pdata.tsv"), size_mb = NA_real_,
          status = "downloaded", path = file.path(dest, paste0(acc, "_pdata.tsv"))
        )
      }
    }, error = function(e) {
      cat(sprintf("    Series matrix fallback failed: %s\n", e$message))
    })

    return(results)
  }

  download_geo_scrna_payload <- function(acc, dest) {
    cat(sprintf("\n  [scRNA] %s — detecting payload type...\n", acc))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    results <- list()

    # Step 1: List supplement files
    supp <- tryCatch({
      GEOquery::getGEOSuppFiles(acc, baseDir = dest, makeDirectory = FALSE,
                                 fetch_files = FALSE)
    }, error = function(e) NULL)

    if (is.null(supp) || nrow(supp) == 0) {
      cat("    No supplement files found. Trying supplement download...\n")
      download_geo_supplement(acc, dest)
      # After supplement download, re-scan for scRNA payloads
      supp <- tryCatch({
        GEOquery::getGEOSuppFiles(acc, baseDir = dest, makeDirectory = FALSE,
                                   fetch_files = FALSE)
      }, error = function(e) NULL)
    }

    # Step 2: Detect payload type from supplement listing
    payload_types <- c("\\.h5ad\\.gz$", "\\.h5\\.gz$", "\\.loom\\.gz$",
                       "\\.h5ad$", "\\.h5$", "\\.loom$",
                       "RAW\\.tar$", "_matrix\\.mtx\\.gz$", "_mtx\\.gz$")

    payload_files <- list()
    seen_fnames  <- character()
    if (!is.null(supp) && nrow(supp) > 0) {
      for (ptype in payload_types) {
        idx <- grep(ptype, supp$fname, ignore.case = TRUE, perl = TRUE)
        if (length(idx) > 0) {
          new_idx <- idx[!supp$fname[idx] %in% seen_fnames]
          if (length(new_idx) > 0) {
            payload_files[[ptype]] <- supp[new_idx, ]
            seen_fnames <- c(seen_fnames, supp$fname[new_idx])
          }
        }
      }
    }

    # Step 3: Also check local directory for existing files
    if (dir.exists(dest)) {
      local_files <- list.files(dest, recursive = TRUE, full.names = TRUE)
      for (ptype in c(".h5", ".h5ad", ".loom")) {
        idx <- grep(ptype, tolower(basename(local_files)), fixed = TRUE)
        hits <- local_files[idx]
        for (h in hits) {
          cat(sprintf("    Found local payload: %s\n", basename(h)))
          results[[length(results) + 1]] <- list(
            file = basename(h), size_mb = file.info(h)$size / 1e6,
            status = "exists_local", path = h
          )
        }
      }
    }

    if (length(payload_files) == 0) {
      cat("    No scRNA payload detected (no .h5/.h5ad/.loom/.tar files).\n")
      cat("    Falling back to standard supplement download...\n")
      return(download_geo_supplement(acc, dest))
    }

    # Step 4: Download detected payloads
    for (ptype in names(payload_files)) {
      for (i in seq_len(nrow(payload_files[[ptype]]))) {
        fname <- payload_files[[ptype]]$fname[i]
        furl  <- payload_files[[ptype]]$url[i]
        fsize_raw <- payload_files[[ptype]]$size[i]
        fsize <- if (length(fsize_raw) > 0) as.numeric(fsize_raw) else NA_real_
        if (length(fsize) == 0) fsize <- NA_real_
        dest_f <- file.path(dest, fname)

        if (skip_raw && is_raw(fname) && !grepl("\\.(h5|h5ad|loom)", fname, ignore.case = TRUE)) {
          next
        }

        if (file.exists(dest_f) && resume) {
          local_sz <- file.info(dest_f)$size
          if (abs(local_sz - fsize) < 1000) {
            cat(sprintf("    SKIP (exists): %s\n", fname))
            results[[length(results) + 1]] <- list(
              file = fname, size_mb = fsize / 1e6, status = "skipped_exists", path = dest_f
            )
            next
          }
        }

        if (download) {
          cat(sprintf("    DOWNLOAD: %s (%s)\n", fname, fmt_bytes(fsize)))
          ok <- download_with_retry(furl, dest_f, fsize, max_retries, timeout_sec)
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6,
            status = if (ok) "downloaded" else "failed", path = dest_f
          )
        } else {
          cat(sprintf("    WOULD: %s (%s)\n", fname, fmt_bytes(fsize)))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = if (is.na(fsize)) NA_real_ else fsize / 1e6,
            status = "dry_run", path = NA_character_
          )
        }
      }
    }

    # Step 5: Also get metadata
    cat("    Fetching series matrix for metadata...\n")
    tryCatch({
      gse <- GEOquery::getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)
      gse <- Filter(function(x) methods::is(x, "ExpressionSet"),
                    if (is.list(gse)) gse else list(gse))
      if (length(gse) > 0) {
        pdata <- Biobase::pData(gse[[1]])
        readr::write_tsv(pdata, file.path(dest, paste0(acc, "_pdata.tsv")))
      }
    }, error = function(e) {
      cat(sprintf("    Series matrix metadata failed: %s\n", e$message))
    })

    return(results)
  }

  # ── PRIDE Proteomics Download ─────────────────────────────────────────────

  download_pride_project <- function(acc, dest) {
    cat(sprintf("\n  [PRIDE] %s — proteomics data...\n", acc))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    results <- list()

    # Step 1: Try PRIDE REST API for project info + file listing
    info <- pride_api_request(paste0("projects/", acc))
    if (!is.null(info)) {
      sub_date <- info$submissionDate %||% NA_character_
      pub_date <- info$publicationDate %||% NA_character_
      refs     <- info$references %||% list()
      pmid     <- if (length(refs) > 0) refs[[1]]$pubmedId %||% "?" else "?"
      cat(sprintf("    Submitted: %s | PMID: %s\n", sub_date, pmid))

      # Save project metadata
      meta_tbl <- tibble::tibble(
        accession      = acc,
        title          = info$title %||% NA_character_,
        submission_date = sub_date,
        publication_date = pub_date,
        pmid           = as.character(pmid),
        num_files      = length(info$fileList %||% list())
      )
      readr::write_csv(meta_tbl, file.path(dest, paste0(acc, "_meta.csv")))
    }

    # Step 2: Try PRIDE REST API file listing
    files <- pride_api_request(paste0("files/byProject?accession=", acc))

    if (!is.null(files) && length(files) > 0) {
      # API route — files have downloadLink + fileSize
      for (f in files) {
        fname <- f$fileName
        furl  <- f$downloadLink
        fsize <- as.numeric(f$fileSize %||% 0)
        dest_f <- file.path(dest, basename(fname))

        if (skip_raw && is_raw(fname)) {
          cat(sprintf("    SKIP (raw): %s (%s)\n", fname, fmt_bytes(fsize)))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6, status = "skipped_raw", path = NA_character_
          )
          next
        }

        if (file.exists(dest_f) && resume) {
          local_sz <- file.info(dest_f)$size
          if (abs(local_sz - fsize) < 1000 || fsize == 0) {
            cat(sprintf("    SKIP (exists): %s\n", fname))
            results[[length(results) + 1]] <- list(
              file = fname, size_mb = fsize / 1e6, status = "skipped_exists", path = dest_f
            )
            next
          }
        }

        if (download) {
          cat(sprintf("    DOWNLOAD: %s (%s)\n", fname, fmt_bytes(fsize)))
          ok <- download_with_retry(furl, dest_f, fsize, max_retries, timeout_sec)
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6,
            status = if (ok) "downloaded" else "failed", path = dest_f
          )
        } else {
          cat(sprintf("    WOULD: %s (%s)\n", fname, fmt_bytes(fsize)))
        }
      }
    } else {
      # Step 3: FTP fallback
      cat("    REST API returned no files. Trying FTP...\n")
      info <- pride_api_request(paste0("projects/", acc))
      if (!is.null(info) && !is.na(info$submissionDate %||% NA_character_)) {
        yr <- substr(info$submissionDate, 1, 4)
        mo <- substr(info$submissionDate, 6, 7)

        # Try reported month + adjacent months
        months_to_try <- unique(c(
          mo,
          sprintf("%02d", as.integer(mo) - 1),
          sprintf("%02d", as.integer(mo) + 1)
        ))
        months_to_try <- months_to_try[as.integer(months_to_try) >= 1 &
                                        as.integer(months_to_try) <= 12]

        for (month in months_to_try) {
          ftp_url <- file.path(PRIDE_FTP, yr, month, acc)
          ftp_results <- download_pride_ftp(acc, ftp_url, dest, skip_raw, resume, download)
          if (length(ftp_results) > 0) {
            results <- c(results, ftp_results)
            break
          }
        }
      }
    }

    return(results)
  }

  pride_api_request <- function(endpoint) {
    url <- paste0(PRIDE_API, "/", endpoint)
    req <- curl::curl_fetch_memory(url, handle = curl::new_handle(
      httpheader = c("Accept" = "application/json"),
      timeout_ms = timeout_sec * 1000
    ))

    if (req$status_code >= 200 && req$status_code < 300) {
      raw <- rawToChar(req$content)
      if (nchar(raw) == 0) return(NULL)
      jsonlite::fromJSON(raw, simplifyVector = FALSE)
    } else {
      NULL
    }
  }

  download_pride_ftp <- function(acc, ftp_url, dest, skip_raw, resume, download_mode) {
    results <- list()

    # List FTP directory with curl
    listing <- tryCatch({
      system2("curl", c("-s", "--ftp-pasv", shQuote(paste0(ftp_url, "/"))),
              stdout = TRUE, stderr = FALSE, timeout = 30)
    }, error = function(e) character())

    if (length(listing) == 0) return(results)

    for (line in listing) {
      parts <- strsplit(line, "\\s+")[[1]]
      if (length(parts) < 9) next
      # FTP ls output: perms links user group size month day time name
      fsize <- tryCatch(as.numeric(parts[5]), warning = function(w) NA_real_)
      fname <- paste(parts[9:length(parts)], collapse = " ")

      if (is.na(fname) || fname %in% c(".", "..")) next

      if (skip_raw && is_raw(fname)) {
        cat(sprintf("    SKIP (raw): %s\n", fname))
        results[[length(results) + 1]] <- list(
          file = fname, size_mb = if (is.na(fsize)) NA else fsize / 1e6,
          status = "skipped_raw", path = NA_character_
        )
        next
      }

      dest_f <- file.path(dest, fname)

      if (file.exists(dest_f) && resume && !is.na(fsize)) {
        local_sz <- file.info(dest_f)$size
        if (abs(local_sz - fsize) < 1000) {
          cat(sprintf("    SKIP (exists): %s\n", fname))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = fsize / 1e6, status = "skipped_exists", path = dest_f
          )
          next
        }
      }

      if (download_mode) {
        file_url <- paste0(ftp_url, "/", fname)
        cat(sprintf("    DOWNLOAD: %s\n", fname))
        # Use curl -C - for resume
        args <- if (resume) c("-C", "-", "--ftp-pasv", "-O", file_url)
                else c("--ftp-pasv", "-O", file_url)
        rc <- system2("curl", args, stdout = FALSE, stderr = FALSE,
                      timeout = timeout_sec, wait = TRUE)
        results[[length(results) + 1]] <- list(
          file = fname, size_mb = if (is.na(fsize)) NA else fsize / 1e6,
          status = if (rc == 0) "downloaded" else "failed", path = dest_f
        )
      } else {
        cat(sprintf("    WOULD: %s (%s)\n", fname,
                    if (is.na(fsize)) "?" else fmt_bytes(fsize)))
      }
    }

    return(results)
  }

  # ── MetaboLights Download ─────────────────────────────────────────────────

  MTBLS_API <- "https://www.ebi.ac.uk/metabolights/ws/v1"
  MTBLS_FTP <- "ftp://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"

  download_metabolights_project <- function(acc, dest) {
    cat(sprintf("\n  [MetaboLights] %s — metabolomics data...\n", acc))
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)
    results <- list()

    # Step 1: Query MetaboLights REST API for study info
    study_url <- paste0(MTBLS_API, "/studies/", acc)
    study_info <- tryCatch({
      req <- curl::curl_fetch_memory(study_url, handle = curl::new_handle(
        httpheader = c("Accept" = "application/json"), timeout_ms = 30000))
      if (req$status_code == 200) jsonlite::fromJSON(rawToChar(req$content), simplifyVector = FALSE)
      else NULL
    }, error = function(e) NULL)

    if (!is.null(study_info)) {
      cat(sprintf("    Title: %s\n", substr(study_info$title %||% "?", 1, 120)))
    }

    # Step 2: Try FTP listing for mzML files
    ftp_url <- file.path(MTBLS_FTP, acc)
    cat(sprintf("    FTP: %s\n", ftp_url))

    listing <- tryCatch({
      system2("curl", c("-s", "--ftp-pasv", shQuote(paste0(ftp_url, "/"))),
              stdout = TRUE, stderr = FALSE, timeout = 30)
    }, error = function(e) character())

    if (length(listing) == 0) {
      cat("    No FTP listing available\n")
      return(results)
    }

    # Parse root FTP listing + check FILES/ subdirectory for mzML data
    parse_ftp_listing <- function(lines) {
      files <- list()
      subdirs <- character()
      for (line in lines) {
        parts <- strsplit(line, "\\s+")[[1]]
        if (length(parts) < 9) next
        fsize <- tryCatch(as.numeric(parts[5]), warning = function(w) NA_real_)
        fname <- paste(parts[9:length(parts)], collapse = " ")
        if (is.na(fname) || fname %in% c(".", "..")) next
        # Detect directories
        if (grepl("^d", line)) {
          subdirs <- c(subdirs, fname)
        } else {
          files[[fname]] <- fsize
        }
      }
      list(files = files, subdirs = subdirs)
    }

    root_info <- parse_ftp_listing(listing)
    all_files <- root_info$files

    # Check FILES/ subdirectory for mzML (standard MetaboLights layout)
    for (sub in intersect(root_info$subdirs, c("FILES", "files", "data", "mzML"))) {
      sub_listing <- tryCatch({
        system2("curl", c("-s", "--ftp-pasv", shQuote(paste0(ftp_url, "/", sub, "/"))),
                stdout = TRUE, stderr = FALSE, timeout = 30)
      }, error = function(e) character())
      if (length(sub_listing) > 0) {
        sub_info <- parse_ftp_listing(sub_listing)
        # Prefix with subdirectory path
        for (nm in names(sub_info$files)) {
          all_files[[file.path(sub, nm)]] <- sub_info$files[[nm]]
        }
      }
    }

    # Priority: mzML files > CDF/netCDF > processed tables > metadata > skip raw
    mzml_files  <- names(all_files)[grepl("\\.mzML$", names(all_files), ignore.case = TRUE)]
    cdf_files   <- names(all_files)[grepl("\\.(cdf|nc)$", names(all_files), ignore.case = TRUE)]
    meta_files  <- names(all_files)[grepl("^[isa]_|ISA|metadata|investigation|assay|sample|METADATA|announcement",
                                          names(all_files), ignore.case = TRUE)]
    table_files <- names(all_files)[grepl("\\.(tsv|csv|xlsx)$", names(all_files), ignore.case = TRUE)]
    # Exclude raw vendor files
    raw_pattern <- "\\.(raw|d|wiff|wiff\\.scan|baf|yep|scan)(\\.zip)?$"
    raw_files <- names(all_files)[grepl(raw_pattern, names(all_files), ignore.case = TRUE)]
    table_files <- setdiff(table_files, c(mzml_files, cdf_files, meta_files, raw_files))

    cat(sprintf("    Found: %d mzML, %d CDF, %d metadata, %d tables, %d raw (skip)\n",
                length(mzml_files), length(cdf_files), length(meta_files),
                length(table_files), length(raw_files)))

    # Step 3: Download metadata first (small, always useful)
    for (fname in meta_files) {
      dest_f <- file.path(dest, fname)
      file_url <- paste0(ftp_url, "/", fname)
      if (download) {
        cat(sprintf("    DOWNLOAD: %s\n", fname))
        system2("curl", c("-s", "--ftp-pasv", "-o", shQuote(dest_f), shQuote(file_url)),
                stdout = FALSE, stderr = FALSE, timeout = 60)
        results[[length(results) + 1]] <- list(
          file = fname, size_mb = 0, status = "downloaded", path = dest_f)
      } else {
        cat(sprintf("    WOULD: %s\n", fname))
      }
    }

    # Step 4: Download mzML files (the main data)
    download_files <- c(mzml_files, cdf_files)
    if (length(download_files) == 0 && length(table_files) > 0) {
      download_files <- table_files  # fallback: download processed tables
    }

    for (fname in download_files) {
      fsize <- all_files[[fname]] %||% NA_real_
      dest_f <- file.path(dest, fname)
      file_url <- paste0(ftp_url, "/", fname)

      if (skip_raw && is_raw(fname)) next

      if (file.exists(dest_f) && resume) {
        local_sz <- file.info(dest_f)$size
        if (!is.na(fsize) && abs(local_sz - fsize) < 1000) {
          cat(sprintf("    SKIP (exists): %s\n", fname))
          results[[length(results) + 1]] <- list(
            file = fname, size_mb = if (is.na(fsize)) NA else fsize / 1e6,
            status = "skipped_exists", path = dest_f)
          next
        }
      }

      if (download) {
        cat(sprintf("    DOWNLOAD: %s (%s)\n", fname,
                    if (is.na(fsize)) "?" else fmt_bytes(fsize)))
        system2("curl", c("-C", "-", "--ftp-pasv", "-o", shQuote(dest_f), shQuote(file_url)),
                stdout = FALSE, stderr = FALSE, timeout = 600)
        results[[length(results) + 1]] <- list(
          file = fname, size_mb = if (is.na(fsize)) NA else fsize / 1e6,
          status = "downloaded", path = dest_f)
      } else {
        cat(sprintf("    WOULD: %s (%s)\n", fname,
                    if (is.na(fsize)) "?" else fmt_bytes(fsize)))
        results[[length(results) + 1]] <- list(
          file = fname, size_mb = if (is.na(fsize)) NA_real_ else fsize / 1e6,
          status = "dry_run", path = NA_character_)
      }
    }

    return(results)
  }

  # ── Generic HTTP Download with Retry ──────────────────────────────────────

  download_with_retry <- function(url, dest, expected_size, max_retries, timeout_sec) {
    for (attempt in seq_len(max_retries)) {
      rc <- tryCatch({
        h <- curl::new_handle(timeout_ms = timeout_sec * 1000,
                               followlocation = TRUE)
        curl::curl_download(url, dest, handle = h)
        0
      }, error = function(e) {
        cat(sprintf("      Attempt %d/%d: %s\n", attempt, max_retries, e$message))
        1
      })

      if (rc == 0 && file.exists(dest)) {
        return(TRUE)
      }

      if (attempt < max_retries) {
        Sys.sleep(min(2^attempt, 60))
      }
    }
    return(FALSE)
  }

  # ── Main Loop ─────────────────────────────────────────────────────────────

  `%||%` <- function(x, y) if (is.null(x)) y else x

  if (!requireNamespace("tibble", quietly = TRUE)) {
    install.packages("tibble")
  }
  if (!requireNamespace("readr", quietly = TRUE)) {
    install.packages("readr")
  }

  all_results <- tibble::tibble(
    accession = character(), data_type = character(), file = character(),
    size_mb = numeric(), status = character(), path = character()
  )

  for (acc in accessions) {
    cat(sprintf("\n%s %s %s\n", strrep("=", 60), acc, strrep("=", 40)))

    dtype <- classify_accession(acc)
    dest  <- file.path(out_dir, acc)
    dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    file_results <- switch(
      dtype,
      metadata_only = download_geo_metadata_only(acc, dest),
      microarray    = download_geo_supplement(acc, dest),  # same pattern — supplement + series matrix
      bulk_rnaseq   = download_geo_supplement(acc, dest),
      scrna_seq     = download_geo_scrna_payload(acc, dest),
      proteomics    = download_pride_project(acc, dest),
      metabolomics  = download_metabolights_project(acc, dest),
      download_geo_supplement(acc, dest)  # fallback
    )

    if (is.list(file_results) && !is.data.frame(file_results)) {
      if (length(file_results) > 0) {
        rows <- lapply(file_results, function(r) {
          tibble::tibble(
            accession = acc,
            data_type = dtype,
            file      = r$file %||% NA_character_,
            size_mb   = r$size_mb %||% NA_real_,
            status    = r$status %||% "unknown",
            path      = r$path %||% NA_character_
          )
        })
        all_results <- rbind(all_results, do.call(rbind, rows))
      }
    }

    cat(sprintf("  Done: %s (%s)\n", acc, dtype))
  }

  # ── Summary ───────────────────────────────────────────────────────────────
  cat(sprintf("\n%s SUMMARY %s\n", strrep("=", 30), strrep("=", 30)))
  status_counts <- table(all_results$status)
  for (s in names(status_counts)) {
    cat(sprintf("  %-20s: %d\n", s, status_counts[s]))
  }
  total_size <- sum(all_results$size_mb, na.rm = TRUE)
  cat(sprintf("  Total size: %.1f GB\n", total_size / 1024))

  invisible(all_results)
}
