download_ms_dataset <- function(source = "metabolights",
                                id = NULL,
                                dest_dir = "ms_data",
                                file_pattern = NULL) {
  if (is.null(id) || !nzchar(id)) {
    stop("Please provide a dataset identifier or DOI in `id`.", call. = FALSE)
  }

  dir.create(dest_dir, showWarnings = FALSE, recursive = TRUE)

  local_lib <- file.path(getwd(), "r_libs")
  if (dir.exists(local_lib)) {
    .libPaths(c(local_lib, .libPaths()))
  }

  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }

  if (source == "metabolights") {
    if (!requireNamespace("MsBackendMetaboLights", quietly = TRUE)) {
      BiocManager::install("MsBackendMetaboLights", ask = FALSE)
    }
    if (!requireNamespace("Spectra", quietly = TRUE)) {
      BiocManager::install("Spectra", ask = FALSE)
    }

    sps <- Spectra::Spectra(
      id,
      source = MsBackendMetaboLights::MsBackendMetaboLights()
    )
    return(sps)
  }

  if (source == "massive") {
    dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
    .libPaths(c(local_lib, .libPaths()))

    if ("MsCoreUtils" %in% loadedNamespaces() &&
        packageVersion("MsCoreUtils") < "1.23.9") {
      stop(
        "MsCoreUtils ",
        as.character(packageVersion("MsCoreUtils")),
        " is already loaded, but MsBackendMassIVE requires >= 1.23.9. ",
        "Restart R, then run source('download_ms_dataset.R') before loading ",
        "other mass spectrometry packages.",
        call. = FALSE
      )
    }

    if (!requireNamespace("MsBackendMassIVE", quietly = TRUE)) {
      install.packages(
        "MsBackendMassIVE",
        lib = local_lib,
        repos = c(
          "https://rformassspectrometry.r-universe.dev",
          "https://cloud.r-project.org"
        )
      )
    }
    if (!requireNamespace("Spectra", quietly = TRUE)) {
      BiocManager::install("Spectra", ask = FALSE)
    }

    if (is.null(file_pattern)) {
      sps <- Spectra::Spectra(
        id,
        source = MsBackendMassIVE::MsBackendMassIVE()
      )
    } else {
      sps <- Spectra::Spectra(
        id,
        filePattern = file_pattern,
        source = MsBackendMassIVE::MsBackendMassIVE()
      )
    }
    return(sps)
  }

  if (source == "zenodo") {
    if (!requireNamespace("zen4R", quietly = TRUE)) {
      install.packages("zen4R", repos = "https://cloud.r-project.org")
    }
    if (!requireNamespace("Spectra", quietly = TRUE)) {
      BiocManager::install("Spectra", ask = FALSE)
    }

    zenodo <- zen4R::ZenodoManager$new()
    record <- zenodo$getRecordByDOI(id)
    record$downloadFiles(path = dest_dir)

    # Unzip any archives that may contain mzML/MGF files
    zip_files <- list.files(dest_dir, pattern = "\\.zip$",
                            full.names = TRUE, ignore.case = TRUE)
    for (zf in zip_files) {
      utils::unzip(zf, exdir = dest_dir)
    }

    supported_pattern <- "\\.(mzml|mzxml|mzdata|mgf)$"
    ms_files <- list.files(
      dest_dir,
      pattern = supported_pattern,
      full.names = TRUE,
      recursive = TRUE,
      ignore.case = TRUE
    )

    if (!length(ms_files)) {
      downloaded_files <- list.files(dest_dir,
                                     full.names = FALSE, recursive = TRUE)
      downloaded_summary <- if (length(downloaded_files)) {
        paste(downloaded_files, collapse = ", ")
      } else {
        "none"
      }
      stop(
        "Zenodo record '", id, "' did not contain supported MS files ",
        "(.mzML, .mzXML, .mzData, .mgf). Files in dest_dir: ",
        downloaded_summary, ".",
        call. = FALSE
      )
    }

    Spectra::Spectra(ms_files)
  } else {
    stop("Source not supported. Choose 'metabolights', 'massive', or 'zenodo'.",
         call. = FALSE)
  }
}

# Usage examples:
# metabolomics_data <- download_ms_dataset("metabolights", "MTBLS39")
# massive_data <- download_ms_dataset(
#   "massive", "MSV000080547", file_pattern = "1.mzML$")
# xcms tutorial LC-MS/MS data (mzML files, ~6 MB):
# zenodo_data <- download_ms_dataset(
#   "zenodo", "10.5281/zenodo.3499650")
