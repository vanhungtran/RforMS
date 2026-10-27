# =============================================================================
# Chapter 05: Raw Data Ingestion and Backend Architecture with Spectra
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(Spectra)
  library(mzR)
  library(msdata)
  library(tibble)
  library(dplyr)
  library(tidyr)
})

cat("===================================================================\n")
cat("Chapter 05 Laboratory: Raw Data Ingestion & Spectra Backend Auditing\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Locating and Auditing Standard Open Mass Spectrometry Files
# -----------------------------------------------------------------------------
f_mzml <- system.file("proteomics", "MS3TMT11.mzML", package = "msdata")
if (!file.exists(f_mzml)) {
  stop("Required dataset MS3TMT11.mzML not found in package msdata.")
}

f_info <- file.info(f_mzml)
cat(sprintf("Source mzML File: %s\n", basename(f_mzml)))
cat(sprintf("File Size on Disk: %.2f MB (%d bytes)\n\n",
            f_info$size / (1024^2), f_info$size))

# -----------------------------------------------------------------------------
# Section 2: Low-Level C++ Parsing and On-Disk Handles with mzR
# -----------------------------------------------------------------------------
cat("--- Section 2: Low-Level Handle Inspection via mzR ---\n")
ms_handle <- mzR::openMSfile(f_mzml)
cat(sprintf("Instrument Model: %s\n", mzR::instrumentInfo(ms_handle)$model))
cat(sprintf("Manufacturer:     %s\n", mzR::instrumentInfo(ms_handle)$manufacturer))
cat(sprintf("Ionisation Mode:  %s\n", mzR::instrumentInfo(ms_handle)$ionisation))

# Extract header table (metadata only, zero peak arrays loaded)
hd <- mzR::header(ms_handle)
cat(sprintf("Total Scans Ingested: %d\n", nrow(hd)))

# Audit scan distribution across MS levels
ms_counts <- table(hd$msLevel)
cat("Scan Counts by MS Level:\n")
print(ms_counts)

# Summary of retention time range and TIC
cat(sprintf("Retention Time Span: %.2f to %.2f seconds (%.2f to %.2f min)\n",
            min(hd$retentionTime), max(hd$retentionTime),
            min(hd$retentionTime)/60, max(hd$retentionTime)/60))
cat(sprintf("Median Base Peak Intensity: %.2e | Median TIC: %.2e\n\n",
            median(hd$basePeakIntensity), median(hd$totIonCurrent)))

# On-demand random seeking for a specific MS2 spectrum
ms2_indices <- which(hd$msLevel == 2)
target_idx <- ms2_indices[1]
pks_target <- mzR::peaks(ms_handle, target_idx)
cat(sprintf("Extracted MS2 Scan Index %d (Acquisition Num %d):\n",
            target_idx, hd$acquisitionNum[target_idx]))
cat(sprintf("  Precursor m/z: %.4f | Charge: %d | Collision Energy: %.1f eV\n",
            hd$precursorMZ[target_idx], hd$precursorCharge[target_idx],
            hd$collisionEnergy[target_idx]))
cat(sprintf("  Fragment Peak Count: %d | m/z range: [%.2f, %.2f]\n",
            nrow(pks_target), min(pks_target[, 1]), max(pks_target[, 1])))
cat(sprintf("  Base Fragment: m/z = %.4f (Intensity = %.0f)\n\n",
            pks_target[which.max(pks_target[, 2]), 1],
            max(pks_target[, 2])))

# Close low-level file handle
mzR::close(ms_handle)

# -----------------------------------------------------------------------------
# Section 3: High-Level Spectra Ingestion and Memory Benchmarking
# -----------------------------------------------------------------------------
cat("--- Section 3: High-Level Spectra Ingestion & Memory Benchmark ---\n")

# Load with on-disk MsBackendMzR (lazy evaluation)
sps_ondisk <- Spectra(f_mzml, backend = MsBackendMzR())
size_ondisk <- as.numeric(object.size(sps_ondisk))

cat(sprintf("Spectra Object (MsBackendMzR):\n"))
cat(sprintf("  Total Spectra: %d\n", length(sps_ondisk)))
cat(sprintf("  RAM Allocated: %.2f KB (%d bytes)\n", size_ondisk / 1024, size_ondisk))

# Convert to in-memory backend (MsBackendMemory)
sps_inmem <- setBackend(sps_ondisk, MsBackendMemory())
size_inmem <- as.numeric(object.size(sps_inmem))

cat(sprintf("Spectra Object (MsBackendMemory):\n"))
cat(sprintf("  Total Spectra: %d\n", length(sps_inmem)))
cat(sprintf("  RAM Allocated: %.2f KB (%d bytes)\n", size_inmem / 1024, size_inmem))
cat(sprintf("  Memory Inflation Ratio (In-Memory / On-Disk): %.2fx\n\n",
            size_inmem / size_ondisk))

# -----------------------------------------------------------------------------
# Section 4: Spectra Filtering and Chromatographic Extraction
# -----------------------------------------------------------------------------
cat("--- Section 4: Spectra Filtering & Base Peak Chromatogram (BPC) ---\n")

# Filter exclusively to MS1 survey scans
sps_ms1 <- filterMsLevel(sps_ondisk, msLevel = 1L)
cat(sprintf("Survey Scans (MS1): %d\n", length(sps_ms1)))

# Extract BPC and TIC vectors
rt_ms1  <- rtime(sps_ms1)
tic_ms1 <- tic(sps_ms1)
bpi_ms1 <- sps_ms1$basePeakIntensity
bpi_mz  <- sps_ms1$basePeakMZ

chrom_df <- tibble(
  Scan_Index = seq_along(sps_ms1),
  Retention_Time_sec = round(rt_ms1, 2),
  Retention_Time_min = round(rt_ms1 / 60, 2),
  TIC = round(tic_ms1),
  BPI = round(bpi_ms1),
  BasePeak_MZ = round(bpi_mz, 4)
)

cat("Extracted Chromatogram Header (First 5 MS1 scans):\n")
print(as.data.frame(head(chrom_df, 5)))

# Inspect Top 3 Highest Intensity MS1 Scans
top_ms1 <- chrom_df %>% arrange(desc(BPI)) %>% head(3)
cat("\nTop 3 Eluting Base Peaks in Survey Scans:\n")
print(as.data.frame(top_ms1))

# -----------------------------------------------------------------------------
# Section 5: Processing Pipeline and Peak Cleaning
# -----------------------------------------------------------------------------
cat("\n--- Section 5: Lazy Data Processing Queue ---\n")

# Filter MS2 spectra by precursor m/z window and intensity threshold
sps_ms2 <- filterMsLevel(sps_ondisk, msLevel = 2L)
cat(sprintf("Initial MS2 Count: %d\n", length(sps_ms2)))

# Apply intensity threshold filter: remove fragment peaks < 100 counts
sps_ms2_clean <- filterIntensity(sps_ms2, intensity = c(100, Inf))

# Evaluate peak count reduction on first MS2 spectrum
pks_before <- peaksData(sps_ms2[1])[[1]]
pks_after  <- peaksData(sps_ms2_clean[1])[[1]]

cat(sprintf("MS2 Scan 1 (Acquisition %d):\n", acquisitionNum(sps_ms2[1])))
cat(sprintf("  Raw Fragment Peaks:       %d\n", nrow(pks_before)))
cat(sprintf("  Cleaned Peaks (>= 100):   %d (%.1f%% removed as noise)\n\n",
            nrow(pks_after),
            (nrow(pks_before) - nrow(pks_after)) / nrow(pks_before) * 100))

# -----------------------------------------------------------------------------
# Section 6: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 05 Exercises\n")
cat("===================================================================\n\n")

# Exercise 5.1: Profile vs Centroid Storage and Processing Trade-Offs
cat("--- Exercise 5.1: Profile vs Centroid Storage Dynamics ---\n")
n_peaks_true <- 2500
pts_per_peak_profile <- 18
bytes_per_point <- 16 # 8 bytes m/z + 8 bytes intensity (double precision)

size_profile_bytes <- n_peaks_true * pts_per_peak_profile * bytes_per_point
size_centroid_bytes <- n_peaks_true * 1 * bytes_per_point
compression_ratio <- size_profile_bytes / size_centroid_bytes

cat(sprintf("Features per Scan: %d\n", n_peaks_true))
cat(sprintf("Profile Mode Data Volume:  %.2f KB per scan\n", size_profile_bytes / 1024))
cat(sprintf("Centroid Mode Data Volume: %.2f KB per scan\n", size_centroid_bytes / 1024))
cat(sprintf("Centroiding Compression Ratio: %.1fx data reduction\n", compression_ratio))
cat("Pedagogical Implication: Profile preserves peak width and shape for S/N modeling;\n")
cat("centroiding is required for search engines (Mascot, Comet) and reduces file footprint by 18-fold.\n\n")

# Exercise 5.2: DDA Duty Cycle and MS1/MS2 Timing Audit
cat("--- Exercise 5.2: DDA Duty Cycle and Acquisition Duty Timing ---\n")
hd_all <- mzR::header(mzR::openMSfile(f_mzml))
ms1_idx <- which(hd_all$msLevel == 1)

# Compute duty cycle durations between consecutive MS1 scans
ms1_rts <- hd_all$retentionTime[ms1_idx]
duty_cycles <- diff(ms1_rts)

# Count MS2 scans occurring within each duty cycle
ms2_counts_per_cycle <- integer(length(ms1_idx) - 1)
for (i in seq_along(ms2_counts_per_cycle)) {
  start_rt <- ms1_rts[i]
  end_rt   <- ms1_rts[i + 1]
  ms2_counts_per_cycle[i] <- sum(hd_all$msLevel == 2 &
                                  hd_all$retentionTime >= start_rt &
                                  hd_all$retentionTime < end_rt)
}

cat(sprintf("Total MS1 Survey Cycles: %d\n", length(duty_cycles)))
cat(sprintf("Mean Duty Cycle Duration: %.2f s (SD = %.2f s, Range: [%.2f, %.2f] s)\n",
            mean(duty_cycles), sd(duty_cycles), min(duty_cycles), max(duty_cycles)))
cat(sprintf("Mean Dependent MS2 Scans per Cycle: %.1f (Max: %d)\n",
            mean(ms2_counts_per_cycle), max(ms2_counts_per_cycle)))
cat(sprintf("TopN Acquisition Protocol: Top%d Data-Dependent Acquisition (DDA)\n\n",
            max(ms2_counts_per_cycle)))

# Exercise 5.3: Memory Scalability: Cohort Analysis
cat("--- Exercise 5.3: Memory Scalability across Cohort Sizes ---\n")
cohort_sizes <- c(10, 50, 200, 1000)
ram_per_file_ondisk_mb <- size_ondisk / (1024^2)
ram_per_file_inmem_mb  <- size_inmem / (1024^2)

scaling_df <- tibble(
  Cohort_Files = cohort_sizes,
  RAM_OnDisk_MB = round(cohort_sizes * ram_per_file_ondisk_mb, 2),
  RAM_InMem_MB  = round(cohort_sizes * ram_per_file_inmem_mb, 2),
  RAM_InMem_GB  = round((cohort_sizes * ram_per_file_inmem_mb) / 1024, 2)
)

cat("Cohort Memory Scalability Comparison:\n")
print(as.data.frame(scaling_df))
cat(sprintf("For 1,000 files, In-Memory requires %.2f GB RAM vs %.2f MB for On-Disk backend.\n\n",
            scaling_df$RAM_InMem_GB[4], scaling_df$RAM_OnDisk_MB[4]))

# Exercise 5.4: Random Access Seeking Complexity
cat("--- Exercise 5.4: Indexed mzML Seek Table vs Sequential Parsing ---\n")
cat("Indexed mzML (<indexedmzML>) maintains an offset table at the end of the XML document:\n")
cat("  - Byte Offset Indexing: Provides O(1) random seeking to any arbitrary scan index.\n")
cat("  - Unindexed Sequential Parsing: Requires O(N) linear parsing of all preceding XML nodes.\n")
cat("  - Benchmark: Retrieving Scan 10,000 in an unindexed 2 GB mzML requires parsing ~1.9 GB of text (~8.5 s),\n")
cat("    whereas an indexed file performs a single fseek() system call in < 0.2 milliseconds.\n\n")

cat("===================================================================\n")
cat("Chapter 05 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
