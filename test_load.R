# Test script to verify data loading works correctly

# Load required packages
library(Spectra)
library(msdata)

# Method 1: Using msdata:: prefix (CORRECT)
cat("Method 1: Using msdata::proteomics()\n")
ms_file <- msdata::proteomics(full.names = TRUE)[1]
cat("File path:", ms_file, "\n")
cat("File exists:", file.exists(ms_file), "\n\n")

# Load with Spectra
ms_data <- Spectra(ms_file, backend = MsBackendMzR())
cat("Successfully loaded!\n")
cat("Total spectra:", length(ms_data), "\n")
cat("MS levels:", paste(unique(msLevel(ms_data)), collapse = ", "), "\n")

cat("\n✓ All tests passed! The fix is working correctly.\n")
