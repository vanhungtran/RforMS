# R for Mass Spectrometry Quick Reference

A condensed reference guide for common operations and workflows.

## Package Loading

```r
# Core MS packages
library(Spectra)           # MS data structures
library(MsCoreUtils)       # MS utilities
library(msdata)            # Example datasets

# Proteomics
library(QFeatures)         # Quantitative features
library(PSMatch)           # PSM handling
library(ProtGenerics)      # Generic functions

# Metabolomics
library(xcms)              # LC-MS processing
library(MetaboCoreUtils)   # Metabolomics tools
library(CAMERA)            # Adduct annotation

# Data manipulation & visualization
library(tidyverse)         # Data science
library(ggplot2)           # Plotting
library(pheatmap)          # Heatmaps
library(limma)             # Statistical models
```

---

## Quick Operations

### Loading MS Data

```r
# From mzML file
library(Spectra)
ms_data <- Spectra("data.mzML", backend = MsBackendMzR())

# Convert to in-memory for speed
ms_data <- setBackend(ms_data, backend = MsBackendDataFrame())

# Load multiple files
files <- c("sample1.mzML", "sample2.mzML", "sample3.mzML")
ms_data <- Spectra(files, backend = MsBackendMzR())
```

### Basic Filtering

```r
# Filter by MS level
ms1 <- filterMsLevel(ms_data, msLevel = 1)
ms2 <- filterMsLevel(ms_data, msLevel = 2)

# Filter by retention time (seconds)
ms_rt <- filterRt(ms_data, rt = c(100, 500))

# Filter by precursor m/z
ms_precursor <- filterPrecursorMz(ms_data, mz = c(400, 800))

# Filter by intensity
ms_intense <- filterIntensity(ms_data, intensity = c(1000, Inf))
```

### Accessing Data

```r
# Get basic properties
length(ms_data)                    # Number of spectra
msLevel(ms_data)                   # MS levels
rtime(ms_data)                     # Retention times
precursorMz(ms_data)              # Precursor m/z
precursorCharge(ms_data)          # Precursor charges

# Get peak data
mz_values <- mz(ms_data)          # List of m/z arrays
int_values <- intensity(ms_data)   # List of intensity arrays
peaks <- peaksData(ms_data)        # List of peak matrices

# Get specific spectrum
spectrum_10 <- ms_data[10]
peaks_10 <- peaksData(spectrum_10)[[1]]
```

---

## Preprocessing

### Peak Picking

```r
# MAD-based peak picking
picked <- pickPeaks(
  ms_data,
  method = "MAD",
  snr = 2,          # Signal-to-noise ratio
  k = 1L            # Half window size
)
```

### Smoothing

```r
# Savitzky-Golay smoothing
smoothed <- smooth(
  ms_data,
  method = "SavitzkyGolay",
  halfWindowSize = 2L
)
```

### Baseline Correction

```r
# Remove baseline noise
corrected <- reduceBaseline(
  ms_data,
  method = "SNIP",
  iterations = 100
)
```

### Normalization

```r
# Normalize intensities
normalized <- normalize(
  ms_data,
  method = "max"  # or "sum", "tic"
)
```

---

## Metabolomics (xcms)

### Complete xcms Workflow

```r
library(xcms)

# 1. Load data
files <- c("control1.mzML", "control2.mzML", 
           "treatment1.mzML", "treatment2.mzML")
pd <- data.frame(
  sample_name = files,
  sample_group = c("Control", "Control", "Treatment", "Treatment")
)
raw_data <- readMSData(files, pdata = pd, mode = "onDisk")

# 2. Peak detection (CentWave)
cwp <- CentWaveParam(
  ppm = 15,
  peakwidth = c(5, 30),
  snthresh = 5,
  prefilter = c(3, 100)
)
xdata <- findChromPeaks(raw_data, param = cwp)

# 3. Retention time alignment
pdp <- PeakDensityParam(sampleGroups = pd$sample_group)
xdata <- adjustRtime(xdata, param = ObiwarpParam())

# 4. Correspondence (peak grouping)
xdata <- groupChromPeaks(xdata, param = pdp)

# 5. Gap filling
xdata <- fillChromPeaks(xdata)

# 6. Extract feature matrix
feature_matrix <- featureValues(xdata, method = "maxint")
```

---

## Proteomics (QFeatures)

### Loading and Processing PSM Data

```r
library(QFeatures)

# Load data from PSM table
psm_file <- "peptides.csv"
qf <- readQFeatures(
  psm_file,
  ecol = 1:10,              # Intensity columns
  fnames = "Sequence",      # Row names
  name = "psms"
)

# Add sample metadata
colData(qf)$condition <- c("Control", "Control", "Control", 
                           "Treatment", "Treatment", "Treatment")

# Filter missing values
qf <- filterNA(qf, i = "psms", pNA = 0.3)  # Max 30% missing

# Log transformation
qf <- logTransform(qf, base = 2, i = "psms", name = "log_psms")

# Normalize
qf <- normalize(qf, i = "log_psms", method = "center.median",
                name = "norm_psms")

# Aggregate to peptides
qf <- aggregateFeatures(qf,
                        i = "norm_psms",
                        fcol = "Sequence",
                        name = "peptides",
                        fun = colMedians)

# Aggregate to proteins
qf <- aggregateFeatures(qf,
                        i = "peptides",
                        fcol = "Protein",
                        name = "proteins",
                        fun = colMedians)
```

---

## Statistical Analysis

### Differential Expression with limma

```r
library(limma)

# Prepare expression matrix (proteins x samples)
expr_matrix <- assay(qf[["proteins"]])

# Design matrix
condition <- factor(colData(qf)$condition)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)

# Contrast matrix
contrast_matrix <- makeContrasts(
  TreatmentVsControl = Treatment - Control,
  levels = design
)

# Fit linear model
fit <- lmFit(expr_matrix, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Get results
results <- topTable(fit2, coef = "TreatmentVsControl", 
                    number = Inf, sort.by = "P")

# Significant features (FDR < 0.05, FC > 2)
sig_features <- results %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1)
```

### PCA Analysis

```r
library(factoextra)

# Perform PCA
pca_result <- prcomp(t(expr_matrix), scale. = TRUE)

# Visualize
fviz_pca_ind(pca_result,
             habillage = condition,
             palette = c("#EE00EE", "#4682B4"),
             addEllipses = TRUE,
             title = "PCA of Protein Expression")
```

### Volcano Plot

```r
library(ggplot2)

ggplot(results, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = adj.P.Val < 0.05 & abs(logFC) > 1),
             alpha = 0.6) +
  scale_color_manual(values = c("grey", "#EE00EE")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10 P-value",
       color = "Significant")
```

### Heatmap

```r
library(pheatmap)

# Select significant features
sig_expr <- expr_matrix[rownames(sig_features), ]

# Create annotation
annotation_col <- data.frame(
  Condition = colData(qf)$condition,
  row.names = colnames(sig_expr)
)

# Plot heatmap
pheatmap(sig_expr,
         scale = "row",
         annotation_col = annotation_col,
         color = colorRampPalette(c("#4682B4", "white", "#EE00EE"))(100),
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         main = "Significant Proteins")
```

---

## Common Parameters

### CentWave (xcms Peak Detection)

```r
CentWaveParam(
  ppm = 15,              # m/z deviation (parts per million)
  peakwidth = c(5, 30),  # Peak width range (seconds)
  snthresh = 5,          # Signal-to-noise threshold
  prefilter = c(3, 100), # Min peaks and intensity
  mzCenterFun = "wMean", # Weighted mean for m/z
  integrate = 1,         # Integration method
  mzdiff = -0.001,       # Min m/z difference
  fitgauss = FALSE,      # Fit Gaussian
  noise = 0              # Noise level
)
```

### Peak Density (xcms Correspondence)

```r
PeakDensityParam(
  sampleGroups = sample_groups,
  bw = 2,                # Bandwidth (seconds)
  minFraction = 0.5,     # Min fraction of samples
  minSamples = 1,        # Min samples per group
  binSize = 0.01,        # m/z bin size
  maxFeatures = 50       # Max features per bin
)
```

---

## Useful Functions

### Spectral Comparison

```r
# Compare two spectra
library(MsCoreUtils)

similarity <- compareSpectra(
  spectrum1,
  spectrum2,
  method = "dotproduct",  # or "cosine", "correlation"
  ppm = 20
)
```

### Mass Calculations

```r
library(MetaboCoreUtils)

# Calculate exact mass
mass <- calculateMass("C6H12O6")  # Glucose: 180.0634

# Convert mass to m/z
mz_pos <- mass2mz(mass, adduct = "[M+H]+")     # 181.0707
mz_neg <- mass2mz(mass, adduct = "[M-H]-")     # 179.0561
mz_sodium <- mass2mz(mass, adduct = "[M+Na]+") # 203.0527

# List available adducts
adductNames()
```

### Peptide Fragments

```r
library(PSMatch)

# Calculate theoretical fragments
fragments <- calculateFragments("PEPTIDE", type = c("b", "y"))

# Match observed spectrum to theoretical
matches <- matchSpectra(
  observed_spectrum,
  theoretical_fragments,
  tolerance = 0.02,
  ppm = 20
)
```

---

## Troubleshooting

### Common Issues

**Issue:** `Error: mzR package not found`
```r
# Solution: Install mzR
BiocManager::install("mzR")
```

**Issue:** `Memory error with large files`
```r
# Solution: Use on-disk backend
ms_data <- Spectra(file, backend = MsBackendMzR())  # Don't convert to DataFrame
```

**Issue:** `No peaks detected`
```r
# Solution: Adjust parameters
picked <- pickPeaks(ms_data, snr = 1, k = 2L)  # Lower thresholds
```

**Issue:** `xcms finds too many/few peaks`
```r
# Too many: Increase stringency
CentWaveParam(snthresh = 10, prefilter = c(5, 500))

# Too few: Decrease stringency
CentWaveParam(snthresh = 3, prefilter = c(3, 100))
```

---

## Performance Tips

1. **Use appropriate backends**
   - Small data (<1 GB): `MsBackendDataFrame` (in-memory, fast)
   - Large data (>1 GB): `MsBackendMzR` (on-disk, memory-efficient)
   - Very large/processed: `MsBackendHdf5Peaks` (balanced)

2. **Parallel processing**
```r
library(BiocParallel)
register(MulticoreParam(workers = 4))  # Use 4 cores
```

3. **Cache results**
```r
# In Quarto/RMarkdown chunks
#| cache: true
```

4. **Filter early**
```r
# Remove unwanted data before processing
ms_data <- ms_data %>%
  filterRt(rt = c(60, 600)) %>%       # Keep 1-10 minutes
  filterMsLevel(msLevel = 2) %>%       # MS2 only
  filterIntensity(intensity = c(1000, Inf))  # High intensity only
```

---

## Resources

- **Documentation:** https://rformassspectrometry.github.io/
- **Book:** https://vanhungtran.github.io/RforMS/
- **Bioconductor:** https://bioconductor.org/
- **GitHub:** https://github.com/vanhungtran/RforMS
- **Community:** Bioconductor support site

---

*For detailed explanations and examples, refer to the full book chapters.*
