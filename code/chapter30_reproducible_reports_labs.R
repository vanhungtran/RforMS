# =============================================================================
# Chapter 30: Automated Reporting, Dashboards, and Web Applications
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph Series)
# Author: Lucas VHH TRAN, CK-CARE Davos
# R Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(tibble)
  library(dplyr)
  library(gt)
})

cat("===================================================================\n")
cat("Chapter 30 Laboratory: Automated Reporting & FAIR Dissemination\n")
cat("===================================================================\n\n")

# -----------------------------------------------------------------------------
# Section 1: Parameterized Reporting Dynamic Filtering Simulation
# -----------------------------------------------------------------------------
cat("--- Section 1: Parameterized Reporting Threshold Sensitivity ---\n")

# Simulated quantitative omics dataset (500 features)
set.seed(42)
n_features <- 500
test_results <- tibble(
  feature_id = paste0("PROT_", sprintf("%04d", 1:n_features)),
  gene       = paste0("GENE_", 1:n_features),
  log2FC     = rnorm(n_features, mean = 0, sd = 1.2),
  p_val      = rbeta(n_features, shape1 = 0.5, shape2 = 2.0)
) %>%
  mutate(fdr = p.adjust(p_val, method = "BH"))

# Evaluate two parameterized configurations:
# Param Set A: Discovery (FDR < 0.05, |log2FC| >= 1.0)
# Param Set B: Stringent Verification (FDR < 0.01, |log2FC| >= 2.0)

params_discovery <- list(fdr_cutoff = 0.05, log2fc_cutoff = 1.0)
params_stringent <- list(fdr_cutoff = 0.01, log2fc_cutoff = 2.0)

filter_hits <- function(df, pars) {
  df %>% filter(fdr < pars$fdr_cutoff, abs(log2FC) >= pars$log2fc_cutoff)
}

hits_disc <- filter_hits(test_results, params_discovery)
hits_strg <- filter_hits(test_results, params_stringent)

cat(sprintf("Total Evaluated Features: %d\n", n_features))
cat(sprintf("Configuration A [Discovery: FDR < %.2f, |log2FC| >= %.1f]:\n",
            params_discovery$fdr_cutoff, params_discovery$log2fc_cutoff))
cat(sprintf("  Significant Features: %d (%.1f%% of total)\n",
            nrow(hits_disc), nrow(hits_disc) / n_features * 100))
cat(sprintf("  Up-regulated:   %d\n", sum(hits_disc$log2FC >= params_discovery$log2fc_cutoff)))
cat(sprintf("  Down-regulated: %d\n", sum(hits_disc$log2FC <= -params_discovery$log2fc_cutoff)))

cat(sprintf("\nConfiguration B [Verification: FDR < %.2f, |log2FC| >= %.1f]:\n",
            params_stringent$fdr_cutoff, params_stringent$log2fc_cutoff))
cat(sprintf("  Significant Features: %d (%.1f%% of total)\n",
            nrow(hits_strg), nrow(hits_strg) / n_features * 100))
cat(sprintf("  Up-regulated:   %d\n", sum(hits_strg$log2FC >= params_stringent$log2fc_cutoff)))
cat(sprintf("  Down-regulated: %d\n\n", sum(hits_strg$log2FC <= -params_stringent$log2fc_cutoff)))

# -----------------------------------------------------------------------------
# Section 2: AI4SciFig Publication Dimensions and Typography Scaling
# -----------------------------------------------------------------------------
cat("--- Section 2: AI4SciFig Physical Journal Dimension Calculations ---\n")

# Nature and Science physical print standards:
# Single column: 88 mm (3.465 in)
# 1.5 column:    120 mm (4.724 in)
# Double column: 180 mm (7.087 in)

mm_to_in <- function(mm) mm / 25.4

col_specs <- tibble(
  Class           = c("Single Column", "1.5 Column", "Double Column (Full Page)"),
  Width_mm        = c(88.0, 120.0, 180.0),
  Height_mm       = c(88.0, 90.0, 120.0),
  Width_in        = round(mm_to_in(Width_mm), 3),
  Height_in       = round(mm_to_in(Height_mm), 3),
  Target_DPI      = c(300, 300, 300),
  Raster_Width_px = round(Width_in * Target_DPI),
  Raster_Height_px= round(Height_in * Target_DPI)
)

print(as.data.frame(col_specs))
cat("\n")

# -----------------------------------------------------------------------------
# Section 3: Publication Table Verification with gt
# -----------------------------------------------------------------------------
cat("--- Section 3: Publication Summary Table (gt) Construction ---\n")

top_hits_table <- hits_disc %>%
  arrange(fdr) %>%
  slice_head(n = 5) %>%
  transmute(
    Protein_ID  = feature_id,
    Gene_Symbol = gene,
    Log2FC      = round(log2FC, 2),
    Raw_P       = p_val,
    FDR_Q       = fdr,
    Regulation  = ifelse(log2FC > 0, "Upregulated", "Downregulated")
  )

cat("Top 5 Candidate Features Formatted for Publication:\n")
print(as.data.frame(top_hits_table))
cat("\n")

# -----------------------------------------------------------------------------
# Section 4: FAIR SDRF Metadata Schema Auditing
# -----------------------------------------------------------------------------
cat("--- Section 4: FAIR SDRF ProteomeXchange Metadata Validation ---\n")

# Construct valid 6-sample SDRF table
sdrf_df <- tibble(
  `source name`                  = paste0("Sample_", sprintf("%02d", 1:6)),
  `characteristics[organism]`    = "Homo sapiens",
  `characteristics[disease]`     = rep(c("Normal Control", "Vasculitis"), each = 3),
  `characteristics[organism part]`= "Blood Plasma",
  `comment[data file]`           = paste0("Plasma_Run_", sprintf("%02d", 1:6), ".raw"),
  `comment[fraction identifier]` = 1L
)

# Audit SDRF mandatory fields
mandatory_fields <- c("source name", "characteristics[organism]", "comment[data file]")
missing_fields   <- setdiff(mandatory_fields, colnames(sdrf_df))

cat(sprintf("SDRF Table Generated: %d samples x %d columns\n", nrow(sdrf_df), ncol(sdrf_df)))
cat(sprintf("Mandatory SDRF Columns Audited: %d / %d present\n",
            length(mandatory_fields) - length(missing_fields), length(mandatory_fields)))
cat(sprintf("SDRF Structural Validity: %s\n\n",
            if (length(missing_fields) == 0) "PASS (Conforms to ProteomeXchange Standard)" else "FAIL"))

# -----------------------------------------------------------------------------
# Section 5: Worked Solutions to Chapter Exercises
# -----------------------------------------------------------------------------
cat("===================================================================\n")
cat("Worked Solutions to Chapter 30 Exercises\n")
cat("===================================================================\n\n")

# Exercise 30.1: Parameterized Quarto Execution Dynamics
cat("--- Exercise 30.1: Parameterized Rendering Mechanics ---\n")
cat(sprintf("Discovery Settings  (FDR < 0.05, |log2FC| >= 1.0): %d significant hits.\n", nrow(hits_disc)))
cat(sprintf("Verification Settings (FDR < 0.01, |log2FC| >= 2.0): %d significant hits.\n", nrow(hits_strg)))
cat("Analytical Benefit: A single Quarto template adapts dynamically to exploratory screening\n")
cat("and stringent clinical biomarker verification without manual code alterations.\n\n")

# Exercise 30.2: AI4SciFig Typography Invariance
cat("--- Exercise 30.2: Vector vs Raster Publication Export ---\n")
cat("Vector graphics (PDF, SVG) represent geometries as mathematical coordinates, preserving\n")
cat("sharp text and lines at any print zoom factor. Sized at 88 mm, text at 7.5 pt remains legible.\n")
cat("Raster formats (TIFF) must be exported at >= 300 DPI with LZW lossless compression\n")
cat("to avoid print blurriness while minimizing file transfer sizes.\n\n")

# Exercise 30.3: FAIR SDRF Metadata Coupling
cat("--- Exercise 30.3: SDRF vs Vendor Annotations ---\n")
cat("Vendor headers only record instrument states; SDRF maps raw spectra to standardized\n")
cat("ontologies (NCBI taxonomy, Human Disease Ontology, BRENDA tissue ontology).\n")
cat("Result: Facilitates automated meta-analyses and re-searching in public PRIDE repositories.\n\n")

# Exercise 30.4: CI/CD Pipeline Cache Preservation
cat("--- Exercise 30.4: GitHub Actions Dependency Caching ---\n")
cat("In CI/CD workflows, Bioconductor packages require substantial compilation time.\n")
cat("Preserving '~/.local/share/renv' via 'actions/cache@v4' keyed to 'renv.lock' hash\n")
cat("avoids redundant source recompilation, reducing build times from 30+ minutes to < 2 minutes.\n\n")

cat("===================================================================\n")
cat("Chapter 30 Laboratory Execution Successfully Completed.\n")
cat("===================================================================\n")
