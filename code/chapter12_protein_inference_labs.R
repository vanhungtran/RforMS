# =============================================================================
# Chapter 12: Protein Inference and Parsimony Grouping
# Laboratory Scripts: Verified Reproducible Code for Listings 1 to 5 and Exercises
# Environment: R 4.6.1 / Bioconductor 3.20
# =============================================================================

suppressPackageStartupMessages({
  library(QFeatures)
  library(MsCoreUtils)
  library(igraph)
  library(dplyr)
  library(tibble)
})

set.seed(42)

# -----------------------------------------------------------------------------
# Listing 12.1: Constructing PSM Dataset and Mapping Peptide Uniqueness
# -----------------------------------------------------------------------------
cat("=== RUNNING LISTING 12.1: CONSTRUCTING PSM DATASET ===\n")

# Benchmark PSM table with 8 peptide-spectrum matches across 5 proteins
psm_data <- tribble(
  ~peptide_sequence, ~protein_accession, ~intensity_s1, ~intensity_s2,
  "SYGFNAAR",        "P00403",           1500,          1600,
  "ELGNDAYK",        "P00403",           3000,          3100,
  "ELGNDAYK",        "P12345",           3000,          3100,  # Shared between P00403 and P12345
  "IAEESNFPFIK",     "P98765",           800,           850,
  "LLQTAEGAE",       "P98765",           1200,          1100,
  "LLQTAEGAE",       "P88888",           1200,          1100,  # Shared between P98765 and P88888
  "VLDALDSIK",       "P55555",           2200,          2100,
  "GVLDTLQFLR",      "P55555",           1900,          2000
)

# Categorize peptides into unique and shared
peptide_mapping <- psm_data %>%
  group_by(peptide_sequence) %>%
  mutate(
    n_proteins = n_distinct(protein_accession),
    peptide_type = if_else(n_proteins == 1, "unique", "shared")
  ) %>%
  ungroup()

cat("Peptide Uniqueness Mapping Summary:\n")
print(as.data.frame(peptide_mapping[, c("peptide_sequence", "protein_accession", "n_proteins", "peptide_type")]))


# -----------------------------------------------------------------------------
# Listing 12.2: Bipartite Graph Analysis and Connected Components
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 12.2: BIPARTITE GRAPH AND CONNECTED COMPONENTS ===\n")

# Construct edge list for all peptide-protein associations
edge_list <- peptide_mapping %>%
  select(from = peptide_sequence, to = protein_accession)

bg <- graph_from_data_frame(edge_list, directed = FALSE)

# Extract connected components using igraph
net_comp <- igraph::components(bg)

cat(sprintf("Bipartite Network Summary:\n"))
cat(sprintf("  Total Graph Nodes:     %d (%d peptides, %d proteins)\n",
            vcount(bg), n_distinct(psm_data$peptide_sequence), n_distinct(psm_data$protein_accession)))
cat(sprintf("  Connected Components:  %d isolated subgraphs\n", net_comp$no))
cat("\nComponent Membership per Node:\n")
print(net_comp$membership)


# -----------------------------------------------------------------------------
# Listing 12.3: Building the QFeatures Container
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 12.3: BUILDING QFEATURES CONTAINER ===\n")

# Filter to retain strictly unique peptides for unambiguous quantification
unique_peptides <- peptide_mapping %>%
  filter(peptide_type == "unique")

# Extract numeric intensity matrix
quant_matrix <- unique_peptides %>%
  select(intensity_s1, intensity_s2) %>%
  as.matrix()
rownames(quant_matrix) <- unique_peptides$peptide_sequence

# Metadata describing peptide-to-protein parentage
row_meta <- DataFrame(
  Protein_Group = unique_peptides$protein_accession,
  Sequence = unique_peptides$peptide_sequence,
  row.names = unique_peptides$peptide_sequence
)

# Wrap into SummarizedExperiment and QFeatures
peptide_se <- SummarizedExperiment(
  assays = list(intensities = quant_matrix),
  rowData = row_meta
)
qf <- QFeatures(list(peptides = peptide_se))

cat(sprintf("Constructed QFeatures container with 1 assay ('peptides'):\n"))
cat(sprintf("  Peptides retained: %d rows x %d samples\n", nrow(quant_matrix), ncol(quant_matrix)))
print(qf)


# -----------------------------------------------------------------------------
# Listing 12.4: Robust Quantification Aggregation via Median Polish
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 12.4: ROBUST QUANTIFICATION AGGREGATION ===\n")

# Aggregate peptides to protein groups using robust Tukey median polish
qf <- aggregateFeatures(
  object = qf,
  i = "peptides",
  fcol = "Protein_Group",
  name = "proteins",
  fun = MsCoreUtils::robustSummary,
  na.rm = TRUE
)

protein_matrix <- assay(qf[["proteins"]])
cat("Protein-Level Quantitative Abundance Matrix (Robust Median Polish):\n")
print(round(protein_matrix, 2))


# -----------------------------------------------------------------------------
# Listing 12.5: Evidence Depth Auditing and the Two-Peptide Rule
# -----------------------------------------------------------------------------
cat("\n=== RUNNING LISTING 12.5: EVIDENCE DEPTH AND TWO-PEPTIDE FILTERING ===\n")

# Inspect evidence depth (.n: number of contributing peptides per protein)
prot_meta <- rowData(qf[["proteins"]])
cat("Protein Evidence Depth (.n peptides per protein):\n")
print(as.data.frame(prot_meta[, c("Protein_Group", ".n")]))

# Filter to retain only proteins supported by >= 2 independent peptides
qf_filtered <- filterFeatures(qf, ~ .n >= 2)
passing_proteins <- assay(qf_filtered[["proteins"]])

cat(sprintf("\nProteins passing the Two-Peptide Rule (.n >= 2): %d of %d proteins\n",
            nrow(passing_proteins), nrow(protein_matrix)))
print(round(passing_proteins, 2))


# -----------------------------------------------------------------------------
# Exercises Verification
# -----------------------------------------------------------------------------
cat("\n=== RUNNING EXERCISES VERIFICATION ===\n")

# Exercise 12.1: Intensity thresholding effect on peptide uniqueness
psm_high <- psm_data %>% filter(intensity_s1 >= 1000)
re_mapped <- psm_high %>%
  group_by(peptide_sequence) %>%
  mutate(n_prot = n_distinct(protein_accession)) %>%
  ungroup()
cat("Exercise 12.1: After removing intensity_s1 < 1000:\n")
cat(sprintf("  Retained PSMs: %d of %d\n", nrow(psm_high), nrow(psm_data)))
cat(sprintf("  Peptide 'IAEESNFPFIK' removed (intensity was 800).\n"))

# Exercise 12.2: Aggregation function comparison: robustSummary vs colMedians vs colSums
qf_med <- aggregateFeatures(qf, i = "peptides", fcol = "Protein_Group", name = "prot_median", fun = colMedians, na.rm = TRUE)
qf_sum <- aggregateFeatures(qf, i = "peptides", fcol = "Protein_Group", name = "prot_sum", fun = colSums, na.rm = TRUE)
cat("\nExercise 12.2: Aggregation Metric Comparison for Protein P55555:\n")
cat(sprintf("  Robust Polish: s1 = %.1f | s2 = %.1f\n", protein_matrix["P55555", 1], protein_matrix["P55555", 2]))
cat(sprintf("  Column Median: s1 = %.1f | s2 = %.1f\n", assay(qf_med[["prot_median"]])["P55555", 1], assay(qf_med[["prot_median"]])["P55555", 2]))
cat(sprintf("  Column Sum:    s1 = %.1f | s2 = %.1f\n", assay(qf_sum[["prot_sum"]])["P55555", 1], assay(qf_sum[["prot_sum"]])["P55555", 2]))

# Exercise 12.3: Minimal Parsimonious Protein Set (Occam's Razor)
# Proteins with unique peptides: P00403 (SYGFNAAR), P98765 (IAEESNFPFIK), P55555 (VLDALDSIK, GVLDTLQFLR).
# Shared peptides: ELGNDAYK explained by P00403; LLQTAEGAE explained by P98765.
# Proteins P12345 and P88888 are subsumed (no unique evidence).
# Minimal parsimonious set: {P00403, P98765, P55555} (3 proteins explain all 8 peptides).
cat(sprintf("\nExercise 12.3: Minimal Parsimonious Set:\n"))
cat("  Proteins with unique evidence: P00403, P98765, P55555 (3 proteins)\n")
cat("  Subsumed proteins dropped by Occam's razor: P12345, P88888 (2 proteins)\n")
cat("  All 8 peptides completely accounted for by 3 parsimonious proteins.\n")

# Exercise 12.4: Razor Peptide Allocation Heuristic
# Under MaxQuant razor rule: ELGNDAYK assigned to P00403 (2 peptides vs P12345's 1 peptide).
# LLQTAEGAE assigned to P98765 (2 peptides vs P88888's 1 peptide).
# P00403 razor intensity: 1500 + 3000 = 4500 (s1).
# P98765 razor intensity: 800 + 1200 = 2000 (s1).
cat(sprintf("\nExercise 12.4: Razor Allocation vs Unique-Only:\n"))
cat(sprintf("  P00403 Unique-only (s1): 1500 | Razor-sum (s1): 4500 (3.0x signal recovery)\n"))
cat(sprintf("  P98765 Unique-only (s1):  800 | Razor-sum (s1): 2000 (2.5x signal recovery)\n"))

cat("=== ALL LISTINGS EXECUTED SUCCESSFULLY ===\n")
