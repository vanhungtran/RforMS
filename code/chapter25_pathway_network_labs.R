# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 25 Laboratory: Pathway Enrichment & Molecular Network Analysis
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(igraph)
  library(ggplot2)
})

message("\n=== Listing 25.1: Over-Representation Analysis with Universe Control ===")
set.seed(42)
ms_universe <- paste0("PROT_", sprintf("%04d", 1:3000))

pathway_db <- list(
  Glycolysis_TCA         = paste0("PROT_", sprintf("%04d", 1:80)),
  Inflammatory_Signaling = paste0("PROT_", sprintf("%04d", 81:120)),
  Olfactory_Receptors    = paste0("PROT_", sprintf("%04d", 121:125))
)

true_hits <- c(
  paste0("PROT_", sprintf("%04d", 1:25)),
  paste0("PROT_", sprintf("%04d", 81:95)),
  sample(ms_universe[126:3000], size = 20)
)

audit_pathway_enrichment <- function(sig_genes, pathway_list, background_universe) {
  N <- length(background_universe)
  n <- length(sig_genes)
  res_list <- list()
  
  for (pname in names(pathway_list)) {
    p_genes <- intersect(pathway_list[[pname]], background_universe)
    M <- length(p_genes)
    overlap_genes <- intersect(sig_genes, p_genes)
    k <- length(overlap_genes)
    
    p_val <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
    odds_ratio <- ((k + 0.5) / (n - k + 0.5)) / ((M - k + 0.5) / (N - M - n + k + 0.5))
    
    res_list[[pname]] <- data.frame(
      Pathway       = pname,
      Universe_Size = M,
      Hit_Count     = k,
      Expected_Hits = round(n * (M / N), 2),
      Odds_Ratio    = round(odds_ratio, 2),
      P_Value       = p_val
    )
  }
  
  out_df <- do.call(rbind, res_list)
  out_df$FDR <- p.adjust(out_df$P_Value, method = "BH")
  return(out_df[order(out_df$P_Value), ])
}

test_res <- audit_pathway_enrichment(true_hits, pathway_db, ms_universe)
cat(sprintf("Audited %d pathways against MS universe of %d proteins.\n", 
            nrow(test_res), length(ms_universe)))
print(test_res[, c("Pathway", "Universe_Size", "Hit_Count", "Expected_Hits", "P_Value", "FDR")])

message("\n=== Exercise 25.1: Exact Hypergeometric Calculation Verification ===")
# N = 4000, M = 100, n = 80, k = 8
N <- 4000
M <- 100
n <- 80
k <- 8

p_exact <- dhyper(k, M, N - M, n)
p_greater_equal <- phyper(k - 1, M, N - M, n, lower.tail = FALSE)
e_k <- n * (M / N)

cat(sprintf("Expected overlap E[K]: %.2f\n", e_k))
cat(sprintf("Exact P(K = 8): %.6f (%.4e)\n", p_exact, p_exact))
cat(sprintf("Enrichment P(K >= 8): %.6f (%.4e)\n", p_greater_equal, p_greater_equal))

message("\n=== Exercise 25.2: Network Centrality Topology in igraph ===")
# Star topology: Node 1 hub connected to nodes 2, 3, 4
adj_mat <- matrix(c(
  0, 1, 1, 1,
  1, 0, 0, 0,
  1, 0, 0, 0,
  1, 0, 0, 0
), nrow = 4, byrow = TRUE)

g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
deg <- degree(g)
bet <- betweenness(g, normalized = FALSE)
clo <- closeness(g)

cat("Graph topology metrics:\n")
cat("Degree centrality:     ", paste(names(deg), deg, sep = "=", collapse = ", "), "\n")
cat("Betweenness centrality:", paste(names(bet), bet, sep = "=", collapse = ", "), "\n")
cat("Closeness centrality:  ", paste(names(clo), round(clo, 3), sep = "=", collapse = ", "), "\n")
message("\n[PASS] Chapter 25 Laboratory successfully executed.")
