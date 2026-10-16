testthat::test_that("adduct math functions operate correctly", {
  glu_mass <- 180.063388
  adducts <- calc_adduct_mz(glu_mass, c("[M+H]+", "[M+Na]+"))

  testthat::expect_equal(round(adducts[["[M+H]+"]], 4), 181.0707)
  testthat::expect_equal(round(adducts[["[M+Na]+"]], 4), 203.0526)

  inferred <- mz_to_neutral_mass(adducts[["[M+Na]+"]], "[M+Na]+")
  testthat::expect_equal(round(inferred, 4), round(glu_mass, 4))
})

testthat::test_that("accession cleaning strips prefixes reliably", {
  testthat::expect_equal(clean_accession("sp|P12345|HBA_HUMAN"), "P12345")
  testthat::expect_equal(clean_accession("tr|Q9Y263|Q9Y263_HUMAN"), "Q9Y263")
  testthat::expect_equal(clean_accession("Cont_P02768"), "P02768")
  testthat::expect_equal(clean_accession("1/P04637;P04638"), "P04637")
})

testthat::test_that("median normalization centers correctly", {
  mat <- matrix(c(10, 12, 11, 20, 22, 21), nrow = 3, ncol = 2)
  norm_mat <- normalize_median(mat)
  medians <- apply(norm_mat, 2, stats::median)
  testthat::expect_equal(medians[1], medians[2])
})

testthat::test_that("volcano data classification functions properly", {
  df <- data.frame(
    logFC = c(2.0, -2.5, 0.1),
    adj.P.Val = c(0.001, 0.0001, 0.5)
  )
  res <- make_volcano_data(df, log2fc_thresh = 1.0, padj_thresh = 0.05)
  testthat::expect_equal(as.character(res$significance), c("Upregulated", "Downregulated", "Not Significant"))
})
