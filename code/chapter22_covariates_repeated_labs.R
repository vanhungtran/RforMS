# =============================================================================
# Mass Spectrometry Data Analysis with R (Springer Nature Monograph)
# Chapter 22 Laboratory: Covariates, Longitudinal Profiling, and Repeated Measures
# Author: Lucas VHH TRAN (CK-CARE Davos)
# Environment: R 4.6.1 / Bioconductor 3.23 (Windows x86_64)
# =============================================================================

suppressPackageStartupMessages({
  library(MSstats)
  library(variancePartition)
  library(limma)
  library(splines)
  library(reshape2)
  library(lme4)
  library(nlme)
  library(BiocParallel)
})

# Register SerialParam for deterministic single-threaded execution on Windows
register(SerialParam(), default = TRUE)

message("\n=== Listing 22.1: Henderson's Mixed Model Equations in Base R ===")
solve_henderson <- function(X, Z, y, sigma_b, sigma_e) {
  n <- length(y); q <- ncol(Z)
  R_inv <- diag(1 / sigma_e^2, n)
  G_inv <- diag(1 / sigma_b^2, q)
  
  C11 <- t(X) %*% R_inv %*% X
  C12 <- t(X) %*% R_inv %*% Z
  C21 <- t(Z) %*% R_inv %*% X
  C22 <- t(Z) %*% R_inv %*% Z + G_inv
  
  C   <- rbind(cbind(C11, C12), cbind(C21, C22))
  rhs <- c(t(X) %*% R_inv %*% y, t(Z) %*% R_inv %*% y)
  sol <- solve(C, rhs)
  
  list(beta = sol[1:ncol(X)], u = sol[(ncol(X) + 1):length(sol)])
}

# Quick simulation verification
set.seed(42)
X_sim <- cbind(1, rnorm(20))
Z_sim <- matrix(0, nrow = 20, ncol = 4)
for(i in 1:4) Z_sim[((i-1)*5 + 1):(i*5), i] <- 1
y_sim <- X_sim %*% c(10, 2) + Z_sim %*% rnorm(4, 0, 1.5) + rnorm(20, 0, 0.5)
h_sol <- solve_henderson(X_sim, Z_sim, y_sim, sigma_b = 1.5, sigma_e = 0.5)
cat(sprintf("Henderson MME solved: beta0 = %.2f, beta1 = %.2f\n", h_sol$beta[1], h_sol$beta[2]))

message("\n=== Listing 22.2: Empirical Verification: REML vs ML Variance Estimates in lme4 ===")
data(SRMRawData)
proc <- dataProcess(SRMRawData, use_log_file = FALSE, verbose = FALSE)
df_idhc <- proc$ProteinLevelData[proc$ProteinLevelData$Protein == "IDHC", ]
df_idhc$Time <- as.numeric(as.character(df_idhc$GROUP))

fit_reml <- lmer(LogIntensities ~ Time + (1 | SUBJECT), data = df_idhc, REML = TRUE)
fit_ml   <- lmer(LogIntensities ~ Time + (1 | SUBJECT), data = df_idhc, REML = FALSE)

cat(sprintf("REML sigma: %.4f | ML sigma: %.4f | Ratio: %.6f\n",
            sigma(fit_reml), sigma(fit_ml), sigma(fit_reml) / sigma(fit_ml)))
cat(sprintf("Theoretical bias adjustment factor sqrt(n / (n - p)): %.6f\n",
            sqrt(30 / (30 - 2))))

message("\n=== Listing 22.3: Hands-On Comparison: CS vs AR(1) Covariance Structures in nlme ===")
m_cs  <- gls(LogIntensities ~ Time, correlation = corCompSymm(form = ~ 1 | SUBJECT), data = df_idhc)
m_ar1 <- gls(LogIntensities ~ Time, correlation = corAR1(form = ~ Time | SUBJECT), data = df_idhc)

cat(sprintf("Compound Symmetry:   AIC = %.2f | Estimated Rho = %.4f\n",
            AIC(m_cs), coef(m_cs$modelStruct$corStruct, unconstrained = FALSE)))
cat(sprintf("Autoregressive AR(1): AIC = %.2f | Estimated Phi = %.4f\n",
            AIC(m_ar1), coef(m_ar1$modelStruct$corStruct, unconstrained = FALSE)))

message("\n=== Listing 22.4: Loading and Preprocessing Longitudinal SRM Proteomics Data ===")
df_prot <- proc$ProteinLevelData
df_wide <- reshape2::dcast(df_prot, Protein ~ SUBJECT + GROUP, value.var = "LogIntensities")
expr_mat <- as.matrix(df_wide[, -1])
rownames(expr_mat) <- df_wide$Protein

meta_df <- data.frame(
  Sample   = colnames(expr_mat),
  SUBJECT  = factor(sub("_.*", "", colnames(expr_mat))),
  Time_fac = factor(sub(".*_", "", colnames(expr_mat)), levels = as.character(1:10)),
  Time_num = as.numeric(sub(".*_", "", colnames(expr_mat)))
)
rownames(meta_df) <- colnames(expr_mat)

cat(sprintf("Loaded real SRM cohort: %d proteins across %d runs (%d timepoints x %d replicates).\n",
            nrow(expr_mat), ncol(expr_mat), length(unique(meta_df$Time_num)), length(unique(meta_df$SUBJECT))))

message("\n=== Listing 22.5: Automated dream Pipeline Execution with Kenward-Roger DDF ===")
form <- ~ 0 + Time_fac + (1 | SUBJECT)

L <- makeContrastsDream(form, meta_df, 
                        contrasts = c(Diauxic_T7_vs_T1     = "Time_fac7 - Time_fac1",
                                      Stationary_T10_vs_T1 = "Time_fac10 - Time_fac1"))

cat("Fitting dream linear mixed model with Kenward-Roger DDF...\n")
fit_mm <- dream(expr_mat, form, meta_df, L = L, ddf = "Kenward-Roger", BPPARAM = SerialParam())
fit_eb <- eBayes(fit_mm)

res_diauxic   <- variancePartition::topTable(fit_eb, coef = "Diauxic_T7_vs_T1", number = Inf)
res_stationary <- variancePartition::topTable(fit_eb, coef = "Stationary_T10_vs_T1", number = Inf)

cat("\nDifferential Abundance: Diauxic Shift (T7 vs T1):\n")
print(head(round(res_diauxic[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 4), 3))

cat("\nDifferential Abundance: Stationary Phase (T10 vs T1):\n")
print(head(round(res_stationary[, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")], 4), 3))

message("\n=== Listing 22.6: Variance Partitioning on Real Mass Spectrometry Time-Course ===")
form_vp <- ~ (1 | Time_fac) + (1 | SUBJECT)
varPart <- fitExtractVarPartModel(expr_mat, form_vp, meta_df, BPPARAM = SerialParam())

cat("Variance Partitioning on Real SRM Time-Course Data (head):\n")
print(head(round(varPart, 3), 4))

message("\n=== Listing 22.7: Continuous Non-Linear Trajectory Modeling via Natural Cubic Splines ===")
form_spline <- ~ ns(Time_num, df = 3) + (1 | SUBJECT)
fit_spline <- dream(expr_mat, form_spline, meta_df, ddf = "Kenward-Roger", BPPARAM = SerialParam())
fit_spline_eb <- eBayes(fit_spline)

spline_terms <- grep("ns\\(Time_num", colnames(coef(fit_spline_eb)), value = TRUE)
res_spline <- variancePartition::topTable(fit_spline_eb, coef = spline_terms)

cat("Joint Moderated F-Test for Non-Linear Spline Trajectories (head):\n")
print(head(round(res_spline[, c("AveExpr", "F", "P.Value", "adj.P.Val")], 4), 3))

m_linear <- lme4::lmer(LogIntensities ~ Time + (1 | SUBJECT), data = df_idhc, REML = FALSE)
m_spline <- lme4::lmer(LogIntensities ~ ns(Time, df = 3) + (1 | SUBJECT), data = df_idhc, REML = FALSE)
cat("\nLikelihood Ratio Test (Linear vs. Natural Cubic Spline):\n")
print(anova(m_linear, m_spline))
message("\n[PASS] Chapter 22 Laboratory successfully executed.")
