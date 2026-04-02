rm(list = ls())
set.seed(1)

# =========================================================
# CCA with permutation + bootstrap
# CSV columns:
#   1-4  : SES scales (Y block)
#   5    : Behavior PC1
#   6-8  : EEG PC1-3
# N ~ 72
# =========================================================

# ---------------------------
# 0) User settings
# ---------------------------
setwd("/Users/Re-Re/Desktop/Researches/SES/R1/")  
dataFile <- "CCAData.csv"                   

B_boot <- 1000   # bootstrap iterations
M_perm <- 1000   # permutation iterations

out_dir <- "CCA_Results"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------
# 1) Read data
# ---------------------------
dat <- read.csv(dataFile, stringsAsFactors = FALSE)

# Optional: if you have SubjectID column, keep it for merging later
has_id <- "SubjectID" %in% colnames(dat)
if (has_id) subj_id <- dat$SubjectID

# Enforce the assumed column structure
Y_raw <- dat[, 1:4]      # SES
X_raw <- dat[, 5:8]      # Beh_PC1 + EEG_PC1-3

# ---------------------------
# 2) Mean imputation (column-wise)
# ---------------------------
mean_impute <- function(M) {
  M <- as.data.frame(M)
  for (j in seq_len(ncol(M))) {
    v <- M[[j]]
    if (anyNA(v)) {
      M[[j]][is.na(v)] <- mean(v, na.rm = TRUE)
    }
  }
  return(M)
}

Y_imp <- mean_impute(Y_raw)
X_imp <- mean_impute(X_raw)

# ---------------------------
# 3) Z-score each column (across subjects)
# ---------------------------
Y <- scale(as.matrix(Y_imp))
X <- scale(as.matrix(X_imp))

# Column names for outputs
colnames(Y) <- paste0("SES_", 1:ncol(Y))
colnames(X) <- c("Beh_PC1", "EEG_PC1", "EEG_PC2", "EEG_PC3")

# ---------------------------
# 4) Fit CCA (classic)
# ---------------------------
fit_cca <- function(X, Y) {
  # cancor expects complete numeric matrices
  cc <- cancor(X, Y)
  
  # canonical correlations
  r <- cc$cor
  
  # canonical variates (scores)
  U <- X %*% cc$xcoef
  V <- Y %*% cc$ycoef
  
  # canonical loadings: corr(original vars, canonical variates)
  x_load <- cor(X, U)
  y_load <- cor(Y, V)
  
  list(
    cc = cc,
    r = r,
    U = U,
    V = V,
    x_load = x_load,
    y_load = y_load
  )
}

obs <- fit_cca(X, Y)

# Save observed canonical correlations
obs_r_df <- data.frame(
  Pair = paste0("CV", seq_along(obs$r)),
  CanonicalCorrelation = obs$r
)
write.csv(obs_r_df, file.path(out_dir, "Observed_CanonicalCorrelations.csv"), row.names = FALSE)

# Save observed weights and loadings
write.csv(obs$cc$xcoef, file.path(out_dir, "Observed_X_CanonicalWeights.csv"))
write.csv(obs$cc$ycoef, file.path(out_dir, "Observed_Y_CanonicalWeights.csv"))
write.csv(obs$x_load,   file.path(out_dir, "Observed_X_CanonicalLoadings.csv"))
write.csv(obs$y_load,   file.path(out_dir, "Observed_Y_CanonicalLoadings.csv"))

# ---------------------------
# 5) Permutation test
#    Null: X-Y association is absent; permute Y rows
#    Test statistic: first canonical correlation (CV1)
#    (Optional) also compute Wilks' Lambda across all dimensions
# ---------------------------
wilks_lambda <- function(r_vec, p, q) {
  # For canonical correlations r1..rk, Wilks' Lambda = Π (1 - r_i^2)
  # where k = min(p,q)
  k <- min(p, q, length(r_vec))
  prod(1 - r_vec[1:k]^2)
}

perm_r1 <- numeric(M_perm)
perm_lambda <- numeric(M_perm)

p <- ncol(X); q <- ncol(Y)

for (m in seq_len(M_perm)) {
  idx <- sample.int(nrow(Y), replace = FALSE)
  Yp <- Y[idx, , drop = FALSE]
  fp <- fit_cca(X, Yp)
  perm_r1[m] <- fp$r[1]
  perm_lambda[m] <- wilks_lambda(fp$r, p, q)
}

# Empirical p-values (one-sided for correlation; smaller lambda indicates stronger association)
p_r1 <- (1 + sum(perm_r1 >= obs$r[1])) / (M_perm + 1)
obs_lambda <- wilks_lambda(obs$r, p, q)
p_lambda <- (1 + sum(perm_lambda <= obs_lambda)) / (M_perm + 1)

perm_summary <- data.frame(
  Statistic = c("CV1_canonical_r", "Wilks_Lambda_all"),
  Observed  = c(obs$r[1], obs_lambda),
  P_value   = c(p_r1, p_lambda),
  Permutations = c(M_perm, M_perm)
)
write.csv(perm_summary, file.path(out_dir, "Permutation_Test_Summary.csv"), row.names = FALSE)

# Save permutation distributions
write.csv(data.frame(perm_r1 = perm_r1), file.path(out_dir, "Permutation_Null_CV1_r.csv"), row.names = FALSE)
write.csv(data.frame(perm_lambda = perm_lambda), file.path(out_dir, "Permutation_Null_WilksLambda.csv"), row.names = FALSE)

# ---------------------------
# 6) Bootstrap
#    Resample subjects with replacement; refit CCA
#    Collect:
#      - CV1 canonical correlation
#      - CV1 weights & loadings (with sign alignment)
# ---------------------------
align_sign <- function(vec, ref) {
  # Align sign to maximize similarity to reference
  if (sum(vec * ref) < 0) -vec else vec
}

boot_r1 <- numeric(B_boot)

# For CV1 only (most common reporting); extend to more pairs if needed
boot_xw1 <- matrix(NA_real_, nrow = B_boot, ncol = ncol(X))
boot_yw1 <- matrix(NA_real_, nrow = B_boot, ncol = ncol(Y))
boot_xl1 <- matrix(NA_real_, nrow = B_boot, ncol = ncol(X))
boot_yl1 <- matrix(NA_real_, nrow = B_boot, ncol = ncol(Y))

colnames(boot_xw1) <- colnames(X)
colnames(boot_yw1) <- colnames(Y)
colnames(boot_xl1) <- colnames(X)
colnames(boot_yl1) <- colnames(Y)

# Reference (observed) for sign alignment: use CV1 X-weights and Y-weights
ref_xw1 <- as.numeric(obs$cc$xcoef[, 1])
ref_yw1 <- as.numeric(obs$cc$ycoef[, 1])

for (b in seq_len(B_boot)) {
  idx <- sample.int(nrow(X), replace = TRUE)
  Xb <- X[idx, , drop = FALSE]
  Yb <- Y[idx, , drop = FALSE]
  
  fb <- fit_cca(Xb, Yb)
  
  # CV1 correlation
  boot_r1[b] <- fb$r[1]
  
  # CV1 weights
  xw1 <- as.numeric(fb$cc$xcoef[, 1])
  yw1 <- as.numeric(fb$cc$ycoef[, 1])
  
  # Align signs (CCA sign is arbitrary)
  xw1 <- align_sign(xw1, ref_xw1)
  yw1 <- align_sign(yw1, ref_yw1)
  
  # Also align loadings consistently with X-weights alignment
  xl1 <- as.numeric(fb$x_load[, 1])
  yl1 <- as.numeric(fb$y_load[, 1])
  
  xl1 <- align_sign(xl1, as.numeric(obs$x_load[, 1]))
  yl1 <- align_sign(yl1, as.numeric(obs$y_load[, 1]))
  
  boot_xw1[b, ] <- xw1
  boot_yw1[b, ] <- yw1
  boot_xl1[b, ] <- xl1
  boot_yl1[b, ] <- yl1
}

# Bootstrap CI for CV1 canonical correlation
ci_r1 <- quantile(boot_r1, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
boot_summary <- data.frame(
  Statistic = "CV1_canonical_r",
  Observed = obs$r[1],
  CI_2.5 = ci_r1[1],
  Median = ci_r1[2],
  CI_97.5 = ci_r1[3],
  Bootstraps = B_boot
)
write.csv(boot_summary, file.path(out_dir, "Bootstrap_CV1_r_CI.csv"), row.names = FALSE)

# Save bootstrap distributions
write.csv(data.frame(boot_r1 = boot_r1), file.path(out_dir, "Bootstrap_CV1_r.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_xw1), file.path(out_dir, "Bootstrap_CV1_X_Weights.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_yw1), file.path(out_dir, "Bootstrap_CV1_Y_Weights.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_xl1), file.path(out_dir, "Bootstrap_CV1_X_Loadings.csv"), row.names = FALSE)
write.csv(as.data.frame(boot_yl1), file.path(out_dir, "Bootstrap_CV1_Y_Loadings.csv"), row.names = FALSE)

# ---------------------------
# 7) (Optional) Quick plots (base R)
# ---------------------------
png(file.path(out_dir, "Permutation_CV1_r_hist.png"), width = 1200, height = 900, res = 150)
hist(perm_r1, breaks = 40, main = "Permutation null: CV1 canonical r", xlab = "CV1 canonical correlation")
abline(v = obs$r[1], lwd = 3)
mtext(paste0("Observed r = ", round(obs$r[1], 3), " | p = ", signif(p_r1, 3)), side = 3, line = 0.5)
dev.off()

png(file.path(out_dir, "Bootstrap_CV1_r_hist.png"), width = 1200, height = 900, res = 150)
hist(boot_r1, breaks = 40, main = "Bootstrap: CV1 canonical r", xlab = "CV1 canonical correlation")
abline(v = obs$r[1], lwd = 3)
abline(v = ci_r1[c(1,3)], lwd = 2, lty = 2)
mtext(paste0("Observed r = ", round(obs$r[1], 3),
             " | 95% CI [", round(ci_r1[1], 3), ", ", round(ci_r1[3], 3), "]"),
      side = 3, line = 0.5)
dev.off()

cat("Done. Results saved to:", normalizePath(out_dir), "\n")