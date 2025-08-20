# ===============================
# Generative Invariance on SPRINT
# ===============================
# Follows the GI setup in "Domain adaptation under hidden confounding"
# Key references in code comments: Eq. (3) generator, Thm 4 identifiability, and the closed-form

# ---- Parameters --------------------------------------------------------------
ONLY_BEST=F
project_dir     <- "~/Projects/dahic"              # change if needed
data_csv        <- "/Users/cgmeixide/Desktop/SPRINT_2021a/SPRINT-POP/data/baseline.csv"
min_per_site    <- 100                             # filter out small sites
cont_vars       <- c("AGE","SBP","BMI")            # covariates
resp_var        <- "RISK10YRS"                     # response
n_envs_combo    <- length(cont_vars) + 1           # k in combn (Z ~ p+1)
n_top_combos    <- 100                             # how many top Gram-determinant env-sets to try
set.seed(1)

# ---- Libraries ---------------------------------------------------------------

library(dplyr)
library(energy)
library(combinat) # provides combn alternative, but we use base::combn


# ---- Sources (helpers & any estimators you already have) ---------------------
setwd(project_dir)
source("aux_gi.R")
# Your high-dimensional estimator providing betahat, khat, vary:
#   hd(Y, X_df) -> list(betahat = ..., khat = ..., vary = ...)
source("hd.R")
# DRIG helpers you already use (for baseline):
source("~/Dropbox/invariant/aux_drig_functions.R")

# ---- Load & basic preprocessing ---------------------------------------------
bas <- na.omit(read.csv(data_csv))
# Sort by site for reproducibility
bas <- bas[order(bas$NEWSITEID), ]
# Keep only sites with enough size
eligible_sites <- names(which(table(bas$NEWSITEID) > min_per_site))
bas <- bas[bas$NEWSITEID %in% as.numeric(eligible_sites), ]
cat("Unique eligible sites:", length(unique(bas$NEWSITEID)), "\n")

# Response & design frame with ENV tag
Yall <- bas[[resp_var]]
Xall <- bas[, cont_vars, drop = FALSE]
Xall$ENV <- bas$NEWSITEID

# ---- Choose TEST environments using Energy Distance -------------------------
# We compute inter-environment energy distances on X (cont. vars only)
X_sorted      <- bas[, cont_vars, drop = FALSE]
envs_sorted   <- bas$NEWSITEID
sizes         <- table(envs_sorted)
env_ids       <- as.numeric(names(sizes))

# Energy distance matrix between environments (energy::edist with sizes)
edist_matrix  <- edist(x = X_sorted, sizes = sizes, distance = FALSE)
dmat          <- as.matrix(edist_matrix)
cat("Energy distance matrix (env-by-env):\n")
print(round(dmat, 3))

# Score environments by median distance to others (high = more "different")
med_scores    <- apply(dmat, 2, median)
ord_by_med    <- order(med_scores, decreasing = TRUE)
# Take the top 10 most different as TEST; the rest are TRAIN
test_envs     <- env_ids[ord_by_med][1:10]
train_envs    <- setdiff(env_ids, test_envs)

cat("Test envs:", paste(test_envs, collapse = ", "), "\n")
cat("Train envs:", paste(train_envs, collapse = ", "), "\n")

# ---- Environment summaries for GI (means/covs and counts) -------------------
# (used later for keys and for forming the GI generator mean-predictions)
means_list <- compute_averages_env(Xall, batch_col = ncol(Xall))     # list per ENV: c(1, mu)
covs_list  <- compute_covariances_env(Xall, batch_col = ncol(Xall))  # list per ENV: cov matrix of X
counts     <- compute_batch_counts_env(Xall, batch_col = ncol(Xall)) # list per ENV: n_z

# Build the "means" matrix per TRAIN env (for determinant selection)
means_df <- Xall %>%
  filter(ENV %in% train_envs) %>%
  group_by(ENV) %>%
  summarise(across(everything(), mean), .groups = "drop")

# Muss = [1, mu_z] by env, stacked as rows
mus      <- as.matrix(means_df[, cont_vars, drop = FALSE])
muss     <- cbind(Intercept = 1, mus)

# ---- Combinatorial selection of env-sets by Gram determinant ----------------
# We want size 'k = n_envs_combo' subsets with large det(t(M)M) where M is rows of muss
k <- n_envs_combo
det_out <- compute_determinants_fast(muss, k = k) # returns $determinants, $combination_indices

ord_desc      <- order(det_out$determinants, decreasing = TRUE)
comb_top_idx  <- det_out$combination_indices[ord_desc, , drop = FALSE]
comb_top_idx  <- comb_top_idx[seq_len(min(n_top_combos, nrow(comb_top_idx))), , drop = FALSE]

B <- nrow(comb_top_idx)
cat("Trying top", B, "env-combinations (size", k, ")\n")

# ---- Fit hd() per combination and cache estimates ---------------------------
p <- length(cont_vars)
betahats <- matrix(NA_real_, nrow = B, ncol = p + 1) # including intercept
khats    <- matrix(NA_real_, nrow = B, ncol = p)     # K dimension = p
sigma2s  <- numeric(B)

for (b in seq_len(B)) {
  env_idx <- comb_top_idx[b, ]
  chosen_envs <- train_envs[env_idx]    # map row indices to ENV ids
  tr_rows     <- which(Xall$ENV %in% chosen_envs)
  
  Y <- Yall[tr_rows]
  X <- Xall[tr_rows, , drop = FALSE]
  
  # hd() expects a data.frame of covariates; drop ENV column inside hd if needed
  Xdf <- data.frame(X)
  fit <- hd(Y, Xdf)
  betahats[b, ] <- fit$betahat   # length p+1 (intercept first)
  khats[b, ]    <- fit$khat      # length p
  sigma2s[b]    <- fit$vary      # scalar
}

# ---- Score each combination with the "causal key" ---------------------------
# Basis from current b (orthogonal to mean betahat), then per-env cov/means/counts
keys <- numeric(B)

beta_bar <- colMeans(betahats)

for (b in seq_len(B)) {
  # orthonormal basis vectors spanning directions orthogonal to beta_bar (drop first/Intercept)
  Bmat  <- orthogonal_basis(betahats[b, ] - beta_bar)[, 2:(p+1), drop = FALSE]
  
  # Recompute stats on current training split only (weighted by in-split counts)
  env_idx    <- comb_top_idx[b, ]
  chosen_envs <- train_envs[env_idx]
  tr_rows    <- which(Xall$ENV %in% chosen_envs)
  
  Xtr <- Xall[tr_rows, , drop = FALSE]
  
  covs_b   <- compute_covariances_env(Xtr, batch_col = ncol(Xtr))
  avs_b    <- compute_averages_env(Xtr, batch_col = ncol(Xtr))
  counts_b <- compute_batch_counts_env(Xtr, batch_col = ncol(Xtr))
  
  keys[b] <- causal_key(
    covs   = covs_b,
    avs    = avs_b,
    counts = counts_b,
    avgbeta = beta_bar,
    betastar = betahats[b, ],
    basis  = Bmat,
    sigma2 = sigma2s[b]
  )
}

# Quick diagnostic plot in khat-space
plot(khats[,1], khats[,2], main = "k-hat scatter (first two coords)", xlab = "k1", ylab = "k2")
points(mean(khats[,1]), mean(khats[,2]), col = "red", pch = 19)
opt <- which.max(keys)
points(khats[opt,1], khats[opt,2], col = "green", pch = 19)
pairs(betahats, main = "betahat pairs")

# TRAIN set for OLS/DRIG baselines on *all* training envs (not just top combos)
train_rows <- which(Xall$ENV %in% train_envs)
Ytr        <- Yall[train_rows]
Xtr        <- Xall[train_rows, , drop = FALSE]
Xtr_df     <- data.frame(Xtr)

# DRIG/OLS fits (gamma controls robustness; gamma=1 ~ OLS baseline)
data_drig1 <- as_drig_list(Ytr, Xtr_df, y_first = TRUE)  # helper to build list per env
set.seed('0')
data_drig <- data_drig1[sample(length(data_drig1))]
drig_coef <- est_drig(data_drig, gamma = 0.3, y_idx = 1, del_idx = NULL, unif_weight = FALSE)
ols_coef  <- est_drig(data_drig, gamma = 1.0, y_idx = 1, del_idx = NULL, unif_weight = FALSE)

# ---- GI mean-predictions on TEST envs (Eq. (3) with ξ=0) --------------------
# For GI mean-prediction in test env e:
#  Ŷ_GI(x) = [1, x]^T beta_hat  +  K_hat^T Σ0^{-1} (x - μ0)

mse_ours <- numeric(length(test_envs))
mse_drig <- numeric(length(test_envs))
mse_ols  <- numeric(length(test_envs))

for (i in seq_along(test_envs)) {
  env_e  <- test_envs[i]
  rows_e <- which(Xall$ENV == env_e)
  
  xnew <- Xall[rows_e, -(p+1), drop = FALSE]  # drop ENV col
  Yte  <- Yall[rows_e]
  
  # Bag the top combos
  gens <- order(keys, decreasing = TRUE)[seq_len(max(1, B - 50))]
  norm <- length(gens)
  
  if(ONLY_BEST) {
    gens <- which.max(keys)  # index of the single best combination
    norm     <- 1
  }
  
  pred_ours <- generator(betahats[gens[1], ], khats[gens[1], ], xnew, 0, F) / norm
  if (norm > 1) {
    for (g in 2:norm) {
      pred_ours <- pred_ours + generator(betahats[gens[g], ], khats[gens[g], ], xnew, 0, F) / norm
    }
  }
  
  x_ <- as.matrix(cbind(1, xnew))
  pred_drig <- drop(x_ %*% drig_coef)
  pred_ols  <- drop(x_ %*% ols_coef)
  
  mse_ours[i] <- median((Yte - pred_ours)^2)
  mse_drig[i] <- median((Yte - pred_drig)^2)
  mse_ols[i]  <- median((Yte - pred_ols)^2)
}

cat("GI (generator) median MSE:", median(mse_ours), " | 30–70%:", quantile(mse_ours, c(0.3, 0.7)), "\n")
cat("DRIG           median MSE:", median(mse_drig), " | 30–70%:", quantile(mse_drig, c(0.3, 0.7)), "\n")
cat("OLS            median MSE:", median(mse_ols),  " | 30–70%:", quantile(mse_ols,  c(0.3, 0.7)), "\n")
# ---- Small visual: env-wise MSE comparison ----------------------------------
boxplot(list(GI = mse_ours, DRIG = mse_drig, OLS = mse_ols),
        main = "Test MSE across held-out environments",
        ylab = "Median squared error (per env)")


