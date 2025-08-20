# ===============================
# Helpers for environments selection
# ===============================

# ---------- Determinants over env-subsets ----------
# Compute Gram(det) for all k-row subsets of 'M' (rows = envs, cols = features incl. intercept)
compute_determinants_fast <- function(M, k) {
  n <- nrow(M)
  combs <- utils::combn(n, k)
  nC <- ncol(combs)
  
  dets <- numeric(nC)
  # simple loop (vectorizing over huge combn objects is memory-heavy)
  pb <- utils::txtProgressBar(min = 0, max = nC, style = 3)
  for (i in seq_len(nC)) {
    R <- M[combs[, i], , drop = FALSE]
    G <- crossprod(R)   # t(R) %*% R
    dets[i] <- det(G)
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  list(determinants = dets, combination_indices = t(combs))
}

# ---------- Orthogonal basis via Gram–Schmidt ----------
orthogonal_basis <- function(v) {
  p <- length(v)
  stopifnot(any(v != 0))
  B <- matrix(0, nrow = p, ncol = p)
  B[,1] <- v / sqrt(sum(v^2))
  for (i in 2:p) {
    e <- rep(0, p); e[i] <- 1
    for (j in 1:(i-1)) {
      e <- e - sum(B[,j] * e) * B[,j]
    }
    e_norm <- sqrt(sum(e^2))
    if (e_norm < .Machine$double.eps) {
      # fall back: pick a random orthogonal direction
      e <- rnorm(p); e <- e - sum(B[,1] * e) * B[,1]
      e <- e / sqrt(sum(e^2))
    } else {
      e <- e / e_norm
    }
    B[,i] <- e
  }
  B
}

# ---------- Per-environment summaries ----------
# batch_col is the column index of ENV in X (data.frame or matrix)
compute_covariances_env <- function(X, batch_col) {
  envs <- sort(unique(X[, batch_col]))
  out  <- vector("list", length(envs)); names(out) <- as.character(envs)
  for (i in seq_along(envs)) {
    e <- envs[i]
    M <- X[X[, batch_col] == e, -batch_col, drop = FALSE]
    out[[i]] <- stats::cov(M)
  }
  out
}

compute_averages_env <- function(X, batch_col) {
  envs <- sort(unique(X[, batch_col]))
  out  <- vector("list", length(envs)); names(out) <- as.character(envs)
  for (i in seq_along(envs)) {
    e <- envs[i]
    M <- X[X[, batch_col] == e, -batch_col, drop = FALSE]
    out[[i]] <- c(1, colMeans(M))
  }
  out
}

compute_batch_counts_env <- function(X, batch_col) {
  envs <- sort(unique(X[, batch_col]))
  out  <- vector("list", length(envs)); names(out) <- as.character(envs)
  for (i in seq_along(envs)) {
    e <- envs[i]
    out[[i]] <- nrow(X[X[, batch_col] == e, , drop = FALSE])
  }
  out
}

# ---------- Causal key score (your original structure) ----------
causal_key <- function(covs, avs, counts, avgbeta, betastar, basis, sigma2) {
  Z <- length(counts)
  N <- sum(as.numeric(counts))
  termz <- numeric(ncol(basis))
  for (b in seq_len(ncol(basis))) {
    v <- basis[, b]
    acc <- 0
    for (z in seq_len(Z)) {
      muz   <- avs[[z]]          # length p+1
      Sigz  <- covs[[z]]         # p x p
      nz    <- as.numeric(counts[[z]])
      wz    <- nz / N
      # NOTE: muz[-1] aligns with cov matrix; v[-1] drops intercept direction
      acc <- acc + ((wz / sqrt(nz)) *
                      ( t(muz) %*% betastar * sqrt(t(v[-1]) %*% Sigz %*% v[-1]) -
                          sqrt(sigma2) * t(muz) %*% v ))^2
    }
    termz[b] <- acc
  }
  min(termz)
}

# ---------- DRIG data helper ----------
# Build a list per environment with [y, 1, X] as required by est_drig()
as_drig_list <- function(Y, X_df, y_first = TRUE) {
  stopifnot("ENV" %in% names(X_df))
  envs <- sort(unique(X_df$ENV))
  out  <- vector("list", length(envs)); names(out) <- as.character(envs)
  for (i in seq_along(envs)) {
    e <- envs[i]
    rows <- which(X_df$ENV == e)
    Xsub <- X_df[rows, , drop = FALSE]
    Xcov <- as.matrix(Xsub[, setdiff(names(Xsub), "ENV"), drop = FALSE])
    if (y_first) {
      out[[i]] <- cbind(Y[rows], cbind(1, Xcov))
    } else {
      out[[i]] <- cbind(cbind(1, Xcov), Y[rows])
    }
  }
  out
}