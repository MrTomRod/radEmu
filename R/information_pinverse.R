# Generalised inverse of the Fisher information on the directions the data identify.
#
# The information matrix `I` of the profiled likelihood (see `f_info`) is singular, or numerically so, in
# directions the data carry no information about: a taxon that is absent from a level of a factor
# (its coefficient runs to -Inf), a taxon seen in a single sample, a coefficient at the `max_abs_B` bound,
# a pair of nearly collinear covariates. `get_score_stat` used to add `lambda * max(diag(I))` to the diagonal
# (lambda = 1e-4, 1e-3, ...), which makes the statistic depend on lambda and on `j_ref`; `micro_wald` stopped
# with "computationally singular".
#
# This function finds those directions and returns the Moore-Penrose inverse of `I` restricted to the others
# (computed in covariate-scaled coordinates, so that it does not depend on the units of the columns of X):
#
#  1. Per taxon `j`, the own block `I[jj]` (p x p) is rescaled covariate-wise (the largest diagonal entry of
#     covariate `a` over all taxa becomes 1) and diagonalised. Eigen-directions with eigenvalue below
#     `block_tol` carry no information and are removed. (`I` is positive semi-definite, so its curvature in a
#     direction supported on one taxon is at most that of the principal block.)
#  2. As everywhere in radEmu the block of the reference taxon `j_ref` is deleted (this removes the row-shift
#     directions `1 (x) c`). What remains is diagonalised once more; eigenvalues below `null_tol` times the
#     largest are directions that involve several taxa (samples that are disconnected from the reference) and
#     are removed too.
#
# Returns NULL when nothing has to be removed: the caller then uses the ordinary solve, unchanged.
#
# @param I full (pJ x pJ) information matrix, as returned by `f_info`.
# @param p,J number of covariates and categories; `j_ref` the reference taxon (its block is deleted).
# @return NULL, or a list with `inverse` (p(J-1) x p(J-1), raw coordinates, rows and columns as `I[-block, -block]`),
#   `unidentified_fraction(H)` (relative norm of the columns of `H` (p(J-1) x m) outside the identified subspace)
#   and `dropped` (data frame: taxon, curvature of each removed single-taxon direction).
information_pinverse <- function(I, p, J, j_ref, block_tol = 1e-10, null_tol = 1e-12) {
  I <- as.matrix(I)
  block <- function(j) (j - 1) * p + seq_len(p)
  # covariate scale: largest own-block diagonal entry of each covariate over all taxa
  own_diag <- matrix(vapply(seq_len(J), function(j) diag(I[block(j), block(j), drop = FALSE]), numeric(p)), nrow = p)
  s_cov <- sqrt(pmax(apply(own_diag, 1, max), .Machine$double.xmin))
  s <- rep(s_cov, times = J)
  I_s <- I / outer(s, s)

  # 1. per-taxon directions without information
  keep_idx <- setdiff(seq_len(p * J), block(j_ref))
  P <- matrix(0, p * J, 0)
  dropped <- data.frame(taxon = integer(0), curvature = numeric(0))
  for (j in setdiff(seq_len(J), j_ref)) {
    idx <- block(j)
    e <- eigen(I_s[idx, idx, drop = FALSE], symmetric = TRUE)
    keep <- e$values >= block_tol
    if (all(keep)) {
      V <- diag(p)    # exactly the identity: no rounding noise when nothing is removed
    } else {
      V <- e$vectors[, keep, drop = FALSE]
      dropped <- rbind(dropped, data.frame(taxon = rep(j, sum(!keep)), curvature = pmax(e$values[!keep], 0)))
    }
    Pj <- matrix(0, p * J, ncol(V))
    Pj[idx, ] <- V
    P <- cbind(P, Pj)
  }

  # 2. cross-taxa directions
  P <- P[keep_idx, , drop = FALSE]
  s <- s[keep_idx]
  I_s <- I_s[keep_idx, keep_idx, drop = FALSE]
  I_P <- crossprod(P, I_s %*% P)
  e2 <- eigen(I_P, symmetric = TRUE)
  keep2 <- e2$values > null_tol * max(e2$values)
  n_null <- sum(!keep2)
  if (nrow(dropped) == 0 && n_null == 0) {
    return(NULL)
  }
  V2 <- e2$vectors[, keep2, drop = FALSE]
  Q <- P %*% V2                                 # orthonormal basis of the identified subspace (scaled)
  inverse_s <- Q %*% (t(Q) / e2$values[keep2])
  inverse <- inverse_s / outer(s, s)
  list(
    inverse = inverse,
    unidentified_fraction = function(H) {
      Hs <- as.matrix(H) / s
      res <- Hs - Q %*% crossprod(Q, Hs)
      sqrt(colSums(res^2) / colSums(Hs^2))
    },
    dropped = dropped,
    n_cross_taxa = n_null
  )
}

# Relative size of the part of H outside the identified subspace above which a contrast is regarded as not
# testable (score test) or not estimable (Wald): it depends on parameters the data carry no information about.
information_unidentified_tol <- 1e-3
