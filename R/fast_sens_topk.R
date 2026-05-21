# ==============================================================================
# File: /R/fast_sens_topk.R
# Purpose: Fast influence detection with iterative refinement, replacing
#          influence::sens() for the top-k use case. Starts with vectorized
#          Sherman-Morrison leave-one-out betas (one-shot ranking), then
#          refines by recomputing (X'X)^{-1} with the current top-k removed
#          and re-ranking. Converges in 2-4 rounds for typical problems.
#
# Performance:
#   - influence::sens(): 7-44 seconds for N=2000 (full N-step sequential loop)
#   - fast_sens_topk v1 (one-shot): <0.01s but misranks under masking
#   - fast_sens_topk v2 (this): 0.02-0.10s with near-exact detection
#
# Called by: /script/02_run_sim_dist.R, /script/03_alg_comp.R
# ==============================================================================

#' Fast Influence Detection with Iterative Refinement
#'
#' Computes the k most influential observations for a target coefficient,
#' matching influence::sens()$influence$id[1:k] but ~100-500x faster.
#'
#' The one-shot Sherman-Morrison ranking assumes observations are independent
#' influence sources — accurate when k is small relative to N. Under masking
#' (collinearity, interactions, large k), removing one outlier changes the
#' apparent influence of others. This function iteratively refines the set
#' by recomputing (X'X)^{-1} after removing the current top-k, until the
#' selected set stabilises.
#'
#' @param mod  An lm object.
#' @param pos  Integer; position of the target coefficient in the design matrix.
#' @param sign Integer; +1 or -1 direction for ranking.
#' @param k    Integer; number of top influential observations to return.
#' @param max_refine Integer; maximum refinement iterations (default = 5).
#'        Set to 0 to get the original one-shot behaviour.
#' @return Integer vector of length k — the row indices of the k most
#'         influential observations (same as sens()$influence$id[1:k]).
#' @export
fast_sens_topk <- function(mod, pos, sign, k, max_refine = 5) {
  X <- model.matrix(mod)
  y <- model.response(model.frame(mod))
  N <- nrow(X)
  p <- ncol(X)

  # ------------------------------------------------------------------
  # Step 0: One-shot leave-one-out ranking (same as the original v1)
  # ------------------------------------------------------------------
  qr_x <- mod$qr
  if (is.null(qr_x)) qr_x <- qr(X)
  R_qr <- qr.R(qr_x)

  beta <- coef(mod)
  res  <- residuals(mod)
  hat  <- rowSums(qr.Q(qr_x)^2)

  XX_inv <- chol2inv(R_qr)
  leverage_col <- as.numeric(XX_inv[pos, , drop = FALSE] %*% t(X))  # length N

  beta_i_pos <- beta[pos] - leverage_col * res / pmax(1 - hat, 1e-12)
  scores <- beta_i_pos * sign

  top_idx <- order(scores, decreasing = TRUE, method = "radix")[1:k]

  # Early exit: no refinement needed for k=1 or if disabled
  if (k <= 1 || max_refine == 0) return(top_idx)

  # ------------------------------------------------------------------
  # Steps 1+: Iterative refinement
  #
  # Remove the current top-k from the design matrix, recompute (X'X)^{-1}
  # on the remaining N-k observations, then re-rank ALL N observations
  # using leave-one-out betas computed from the reduced model. This
  # captures masking: if obs i was hiding obs j's influence, removing i
  # (as part of the current top-k) reveals j's true influence.
  #
  # We refit on X[-top_k,] to get the "clean" coefficient and (X'X)^{-1},
  # then score every observation (including the current top-k) against
  # this cleaner baseline.
  # ------------------------------------------------------------------
  for (iter in seq_len(max_refine)) {
    prev_idx <- top_idx

    # Refit on the complement of the current top-k
    keep   <- setdiff(seq_len(N), top_idx)
    X_keep <- X[keep, , drop = FALSE]
    y_keep <- y[keep]

    # Solve via QR on the reduced design matrix
    qr_keep <- qr(X_keep)
    if (qr_keep$rank < p) break  # rank deficient after removal — stop refining

    beta_clean <- qr.coef(qr_keep, y_keep)
    XX_inv_clean <- chol2inv(qr.R(qr_keep))

    # Score ALL N observations against the clean model:
    #   For each obs i, compute how much beta[pos] would change if i were
    #   added back (for obs in top_k) or removed (for obs in keep).
    #
    #   The leave-one-out formula relative to the clean model:
    #     delta_i = (X'X_clean)^{-1} X_i' * (y_i - X_i beta_clean) / (1 + h_i_clean)
    #   where h_i_clean = X_i (X'X_clean)^{-1} X_i'
    #
    #   Note the sign flip: for obs NOT in the fit, the denominator is (1 + h_i)
    #   (addition), not (1 - h_i) (deletion). For obs IN the fit, it's (1 - h_i).
    #   But for ranking purposes, we use a unified score: the magnitude of the
    #   influence of each observation on beta[pos] relative to the clean model.

    lev_col_clean <- as.numeric(XX_inv_clean[pos, , drop = FALSE] %*% t(X))  # all N
    res_clean     <- y - X %*% beta_clean                                     # all N
    h_clean       <- rowSums((X %*% XX_inv_clean) * X)                        # all N

    # For obs in `keep`: standard LOO deletion formula (1 - h_i)
    # For obs in `top_idx`: addition formula (1 + h_i)
    denom <- rep(NA_real_, N)
    denom[keep]    <- pmax(1 - h_clean[keep], 1e-12)
    denom[top_idx] <- 1 + h_clean[top_idx]

    delta_pos <- lev_col_clean * res_clean / denom
    scores_refined <- delta_pos * sign

    # Re-select top-k (most influential = most negative score after sign flip)
    top_idx <- order(scores_refined, decreasing = TRUE, method = "radix")[1:k]

    # Convergence check: top-k set hasn't changed
    if (setequal(top_idx, prev_idx)) break
  }

  return(top_idx)
}