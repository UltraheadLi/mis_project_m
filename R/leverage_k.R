# ==============================================================================
# File: /R/leverage_k.R
# Purpose: Leverage-aware k-selection for MIS. Replaces pure residual-based
#          k-selectors (gap_k, alpha_k) with an influence-based score that
#          captures both residual outlyingness AND leverage simultaneously.
#
# Root cause addressed:
#   alpha_k uses |r_i / scale| from the MM fit. Bad-leverage points have
#   extreme X-values, so the MM hat matrix partially absorbs them —
#   their robust residuals are SMALLER than their true influence on beta.
#   Result: alpha_k undercounts k in bad_leverage + heavy-tail scenarios.
#
#   The fix: score each observation by its approximate influence on beta,
#   not just by its residual. The Sherman-Morrison LOO formula gives:
#
#     influence_i ≈ |e_j' (X'X)^{-1} x_i * r_i / (1 - h_i)|
#
#   which is proportional to |DFBETAS_i| — i.e., how much beta_j changes
#   when obs i is deleted. Points with high leverage AND moderate residuals
#   score high, exactly the cases alpha_k misses.
#
#   We compute this from the ROBUST fit (lmrob), so the residuals and
#   (X'X)^{-1} are resistant to contamination. Then apply the same
#   ratio-gap detection logic that works well in alpha_k.
#
# Performance: O(n log n) — one matrix multiply + one sort. Negligible.
# ==============================================================================

#' Leverage-Aware k-Selection via Robust Influence Scores
#'
#' Computes a robust influence score for each observation that combines
#' residual magnitude with leverage, then uses ratio-gap detection on the
#' sorted scores to find the number of outliers.
#'
#' The influence score is derived from the MM-estimator's fit:
#'   score_i = |lev_col_j * r_i / (1 - h_i)|
#' where:
#'   - lev_col_j = e_j' (X'X)^{-1} x_i  (leverage on target coefficient j)
#'   - r_i = robust residual from lmrob
#'   - h_i = hat value from the robust fit
#'
#' This is proportional to the LOO change in beta_j — it captures BOTH
#' residual-outliers (large r_i, normal h_i) AND leverage-outliers
#' (moderate r_i, extreme h_i), which pure residual methods miss.
#'
#' @param mod_mm    A fitted lmrob object (or NULL if fitting failed).
#' @param target_pos Integer; position of the target coefficient in the
#'        design matrix. Default = 2 (first slope with intercept).
#' @param max_k_frac Numeric; cap k at this fraction of n.
#'        Default = 0.05 (5%).
#' @param ratio_threshold Numeric; minimum ratio between consecutive
#'        sorted scores to detect a gap. Default = 1.5.
#' @param gap_threshold Numeric; minimum absolute gap between consecutive
#'        sorted scores. Default = 0 (disabled — rely on ratio only).
#'        Set > 0 to require both ratio AND absolute gap.
#'
#' @return Integer; the estimated number of outliers (0 if none detected).
#' @export
leverage_k <- function(mod_mm,
                       target_pos   = 2L,
                       max_k_frac   = 0.05,
                       ratio_threshold = 1.5,
                       gap_threshold   = 0) {

  if (is.null(mod_mm)) return(0L)

  # ------------------------------------------------------------------
  # 1. Extract components from the robust fit
  # ------------------------------------------------------------------
  X   <- stats::model.matrix(mod_mm)
  n   <- nrow(X)
  p   <- ncol(X)
  max_k <- max(1L, floor(n * max_k_frac))

  r_robust <- mod_mm$residuals           # robust residuals
  h_robust <- stats::hatvalues(mod_mm)   # robust hat values

  # ------------------------------------------------------------------
  # 2. Compute per-observation influence on target coefficient
  #
  # The LOO change in beta[target_pos] when obs i is deleted:
  #   delta_i = (X'X)^{-1}[target_pos, ] %*% x_i * r_i / (1 - h_i)
  #
  # We need (X'X)^{-1} from the robust fit. lmrob stores the QR
  # decomposition; we extract it. If unavailable, compute from X
  # weighted by the final IRWLS weights.
  # ------------------------------------------------------------------
  w <- stats::weights(mod_mm, type = "robustness")
  if (is.null(w)) w <- rep(1, n)

  # Weighted (X'WX)^{-1} — this is what lmrob's coefficient covariance uses
  Xw <- X * sqrt(pmax(w, 0))
  XX_inv <- tryCatch(
    chol2inv(qr.R(qr(Xw))),
    error = function(e) {
      # Fallback: use unweighted if weighted QR fails
      tryCatch(chol2inv(qr.R(qr(X))), error = function(e2) NULL)
    }
  )

  if (is.null(XX_inv)) return(0L)

  # Leverage of each obs on the target coefficient
  lev_col <- as.numeric(XX_inv[target_pos, , drop = FALSE] %*% t(X))

  # ------------------------------------------------------------------
  # 3. Influence scores (absolute LOO DFBETA approximation)
  # ------------------------------------------------------------------
  denom <- pmax(1 - h_robust, 1e-12)
  influence_scores <- abs(lev_col * r_robust / denom)

  # ------------------------------------------------------------------
  # 4. Ratio-gap detection on sorted influence scores
  #
  # Same logic as alpha_k but on influence scores instead of residuals.
  # Sort descending, look at ratios between consecutive scores.
  # A large ratio indicates a natural break: observations above it
  # are qualitatively more influential than those below.
  # ------------------------------------------------------------------
  sorted_scores <- sort(influence_scores, decreasing = TRUE)
  top_vals <- sorted_scores[1:min(max_k + 1L, n)]

  if (length(top_vals) < 2L) return(0L)

  m <- length(top_vals) - 1L
  gaps   <- top_vals[1:m] - top_vals[2:(m + 1L)]
  ratios <- top_vals[1:m] / pmax(top_vals[2:(m + 1L)], 1e-15)

  # Find candidate break points
  if (gap_threshold > 0) {
    candidates <- which(ratios > ratio_threshold & gaps > gap_threshold)
  } else {
    candidates <- which(ratios > ratio_threshold)
  }

  if (length(candidates) == 0L) return(0L)

  # Select the break with the largest gap (most decisive separation)
  best <- candidates[which.max(gaps[candidates])]
  return(as.integer(min(best, max_k)))
}


#' Combined Leverage-Aware Adaptive k with Fallback
#'
#' Primary: leverage_k (captures both residual and leverage outliers).
#' Fallback: alpha_k (catches cases where leverage is normal but
#' residuals are extreme — pure vertical outliers).
#'
#' This mirrors the structure of adaptive_k but replaces gap_k
#' (which always saturates at max_k) with leverage_k as the primary.
#'
#' @param mod_mm         A fitted lmrob object (or NULL).
#' @param target_pos     Integer; coefficient position (default = 2).
#' @param max_k_frac     Numeric; max fraction of n to remove (default = 0.05).
#' @param ratio_threshold Numeric; passed to leverage_k (default = 1.5).
#' @param alpha          Numeric; passed to alpha_k fallback (default = 1).
#'
#' @return Integer; the estimated number of outliers.
#' @export
adaptive_k_v2 <- function(mod_mm,
                          target_pos      = 2L,
                          max_k_frac      = 0.05,
                          ratio_threshold = 1.5,
                          alpha           = 1) {

  # Primary: leverage-aware detection
  k_lev <- leverage_k(mod_mm,
                      target_pos      = target_pos,
                      max_k_frac      = max_k_frac,
                      ratio_threshold = ratio_threshold)

  if (k_lev > 0L) return(k_lev)

  # Fallback: residual-only detection (for pure vertical outliers
  # where leverage is uninformative)
  k_res <- alpha_k(mod_mm, alpha = alpha, max_k_frac = max_k_frac)
  return(k_res)
}
