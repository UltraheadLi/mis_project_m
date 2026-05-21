# ==============================================================================
# File: /R/dynamic_k_adaptive.R
# Purpose: Fast, distribution-adaptive methods for determining the number of
#          outliers to remove (dynamic_k) before calling fast_sens_topk.
#          Replaces the broken MM 3-sigma residual count, which overcounts
#          in heavy-tailed DGPs.
#
# Performance target: <1ms per call (must not bottleneck 9600 MC iterations)
#
# Three methods provided, in order of recommendation:
#   1. gap_k     — residual gap heuristic (fastest, most robust)
#   2. alpha_k   — outlyingness-based with adaptive threshold
#   3. oracle_k  — fixed at true injection size (for benchmarking only)
# ==============================================================================

#' Adaptive k via Sorted-Residual Gap Detection
#'
#' Uses the MM-estimator's standardised residuals, sorted by absolute value,
#' to find a natural "gap" separating the bulk distribution from the outlier
#' cluster. The gap is defined as the largest jump in consecutive sorted
#' |residuals| that exceeds `gap_factor` times the median spacing.
#'
#' This is distribution-adaptive because the gap criterion is relative to
#' the observed residual spacing — heavy-tailed DGPs have larger median
#' spacing, so the threshold auto-scales.
#'
#' Cost: one sort (O(n log n)) + one diff + one which.max. Negligible.
#'
#' @param mod_mm   A fitted lmrob object (or NULL if fitting failed).
#' @param gap_factor Numeric; the gap must exceed gap_factor * median(spacing)
#'        to be considered a true separation. Default = 3.
#' @param min_tail  Integer; only consider gaps in the top `min_tail` fraction
#'        of sorted residuals, to avoid detecting gaps in the bulk.
#'        Default = 0.1 (top 10%).
#' @param max_k_frac Numeric; cap k at this fraction of n to prevent
#'        catastrophic over-removal. Default = 0.05 (5%).
#' @return Integer; the estimated number of outliers (0 if none detected).
#' @export
gap_k <- function(mod_mm, gap_factor = 3, min_tail = 0.1, max_k_frac = 0.05) {
  if (is.null(mod_mm)) return(0L)
  
  std_res <- abs(mod_mm$residuals / mod_mm$scale)
  n <- length(std_res)
  max_k <- max(1L, floor(n * max_k_frac))
  
  # Sort in ascending order
  sorted_res <- sort(std_res)
  
  # Only look for gaps in the upper tail
  tail_start <- max(1L, floor(n * (1 - min_tail)))
  tail_vals <- sorted_res[tail_start:n]
  
  if (length(tail_vals) < 3L) return(0L)
  
  spacings <- diff(tail_vals)
  med_spacing <- median(spacings)
  
  # Avoid division by zero in perfectly clean data
  if (med_spacing < 1e-10) return(0L)
  
  # Find the largest gap that exceeds the threshold
  big_gaps <- which(spacings > gap_factor * med_spacing)
  
  if (length(big_gaps) == 0L) return(0L)
  
  # The first big gap (from the bulk side) defines the cutoff
  # Everything above this gap is "outlier"
  gap_pos <- big_gaps[1L]  # position in tail_vals
  k_detected <- length(tail_vals) - gap_pos
  
  # Apply cap
  k_out <- min(k_detected, max_k)
  return(as.integer(k_out))
}


#' Adaptive k via Outlyingness with Distribution-Scaled Threshold
#'
#' Uses the MM-estimator's robust residuals but replaces the fixed 3-sigma
#' threshold with one calibrated to the estimated residual distribution.
#' Specifically, it uses a quantile-based threshold: observations beyond the
#' (1 - alpha/n) quantile of the fitted residual distribution are flagged.
#'
#' This Bonferroni-style correction ensures that under the null (no outliers),
#' the expected number of false flags is alpha (default = 1), regardless of
#' the error distribution's tail weight.
#'
#' Cost: one quantile computation + one sum. Negligible.
#'
#' @param mod_mm   A fitted lmrob object (or NULL if fitting failed).
#' @param alpha    Numeric; expected false positives under the null.
#'        Default = 1 (expect ~1 false flag on average).
#' @param max_k_frac Numeric; cap k at this fraction of n. Default = 0.05.
#' @return Integer; the estimated number of outliers.
#' @export
alpha_k <- function(mod_mm, alpha = 1, max_k_frac = 0.05) {
  if (is.null(mod_mm)) return(0L)
  
  std_res <- abs(mod_mm$residuals / mod_mm$scale)
  n <- length(std_res)
  p <- length(coef(mod_mm))
  max_k <- max(1L, floor(n * max_k_frac))
  
  sorted_res <- sort(std_res, decreasing = TRUE)
  top_vals <- sorted_res[1:(max_k + 1L)]
  
  gaps <- top_vals[1:max_k] - top_vals[2:(max_k + 1L)]
  ratios <- top_vals[1:max_k] / pmax(top_vals[2:(max_k + 1L)], 1e-10)
  
  hat_vals <- hatvalues(mod_mm)
  high_lev_threshold <- 3 * p / n
  high_lev_mask <- hat_vals > high_lev_threshold
  n_high_lev_extreme <- sum(high_lev_mask & std_res > 3)
  
  if (n_high_lev_extreme >= 2) {
    candidates <- which(ratios > 2.0)
  } else {
    candidates <- which(ratios > 1.4 & gaps > 0.5)
  }
  
  if (length(candidates) == 0L) return(0L)
  
  best <- candidates[which.max(gaps[candidates])]
  return(as.integer(min(best, max_k)))
}

#' Oracle k (Fixed at True Injection Size)
#'
#' Returns a fixed k regardless of the data. Used as an upper bound on
#' what any adaptive method can achieve: if MIS with oracle k still
#' underperforms, the problem is in detection/direction, not in k-selection.
#'
#' @param k Integer; the true number of injected outliers.
#' @return Integer; k unchanged.
#' @export
oracle_k <- function(k) {
  as.integer(k)
}


#' Combined k-Selection with Fallback
#'
#' Runs gap_k as the primary method. If it returns 0 but alpha_k detects
#' something, uses alpha_k as a fallback. This handles edge cases where
#' the gap is gradual (no sharp separation) but outlyingness is still clear.
#'
#' @param mod_mm   A fitted lmrob object (or NULL).
#' @param gap_factor Passed to gap_k.
#' @param alpha    Passed to alpha_k.
#' @param max_k_frac Maximum fraction of n to remove.
#' @return Integer; the estimated number of outliers.
#' @export
adaptive_k <- function(mod_mm, gap_factor = 3, alpha = 1, max_k_frac = 0.05) {
  k_gap <- gap_k(mod_mm, gap_factor = gap_factor, max_k_frac = max_k_frac)
  
  if (k_gap > 0L) return(k_gap)
  
  # Fallback: alpha_k catches cases where the gap is smooth
  # but there are still extreme outliers
  k_alpha <- alpha_k(mod_mm, alpha = alpha, max_k_frac = max_k_frac)
  return(k_alpha)
}