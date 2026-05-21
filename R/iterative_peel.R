# ==============================================================================
# File: /R/iterative_peel.R
# Purpose: Iterative peel-and-refit MIS detection. Bypasses the need to guess
#          k upfront by removing influential points one (or a small batch) at
#          a time, refitting the robust model after each removal, and stopping
#          when no more outliers are detected or when removal harms the fit.
#
# Why this works:
#   The fundamental problem with single-shot k-selection is masking: in
#   bad-leverage contamination, outlier i masks outlier j so that neither
#   appears extreme in the robust residuals. No single residual/influence
#   threshold can reliably count all k outliers at once.
#
#   Iterative peeling sidesteps this entirely. After removing the most
#   obvious outlier(s), the MM refit unmasks the next layer. The k-selector
#   only needs to answer "is there still at least one outlier?" — a much
#   easier binary question than "how many total?"
#
# Two stopping rules (whichever triggers first):
#   1. No-detection stop: the k-selector returns 0 on the cleaned data
#   2. Direction stop: removing the next batch moves the OLS coefficient
#      AWAY from the robust anchor — we're removing signal, not contamination
#   3. Max-iterations cap: hard ceiling to prevent runaway peeling
#
# Performance:
#   Each iteration: 1 lmrob + 1 dinkelbach + 1 OLS  ≈  3-5ms
#   Typical convergence: 5-25 iterations for k_true = 20
#   Worst case at batch_size = 1, max_iter = 50:  ~200ms
#   At batch_size = 5, max_iter = 15:              ~60ms
# ==============================================================================

#' Iterative Peel MIS Detection
#'
#' Removes influential observations iteratively, refitting the robust model
#' after each step, until no more outliers are detected or removal begins
#' to harm the estimate.
#'
#' @param formula  A formula (e.g., y ~ x).
#' @param data     A data.frame containing the variables.
#' @param target_var Character; name of the target coefficient (e.g., "x").
#' @param target_pos Integer; position of the target coefficient in the
#'        design matrix (default = 2, first slope with intercept).
#' @param batch_size Integer; number of points to peel per iteration.
#'        1 = safest (no over-removal risk), 5 = faster.
#'        Default = 1.
#' @param max_iter Integer; maximum peel iterations. Default = 50.
#'        Total observations removed is capped at batch_size * max_iter.
#' @param max_k_frac Numeric; hard cap on total fraction of n removed.
#'        Default = 0.06 (6%, slightly above the 5% injection rate).
#' @param detector Character; which MIS detector to use.
#'        "dinkelbach" (exact, default) or "greedy" (fast_sens_topk).
#' @param k_method Character; which k-selector to use at each refit step.
#'        "leverage" (default), "alpha", or "any" (either detects > 0).
#' @param verbose Logical; print iteration trace. Default = FALSE.
#'
#' @return A list with components:
#'   \item{excluded}{Integer vector of all removed observation indices
#'                   (in original data row numbering).}
#'   \item{k_total}{Integer; total number removed.}
#'   \item{n_iters}{Integer; number of peel iterations performed.}
#'   \item{stop_reason}{Character; why peeling stopped.}
#'   \item{beta_trajectory}{Numeric vector; coefficient after each peel step.}
#' @export
iterative_peel <- function(formula, data,
                           target_var  = "x",
                           target_pos  = 2L,
                           batch_size  = 1L,
                           max_iter    = 50L,
                           max_k_frac  = 0.06,
                           detector    = "dinkelbach",
                           k_method    = "leverage",
                           verbose     = FALSE) {

  n_total   <- nrow(data)
  max_k_abs <- floor(n_total * max_k_frac)

  # ------------------------------------------------------------------
  # 0. Initial robust fit → anchor coefficient
  # ------------------------------------------------------------------
  mod_mm_init <- tryCatch(
    robustbase::lmrob(formula, data = data, setting = "KS2014"),
    error = function(e) NULL
  )
  if (is.null(mod_mm_init)) {
    return(list(excluded = integer(0), k_total = 0L, n_iters = 0L,
                stop_reason = "mm_failed", beta_trajectory = numeric(0)))
  }

  beta_anchor <- tryCatch(
    unname(coef(mod_mm_init)[target_var]),
    error = function(e) NA_real_
  )
  if (is.na(beta_anchor)) {
    return(list(excluded = integer(0), k_total = 0L, n_iters = 0L,
                stop_reason = "anchor_na", beta_trajectory = numeric(0)))
  }

  # ------------------------------------------------------------------
  # State: track which original rows have been excluded
  # ------------------------------------------------------------------
  excluded_all <- integer(0)
  beta_traj    <- numeric(0)
  stop_reason  <- "max_iter"

  # Current OLS on full data → baseline distance from anchor
  mod_ols_current <- stats::lm(formula, data = data)
  beta_current    <- unname(coef(mod_ols_current)[target_var])
  dist_current    <- abs(beta_current - beta_anchor)

  for (it in seq_len(max_iter)) {

    # Check hard cap
    if (length(excluded_all) >= max_k_abs) {
      stop_reason <- "max_k_reached"
      break
    }

    # ----------------------------------------------------------------
    # A. Refit robust model on current clean data
    # ----------------------------------------------------------------
    if (length(excluded_all) > 0) {
      clean_data <- data[-excluded_all, , drop = FALSE]
    } else {
      clean_data <- data
    }

    mod_mm <- tryCatch(
      robustbase::lmrob(formula, data = clean_data, setting = "KS2014"),
      error = function(e) NULL
    )

    # ----------------------------------------------------------------
    # B. Check if k-selector still detects outliers
    # ----------------------------------------------------------------
    k_remaining <- 0L
    if (!is.null(mod_mm)) {
      k_remaining <- switch(k_method,
        "leverage" = leverage_k(mod_mm, target_pos = target_pos, max_k_frac = max_k_frac),
        "alpha"    = alpha_k(mod_mm, max_k_frac = max_k_frac),
        "any"      = max(leverage_k(mod_mm, target_pos = target_pos, max_k_frac = max_k_frac),
                         alpha_k(mod_mm, max_k_frac = max_k_frac)),
        0L
      )
    }

    if (k_remaining == 0L) {
      stop_reason <- "no_outliers_detected"
      break
    }

    # ----------------------------------------------------------------
    # C. Find the next batch of influential points
    #    (operate on full-data OLS, but using only current clean indices)
    # ----------------------------------------------------------------
    mod_ols_clean <- stats::lm(formula, data = clean_data)
    n_clean <- nrow(clean_data)

    # How many to peel this iteration
    this_batch <- min(batch_size, k_remaining,
                      max_k_abs - length(excluded_all))
    if (this_batch < 1L) {
      stop_reason <- "max_k_reached"
      break
    }

    # Detect on the clean-data OLS
    if (detector == "dinkelbach") {
      idx_pos <- dinkelbach_topk_lm(mod_ols_clean, pos = target_pos,
                                     sign =  1, k = this_batch)
      idx_neg <- dinkelbach_topk_lm(mod_ols_clean, pos = target_pos,
                                     sign = -1, k = this_batch)
    } else {
      idx_pos <- fast_sens_topk(mod_ols_clean, pos = target_pos,
                                sign =  1, k = this_batch)
      idx_neg <- fast_sens_topk(mod_ols_clean, pos = target_pos,
                                sign = -1, k = this_batch)
    }

    # Map clean-data row indices back to original-data row indices
    if (length(excluded_all) > 0) {
      original_rows <- seq_len(n_total)[-excluded_all]
    } else {
      original_rows <- seq_len(n_total)
    }
    orig_idx_pos <- original_rows[idx_pos]
    orig_idx_neg <- original_rows[idx_neg]

    # ----------------------------------------------------------------
    # D. Direction selection: which direction moves toward anchor?
    # ----------------------------------------------------------------
    r_pos <- fit_clean_ols(formula, data = data,
                           exclude_idx = c(excluded_all, orig_idx_pos))
    r_neg <- fit_clean_ols(formula, data = data,
                           exclude_idx = c(excluded_all, orig_idx_neg))

    d_pos <- abs(unname(r_pos["coef"]) - beta_anchor)
    d_neg <- abs(unname(r_neg["coef"]) - beta_anchor)

    # Pick direction closer to anchor
    if (is.na(d_pos) && is.na(d_neg)) {
      stop_reason <- "fit_failed"
      break
    }
    if (is.na(d_pos)) {
      best_idx <- orig_idx_neg; best_res <- r_neg; d_new <- d_neg
    } else if (is.na(d_neg)) {
      best_idx <- orig_idx_pos; best_res <- r_pos; d_new <- d_pos
    } else if (d_pos <= d_neg) {
      best_idx <- orig_idx_pos; best_res <- r_pos; d_new <- d_pos
    } else {
      best_idx <- orig_idx_neg; best_res <- r_neg; d_new <- d_neg
    }

    # ----------------------------------------------------------------
    # E. Direction stop: if removal moves AWAY from anchor, stop
    # ----------------------------------------------------------------
    if (d_new > dist_current + 1e-10) {
      stop_reason <- "direction_reversal"
      break
    }

    # ----------------------------------------------------------------
    # F. Accept this peel step
    # ----------------------------------------------------------------
    excluded_all <- c(excluded_all, best_idx)
    beta_new     <- unname(best_res["coef"])
    beta_traj    <- c(beta_traj, beta_new)
    dist_current <- d_new

    if (verbose) {
      cat(sprintf("  Peel %02d: removed %d obs, k_total=%d, beta=%.4f, dist=%.4f\n",
                  it, this_batch, length(excluded_all), beta_new, d_new))
    }
  }

  return(list(
    excluded       = excluded_all,
    k_total        = length(excluded_all),
    n_iters        = min(it, max_iter),
    stop_reason    = stop_reason,
    beta_trajectory = beta_traj
  ))
}


#' Convenience Wrapper: Iterative Peel → Clean OLS Result
#'
#' Runs iterative_peel and returns the clean OLS coefficient + SE in the
#' same format as fit_clean_ols, for direct comparison in simulations.
#'
#' @param formula   A formula.
#' @param data      A data.frame.
#' @param target_var Character; target coefficient name.
#' @param target_pos Integer; target coefficient position.
#' @param batch_size Integer; peel batch size (default = 1).
#' @param max_iter   Integer; max peel iterations (default = 50).
#' @param max_k_frac Numeric; max removal fraction (default = 0.06).
#' @param detector   Character; "dinkelbach" or "greedy".
#' @param k_method   Character; "leverage", "alpha", or "any".
#'
#' @return A named numeric vector c(coef = ..., se = ...) matching
#'         the output format of fit_clean_ols.
#' @export
peel_and_fit <- function(formula, data,
                         target_var  = "x",
                         target_pos  = 2L,
                         batch_size  = 1L,
                         max_iter    = 50L,
                         max_k_frac  = 0.06,
                         detector    = "dinkelbach",
                         k_method    = "leverage") {

  peel_result <- iterative_peel(
    formula     = formula,
    data        = data,
    target_var  = target_var,
    target_pos  = target_pos,
    batch_size  = batch_size,
    max_iter    = max_iter,
    max_k_frac  = max_k_frac,
    detector    = detector,
    k_method    = k_method
  )

  fit_clean_ols(formula, data = data, exclude_idx = peel_result$excluded)
}
