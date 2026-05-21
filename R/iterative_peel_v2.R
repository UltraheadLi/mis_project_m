# ==============================================================================
# File: /R/iterative_peel_v2.R
# Purpose: Anchor-free iterative peel MIS detection. Replaces the fixed
#          MM-estimator anchor with a self-contained direction rule based on
#          residual standard error (sigma) reduction.
#
# Key change from v1:
#   v1: Direction is chosen by proximity to a pre-computed MM beta_anchor.
#       This creates a circular dependency — MIS's quality depends on the
#       MM estimator's quality, which can itself be compromised by masking.
#
#   v2: Direction is chosen by which removal most REDUCES sigma_hat (the
#       residual standard error of the clean OLS). No external anchor needed.
#       The logic: genuine outlier removal tightens the fit (sigma drops);
#       signal removal loosens it (sigma rises or stagnates). This provides
#       both a direction rule AND a natural stopping criterion.
#
# Stopping rules (whichever triggers first):
#   1. Sigma-stagnation stop: neither direction reduces sigma by more than
#      `sigma_tol` fraction — we've run out of outliers to peel.
#   2. No-detection stop: the k-selector returns 0 on the cleaned data.
#   3. Max-iterations cap: hard ceiling to prevent runaway peeling.
#   4. Max-k-fraction cap: total removals capped at max_k_frac of n.
#
# Performance: identical to v1 per iteration (1 lmrob + 1 dinkelbach + 1 OLS)
# ==============================================================================

#' Anchor-Free Iterative Peel MIS Detection (v2)
#'
#' Removes influential observations iteratively using residual standard error
#' reduction as the direction criterion, rather than an external robust anchor.
#'
#' @param formula    A formula (e.g., y ~ x).
#' @param data       A data.frame containing the variables.
#' @param target_var Character; name of the target coefficient (e.g., "x").
#' @param target_pos Integer; position of the target coefficient in the
#'        design matrix (default = 2, first slope with intercept).
#' @param batch_size Integer; number of points to peel per iteration.
#'        Default = 1.
#' @param max_iter   Integer; maximum peel iterations. Default = 50.
#' @param max_k_frac Numeric; hard cap on total fraction of n removed.
#'        Default = 0.06 (6%).
#' @param sigma_tol  Numeric; minimum fractional reduction in sigma to
#'        accept a peel step. Default = 1e-4 (0.01% improvement required).
#'        Set to 0 to disable (rely only on k-selector stopping).
#' @param detector   Character; which MIS detector to use.
#'        "dinkelbach" (exact, default) or "greedy" (fast_sens_topk).
#' @param k_method   Character; which k-selector to use at each refit step.
#'        "leverage" (default), "alpha", or "any".
#' @param verbose    Logical; print iteration trace. Default = FALSE.
#'
#' @return A list with components:
#'   \item{excluded}{Integer vector of all removed observation indices.}
#'   \item{k_total}{Integer; total number removed.}
#'   \item{n_iters}{Integer; number of peel iterations performed.}
#'   \item{stop_reason}{Character; why peeling stopped.}
#'   \item{beta_trajectory}{Numeric vector; coefficient after each peel step.}
#'   \item{sigma_trajectory}{Numeric vector; sigma_hat after each peel step.}
#' @export
iterative_peel_v2 <- function(formula, data,
                              target_var  = "x",
                              target_pos  = 2L,
                              batch_size  = 1L,
                              max_iter    = 50L,
                              max_k_frac  = 0.06,
                              sigma_tol   = 1e-4,
                              detector    = "dinkelbach",
                              k_method    = "leverage",
                              verbose     = FALSE) {

  n_total   <- nrow(data)
  max_k_abs <- floor(n_total * max_k_frac)

  # ------------------------------------------------------------------
  # 0. Initial state: OLS on full data → baseline sigma
  #    No MM anchor needed.
  # ------------------------------------------------------------------
  mod_ols_init   <- stats::lm(formula, data = data)
  sigma_current  <- summary(mod_ols_init)$sigma

  # We still need lmrob for the k-selector (leverage_k / alpha_k),
  # but NOT for direction selection.
  mod_mm_init <- tryCatch(
    robustbase::lmrob(formula, data = data, setting = "KS2014"),
    error = function(e) NULL
  )

  # ------------------------------------------------------------------
  # State tracking
  # ------------------------------------------------------------------
  excluded_all    <- integer(0)
  beta_traj       <- numeric(0)
  sigma_traj      <- numeric(0)
  stop_reason     <- "max_iter"
  it              <- 0L

  for (it in seq_len(max_iter)) {

    # Check hard cap
    if (length(excluded_all) >= max_k_abs) {
      stop_reason <- "max_k_reached"
      break
    }

    # ----------------------------------------------------------------
    # A. Refit robust model on current clean data (for k-selector only)
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
        "leverage" = leverage_k(mod_mm, target_pos = target_pos,
                                max_k_frac = max_k_frac),
        "alpha"    = alpha_k(mod_mm, max_k_frac = max_k_frac),
        "any"      = max(
          leverage_k(mod_mm, target_pos = target_pos, max_k_frac = max_k_frac),
          alpha_k(mod_mm, max_k_frac = max_k_frac)
        ),
        0L
      )
    }

    if (k_remaining == 0L) {
      stop_reason <- "no_outliers_detected"
      break
    }

    # ----------------------------------------------------------------
    # C. Find the next batch of influential points
    # ----------------------------------------------------------------
    mod_ols_clean <- stats::lm(formula, data = clean_data)
    n_clean <- nrow(clean_data)

    this_batch <- min(batch_size, k_remaining,
                      max_k_abs - length(excluded_all))
    if (this_batch < 1L) {
      stop_reason <- "max_k_reached"
      break
    }

    # Detect in both directions on the clean-data OLS
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

    # Map clean-data indices to original-data indices
    if (length(excluded_all) > 0) {
      original_rows <- seq_len(n_total)[-excluded_all]
    } else {
      original_rows <- seq_len(n_total)
    }
    orig_idx_pos <- original_rows[idx_pos]
    orig_idx_neg <- original_rows[idx_neg]

    # ----------------------------------------------------------------
    # D. Direction selection via SIGMA REDUCTION (no anchor)
    #
    #    Fit OLS excluding each candidate set. Whichever produces the
    #    lower sigma_hat is the correct direction: genuine outlier
    #    removal tightens the fit, signal removal doesn't.
    # ----------------------------------------------------------------
    r_pos <- fit_clean_ols(formula, data = data,
                           exclude_idx = c(excluded_all, orig_idx_pos))
    r_neg <- fit_clean_ols(formula, data = data,
                           exclude_idx = c(excluded_all, orig_idx_neg))

    # Compute sigma for each direction
    excl_pos_full <- c(excluded_all, orig_idx_pos)
    excl_neg_full <- c(excluded_all, orig_idx_neg)

    mod_pos <- stats::lm(formula, data = data[-excl_pos_full, , drop = FALSE])
    mod_neg <- stats::lm(formula, data = data[-excl_neg_full, , drop = FALSE])

    sigma_pos <- tryCatch(summary(mod_pos)$sigma, error = function(e) Inf)
    sigma_neg <- tryCatch(summary(mod_neg)$sigma, error = function(e) Inf)

    # Pick the direction that reduces sigma most
    if (is.infinite(sigma_pos) && is.infinite(sigma_neg)) {
      stop_reason <- "fit_failed"
      break
    }

    if (sigma_pos <= sigma_neg) {
      best_idx   <- orig_idx_pos
      best_res   <- r_pos
      sigma_new  <- sigma_pos
    } else {
      best_idx   <- orig_idx_neg
      best_res   <- r_neg
      sigma_new  <- sigma_neg
    }

    # ----------------------------------------------------------------
    # E. Sigma-stagnation stop: if sigma didn't improve meaningfully,
    #    we're removing signal, not contamination.
    # ----------------------------------------------------------------
    sigma_improvement <- (sigma_current - sigma_new) / sigma_current

    if (sigma_improvement < sigma_tol) {
      stop_reason <- "sigma_stagnation"
      break
    }

    # ----------------------------------------------------------------
    # F. Accept this peel step
    # ----------------------------------------------------------------
    excluded_all  <- c(excluded_all, best_idx)
    beta_new      <- unname(best_res["coef"])
    beta_traj     <- c(beta_traj, beta_new)
    sigma_traj    <- c(sigma_traj, sigma_new)
    sigma_current <- sigma_new

    if (verbose) {
      cat(sprintf(
        "  Peel %02d: removed %d obs, k_total=%d, beta=%.4f, sigma=%.6f (improvement=%.4f%%)\n",
        it, this_batch, length(excluded_all), beta_new, sigma_new,
        sigma_improvement * 100
      ))
    }
  }

  return(list(
    excluded         = excluded_all,
    k_total          = length(excluded_all),
    n_iters          = min(it, max_iter),
    stop_reason      = stop_reason,
    beta_trajectory  = beta_traj,
    sigma_trajectory = sigma_traj
  ))
}
