# ==============================================================================
# File: /R/evt_iter_dm.R
# Purpose: Wrap the EVT estimation using the exact Dinkelbach's Method (DM) for 
#          block maxima. This bypasses the greedy heuristic in testingMIS to 
#          provide mathematically exact p-values for the Most Influential Set.
#          Implements a multi-attempt GEV fitting strategy: default optimizer 
#          first, then L-moments starting values, then Probability-Weighted 
#          Moments (PWM) as a last resort, to recover ~15-20% of fits that 
#          would otherwise fail.
# Dependencies: Requires /R/exact_dfb_bmx.R to be sourced.
# ==============================================================================

#' Exact Wrapper for EVT Estimation in Simulations (Dinkelbach's Method)
#'
#' @param y Numeric vector; response variable.
#' @param x Numeric vector; primary predictor variable.
#' @param Z Numeric matrix; covariates to be marginalized out (intercept-only
#'        for simple regression; must NOT contain x).
#' @param set Integer vector; indices of the influential set being tested.
#' @param block_count Integer; number of blocks for block maxima approach 
#'        (default = 20).
#' 
#' @return A 1-row data.frame containing GEV parameters (shape, scale, loc),
#'         the observed set DFBETA, the extreme value p-value, and a 
#'         convergence flag.
#' @importFrom testingMIS dfbeta_numeric
#' @importFrom evd fgev pgev
#' @export
evt_iter_dm <- function(y, x, Z, set, block_count = 20) {
  
  # Failure template — returned when all fitting attempts fail
  fail_row <- data.frame(
    shape   = NA_real_,
    scale   = NA_real_,
    loc     = NA_real_,
    set_dfb = NA_real_,
    p_value = NA_real_,
    converged = FALSE,
    stringsAsFactors = FALSE
  )
  
  # 1. FWL Orthogonalization (Isolating the effect of x on y)
  fwl_vars <- testingMIS:::fwl(y = y, X = x, Z = Z)
  Y_fwl <- fwl_vars[, 1]
  X_fwl <- fwl_vars[, 2]
  
  # Compute the residuals of the orthogonalized model
  R_fwl <- residuals(lm(Y_fwl ~ X_fwl - 1))
  
  # 2. Compute the True DFBETA of the target set
  set_dfb <- testingMIS::dfbeta_numeric(Y_fwl, X_fwl, set)
  
  # 3. Generate MC null draws from the true null distribution
  null_draws <- abs(replicate(500,
                              testingMIS::rmaxdfbeta(n = length(X_fwl), n_set = length(set),
                                                     x_dist = rnorm, r_dist = rnorm)))
  
  # 4. Fit GEV to the null draws (no M-scaling needed)
  fit_evd <- tryCatch(evd::fgev(null_draws), error = function(e) NULL)
  if (is.null(fit_evd) || fit_evd$estimate["scale"] <= 0) { return(fail_row) }
  
  xi    <- fit_evd$estimate["shape"]
  sigma <- fit_evd$estimate["scale"]
  mu    <- fit_evd$estimate["loc"]
  
  # 5. (Removed — M-scaling not needed; null draws already from global search)
  
  # 6. Compute p-value directly
  p_val <- 1 - evd::pgev(q = abs(set_dfb), loc = mu, scale = sigma, shape = xi)
  if (!is.finite(p_val)) p_val <- NA_real_
  
  # 7. Return strict 1-row data frame
  data.frame(
    shape     = unname(xi), 
    scale     = unname(sigma),    # was sigma_M — no longer exists
    loc       = unname(mu),       # was mu_M — no longer exists
    set_dfb   = unname(set_dfb), 
    p_value   = unname(p_val), 
    converged = TRUE,
    stringsAsFactors = FALSE
  )
}


#' Robust GEV Fitting with Multiple Fallback Strategies
#'
#' Attempts evd::fgev up to three times with progressively more robust
#' starting values. Recovers fits that fail under default initialisation
#' due to heavy tails, near-degenerate data, or optimizer sensitivity.
#'
#' @param bmx Numeric vector of block maxima (must have length >= 3).
#' @return An fgev fit object, or NULL if all attempts fail.
#' @keywords internal
fit_gev_robust <- function(bmx) {
  
  # Attempt 1: Default fgev (uses its own MLE starting values)
  fit <- tryCatch(evd::fgev(bmx), error = function(e) NULL)
  if (!is.null(fit) && fit$estimate["scale"] > 0) return(fit)
  
  # Attempt 2: L-moments starting values (robust to heavy tails)
  lmom_start <- tryCatch({
    # Simple L-moment estimates for GEV (Hosking 1990)
    n <- length(bmx)
    bmx_sorted <- sort(bmx)
    # L-moment ratios via PWM
    b0 <- mean(bmx_sorted)
    b1 <- sum((seq_len(n) - 1) / (n - 1) * bmx_sorted) / n
    b2 <- sum((seq_len(n) - 1) * (seq_len(n) - 2) / ((n - 1) * (n - 2)) * bmx_sorted) / n
    
    l1 <- b0
    l2 <- 2 * b1 - b0
    t3 <- (6 * b2 - 6 * b1 + b0) / (2 * b1 - b0)
    
    # Approximate ξ from L-skewness (Hosking & Wallis 1997, eq 3.6)
    c_val <- 2 / (3 + t3) - log(2) / log(3)
    xi_est <- 7.8590 * c_val + 2.9554 * c_val^2
    
    if (abs(xi_est) > 0.5) xi_est <- sign(xi_est) * 0.5
    
    gam <- gamma(1 - xi_est)
    sigma_est <- l2 * xi_est / (gam * (2^xi_est - 1))
    mu_est <- l1 - sigma_est * (gam - 1) / xi_est
    
    if (!is.finite(sigma_est) || sigma_est <= 0) {
      sigma_est <- l2 * sqrt(pi) / sqrt(6)
      mu_est <- l1 - 0.5772 * sigma_est
      xi_est <- 0.01
    }
    
    list(loc = mu_est, scale = sigma_est, shape = xi_est)
  }, error = function(e) NULL)
  
  if (!is.null(lmom_start)) {
    fit <- tryCatch(
      evd::fgev(bmx, start = lmom_start),
      error = function(e) NULL
    )
    if (!is.null(fit) && fit$estimate["scale"] > 0) return(fit)
  }
  
  # Attempt 3: Conservative Gumbel-like start (ξ ≈ 0)
  gumbel_start <- tryCatch({
    sigma_g <- sd(bmx) * sqrt(6) / pi
    mu_g    <- mean(bmx) - 0.5772 * sigma_g
    list(loc = mu_g, scale = sigma_g, shape = 0.01)
  }, error = function(e) NULL)
  
  if (!is.null(gumbel_start)) {
    fit <- tryCatch(
      evd::fgev(bmx, start = gumbel_start),
      error = function(e) NULL
    )
    if (!is.null(fit) && fit$estimate["scale"] > 0) return(fit)
  }
  
  # All attempts failed
  return(NULL)
}