# ==============================================================================
# File: /R/evt_iter_dm.R
# Purpose: Wrap the EVT estimation using the exact Dinkelbach's Method (DM) for 
#          block maxima. This bypasses the greedy heuristic in testingMIS to 
#          provide mathematically exact p-values for the Most Influential Set.
# Dependencies: Requires /R/exact_dfb_bmx.R to be sourced.
# ==============================================================================

#' Exact Wrapper for EVT Estimation in Simulations (Dinkelbach's Method)
#'
#' @param y Numeric vector; response variable.
#' @param x Numeric vector; primary predictor variable.
#' @param Z Numeric matrix; covariates to be marginalized out.
#' @param set Integer vector; indices of the influential set being tested.
#' @param block_count Integer; number of blocks for block maxima approach.
#' 
#' @return A 1-row data.frame containing GEV parameters, the observed set DFBETA,
#'         the extreme value p-value, and a convergence flag.
#' @importFrom testingMIS dfbeta_numeric
#' @importFrom evd fgev pgev
#' @export
evt_iter_dm <- function(y, x, Z, set, block_count = 20) {
  
  # 1. FWL Orthogonalization (Isolating the effect of x on y)
  # We use the internal fwl function from testingMIS to ensure exact data matching
  fwl_vars <- testingMIS:::fwl(y = y, X = x, Z = Z)
  Y_fwl <- fwl_vars[, 1]
  X_fwl <- fwl_vars[, 2]
  
  # Compute the residuals of the orthogonalized model
  R_fwl <- residuals(lm(Y_fwl ~ X_fwl - 1))
  
  # 2. Compute the True DFBETA of the target set
  set_dfb <- testingMIS::dfbeta_numeric(Y_fwl, X_fwl, set)
  
  # 3. Generate mathematically EXACT block maxima
  # This uses your standalone Dinkelbach function rather than the package's greedy one
  Delta_bmx <- abs(exact_dfb_bmx(X_fwl, R_fwl, set, block_count))
  
  # 4. Fit the Extreme Value Distribution (EVD) safely
  fit_evd <- tryCatch({
    evd::fgev(Delta_bmx)
  }, error = function(e) return(NULL))
  
  # Handle failure states (optimizer crash or degenerate scale)
  if (is.null(fit_evd) || fit_evd$estimate["scale"] <= 0) {
    return(data.frame(
      shape   = NA_real_,
      scale   = NA_real_,
      loc     = NA_real_,
      set_dfb = NA_real_,
      p_value = NA_real_,
      converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  xi    <- fit_evd$estimate["shape"]
  sigma <- fit_evd$estimate["scale"]
  mu    <- fit_evd$estimate["loc"]
  
  # 5. EVT Scaling Correction 
  # M is simply the number of independent blocks used to build the baseline
  M <- block_count
  
  if (abs(xi) > 1e-6) {
    # Fréchet / Weibull domain
    mu_M    <- mu + sigma * (M^xi - 1) / xi
    sigma_M <- sigma * M^xi
  } else {
    # Gumbel limit
    mu_M    <- mu + sigma * log(M)
    sigma_M <- sigma
  }
  
  # 6. Compute exact P-Value
  p_val <- 1 - evd::pgev(
    q     = abs(set_dfb), 
    loc   = mu_M, 
    scale = sigma_M, 
    shape = xi
  )
  
  # 7. Return strict 1-row data frame
  data.frame(
    shape     = unname(xi), 
    scale     = unname(sigma_M), 
    loc       = unname(mu_M),
    set_dfb   = unname(set_dfb), 
    p_value   = unname(p_val), 
    converged = TRUE,
    stringsAsFactors = FALSE
  )
}