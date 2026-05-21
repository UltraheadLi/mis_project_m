# ==============================================================================
# File: /R/evt_iter.R
# Purpose: Wrap the EVT estimation function (estimate_dfb_evd) to safely compute 
#          the p-value of observed set influence under the fitted GEV distribution. 
#          Formats output as a flat 1-row data.frame for robust row-binding in 
#          large-scale simulation loops, catching optimization failures cleanly.
# Dependencies: Requires helpers_local.R (provides estimate_dfb_evd)
# ==============================================================================

#' Wrapper for EVT Estimation in Simulations
#'
#' @param y Numeric vector; response variable.
#' @param x Numeric vector; primary predictor variable.
#' @param Z Numeric matrix; covariates to be marginalized out.
#' @param set Integer vector; indices of the influential set being tested.
#' @param block_count Integer; number of blocks for block maxima approach.
#' 
#' @return A 1-row data.frame containing GEV parameters, the observed set DFBETA,
#'         the extreme value p-value, and a convergence flag.
#' @importFrom evd pgev
#' @export
evt_iter <- function(y, x, Z, set, block_count = 20) {
  
  # 1. EVD
  res <- tryCatch({
    estimate_dfb_evd(
      y = y,
      x = x,
      Z = Z,
      set = set,
      block_count = block_count,
      verbose = FALSE 
    )
  }, error = function(e) {
    return(NULL)
  })
  
  # 2. Failure state
  if (is.null(res) || 
      !all(c("loc", "scale", "shape") %in% names(res$params)) || 
      res$params["scale"] <= 0) {
    
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
  
  # 3. Scale GEV parameters for combinatorial search space
  xi    <- res$params["shape"]
  sigma <- res$params["scale"]
  mu    <- res$params["loc"]
  
  M <- block_count
  
  # EVT location-scale scaling: max of M blocks
  if (abs(xi) > 1e-6) {
    mu_M    <- mu + sigma * (M^xi - 1) / xi
    sigma_M <- sigma * M^xi
  } else {
    mu_M    <- mu + sigma * log(M)
    sigma_M <- sigma
  }
  
  # 4. Compute p-value under the corrected distribution
  p_val <- 1 - evd::pgev(
    q     = abs(res$set_dfb),
    loc   = mu_M,
    scale = sigma_M,
    shape = xi
  )
  # 5. Return dataframe
  data.frame(
    shape     = unname(res$params["shape"]),
    scale     = unname(res$params["scale"]),
    loc       = unname(res$params["loc"]),
    set_dfb   = unname(res$set_dfb),
    p_value   = unname(p_val),
    converged = TRUE,
    stringsAsFactors = FALSE
  )
}