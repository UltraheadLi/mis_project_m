# ==============================================================================
# File: /R/exact_dfb_bmx.R
# Purpose: Exact implementation of the block maxima DFBETA search using 
#          Dinkelbach's method for linear-fractional programming. Intended for 
#          direct comparison against the greedy testingMIS heuristic.
# ==============================================================================

#' Exact Set Influence Block Maxima
#'
#' @param X Vector of predictor values (full sample)
#' @param R Vector of residual values (full sample)
#' @param set Vector of indices for the influential observation set
#' @param block_count Number of blocks to divide the data into
#' @return Numeric vector of mathematically exact block maxima DFBETA values
exact_dfb_bmx <- function(X, R, set, block_count) {
  sgn <- sign(sum(X[set] * R[set]))
  if (sgn == 0) stop("dfbeta of set is exactly zero")
  
  # Full sum of squares (constant for denominator)
  sumX2 <- sum(X[-set]^2) 
  
  # Remove the target set from the search pool
  X_inf <- X[-set]
  R_inf <- R[-set]
  
  nS <- length(set)
  block_size <- length(X_inf) %/% block_count
  
  # Access the internal make_blocks function from testingMIS
  Xbl <- testingMIS:::make_blocks(X_inf, block_size)
  Rbl <- testingMIS:::make_blocks(R_inf, block_size)
  
  res <- numeric(block_count)
  
  # Dinkelbach's Method to find exact block maxima
  for (i in seq_len(block_count)) {
    x_bl <- Xbl[, i]
    r_bl <- Rbl[, i]
    
    # Formulate fractional program: sgn * sum(xr) / (sumX2 - sum(x^2))
    n_val <- sgn * (x_bl * r_bl)
    d_val <- - (x_bl^2)
    
    lambda <- 0
    
    for (iter in 1:50) {
      # Linearized objective
      w <- n_val - lambda * d_val
      
      # Select top k elements
      idx <- order(w, decreasing = TRUE)[seq_len(nS)]
      
      # Update lambda
      new_lambda <- sum(n_val[idx]) / (sumX2 + sum(d_val[idx]))
      
      if (abs(new_lambda - lambda) < 1e-9) break
      lambda <- new_lambda
    }
    
    # Reapply sign to match true DFBETA
    res[i] <- sgn * lambda
  }
  
  return(res)
}