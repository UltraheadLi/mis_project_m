# ==============================================================================
# File: /R/influence_injector.R
# Purpose: Injects controlled adversarial influence into targeted datasets.
#          Implements scale-aware transformations (+/- mu, * mu, / phi) 
#          to avoid "shrinking-to-zero" bugs when E[y] = 0, and tags the 
#          ground-truth indices for exact MIS algorithm validation.
# ==============================================================================
#'
#' @param data A list containing `X`, `y`, `error`, and `true_beta` (from DGP)
#' @param method Character: "vertical_outlier", "good_leverage", "bad_leverage"
#' @param k Integer: The size of the influential set to inject (default = 1)
#' @param magnitude Numeric: The multiplier for the shift (must be > 0, default = 5)
#'
#' @return The modified data list, including exact injection metadata.
#' @export
apply_influence_shift <- function(data, method = "vertical_outlier", k = 1, magnitude = 5) {
  
  # Defensive Programming Checks
  if (is.null(data$true_beta)) stop("`data` must contain `true_beta` from the DGP.")
  if (magnitude <= 0) stop("Magnitude must be strictly positive.")
  if (k < 1 || k > length(data$y)) stop("k must be between 1 and n.")
  
  y <- data$y
  X <- data$X
  n <- length(y)
  
  outlier_idx <- sample(seq_len(n), size = k)
  
  # Robust scaling for Y
  scale_y <- max(mad(y, constant = 1.4826), 1e-4)
  
  y_new <- y
  X_new <- X
  
  # Initialize tracking directions to avoid NA scope issues
  y_shift_dir <- NA
  x_shift_dir <- NA
  
  switch(method,
         
         "vertical_outlier" = {
           y_shift_dir <- sample(c(1, -1), size = 1)
           y_new[outlier_idx] <- y[outlier_idx] + (y_shift_dir * magnitude * scale_y)
         },
         
         "good_leverage" = {
           # Robust scaling for X
           scale_x <- apply(X, 2, function(col) max(mad(col, constant = 1.4826), 1e-4))
           shift_mat <- matrix(magnitude * scale_x, nrow = k, ncol = ncol(X), byrow = TRUE)
           
           # Force the good leverage coalition to move together on the X-axis
           x_shift_dir <- sample(c(1, -1), size = 1)
           sign_mat <- matrix(x_shift_dir, nrow = k, ncol = ncol(X))
           
           X_new[outlier_idx, ] <- X[outlier_idx, , drop = FALSE] + (shift_mat * sign_mat)
           
           # Realistic good leverage: exact plane + small noise
           true_y <- as.vector(X_new[outlier_idx, , drop = FALSE] %*% data$true_beta)
           y_new[outlier_idx] <- true_y + rnorm(k, mean = 0, sd = scale_y * 0.1)
         },
         
         "bad_leverage" = {
           # Robust scaling for X
           scale_x <- apply(X, 2, function(col) max(mad(col, constant = 1.4826), 1e-4))
           shift_mat <- matrix(magnitude * scale_x, nrow = k, ncol = ncol(X), byrow = TRUE)
           
           # FUNDAMENTAL FIX: Force the entire coalition to shift X in the SAME direction
           x_shift_dir <- sample(c(1, -1), size = 1)
           sign_mat <- matrix(x_shift_dir, nrow = k, ncol = ncol(X))
           
           X_new[outlier_idx, ] <- X[outlier_idx, , drop = FALSE] + (shift_mat * sign_mat)
           
           true_y <- as.vector(X_new[outlier_idx, , drop = FALSE] %*% data$true_beta)
           
           # Shift Y uniformly to twist the slope maximally
           y_shift_dir <- sample(c(1, -1), size = 1)
           severe_shift <- max(scale_y, abs(true_y) * 0.5) * magnitude
           y_new[outlier_idx] <- true_y + (y_shift_dir * severe_shift)
         },
         
         stop(sprintf("Influence method '%s' not recognized.", method))
  )
  
  data$y <- y_new
  data$X <- X_new
  data$outlier_indices <- outlier_idx
  data$injection_method <- method
  data$injection_magnitude <- magnitude
  data$injection_k <- k
  
  # Clean, scalar metadata tracking
  data$injection_x_direction <- x_shift_dir
  data$injection_y_direction <- y_shift_dir
  
  return(data)
}