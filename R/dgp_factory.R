# ==============================================================================
# File: /R/dgp_factory.R
# Purpose: Generates synthetic data (X, y, error) under various complex 
#          distributions (e.g., Pareto, GOLM, Skewed-t). Separates X and 
#          error generation to isolate leverage vs. variance effects, and 
#          enforces zero-expectation for errors to preserve OLS assumptions.
# ==============================================================================
#'
#' @param n Sample size (integer)
#' @param p Number of predictors (integer, default = 1)
#' @param x_type Distribution for the design matrix X
#' @param error_type Distribution for the error term
#' @param beta True coefficient vector (default = rep(1, p))
#' @param ... Additional arguments
#'
#' @return A list containing X, y, and error
#' @export
generate_complex_data <- function(n, p = 1, x_type = "normal", error_type = "normal", 
                                  beta = rep(1, p), dist_param = 3, mix_prop = 0.1, ...) {
  
  generate_vector <- function(n_vals, dist_type, center = FALSE) {
    vec <- switch(dist_type,
                  "normal" = rnorm(n_vals, 0, 1),
                  # Symmetric, thin-tailed baseline. No extremes, no skew.
                  
                  "mixed_normal" = {
                    out <- numeric(n_vals)
                    is_comp1 <- sample(c(TRUE, FALSE), n_vals, replace = TRUE, 
                                       prob = c(1 - mix_prop, mix_prop))
                    out[is_comp1]  <- rnorm(sum(is_comp1), 0, 1)
                    out[!is_comp1] <- rnorm(sum(!is_comp1), 0, 10)
                    out# 90% standard normal, 10% high-variance normal. Mimics occasional large shocks.
                  },
                  
                  "skewed_t" = {
                    if (!requireNamespace("sn", quietly = TRUE)) stop("Install 'sn' package.")
                    sn::rst(n_vals, xi = 0, omega = 1, alpha = 5, nu = dist_param)
                    # Right-skewed, heavy-tailed. Combines asymmetry with tail weight.
                  },
                  
                  "golm" = {
                    out <- numeric(n_vals)
                    is_comp1 <- sample(c(TRUE, FALSE), n_vals, replace = TRUE, prob = c(0.7, 0.3))
                    out[is_comp1]  <- rlnorm(sum(is_comp1), meanlog = 0, sdlog = 0.5)
                    out[!is_comp1] <- rlnorm(sum(!is_comp1), meanlog = 1, sdlog = 1.5)
                    out
                    # Mixture of two log-normals. Multi-modal, heavily right-skewed.
                  },
                  
                  "beta_logistic" = {
                    rbeta(n_vals, shape1 = 2, shape2 = 5)
                    # Bounded on [0, 1], right-skewed. Models proportions or probabilities.
                  },
                  
                  "gpd" = {
                    if (!requireNamespace("evd", quietly = TRUE)) stop("Install 'evd' package.")
                    evd::rgpd(n_vals, loc = 0, scale = 1, shape = dist_param)
                    # Generalized Pareto with shape > 0. Fréchet domain; unbounded heavy tail.
                  },
                  
                  "contaminated" = {
                    eps <- rnorm(n_vals, 0, 1)
                    outlier_idx <- sample(1:n_vals, size = floor(mix_prop * n_vals))
                    eps[outlier_idx] <- eps[outlier_idx] + rnorm(length(outlier_idx), 0, 50)
                    eps
                    # Standard normal with 5% gross outliers. Classic contamination model.
                  },
                  
                  "pareto" = {
                    if (!requireNamespace("actuar", quietly = TRUE)) stop("Install 'actuar' package.")
                    actuar::rpareto(n_vals, shape = dist_param, scale = 1)
                    # Classical Pareto with shape = 3. Finite variance, but heavier tail than normal.
                  },
                  
                  stop(sprintf("Distribution type '%s' is not supported.", dist_type))
    )
    
    # Apply centering ONLY if explicitly requested
    if (center) {
      vec <- vec - mean(vec)
    }
    
    return(vec)
  }
  
  # 1. Generate X (Never center the design matrix by default)
  X_vec <- generate_vector(n * p, x_type, center = FALSE)
  X <- matrix(X_vec, nrow = n, ncol = p)
  
  # 2. Generate Errors (ALWAYS center to maintain E[eps] = 0 for OLS)
  error <- generate_vector(n, error_type, center = TRUE)
  
  # 3. Construct y using the provided beta
  y <- as.vector(X %*% beta + error)
  
  return(list(X = X, y = y, error = error, true_beta = beta))
}