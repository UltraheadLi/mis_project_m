# ==============================================================================
# File: /R/diagnostics_classical.R
# Purpose: Compute and standardize classical leave-one-out (LOO) influence 
#          diagnostics (Leverage, Cook's D, DFBETAS) to serve as a baseline 
#          comparison against exact Most Influential Sets (MIS).
# ==============================================================================

#' 1. Compute Leverage (Diagonal Values in the Hat)
#'
#' @param model A fitted linear model object (e.g., from `lm()`).
#' @return A numeric vector of leverage values for each observation.
get_leverage <- function(model) {
  if (!inherits(model, "lm")) stop("Input must be a fitted 'lm' object.")
  stats::hatvalues(model)
}

#' 2. Compute Cook's Distance
#'
#' @param model A fitted linear model object.
#' @return A numeric vector of Cook's distances.
get_cooks_d <- function(model) {
  if (!inherits(model, "lm")) stop("Input must be a fitted 'lm' object.")
  stats::cooks.distance(model)
}

#' 3. Compute DFBETAS for a specific target variable
#'
#' @param model A fitted linear model object.
#' @param target_var Character string specifying the name of the coefficient to evaluate.
#' @return A numeric vector of DFBETAS for the target variable.
get_dfbetas <- function(model, target_var) {
  if (!inherits(model, "lm")) stop("Input must be a fitted 'lm' object.")
  dfb <- stats::dfbetas(model)
  if (!target_var %in% colnames(dfb)) {
    stop(sprintf("Target variable '%s' not found in model coefficients. Available are: %s", 
                 target_var, paste(colnames(dfb), collapse = ", ")))
  }
  dfb[, target_var]
}

#' 4. Aggregate the 3 classical diagnostics into a single data frame
#'
#' @param model A fitted linear model object.
#' @param target_var Character string specifying the coefficient of interest.
#' @return A data.frame containing columns: id, leverage, cooks_d, dfbetas.
get_all_classical <- function(model, target_var) {
  lev    <- get_leverage(model)
  cooksd <- get_cooks_d(model)
  dfb    <- get_dfbetas(model, target_var)
  
  out <- data.frame(
    id = seq_along(lev),
    leverage = lev,
    cooks_d = cooksd,
    dfbetas_target = dfb
  )
  
  # Strip rownames to prevent visual clutter in the output table
  rownames(out) <- NULL
  return(out)
}

#' 5. Extract classical influential set (Top K or Threshold-based)
#'
#' Extracts an influential set based on classical leave-one-out diagnostics.
#' If 'k' is provided, it returns the top K most influential points based on 
#' absolute magnitude. If 'k' is NULL, it dynamically identifies points that 
#' exceed standard statistical thresholds for the chosen metric.
#'
#' @param model A fitted linear model object.
#' @param target_var Character string specifying the coefficient of interest.
#' @param k Integer or NULL. The size of the influential set to extract. 
#'          If NULL (the default), statistical thresholds are used instead.
#' @param metric Character string. Which classical metric to evaluate.
#'               Must be one of "cooks_d", "leverage", or "dfbetas_target".
#' @return A numeric vector of the original observation IDs belonging to the set.
get_classical_set <- function(model, target_var, k = NULL, metric = "cooks_d") {
  
  ## 1. Get raw metrics and validate the chosen metric
  raw_metrics <- get_all_classical(model, target_var)
  available <- setdiff(colnames(raw_metrics), "id")
  
  if (!metric %in% available) {
    stop(sprintf("Metric '%s' not found. Choose from: %s", 
                 metric, paste(available, collapse = ", ")))
  }
  
  # =====================================================================
  # PATH A: Top 'K' Heuristic (If user specifies k)
  # =====================================================================
  if (!is.null(k)) {
    # Sort by absolute magnitude to account for negative DFBETAS
    sorted_data <- raw_metrics[order(abs(raw_metrics[[metric]]), decreasing = TRUE), ]
    
    # Return top K IDs safely
    return(sorted_data$id[1:min(k, nrow(sorted_data))])
  }
  
  # =====================================================================
  # PATH B: Statistical Thresholds (Theoritical Default Setting)
  # =====================================================================
  n <- stats::nobs(model)
  p <- length(stats::coef(model))
  
  if (metric == "cooks_d") {
    # Threshold for Cook's Distance: 4 / n
    threshold <- 4 / n
    flagged_ids <- raw_metrics$id[raw_metrics$cooks_d > threshold]
    
  } else if (metric == "leverage") {
    # Threshold for Leverage (Hat Values): 2p / n
    threshold <- 2 * p / n
    flagged_ids <- raw_metrics$id[raw_metrics$leverage > threshold]
    
  } else if (metric == "dfbetas_target") {
    # Threshold for DFBETAS: 2 / sqrt(n)
    threshold <- 2 / sqrt(n)
    flagged_ids <- raw_metrics$id[abs(raw_metrics$dfbetas_target) > threshold]
  }
  
  return(flagged_ids)
}