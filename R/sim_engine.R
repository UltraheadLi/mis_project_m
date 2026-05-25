# ==============================================================================
# File: /R/sim_engine.R
# Purpose: Core simulation wrapper combining DGP, influence injection, MIS 
#          detection, classical LOO diagnostics, coverage comparison, and EVD 
#          parameter estimation. Orchestrates one complete Monte Carlo iteration.
#
# Dependencies: Requires /R/dgp_factory.R, /R/influence_injector.R, 
#               /R/evt_iter_dm.R, /R/exact_dfb_bmx.R, /R/dinkelbach_topk.R
# ==============================================================================
#'
#' @param iter Integer: Iteration ID for tracking across simulation grid
#' @param n Integer: Sample size for the generated dataset
#' @param p Integer: Number of predictors (currently only p = 1 supported)
#' @param x_type Character: Distribution type for the design matrix X 
#'        (from dgp_factory: "normal", "skewed_t", "pareto")
#' @param error_type Character: Distribution type for the error term 
#'        (from dgp_factory: "normal", "mixed_normal", "skewed_t", "golm",
#'        "beta_logistic", "gpd", "contaminated", "pareto")
#' @param outlier_method Character: Injection strategy from influence_injector
#'        ("none", "vertical_outlier", "good_leverage", "bad_leverage")
#' @param k Integer: Size of the influential set to inject and detect (default = 1)
#' @param magnitude Numeric: Scale multiplier for the adversarial shift (default = 5)
#' @param block_count Integer: REQUESTED number of blocks for block-maxima EVD.
#'        Will be adaptively reduced if k is too large relative to n. (default = 50)
#' @param dist_param Numeric: Shape/tail parameter passed to dgp_factory (default = 3)
#' @param mix_prop Numeric: Mixture proportion passed to dgp_factory (default = 0.1)
#' @param alpha Numeric: Significance level for coverage calculations (default = 0.05)
#'
#' @return A flat, 1-row data.frame containing all simulation inputs, detection
#'         metrics, coverage flags, EVD parameters, and wall-clock time.
#' @export
run_mis_iteration <- function(iter = 1, n = 1000, p = 1,
                              x_type = "normal", error_type = "normal",
                              outlier_method = "none", k = 1, magnitude = 5,
                              block_count = 50, dist_param = 3, mix_prop = 0.1,
                              alpha = 0.05) {
  
  if (p != 1) stop("sim_engine currently only supports p = 1.")
  start_time <- Sys.time()
  
  # -------------------------------------------------------------------------
  # 0. Adaptive block_count — prevent k vs block_size collapse
  # -------------------------------------------------------------------------
  # In exact_dfb_bmx, the clean sample (n - k) is split into M blocks.
  # Each block has block_size = (n - k) / M observations.
  # Dinkelbach selects k observations from each block, so we need:
  #   block_size >= k  →  (n - k) / M >= k  →  M <= (n - k) / k
  # GEV fitting requires >= 3 block maxima, so M >= 3.
  
  n_clean <- n - k
  max_feasible_blocks <- floor(n_clean / k)
  effective_block_count <- min(block_count, max_feasible_blocks)
  
  # Quality flag for the EVD estimate based on block count
  evd_quality <- if (effective_block_count < 3) {
    "infeasible"
  } else if (effective_block_count < 10) {
    "low"
  } else if (effective_block_count < 20) {
    "moderate"
  } else {
    "good"
  }
  
  if (effective_block_count < 3) {
    # Cannot fit GEV with fewer than 3 block maxima — return failure row
    warning(sprintf(
      "Infeasible (n=%d, k=%d): max %d blocks (need >=3). Skipping EVD.",
      n, k, max_feasible_blocks
    ))
    
    end_time <- Sys.time()
    compute_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    return(data.frame(
      iter = iter, n_obs = n, x_type = x_type, error_type = error_type,
      dist_param = dist_param, outlier_method = outlier_method,
      magnitude = magnitude, set_size = k, contam_prop = k / n,
      block_count = effective_block_count, evd_quality = evd_quality,
      alpha = alpha,
      detection_success = NA, detect_cooks = NA, detect_lev = NA,
      detect_dfbetas = NA, overlap_mis = NA, overlap_cooks = NA,
      overlap_lev = NA, overlap_dfbetas = NA,
      dfb_detected = NA_real_, dfb_injected = NA_real_,
      mean_lev_detected = NA_real_, mean_lev_injected = NA_real_,
      mean_res_detected = NA_real_, mean_res_injected = NA_real_,
      cover_evd = NA, cover_cooks = NA, cover_lev = NA, cover_dfbetas = NA,
      evd_test_type = if (outlier_method != "none") "injected" else "random",
      compute_time = compute_time,
      shape = NA_real_, scale = NA_real_, loc = NA_real_,
      set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
      stringsAsFactors = FALSE
    ))
  }
  
  # -------------------------------------------------------------------------
  # 1. Data Generation & Injection
  # -------------------------------------------------------------------------
  dat_clean <- generate_complex_data(n = n, p = p, x_type = x_type,
                                     error_type = error_type,
                                     dist_param = dist_param,
                                     mix_prop = mix_prop)
  
  if (outlier_method != "none") {
    dat <- apply_influence_shift(dat_clean, method = outlier_method,
                                 k = k, magnitude = magnitude)
    true_injected_indices <- dat$outlier_indices
  } else {
    dat <- dat_clean
    true_injected_indices <- NA
  }
  
  df_sim <- data.frame(y = dat$y, x = dat$X[, 1])
  
  # -------------------------------------------------------------------------
  # 2. Fit Base Model & Detect MIS
  # -------------------------------------------------------------------------
  base_model <- lm(y ~ x, data = df_sim)
  
  # Exact influence detection via Dinkelbach
  empirical_mis_pos <- dinkelbach_topk_lm(base_model, pos = 2, sign = 1, k = k)
  empirical_mis_neg <- dinkelbach_topk_lm(base_model, pos = 2, sign = -1, k = k)
  
  if (outlier_method != "none") {
    # Pick the MIS direction with better overlap for detection metrics
    overlap_pos <- length(intersect(true_injected_indices, empirical_mis_pos))
    overlap_neg <- length(intersect(true_injected_indices, empirical_mis_neg))
    if (overlap_neg > overlap_pos) {
      empirical_mis <- empirical_mis_neg
    } else {
      empirical_mis <- empirical_mis_pos
    }
    detection_success <- (length(intersect(true_injected_indices, 
                                           empirical_mis)) / k) >= 0.90
  } else {
    detection_success <- NA
    empirical_mis <- empirical_mis_pos
  }
  
  # -------------------------------------------------------------------------
  # 2b. Classical LOO Diagnostics — Detection (Top-K overlap)
  # -------------------------------------------------------------------------
  cooks_vals   <- cooks.distance(base_model)
  lev_vals     <- hatvalues(base_model)
  dfbetas_vals <- dfbetas(base_model)[, "x"]  # slope coefficient
  
  top_k_cooks   <- order(abs(cooks_vals),   decreasing = TRUE)[1:k]
  top_k_lev     <- order(abs(lev_vals),     decreasing = TRUE)[1:k]
  top_k_dfbetas <- order(abs(dfbetas_vals), decreasing = TRUE)[1:k]
  
  if (outlier_method != "none") {
    detect_cooks   <- (length(intersect(true_injected_indices, 
                                        top_k_cooks)) / k) >= 0.80
    detect_lev     <- (length(intersect(true_injected_indices, 
                                        top_k_lev)) / k) >= 0.80
    detect_dfbetas <- (length(intersect(true_injected_indices, 
                                        top_k_dfbetas)) / k) >= 0.80
    
    overlap_mis     <- length(intersect(true_injected_indices, empirical_mis)) / k
    overlap_cooks   <- length(intersect(true_injected_indices, top_k_cooks)) / k
    overlap_lev     <- length(intersect(true_injected_indices, top_k_lev)) / k
    overlap_dfbetas <- length(intersect(true_injected_indices, top_k_dfbetas)) / k
  } else {
    detect_cooks   <- NA
    detect_lev     <- NA
    detect_dfbetas <- NA
    overlap_mis     <- NA
    overlap_cooks   <- NA
    overlap_lev     <- NA
    overlap_dfbetas <- NA
  }
  
  # -------------------------------------------------------------------------
  # 2c. FWL Decomposition Diagnostics
  # -------------------------------------------------------------------------
  if (outlier_method != "none") {
    Z_diag <- matrix(1, nrow = n, ncol = 1)
    fwl_vars <- fwl(y = dat$y, X = dat$X[, 1], Z = Z_diag)
    X_fwl <- fwl_vars[, 2]
    Y_fwl <- fwl_vars[, 1]

    dfb_detected <- dfbeta_numeric(Y_fwl, cbind(X_fwl), empirical_mis, col_X = 1L)
    dfb_injected <- dfbeta_numeric(Y_fwl, cbind(X_fwl), true_injected_indices, col_X = 1L)

    h_vals  <- hatvalues(base_model)
    r_vals  <- abs(residuals(base_model))

    mean_lev_detected <- mean(h_vals[empirical_mis])
    mean_lev_injected <- mean(h_vals[true_injected_indices])
    mean_res_detected <- mean(r_vals[empirical_mis])
    mean_res_injected <- mean(r_vals[true_injected_indices])
  } else {
    dfb_detected      <- NA_real_
    dfb_injected      <- NA_real_
    mean_lev_detected <- NA_real_
    mean_lev_injected <- NA_real_
    mean_res_detected <- NA_real_
    mean_res_injected <- NA_real_
  }

  # -------------------------------------------------------------------------
  # 3. Determine the EVD test set
  # -------------------------------------------------------------------------
  if (outlier_method != "none") {
    evd_test_set <- true_injected_indices
  } else {
    evd_test_set <- sample(seq_len(n), k)
  }
  
  # -------------------------------------------------------------------------
  # 2c. Coverage — Threshold-based flagging (apples-to-apples comparison)
  #     All methods are evaluated on the SAME set: evd_test_set.
  #     For outlier scenarios → power. For "none" → empirical size.
  # -------------------------------------------------------------------------
  n_model <- stats::nobs(base_model)
  p_model <- length(stats::coef(base_model))  # includes intercept
  
  thresh_cooks   <- 4 / n_model
  thresh_lev     <- 2 * p_model / n_model
  thresh_dfbetas <- 2 / sqrt(n_model)
  
  # Evaluate: does ANY point in evd_test_set exceed the threshold?
  cover_cooks   <- any(cooks_vals[evd_test_set]          > thresh_cooks)
  cover_lev     <- any(lev_vals[evd_test_set]            > thresh_lev)
  cover_dfbetas <- any(abs(dfbetas_vals[evd_test_set])   > thresh_dfbetas)
  
  # -------------------------------------------------------------------------
  # 4. Estimate Extreme Value Distribution (EVD) Parameters
  #    Uses the ADAPTIVE effective_block_count computed in Section 0.
  # -------------------------------------------------------------------------
  Z_matrix <- matrix(1, nrow = length(dat$y), ncol = 1)
  
  res_exact <- tryCatch({
    evt_iter_dm(
      y = dat$y,
      x = dat$X[, 1],
      Z = Z_matrix,
      set = evd_test_set,
      block_count = effective_block_count
    )
  }, error = function(e) {
    warning(sprintf("evt_iter_dm failed (iter=%d, x=%s, err=%s, outlier=%s): %s",
                    iter, x_type, error_type, outlier_method,
                    conditionMessage(e)))
    data.frame(shape = NA_real_, scale = NA_real_, loc = NA_real_,
               set_dfb = NA_real_, p_value = NA_real_, converged = FALSE,
               stringsAsFactors = FALSE)
  })
  
  # If evt_iter_dm returned a list, coerce to 1-row data.frame
  expected_names <- c("shape", "scale", "loc", "set_dfb", "p_value", "converged")
  if (is.list(res_exact) && !is.data.frame(res_exact)) {
    missing <- setdiff(expected_names, names(res_exact))
    if (length(missing) > 0) {
      stop(sprintf("evt_iter_dm list is missing fields: %s  (got: %s)",
                   paste(missing, collapse = ", "),
                   paste(names(res_exact), collapse = ", ")))
    }
    res_exact <- as.data.frame(res_exact[expected_names],
                               stringsAsFactors = FALSE)
  }
  
  # EVD coverage: is p_value < alpha?
  cover_evd <- if (!is.na(res_exact$p_value)) {
    res_exact$p_value < alpha
  } else {
    NA
  }
  
  end_time <- Sys.time()
  compute_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # -------------------------------------------------------------------------
  # 5. Assemble Tidy Output
  # -------------------------------------------------------------------------
  res_df <- data.frame(
    iter = iter,
    n_obs = n,
    x_type = x_type,
    error_type = error_type,
    dist_param = dist_param,
    outlier_method = outlier_method,
    magnitude = magnitude,
    set_size = k,
    contam_prop = k / n,
    block_count = effective_block_count,
    evd_quality = evd_quality,
    alpha = alpha,
    # MIS detection
    detection_success = detection_success,
    # Classical detection (exact set match)
    detect_cooks = detect_cooks,
    detect_lev = detect_lev,
    detect_dfbetas = detect_dfbetas,
    # Partial overlap (fraction of injected points found in top-k)
    overlap_mis = overlap_mis,
    overlap_cooks = overlap_cooks,
    overlap_lev = overlap_lev,
    overlap_dfbetas = overlap_dfbetas,
    # FWL decomposition diagnostics
    dfb_detected = as.numeric(dfb_detected),
    dfb_injected = as.numeric(dfb_injected),
    mean_lev_detected = mean_lev_detected,
    mean_lev_injected = mean_lev_injected,
    mean_res_detected = mean_res_detected,
    mean_res_injected = mean_res_injected,
    # Coverage — threshold-based flagging (power or size)
    cover_evd = cover_evd,
    cover_cooks = cover_cooks,
    cover_lev = cover_lev,
    cover_dfbetas = cover_dfbetas,
    # Which set was tested by EVD
    evd_test_type = if (outlier_method != "none") "injected" else "random",
    # Timing
    compute_time = compute_time,
    stringsAsFactors = FALSE
  )
  
  if (!is.data.frame(res_exact) || nrow(res_exact) != 1) {
    stop("evt_iter_dm returned unexpected output structure.")
  }
  
  res_final <- cbind(res_df, res_exact)
  return(res_final)
}