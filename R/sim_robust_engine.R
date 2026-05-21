# ==============================================================================
# File: /R/sim_robust_engine.R
# Purpose: Core iteration engine for the robust comparison simulation. Generates
#          one Monte Carlo draw: creates clean data, injects contamination,
#          fits all estimators, and returns a flat 1-row data.frame of bias,
#          raw coefficients, and 95% CI coverage for robust row-binding across
#          parallel simulation loops.
#
# Revision: v3 — Integrates Dinkelbach detector, leverage_k, iterative peel.
#           Replaces fast_sens_topk with dinkelbach_topk_lm throughout.
#           MIS variants reduced to: alpha_k, peel (leverage, batch=1), oracle_k.
#           Dropped: old_k (3-sigma), gap_k, adaptive_k (dominated by alpha_k).
# ==============================================================================

#' Check 95% Wald Interval Coverage
#'
#' @param coef Numeric; the point estimate of the coefficient.
#' @param se Numeric; the standard error of the coefficient estimate.
#' @param true_b Numeric; the true population parameter value.
#'
#' @return Integer: 1L if the true value falls within the 95% CI, 0L if not,
#'         NA_integer_ if either input is NA.
check_coverage <- function(coef, se, true_b) {
  if (is.na(coef) || is.na(se)) return(NA_integer_)
  lo <- coef - 1.96 * se
  hi <- coef + 1.96 * se
  as.integer(true_b >= lo & true_b <= hi)
}

#' Single Iteration of the Robust Comparison Simulation
#'
#' @param iter Integer; the current iteration index.
#' @param n Integer; sample size (default = 1000).
#' @param p Integer; number of predictors (default = 1).
#' @param x_type Character; distribution for the design matrix X. One of
#'        "normal", "mixed_normal", "contaminated".
#' @param error_type Character; distribution for the error term. One of
#'        "normal", "mixed_normal", "skewed_t", "golm", "beta_logistic",
#'        "gpd", "contaminated", "pareto".
#' @param outlier_method Character; contamination topology. One of "none",
#'        "vertical_outlier", "good_leverage", "bad_leverage".
#' @param k Integer; number of observations to contaminate.
#' @param magnitude Numeric; severity multiplier for the injected shift.
#'
#' @return A 1-row data.frame containing iteration metadata, absolute bias,
#'         raw coefficients, and CI coverage flags for all estimators.
run_robust_comparison_iter <- function(iter,
                                       n = 1000,
                                       p = 1,
                                       x_type = "normal",
                                       error_type = "normal",
                                       outlier_method,
                                       k,
                                       magnitude) {
  
  # ---------------------------------------------------------
  # 1. Data Generation & Injection
  # ---------------------------------------------------------
  dat_clean <- generate_complex_data(
    n = n, p = p,
    x_type = x_type,
    error_type = error_type
  )
  true_b <- dat_clean$true_beta[1]
  
  if (outlier_method != "none") {
    dat <- apply_influence_shift(
      dat_clean,
      method = outlier_method,
      k = k,
      magnitude = magnitude
    )
  } else {
    dat <- dat_clean
  }
  df <- data.frame(y = dat$y, x = dat$X[, 1])
  
  # ---------------------------------------------------------
  # 2. Baseline Contaminated OLS
  # ---------------------------------------------------------
  res_full <- fit_clean_ols(y ~ x, data = df, exclude_idx = integer(0))
  mod_full <- stats::lm(y ~ x, data = df)
  beta_full <- unname(stats::coef(mod_full)["x"])
  
  # ---------------------------------------------------------
  # 3. Classical Diagnostics (Using Default Statistical Thresholds)
  # Passing k = NULL forces get_classical_set to use theoretical defaults
  # (4/n for Cook's D, 2p/n for Leverage, 2/sqrt(n) for DFBETAS)
  # ---------------------------------------------------------
  cd_idx  <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "cooks_d")
  lev_idx <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "leverage")
  dfb_idx <- get_classical_set(mod_full, target_var = "x", k = NULL, metric = "dfbetas_target")
  
  res_cd  <- fit_clean_ols(y ~ x, data = df, exclude_idx = cd_idx)
  res_lev <- fit_clean_ols(y ~ x, data = df, exclude_idx = lev_idx)
  res_dfb <- fit_clean_ols(y ~ x, data = df, exclude_idx = dfb_idx)
  
  # ---------------------------------------------------------
  # 4. Direct Robust Estimation (MM and LTS)
  # ---------------------------------------------------------
  res_mm  <- fit_mm_estimator(y ~ x, data = df)
  res_lts <- fit_lts_estimator(y ~ x, data = df)
  
  # ---------------------------------------------------------
  # 5. Shared robust fit for k-selection
  # ---------------------------------------------------------
  mod_mm_obj <- tryCatch(
    robustbase::lmrob(y ~ x, data = df, setting = "KS2014"),
    error = function(e) NULL
  )
  
  beta_mm <- tryCatch(unname(coef(mod_mm_obj)["x"]), error = function(e) NA_real_)
  if (is.null(beta_mm) || is.na(beta_mm)) beta_mm <- beta_full
  
  # ---------------------------------------------------------
  # 6. MIS via Dinkelbach: shared direction-selection helper
  #    (used by alpha_k and oracle_k single-shot methods)
  # ---------------------------------------------------------
  run_mis_for_k <- function(k_val) {
    if (k_val == 0L) return(res_full)
    idx_pos <- dinkelbach_topk_lm(mod_full, pos = 2, sign =  1, k = k_val)
    idx_neg <- dinkelbach_topk_lm(mod_full, pos = 2, sign = -1, k = k_val)
    r_pos   <- fit_clean_ols(y ~ x, data = df, exclude_idx = idx_pos)
    r_neg   <- fit_clean_ols(y ~ x, data = df, exclude_idx = idx_neg)
    d_pos   <- abs(r_pos["coef"] - beta_mm)
    d_neg   <- abs(r_neg["coef"] - beta_mm)
    if (is.na(d_pos) && is.na(d_neg)) return(res_full)
    if (is.na(d_pos)) return(r_neg)
    if (is.na(d_neg)) return(r_pos)
    if (d_pos <= d_neg) r_pos else r_neg
  }
  
  # ---------------------------------------------------------
  # 7a. MIS — alpha_k (best single-shot baseline)
  # ---------------------------------------------------------
  k_alpha_val <- alpha_k(mod_mm_obj)
  res_mis_alpha <- run_mis_for_k(k_alpha_val)
  
  # ---------------------------------------------------------
  # 7b. MIS — oracle_k (theoretical ceiling)
  # ---------------------------------------------------------
  k_oracle_val <- oracle_k(if (outlier_method == "none") 0L else k)
  res_mis_oracle <- run_mis_for_k(k_oracle_val)
  
  # ---------------------------------------------------------
  # 7c. MIS — Iterative Peel (leverage, batch=1, Dinkelbach)
  # ---------------------------------------------------------
  peel_result <- tryCatch(
    iterative_peel_v2(
      formula    = y ~ x,
      data       = df,
      target_var = "x",
      target_pos = 2L,
      batch_size = 1L,
      max_iter   = 50L,
      max_k_frac = 0.06,
      detector   = "dinkelbach",
      k_method   = "leverage"
    ),
    error = function(e) {
      list(excluded = integer(0), k_total = 0L, n_iters = 0L,
           stop_reason = "error", beta_trajectory = numeric(0))
    }
  )
  
  res_mis_peel <- fit_clean_ols(y ~ x, data = df,
                                exclude_idx = peel_result$excluded)
  k_peel_val   <- peel_result$k_total
  
  # ---------------------------------------------------------
  # 8. Compile Bias and Coverage Metrics
  # ---------------------------------------------------------
  res <- data.frame(
    iter           = iter,
    x_type         = x_type,
    error_type     = error_type,
    outlier_method = outlier_method,
    set_size       = if (outlier_method == "none") 0L else k,
    
    # k counts
    k_alpha    = k_alpha_val,
    k_oracle   = k_oracle_val,
    k_peel     = k_peel_val,
    peel_stop  = peel_result$stop_reason,
    peel_iters = peel_result$n_iters,
    
    # Absolute Bias
    bias_full       = unname(abs(res_full["coef"]        - true_b)),
    bias_cd         = unname(abs(res_cd["coef"]          - true_b)),
    bias_lev        = unname(abs(res_lev["coef"]         - true_b)),
    bias_dfb        = unname(abs(res_dfb["coef"]         - true_b)),
    bias_mis_alpha  = unname(abs(res_mis_alpha["coef"]   - true_b)),
    bias_mis_oracle = unname(abs(res_mis_oracle["coef"]  - true_b)),
    bias_mis_peel   = unname(abs(res_mis_peel["coef"]    - true_b)),
    bias_mm         = unname(abs(res_mm["coef"]           - true_b)),
    bias_lts        = unname(abs(res_lts["coef"]          - true_b)),
    
    # 95% CI Coverage
    cov_full        = check_coverage(res_full["coef"],        res_full["se"],        true_b),
    cov_cd          = check_coverage(res_cd["coef"],          res_cd["se"],          true_b),
    cov_lev         = check_coverage(res_lev["coef"],         res_lev["se"],         true_b),
    cov_dfb         = check_coverage(res_dfb["coef"],         res_dfb["se"],         true_b),
    cov_mis_alpha   = check_coverage(res_mis_alpha["coef"],   res_mis_alpha["se"],   true_b),
    cov_mis_oracle  = check_coverage(res_mis_oracle["coef"],  res_mis_oracle["se"],  true_b),
    cov_mis_peel    = check_coverage(res_mis_peel["coef"],    res_mis_peel["se"],    true_b),
    cov_mm          = check_coverage(res_mm["coef"],          res_mm["se"],          true_b),
    cov_lts         = check_coverage(res_lts["coef"],         res_lts["se"],         true_b),
    
    # Raw Coefficients
    coef_full        = unname(res_full["coef"]),
    coef_cd          = unname(res_cd["coef"]),
    coef_lev         = unname(res_lev["coef"]),
    coef_dfb         = unname(res_dfb["coef"]),
    coef_mis_alpha   = unname(res_mis_alpha["coef"]),
    coef_mis_oracle  = unname(res_mis_oracle["coef"]),
    coef_mis_peel    = unname(res_mis_peel["coef"]),
    coef_mm          = unname(res_mm["coef"]),
    coef_lts         = unname(res_lts["coef"]),
    
    stringsAsFactors = FALSE
  )
  
  return(res)
}