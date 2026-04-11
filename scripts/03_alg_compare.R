# ==============================================================================
# File: /script/03_algorithmic_comparison.R
# Purpose: Systematically compare the literature's standard Greedy heuristic 
#          against our Exact Dinkelbach method across varying set sizes.
#          Quantifies the algorithmic bias and failure rate of heuristics.
# Inputs: None (Synthetic DGP)
# Outputs: ../output/algorithmic_bias_results.rds
# Paper Section: Methods / Algorithmic Exactness
# ==============================================================================

library(dplyr)
library(purrr)
library(tidyr)
library(testingMIS)

# Source our exact mathematical implementation
source("../R/exact_dfb_bmx.R")

# ------------------------------------------------------------------------------
# 1. Configuration & Parameter Grid
# ------------------------------------------------------------------------------
# We will test how the algorithms perform as the size of the set (k) increases.
# Hypothesis: Greedy fails worse as k grows.
SIM_GRID <- expand.grid(
  set_size = c(2, 5, 10),
  iter = 1:50 # 50 iterations per set size for robust averaging
)

n_obs <- 500
block_count <- 20
set.seed(20260409)

# ------------------------------------------------------------------------------
# 2. Data Generating Process
# ------------------------------------------------------------------------------
generate_sim_data <- function(n) {
  x <- rnorm(n)
  Z <- matrix(rnorm(n * 2), ncol = 2)
  error <- rnorm(n)
  y <- 2 * x + Z[, 1] - 0.5 * Z[, 2] + error
  list(y = y, x = x, Z = Z)
}

# ------------------------------------------------------------------------------
# 3. Main Head-to-Head Loop
# ------------------------------------------------------------------------------
cat("Starting Systematic Algorithmic Comparison...\n")

results <- purrr::map_dfr(1:nrow(SIM_GRID), function(row_idx) {
  
  current_k <- SIM_GRID$set_size[row_idx]
  current_iter <- SIM_GRID$iter[row_idx]
  
  if (row_idx %% 10 == 0) cat(sprintf("  Processing run %d / %d (k=%d)...\n", 
                                      row_idx, nrow(SIM_GRID), current_k))
  
  # 3a. Generate Data & FWL components
  dat <- generate_sim_data(n_obs)
  fwl_vars <- testingMIS:::fwl(y = dat$y, X = dat$x, Z = dat$Z)
  X_fwl <- fwl_vars[, 2]
  fwl_lm <- lm(fwl_vars[, 1] ~ X_fwl - 1)
  R_fwl <- residuals(fwl_lm)
  
  # 3b. Find an influential set to test using the exact method
  df_sim <- data.frame(y = dat$y, x = dat$x, Z1 = dat$Z[,1], Z2 = dat$Z[,2])
  base_model <- lm(y ~ x + Z1 + Z2, data = df_sim)
  
  sens_obj <- influence::sens(
    base_model,
    lambda = influence::set_lambda("beta_i", pos = 2, sign = sign(coef(base_model)[2]))
  )
  target_set <- sens_obj$influence$id[1:current_k]
  
  # 3c. Calculate Block Maxima using BOTH methods
  bmx_greedy <- abs(testingMIS::dfb_bmx(X_fwl, R_fwl, set = target_set, block_count = block_count))
  bmx_exact  <- abs(exact_dfb_bmx(X_fwl, R_fwl, set = target_set, block_count = block_count))
  
  # 3d. Return comparison metrics
  data.frame(
    iter = current_iter,
    k_size = current_k,
    mean_greedy = mean(bmx_greedy),
    mean_exact  = mean(bmx_exact),
    max_greedy  = max(bmx_greedy),
    max_exact   = max(bmx_exact)
  )
})

# ------------------------------------------------------------------------------
# 4. Analysis and Output
# ------------------------------------------------------------------------------
# Calculate the failure rate and magnitude of the greedy heuristic
summary_table <- results %>%
  mutate(
    # Did Exact find a mathematically larger influence? (Adding 1e-8 for float precision)
    exact_wins = mean_exact > (mean_greedy + 1e-8),
    # How much did Greedy underestimate the true maximum?
    percent_underestimation = ((mean_exact - mean_greedy) / mean_exact) * 100
  ) %>%
  group_by(k_size) %>%
  summarize(
    greedy_failure_rate = mean(exact_wins) * 100,
    avg_underestimation_pct = mean(percent_underestimation),
    max_underestimation_pct = max(percent_underestimation),
    .groups = "drop"
  )

print(summary_table)

if (!dir.exists("../output")) dir.create("../output")
saveRDS(results, "../output/algorithmic_bias_results.rds")