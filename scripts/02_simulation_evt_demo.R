# ==============================================================================
# File: /script/02_simulation_evt_demo.R
# Purpose: Execute simulation demo to evaluate the statistical properties 
#          (size and power) of the EVT-based Most Influential Set hypothesis test.
#          Compares Greedy heuristic vs. Exact Dinkelbach on the same datasets.
# Inputs: None (Synthetic Data Generation Process defined internally)
# Outputs: ../output/sim_evt_null_distribution.rds
# ==============================================================================

library(dplyr)
library(purrr)
library(testingMIS)
library(influence)

source("../R/evt_iter.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")

# 1. Configuration & Parameters
sim_params <- list(
  n_iters     = 100,
  n_obs       = 2000,
  set_size    = 1,
  block_count = 40,
  seed        = 20260408
)

set.seed(sim_params$seed)

# 2. DGP
#' 
#' @param n Number of observations
#' @param set_size Number of observations to randomly select as the "test set"
#' @param inject_influence Logical; if TRUE, forces the selected set to be highly influential
#' @return A list containing y, x, Z, and the indices of the test set
generate_sim_data <- function(n, set_size, inject_influence = FALSE) {
  # Base covariates
  x <- rnorm(n)
  Z <- matrix(rnorm(n * 5), ncol = 5)
  
  # DGP: y = 2x + z1 - 0.5z2 + 0.3z3 - z4 + 0.7z5 + error
  error <- rnorm(n)
  y <- 2*x + Z[,1] - 0.5*Z[,2] + 0.3*Z[,3] - Z[,4] + 0.7*Z[,5] +
     0.4 * x * Z[,1] +
     0.3 * Z[,2] * Z[,3] +
     error
  
  # Randomly pick a set of indices to test
  test_set <- sample(1:n, size = set_size)
  
  if (isTRUE(inject_influence)) {
    # To engineer influence, we push the test set to extreme leverage (x) 
    # and extreme residual (y) space.
    x[test_set] <- x[test_set] + 5 
    y[test_set] <- y[test_set] - 10 
  }
  
  list(y = y, x = x, Z = Z, set = test_set)
}

# 3. Simulation Loop
cat(sprintf("Starting simulation with %d iterations...\n", sim_params$n_iters))

# purrr::map_dfr iterates from 1 to n_iters, runs the block, and row-binds the results
sim_results <- purrr::map_dfr(1:sim_params$n_iters, function(i) {
  
  # Progress tracker
  if (i %% 50 == 0) cat(sprintf("  Completed %d / %d...\n", i, sim_params$n_iters))
  
  # 3a. Generate Data (Testing the Null: inject_influence = FALSE)
  dat <- generate_sim_data(
    n = sim_params$n_obs, 
    set_size = sim_params$set_size, 
    inject_influence = FALSE
  )
  
  # Find MIS
  # 1). Fit the OLS model
  # NEW — all 5 Z columns
  df_sim <- data.frame(y = dat$y, x = dat$x, 
                       Z1 = dat$Z[,1], Z2 = dat$Z[,2], Z3 = dat$Z[,3],
                       Z4 = dat$Z[,4], Z5 = dat$Z[,5])
  base_model <- lm(y ~ x * Z1 + Z2 * Z3 + Z4 + Z5, data = df_sim)
  
  # 2). Use 'sens' generic.
  sens_obj <- influence::sens(
    base_model,
    lambda = influence::set_lambda("beta_i", pos = 2, sign = sign(coef(base_model)[2]))
  )
  
  # Extract the exact top-k indices of the MIS
  true_mis_indices <- sens_obj$influence$id[1:sim_params$set_size]

  non_mis  <- setdiff(1:sim_params$n_obs, true_mis_indices)
  null_set <- sample(non_mis, size = sim_params$set_size)
  
  
  
  # 3b. Run the EVT wrapper on the TRUE MIS — GREEDY
  res_greedy <- evt_iter(
    y = dat$y, 
    x = dat$x, 
    Z = dat$Z, 
    set = true_mis_indices, 
    block_count = sim_params$block_count
  )
  res_greedy$algorithm <- "Greedy"
  
  # 3c. Run the EVT wrapper on the TRUE MIS — EXACT Dinkelbach
  res_exact <- evt_iter_dm(
    y = dat$y,
    x = dat$x,
    Z = dat$Z,
    set = true_mis_indices,
    block_count = sim_params$block_count
  )
  res_exact$algorithm <- "Exact" 
  
  # 3d. Append iteration metadata and combine 
  for (res in list(res_greedy, res_exact)) {
    res$iter             <- i
    res$n_obs            <- sim_params$n_obs
    res$set_size         <- sim_params$set_size
    res$inject_influence <- FALSE
  }
  
  bind_rows(res_greedy, res_exact)
})

# 4. Save and Quick Diagnostic
if (!dir.exists("../output")) dir.create("../output")

output_file <- "../output/sim_evt_null_distribution.rds"
saveRDS(sim_results, output_file)

cat(sprintf("\nSimulation complete. Saved %d rows to %s\n", nrow(sim_results), output_file))

# Sanity check
summary_stats <- sim_results %>%
  group_by(algorithm) %>%
  summarize(
    convergence_rate = mean(converged, na.rm = TRUE) * 100,
    empirical_size   = mean(p_value[converged == TRUE] < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_stats)