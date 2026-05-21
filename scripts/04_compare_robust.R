# ==============================================================================
# File: /script/04_compare_robust.R
# Purpose: Orchestrates the cross-comparison with Dinkelbach detector and
#          iterative peel MIS. Compares 9 estimators:
#            Full OLS, CD, Leverage, DFBETAS,
#            MIS-alpha, MIS-oracle, MIS-peel,
#            MM, LTS
#          Produces k-diagnostic tables alongside bias/coverage/RMSE.
#
# Outputs: ../output/04_robust_comparison_results.rds
# ==============================================================================

# 1. Load Dependencies
library(dplyr)
library(purrr)
library(future)
library(furrr)
library(robustbase)

source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/diagnostics_classical.R")
source("../R/estimators_robust.R")
source("../R/dynamic_k_adaptive.R")
source("../R/dinkelbach_topk.R")
source("../R/leverage_k.R")
source("../R/iterative_peel_v2.R")
source("../R/sim_robust_engine.R")
source("../R/utils_checkpoint.R")

# 2. Global Configuration
sim_params <- list(
  n_iters   = 100,
  n_obs     = 5000,
  set_size  = 50,
  magnitude = 10,
  seed      = 20260503
)

set.seed(sim_params$seed)

num_workers <- max(1, future::availableCores() - 2)
cat(sprintf("Local environment, using %d workers.\n", num_workers))
plan(multisession, workers = num_workers)

# 3. Define the Full Reality Grid
param_grid <- expand.grid(
  x_type         = c("normal", "mixed_normal", "contaminated"),
  error_type     = c("normal", "mixed_normal", "skewed_t", "golm",
                     "beta_logistic", "gpd", "contaminated", "pareto"),
  outlier_method = c("none", "vertical_outlier", "good_leverage", "bad_leverage"),
  stringsAsFactors = FALSE
)

cat(sprintf(
  paste0("Starting 04 Robust Comparison Suite (Dinkelbach + Peel).\n",
         "  Total scenarios:     %d\n",
         "  Iterations each:     %d\n",
         "  Total MC draws:      %d\n\n"),
  nrow(param_grid),
  sim_params$n_iters,
  nrow(param_grid) * sim_params$n_iters
))

# 4. The Orchestrator Loop
for (i in seq_len(nrow(param_grid))) {
  
  p_current <- param_grid[i, , drop = FALSE]
  chunk_file <- sprintf("../output/temp_04/04_chunk_%03d.rds", i)
  
  if (is_computed(chunk_file)) {
    cat(sprintf("[%03d/%03d] Skipping (cached): x=%s | err=%s | out=%s\n",
                i, nrow(param_grid),
                p_current$x_type, p_current$error_type, p_current$outlier_method))
    next
  }
  
  cat(sprintf("[%03d/%03d] Computing: x=%s | err=%s | out=%s ... ",
              i, nrow(param_grid),
              p_current$x_type, p_current$error_type, p_current$outlier_method))
  
  scenario_results <- furrr::future_map_dfr(
    seq_len(sim_params$n_iters), function(iter_id) {
      
      tryCatch({
        run_robust_comparison_iter(
          iter           = iter_id,
          n              = sim_params$n_obs,
          p              = 1,
          x_type         = p_current$x_type,
          error_type     = p_current$error_type,
          outlier_method = p_current$outlier_method,
          k              = sim_params$set_size,
          magnitude      = sim_params$magnitude
        )
      }, error = function(e) {
        warning(sprintf("Iter %d failed for x=%s|err=%s|out=%s: %s",
                        iter_id, p_current$x_type, p_current$error_type,
                        p_current$outlier_method, e$message))
        return(NULL)
      })
      
    }, .options = furrr_options(seed = TRUE))
  
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# 5. Final Assembly
cat("\nAll scenarios completed. Assembling final dataset...\n")
final_dataset <- compile_checkpoints(
  temp_dir          = "../output/temp_04",
  pattern           = "^04_chunk_.*\\.rds$",
  final_output_path = "../output/04_robust_comparison_results.rds",
  clear_temp        = FALSE
)

cat("\nScript 04 execution finished successfully.\n")

# ==============================================================================
# 6. Diagnostics & Sanity Check
# ==============================================================================
library(tidyr)

results <- readRDS("../output/04_robust_comparison_results.rds")

# 6a. k-Selection Diagnostic: How does each method count?
cat("\n--- Mean k Selected by Each Method ---\n")
k_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    true_k   = mean(set_size, na.rm = TRUE),
    k_alpha  = mean(k_alpha, na.rm = TRUE),
    k_oracle = mean(k_oracle, na.rm = TRUE),
    k_peel   = mean(k_peel, na.rm = TRUE),
    .groups = "drop"
  )
print(k_table, n = Inf)

# 6b. Peel stopping diagnostics
cat("\n--- Peel Stop Reason Distribution ---\n")
peel_stops <- results %>%
  group_by(outlier_method, peel_stop) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(outlier_method) %>%
  mutate(pct = round(100 * n / sum(n), 1)) %>%
  arrange(outlier_method, desc(n))
print(peel_stops, n = Inf)

# 6c. Coverage comparison
cat("\n--- 95% CI Coverage (%) ---\n")
cov_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    full      = mean(cov_full, na.rm = TRUE) * 100,
    cd        = mean(cov_cd, na.rm = TRUE) * 100,
    lev       = mean(cov_lev, na.rm = TRUE) * 100,
    dfb       = mean(cov_dfb, na.rm = TRUE) * 100,
    mis_alpha = mean(cov_mis_alpha, na.rm = TRUE) * 100,
    mis_oracle= mean(cov_mis_oracle, na.rm = TRUE) * 100,
    mis_peel  = mean(cov_mis_peel, na.rm = TRUE) * 100,
    mm        = mean(cov_mm, na.rm = TRUE) * 100,
    lts       = mean(cov_lts, na.rm = TRUE) * 100,
    .groups = "drop"
  )
print(cov_table, n = Inf)

# 6d. Bias comparison
cat("\n--- Mean Absolute Bias ---\n")
bias_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    full      = mean(bias_full, na.rm = TRUE),
    cd        = mean(bias_cd, na.rm = TRUE),
    lev       = mean(bias_lev, na.rm = TRUE),
    dfb       = mean(bias_dfb, na.rm = TRUE),
    mis_alpha = mean(bias_mis_alpha, na.rm = TRUE),
    mis_oracle= mean(bias_mis_oracle, na.rm = TRUE),
    mis_peel  = mean(bias_mis_peel, na.rm = TRUE),
    mm        = mean(bias_mm, na.rm = TRUE),
    lts       = mean(bias_lts, na.rm = TRUE),
    .groups = "drop"
  )
print(bias_table, n = Inf)

# 6e. RMSE comparison
cat("\n--- RMSE ---\n")
rmse_table <- results %>%
  group_by(x_type, error_type, outlier_method) %>%
  summarise(
    full      = sqrt(mean((coef_full - 1)^2, na.rm = TRUE)),
    cd        = sqrt(mean((coef_cd - 1)^2, na.rm = TRUE)),
    lev       = sqrt(mean((coef_lev - 1)^2, na.rm = TRUE)),
    dfb       = sqrt(mean((coef_dfb - 1)^2, na.rm = TRUE)),
    mis_alpha = sqrt(mean((coef_mis_alpha - 1)^2, na.rm = TRUE)),
    mis_oracle= sqrt(mean((coef_mis_oracle - 1)^2, na.rm = TRUE)),
    mis_peel  = sqrt(mean((coef_mis_peel - 1)^2, na.rm = TRUE)),
    mm        = sqrt(mean((coef_mm - 1)^2, na.rm = TRUE)),
    lts       = sqrt(mean((coef_lts - 1)^2, na.rm = TRUE)),
    .groups = "drop"
  )
print(rmse_table, n = Inf)