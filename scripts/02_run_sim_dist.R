# ==============================================================================
# File: /script/02_run_sim_dist.R
# Purpose: Execute simulation mapping MIS EVD convergence across complex 
#          error distributions and adversarial leverage architectures.
# Inputs: ../R/ (All fat engine scripts)
# Outputs: ../output/02_final_distributions.rds
# ==============================================================================

# 1. Load Dependencies & Engines
library(dplyr)
library(purrr)
library(testingMIS)
library(influence)

source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/sim_engine.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")

unlink("../output/temp", recursive = TRUE)
unlink("../output/02_final_distributions.rds")

# 2. Global Configuration
sim_params <- list(
  n_iters = 100,      # Number of simulations per grid combination
  n_obs = 1000,     # Standard sample size for 02_script
  set_size = 3,        # k: Size of the influential set to inject
  magnitude = 5,        # 5-sigma shift for the outliers
  block_count = 25,       # Blocks for EVD estimation
  seed = 20260415
)

set.seed(sim_params$seed)

# Set up Parallel Processing
# Using 14 workers (threads) to leave 2 threads for OS stability on a 16-thread CPU
#install.packages("future")
#install.packages("furrr")
library(future)
library(furrr)
plan(multisession, workers = 14)
cat("Parallel processing initialized with 14 workers.\n\n")

# 3. Define the Reality Grid (The Parameter Space)
# We expand all combinations of X distributions, Error distributions, and Outlier types
param_grid <- expand.grid(
  x_type = c("normal", "skewed_t", "pareto"),
  error_type = c("normal", "mixed_normal", "skewed_t", "golm", 
                 "beta_logistic", "gpd", "contaminated", "pareto"),
  outlier_method = c("none", "vertical_outlier", "good_leverage", "bad_leverage"),
  
  dist_param = c(1.5, 2.1, 3.0, 5.0, 10.0),
  mix_prop   = c(0.05, 0.15, 0.30),
  
  stringsAsFactors = FALSE
)

distributions_using_shape <- c("skewed_t", "pareto", "gpd")
distributions_using_mix   <- c("mixed_normal", "contaminated")

param_grid <- param_grid %>%
  filter(
    ((x_type %in% distributions_using_shape | error_type %in% distributions_using_shape) | dist_param == 3.0) &
      
    ((x_type %in% distributions_using_mix | error_type %in% distributions_using_mix) | mix_prop == 0.05)
  )

cat(sprintf("Starting 02 Simulation Suite.\nTotal Scenarios: %d\nIterations per Scenario: %d\n\n", 
            nrow(param_grid), sim_params$n_iters))

# 4. The Orchestrator Loop
for (i in seq_len(nrow(param_grid))) {
  
  # Extract current scenario parameters
  p_current <- param_grid[i, ]
  
  # Define safe chunk path
  chunk_file <- sprintf("../output/temp/02_chunk_%04d.rds", i)
  
  if (is_computed(chunk_file)) {
    cat(sprintf("[%04d/%04d] Skipping (Already Computed): X=%s, Error=%s, Outlier=%s\n", 
                i, nrow(param_grid), p_current$x_type, p_current$error_type, p_current$outlier_method))
    next
  }
  
  cat(sprintf("[%04d/%04d] Computing: X=%s, Error=%s, Outlier=%s ... ", 
              i, nrow(param_grid), p_current$x_type, p_current$error_type, p_current$outlier_method))
  
  # Run the n_iters for this specific combination
  # Run the n_iters for this specific combination
  scenario_results <- furrr::future_map_dfr(
    1:sim_params$n_iters, function(iter_id) {
    
    # FIX: Removed the double assignment
    iteration_seed <- sim_params$seed + (i - 1) * sim_params$n_iters + iter_id
    set.seed(iteration_seed)
    
    # FIX: Added curly braces {} and the crucial comma before error = ...
    tryCatch({
      run_mis_iteration(
        iter = iter_id,
        n = sim_params$n_obs,
        p = 1,
        x_type = p_current$x_type,
        error_type = p_current$error_type,
        outlier_method = p_current$outlier_method,
        k = sim_params$set_size,
        magnitude = sim_params$magnitude,
        block_count = sim_params$block_count,
        dist_param = p_current$dist_param, 
        mix_prop = p_current$mix_prop
      )
    }, error = function(e) {
      warning(sprintf("Iter %d failed: %s", iter_id, e$message))
      return(NULL) # This will silently drop the failed iteration from the final data.frame
    }, .options = furrr_options(seed = TRUE))
    
  })
  
  # Safely save the chunk
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# 5. Final Assembly
cat("\nAll scenarios completed. Assembling final dataset...\n")


final_dataset <- compile_checkpoints(
  temp_dir = "../output/temp", 
  pattern = "^02_chunk_.*\\.rds$", 
  final_output_path = "../output/02_final_distributions.rds",
  clear_temp = FALSE  # Keep chunks until you are 100% sure the final file is perfect
)
results <- final_dataset
cat("\nScript 02 execution finished successfully.\n")

# Quick Sanity Check for 02_run_sim_dist.R Output
library(tidyverse)

# Uncomment and run if you haven't loaded the data into your environment yet:
# results <- readRDS("../output/02_final_distributions.rds")

cat("\n=== 1. QUICK STRUCTURAL AUDIT ===\n")
expected_rows <- 57600
cat(sprintf("Rows: %d / %d (Missing: %d)\n", nrow(results), expected_rows, expected_rows - nrow(results)))
cat(sprintf("Overall Convergence: %.1f%%\n", mean(results$converged, na.rm = TRUE) * 100))

# Check for sketchy NAs (should only be NA if converged == FALSE)
bad_na <- results %>% filter(converged == TRUE & (is.na(shape) | is.na(scale) | is.na(loc)))
if(nrow(bad_na) > 0) cat("⚠️ WARNING: Found", nrow(bad_na), "rows where converged=TRUE but EVD params are NA!\n") else cat("✅ NA logic holds up.\n")


cat("\n=== 2. CORE METRICS ===\n")
cat("\n-- Detection Success (%) --\n")
results %>%
  filter(outlier_method != "none") %>%
  group_by(error_type, outlier_method) %>%
  summarise(det_pct = round(mean(detection_success, na.rm = TRUE) * 100, 1), .groups = "drop") %>%
  pivot_wider(names_from = outlier_method, values_from = det_pct) %>%
  print()

cat("\n-- EVD Shape (ξ) Summaries (Converged Only) --\n")
results %>%
  filter(converged == TRUE) %>%
  group_by(error_type) %>%
  summarise(
    n = n(),
    shape_mean = round(mean(shape, na.rm = TRUE), 3),
    shape_med = round(median(shape, na.rm = TRUE), 3)
  ) %>%
  arrange(desc(shape_mean)) %>%
  print()


# === 3. QUICK PLOTS (Will render in your plot pane) ===
cat("\n=== Rendering Plots to Viewer... ===\n")

# Plot 1: Convergence Heatmap
results %>%
  group_by(error_type, x_type) %>%
  summarise(conv_rate = mean(converged, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = x_type, y = error_type, fill = conv_rate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(conv_rate, accuracy = 1)), size = 3) +
  scale_fill_gradient2(low = "#A32D2D", mid = "#FAC775", high = "#0F6E56", midpoint = 0.5, labels = scales::percent) +
  labs(title = "EVD Convergence Rate", x = "X Dist", y = "Error Dist") +
  theme_minimal() -> p1; print(p1)

# Plot 2: Detection Heatmap
results %>%
  filter(outlier_method != "none") %>%
  group_by(error_type, outlier_method) %>%
  summarise(det = mean(detection_success, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = outlier_method, y = error_type, fill = det)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(det, accuracy = 1)), size = 3) +
  scale_fill_gradient2(low = "#A32D2D", mid = "#FAC775", high = "#0F6E56", midpoint = 0.5, labels = scales::percent) +
  labs(title = "MIS Detection Success Rate", x = "Injection", y = "Error Dist") +
  theme_minimal() -> p2; print(p2)

# Plot 3: Shape Parameter Distributions
results %>%
  filter(converged == TRUE) %>%
  ggplot(aes(x = reorder(error_type, shape, median), y = shape, fill = error_type)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  labs(title = "EVD Shape Parameter (ξ)", x = NULL, y = "Shape") +
  theme_minimal() -> p3; print(p3)
