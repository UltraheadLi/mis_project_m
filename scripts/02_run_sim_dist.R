# ==============================================================================
# File: /script/02_run_sim_dist.R
# Purpose: Execute simulation mapping MIS EVD convergence across complex 
#          error distributions, adversarial leverage architectures, and a 
#          grid of sample sizes × contamination proportions.
#
# Revision notes (v2):
#   - Replaced fixed n_obs/set_size with n_obs_grid × contam_prop_grid
#   - k is computed as floor(contam_prop * n) with min-k guard
#   - Trimmed dist_param and mix_prop grids for tractability
#   - Reduced n_iters from 100 to 50
#   - CRITICAL FIX: feasibility guard uses correct constraint:
#       floor((n - k) / k) >= 3, i.e., n >= 4k
#     This ensures each block in exact_dfb_bmx has at least k observations
#     for the Dinkelbach solver, AND that we have >= 3 blocks for GEV fitting.
#   - Removed dead iteration_seed code inside furrr lambda
#   - Added coverage analysis section comparing EVD vs classical tools
#   - Checkpoint files now encode (n, contam_prop, scenario_id)
#
# Inputs: ../R/ (All engine scripts)
# Outputs: ../output/02_final_distributions.rds
# ==============================================================================

# 1. Load Dependencies & Engines
library(dplyr)
library(purrr)
library(evd)

source("../R/helpers_local.R")
source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")
source("../R/dinkelbach_topk.R")
source("../R/sim_engine.R")

# 2. Global Configuration
sim_params <- list(
  n_iters     = 50,         # Iterations per grid cell (reduced for tractability)
  magnitude   = 10,         # Sigma shift for outliers
  block_count = 50,         # REQUESTED blocks (sim_engine will adaptively cap)
  alpha       = 0.05,       # Significance level for coverage
  seed        = 20260415
)

set.seed(sim_params$seed)

# --- Sample size × Contamination proportion grid ---
n_obs_grid       <- c(500, 1000, 2500, 5000)
contam_prop_grid <- c(0.01, 0.05, 0.10, 0.20)

# Build the (n, contam_prop) design with feasibility guard
nk_grid <- expand.grid(
  n_obs       = n_obs_grid,
  contam_prop = contam_prop_grid,
  stringsAsFactors = FALSE
)
nk_grid$set_size <- pmax(floor(nk_grid$contam_prop * nk_grid$n_obs), 2L)

# CORRECT feasibility constraint:
#   block_size = (n - k) / M  must be >= k for Dinkelbach
#   → M <= (n - k) / k = floor((n - k) / k)
#   GEV needs M >= 3, so: floor((n - k) / k) >= 3 → n >= 4k
#   (sim_engine also handles this adaptively, but we filter obviously
#    infeasible cells here to avoid wasting orchestrator time)
nk_grid$max_blocks <- floor((nk_grid$n_obs - nk_grid$set_size) / nk_grid$set_size)
nk_grid$effective_blocks <- pmin(sim_params$block_count, nk_grid$max_blocks)

nk_grid <- nk_grid %>%
  filter(max_blocks >= 3)

cat(sprintf("Sample size × contamination grid: %d feasible cells\n", nrow(nk_grid)))
cat("  (Constraint: each block must hold >= k obs for Dinkelbach, need >= 3 blocks for GEV)\n")
print(nk_grid)

# --- Parallelism ---
library(future)
library(furrr)
slurm_cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
if (slurm_cpus != "") {
  num_workers <- as.numeric(slurm_cpus)
  cat(sprintf("Detected SLURM environment, using %d workers.\n", num_workers))
} else {
  num_workers <- max(1, parallel::detectCores() - 2)
  cat(sprintf("Local environment, using %d workers.\n", num_workers))
}
plan(multisession, workers = num_workers)

# 3. Define the Distribution × Outlier Grid
# Trimmed for tractability with the (n, k) grid
param_grid <- expand.grid(
  x_type = c("normal", "skewed_t", "pareto"),
  error_type = c("normal", "mixed_normal", "skewed_t", "golm", 
                 "beta_logistic", "gpd", "contaminated", "pareto"),
  outlier_method = c("none", "vertical_outlier", "good_leverage", "bad_leverage"),
  
  # Trimmed: 3 shape values instead of 5
  dist_param = c(1.5, 3.0, 10.0),
  # Trimmed: 2 mixture values instead of 3
  mix_prop   = c(0.05, 0.30),
  
  stringsAsFactors = FALSE
)

# Filter: dist_param only matters for distributions with shape parameters
distributions_using_shape <- c("skewed_t", "pareto", "gpd")
distributions_using_mix   <- c("mixed_normal", "contaminated")

param_grid <- param_grid %>%
  filter(
    ((x_type %in% distributions_using_shape | error_type %in% distributions_using_shape) | dist_param == 3.0) &
    ((x_type %in% distributions_using_mix | error_type %in% distributions_using_mix) | mix_prop == 0.05)
  )

total_cells <- nrow(nk_grid) * nrow(param_grid)
cat(sprintf(
  paste0("\nStarting 02 Simulation Suite.\n",
         "  (n, contam) cells: %d\n",
         "  Distribution scenarios: %d\n",
         "  Total grid cells: %d\n",
         "  Iterations per cell: %d\n",
         "  Total iterations: %d\n\n"),
  nrow(nk_grid), nrow(param_grid), total_cells, 
  sim_params$n_iters, total_cells * sim_params$n_iters
))

# 4. The Orchestrator Loop — outer over (n, contam_prop), inner over scenarios
global_counter <- 0L

for (g in seq_len(nrow(nk_grid))) {
  
  current_n   <- nk_grid$n_obs[g]
  current_cp  <- nk_grid$contam_prop[g]
  current_k   <- nk_grid$set_size[g]
  eff_blocks  <- nk_grid$effective_blocks[g]
  
  cat(sprintf(
    "\n========== n=%d, contam_prop=%.2f, k=%d, effective_blocks=%d ==========\n",
    current_n, current_cp, current_k, eff_blocks
  ))
  
  for (i in seq_len(nrow(param_grid))) {
    
    global_counter <- global_counter + 1L
    p_current <- param_grid[i, ]
    
    # Checkpoint file encodes (n, contam_prop, scenario_id)
    chunk_file <- sprintf(
      "../output/temp_02/02_chunk_n%d_cp%03d_%04d.rds",
      current_n, round(current_cp * 100), i
    )
    
    if (is_computed(chunk_file)) {
      cat(sprintf("[%05d/%05d] Skipping (Cached): n=%d, cp=%.2f, X=%s, Err=%s, Out=%s\n",
                  global_counter, total_cells,
                  current_n, current_cp,
                  p_current$x_type, p_current$error_type, p_current$outlier_method))
      next
    }
    
    cat(sprintf("[%05d/%05d] Computing: n=%d, cp=%.2f, k=%d, X=%s, Err=%s, Out=%s ... ",
                global_counter, total_cells,
                current_n, current_cp, current_k,
                p_current$x_type, p_current$error_type, p_current$outlier_method))
    
    # Capture loop variables for safe furrr closure
    local_n     <- current_n
    local_k     <- current_k
    local_bc    <- sim_params$block_count   # requested; sim_engine caps adaptively
    local_mag   <- sim_params$magnitude
    local_alpha <- sim_params$alpha
    local_x     <- p_current$x_type
    local_err   <- p_current$error_type
    local_out   <- p_current$outlier_method
    local_dp    <- p_current$dist_param
    local_mp    <- p_current$mix_prop
    
    # Run n_iters for this specific combination
    scenario_results <- furrr::future_map_dfr(
      1:sim_params$n_iters, function(iter_id) {
        
        tryCatch({
          run_mis_iteration(
            iter           = iter_id,
            n              = local_n,
            p              = 1,
            x_type         = local_x,
            error_type     = local_err,
            outlier_method = local_out,
            k              = local_k,
            magnitude      = local_mag,
            block_count    = local_bc,
            dist_param     = local_dp,
            mix_prop       = local_mp,
            alpha          = local_alpha
          )
        }, error = function(e) {
          warning(sprintf("Iter %d failed: %s", iter_id, e$message))
          return(NULL)
        })
        
      }, .options = furrr_options(seed = TRUE))
    
    safe_save_rds(scenario_results, chunk_file)
    cat("Done.\n")
  }
}

# 5. Final Assembly
cat("\nAll scenarios completed. Assembling final dataset...\n")

final_dataset <- compile_checkpoints(
  temp_dir = "../output/temp_02", 
  pattern = "^02_chunk_.*\\.rds$", 
  final_output_path = "../output/02_final_distributions.rds",
  clear_temp = FALSE
)
results <- final_dataset
cat("\nScript 02 execution finished successfully.\n")

# ==============================================================================
# 6. Sanity Check & Analysis
# ==============================================================================
library(tidyverse)

# Uncomment if loading from saved file:
# results <- readRDS("../output/02_final_distributions.rds")

cat("\n=== 1. STRUCTURAL AUDIT ===\n")
cat(sprintf("Total rows: %d\n", nrow(results)))
cat(sprintf("Unique (n, contam_prop) cells: %d\n", 
            nrow(distinct(results, n_obs, contam_prop))))
cat(sprintf("Unique block_count values: %s\n",
            paste(sort(unique(results$block_count)), collapse = ", ")))
cat(sprintf("Overall EVD convergence: %.1f%%\n", 
            mean(results$converged, na.rm = TRUE) * 100))

bad_na <- results %>% filter(converged == TRUE & (is.na(shape) | is.na(scale) | is.na(loc)))
if (nrow(bad_na) > 0) {
  cat("WARNING: Found", nrow(bad_na), "rows with converged=TRUE but NA EVD params!\n")
} else {
  cat("NA logic OK.\n")
}

# Quick check: which (n, k) cells had block_count reduced?
cat("\n-- Adaptive block_count summary --\n")
results %>%
  group_by(n_obs, set_size, contam_prop) %>%
  summarise(
    block_min = min(block_count),
    block_max = max(block_count),
    block_med = median(block_count),
    .groups = "drop"
  ) %>%
  print(n = 50)


# === 2. DETECTION SUCCESS (grouped by n × contam_prop) ===
cat("\n=== 2. DETECTION SUCCESS (%) — by n, contam_prop, method ===\n")
results %>%
  filter(outlier_method != "none") %>%
  group_by(n_obs, contam_prop, outlier_method) %>%
  summarise(
    MIS    = round(mean(detection_success, na.rm = TRUE) * 100, 1),
    Cooks  = round(mean(detect_cooks, na.rm = TRUE) * 100, 1),
    Lev    = round(mean(detect_lev, na.rm = TRUE) * 100, 1),
    DFB    = round(mean(detect_dfbetas, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  print(n = 100)


# === 3. COVERAGE — THE KEY COMPARISON ===
cat("\n=== 3. COVERAGE COMPARISON (Power & Size) ===\n")

# 3a. Empirical SIZE (false positive rate) — outlier_method == "none"
cat("\n-- 3a. Empirical Size (Should be ~5%) --\n")
results %>%
  filter(outlier_method == "none") %>%
  group_by(n_obs, contam_prop) %>%
  summarise(
    MIS     = round(mean(cover_evd, na.rm = TRUE) * 100, 1),
    Cooks   = round(mean(cover_cooks, na.rm = TRUE) * 100, 1),
    Lev     = round(mean(cover_lev, na.rm = TRUE) * 100, 1),
    DFBETAS = round(mean(cover_dfbetas, na.rm = TRUE) * 100, 1),
    n_conv  = sum(converged, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print(n = 100)

# 3b. POWER — outlier scenarios
cat("\n-- 3b. Power (%) by outlier method --\n")
results %>%
  filter(outlier_method != "none") %>%
  group_by(n_obs, contam_prop, outlier_method) %>%
  summarise(
    MIS     = round(mean(cover_evd, na.rm = TRUE) * 100, 1),
    Cooks   = round(mean(cover_cooks, na.rm = TRUE) * 100, 1),
    Lev     = round(mean(cover_lev, na.rm = TRUE) * 100, 1),
    DFBETAS = round(mean(cover_dfbetas, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  print(n = 100)

# 3c. Power broken down by error distribution (for bad_leverage, most informative)
cat("\n-- 3c. Power by error distribution (bad_leverage only) --\n")
results %>%
  filter(outlier_method == "bad_leverage") %>%
  group_by(n_obs, contam_prop, error_type) %>%
  summarise(
    MIS     = round(mean(cover_evd, na.rm = TRUE) * 100, 1),
    Cooks   = round(mean(cover_cooks, na.rm = TRUE) * 100, 1),
    Lev     = round(mean(cover_lev, na.rm = TRUE) * 100, 1),
    DFBETAS = round(mean(cover_dfbetas, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  print(n = 200)

# 3d. Coverage stratified by EVD quality flag
cat("\n-- 3d. Coverage by EVD Quality (all outlier methods) --\n")
results %>%
  filter(outlier_method != "none") %>%
  group_by(evd_quality, outlier_method) %>%
  summarise(
    n_rows  = n(),
    MIS     = round(mean(cover_evd, na.rm = TRUE) * 100, 1),
    Cooks   = round(mean(cover_cooks, na.rm = TRUE) * 100, 1),
    Lev     = round(mean(cover_lev, na.rm = TRUE) * 100, 1),
    DFBETAS = round(mean(cover_dfbetas, na.rm = TRUE) * 100, 1),
    conv    = round(mean(converged, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(factor(evd_quality, levels = c("good", "moderate", "low", "infeasible"))) %>%
  print(n = 50)

cat("\n-- 3e. Size control stratified by EVD quality (none scenarios) --\n")
results %>%
  filter(outlier_method == "none") %>%
  group_by(evd_quality) %>%
  summarise(
    n_rows  = n(),
    MIS     = round(mean(cover_evd, na.rm = TRUE) * 100, 1),
    Cooks   = round(mean(cover_cooks, na.rm = TRUE) * 100, 1),
    Lev     = round(mean(cover_lev, na.rm = TRUE) * 100, 1),
    DFBETAS = round(mean(cover_dfbetas, na.rm = TRUE) * 100, 1),
    conv    = round(mean(converged, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  ) %>%
  print(n = 20)


# === 4. EVD SHAPE SUMMARIES ===
cat("\n=== 4. EVD Shape (xi) by n and error type ===\n")
results %>%
  filter(converged == TRUE) %>%
  group_by(n_obs, error_type) %>%
  summarise(
    n = n(),
    shape_mean = round(mean(shape, na.rm = TRUE), 3),
    shape_med  = round(median(shape, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(n_obs, desc(shape_mean)) %>%
  print(n = 100)


# ==============================================================================
# === 5. PLOTS ===
# ==============================================================================
cat("\n=== Rendering Plots... ===\n")

# ---------- Plot 1: Convergence Heatmap (faceted by n_obs) ----------
results %>%
  group_by(n_obs, error_type, x_type) %>%
  summarise(conv_rate = mean(converged, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = x_type, y = error_type, fill = conv_rate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(conv_rate, accuracy = 1)), size = 2.5) +
  scale_fill_gradient2(low = "#A32D2D", mid = "#FAC775", high = "#0F6E56",
                       midpoint = 0.5, labels = scales::percent) +
  facet_wrap(~n_obs, labeller = label_both) +
  labs(title = "EVD Convergence Rate by Sample Size",
       x = "X Distribution", y = "Error Distribution") +
  theme_minimal(base_size = 10) -> p1
print(p1)


# ---------- Plot 2: Size Control Panel ----------
# Empirical rejection rate under the null — should be ~alpha
results %>%
  filter(outlier_method == "none") %>%
  group_by(n_obs, contam_prop) %>%
  summarise(
    MIS      = mean(cover_evd, na.rm = TRUE),
    `Cook's D` = mean(cover_cooks, na.rm = TRUE),
    Leverage = mean(cover_lev, na.rm = TRUE),
    DFBETAS  = mean(cover_dfbetas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(MIS, `Cook's D`, Leverage, DFBETAS),
               names_to = "Method", values_to = "rej_rate") %>%
  ggplot(aes(x = factor(n_obs), y = rej_rate, 
             color = factor(contam_prop), group = factor(contam_prop))) +
  geom_point(size = 2) +
  geom_line() +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "grey50") +
  facet_wrap(~Method) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "Size Control: Empirical Rejection Rate Under Null",
       subtitle = "Dashed line = nominal alpha = 5%",
       x = "Sample Size (n)", y = "Rejection Rate",
       color = "Contam.\nProportion") +
  theme_minimal(base_size = 10) -> p2
print(p2)


# ---------- Plot 3: Power Heatmap — 4-panel (bad_leverage) ----------
results %>%
  filter(outlier_method == "bad_leverage") %>%
  group_by(n_obs, contam_prop) %>%
  summarise(
    MIS      = mean(cover_evd, na.rm = TRUE),
    `Cook's D` = mean(cover_cooks, na.rm = TRUE),
    Leverage = mean(cover_lev, na.rm = TRUE),
    DFBETAS  = mean(cover_dfbetas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(MIS, `Cook's D`, Leverage, DFBETAS),
               names_to = "Method", values_to = "power") %>%
  ggplot(aes(x = factor(contam_prop), y = factor(n_obs), fill = power)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(power, accuracy = 1)), size = 3) +
  scale_fill_gradient2(low = "#A32D2D", mid = "#FAC775", high = "#0F6E56",
                       midpoint = 0.5, labels = scales::percent) +
  facet_wrap(~Method) +
  labs(title = "Power Comparison: Bad Leverage Scenarios",
       x = "Contamination Proportion", y = "Sample Size (n)") +
  theme_minimal(base_size = 10) -> p3
print(p3)


# ---------- Plot 4: Power Curves by Contamination Proportion ----------
results %>%
  filter(outlier_method != "none") %>%
  group_by(n_obs, contam_prop, outlier_method) %>%
  summarise(
    MIS      = mean(cover_evd, na.rm = TRUE),
    `Cook's D` = mean(cover_cooks, na.rm = TRUE),
    Leverage = mean(cover_lev, na.rm = TRUE),
    DFBETAS  = mean(cover_dfbetas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(MIS, `Cook's D`, Leverage, DFBETAS),
               names_to = "Method", values_to = "power") %>%
  ggplot(aes(x = contam_prop, y = power, color = Method, 
             linetype = factor(n_obs), group = interaction(Method, n_obs))) +
  geom_point(size = 1.5) +
  geom_line() +
  facet_wrap(~outlier_method) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  scale_x_continuous(labels = scales::percent) +
  labs(title = "Power vs Contamination Proportion",
       x = "Contamination Proportion (k/n)", y = "Power (Coverage Rate)",
       color = "Method", linetype = "n") +
  theme_minimal(base_size = 10) -> p4
print(p4)


# ---------- Plot 5: Detection Heatmap — MIS vs Classical (bad_leverage) ----------
results %>%
  filter(outlier_method == "bad_leverage") %>%
  group_by(n_obs, contam_prop) %>%
  summarise(
    MIS      = mean(detection_success, na.rm = TRUE),
    `Cook's D` = mean(detect_cooks, na.rm = TRUE),
    Leverage = mean(detect_lev, na.rm = TRUE),
    DFBETAS  = mean(detect_dfbetas, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(c(MIS, `Cook's D`, Leverage, DFBETAS),
               names_to = "Method", values_to = "det") %>%
  ggplot(aes(x = factor(contam_prop), y = factor(n_obs), fill = det)) +
  geom_tile(color = "white") +
  geom_text(aes(label = scales::percent(det, accuracy = 1)), size = 3) +
  scale_fill_gradient2(low = "#A32D2D", mid = "#FAC775", high = "#0F6E56",
                       midpoint = 0.5, labels = scales::percent) +
  facet_wrap(~Method) +
  labs(title = "Detection Success: Bad Leverage (Top-K Overlap >= 80%)",
       x = "Contamination Proportion", y = "Sample Size (n)") +
  theme_minimal(base_size = 10) -> p5
print(p5)


# ---------- Plot 6: Shape Parameter by n_obs ----------
results %>%
  filter(converged == TRUE) %>%
  ggplot(aes(x = reorder(error_type, shape, median), y = shape, 
             fill = factor(n_obs))) +
  geom_boxplot(outlier.size = 0.3, width = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  coord_flip() +
  labs(title = "MIS Shape Parameter (xi) by Error Distribution and n",
       x = NULL, y = "Shape (xi)", fill = "n") +
  theme_minimal(base_size = 10) -> p6
print(p6)


# ---------- Plot 7: Effective block_count across the grid ----------
results %>%
  group_by(n_obs, contam_prop) %>%
  summarise(
    med_blocks = median(block_count),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = factor(contam_prop), y = factor(n_obs), fill = med_blocks)) +
  geom_tile(color = "white") +
  geom_text(aes(label = round(med_blocks)), size = 4) +
  scale_fill_gradient(low = "#FAC775", high = "#0F6E56") +
  labs(title = "Effective Block Count (Median) Across Grid",
       subtitle = "Requested = 50; capped at floor((n-k)/k)",
       x = "Contamination Proportion", y = "Sample Size (n)",
       fill = "Blocks") +
  theme_minimal(base_size = 10) -> p7
print(p7)

cat("\n=== All plots rendered. Script 02 complete. ===\n")
