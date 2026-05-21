# ==============================================================================
# File: script/03_alg_comp.R
# Purpose: Execute dimensional scaling and structural robustness simulations.
#          Compares the Greedy heuristic (testingMIS::dfb_bmx) against
#          Exact Dinkelbach (exact_dfb_bmx) across sample sizes, set sizes,
#          block-count strategies, and five model architectures.
#          Quantifies detection accuracy, EVD convergence, and runtime scaling.
# Inputs:  ../R/03_scaling_dgp.R, ../R/exact_dfb_bmx.R, ../R/evt_iter.R,
#          ../R/evt_iter_dm.R, ../R/utils_checkpoint.R
# Outputs: ../output/temp_03/03_chunk_*.rds -> ../output/03_scaling_results_master.rds
# ==============================================================================

# 1. Load Dependencies & Source Engines
library(dplyr)
library(purrr)
library(evd)
library(future)
library(furrr)

# Source order matters:
#   1. helpers_local.R FIRST (provides fwl, make_blocks, dfbeta_numeric,
#      dfb_bmx, estimate_dfb_evd — replaces testingMIS)
#   2. exact_dfb_bmx BEFORE evt_iter_dm (handoff rule #4)
#   3. fast_sens_topk replaces influence::sens()
source("../R/helpers_local.R")
source("../R/03_scaling_dgp.R")
source("../R/utils_checkpoint.R")
source("../R/fast_sens_topk.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")
source("../R/evt_iter.R")



# Set up parallel processing (leave 2 threads for OS stability)
plan(multisession, workers = 14)
cat("Parallel processing initialized with 14 workers.\n\n")

# ------------------------------------------------------------------------------
# 2. Global Configuration
# ------------------------------------------------------------------------------
sim_params <- list(
  n_iters   = 100,
  magnitude = 10,
  seed      = 20260421
)

set.seed(sim_params$seed)

# Minimum non-set observations per block.
# When B is too large relative to (N - k), each block has too few observations
# and exact Dinkelbach produces near-constant block maxima that cause GEV MLE
# to fail (flat likelihood surface). Diagnostic 04 showed:
#   - B=100, N=500, k>=5:  0% exact convergence  (5 obs/block)
#   - B=100, N=500, k=1:  77% convergence        (5 obs/block)
#   - B=20,  N=500, k=5:  86% convergence        (25 obs/block)
#   - B=sqrt, N=5000:      stable until k=50
# Floor of 30 obs/block recovers convergence without altering the EVT logic.
# This does NOT modify evt_iter_dm.R (shared with Script 02).
MIN_OBS_PER_BLOCK <- 30

resolve_block_count <- function(B_type, N, k) {
  B_raw <- if (B_type == "sqrt") floor(sqrt(N)) else as.numeric(B_type)
  B_max <- floor((N - k) / MIN_OBS_PER_BLOCK)
  B_max <- max(B_max, 3L)  # absolute floor: need >=3 blocks for GEV
  min(B_raw, B_max)
}

# ------------------------------------------------------------------------------
# 3. The Safe Injection Adapter (Preserves Interaction Geometry)
# ------------------------------------------------------------------------------
inject_safe_outliers <- function(dgp_res, k, magnitude = 4) {
  df <- dgp_res$data
  N  <- nrow(df)
  outlier_idx <- sample(1:N, k, replace = FALSE)
  
  # 3a. Shift main covariates robustly
  shift_x <- max(mad(df$X, constant = 1.4826), 1e-4) * magnitude
  df$X[outlier_idx] <- df$X[outlier_idx] + shift_x
  
  if ("M" %in% names(df)) {
    if (length(unique(df$M)) <= 2) {
      # Binary M: force outlier rows into the active category (M=1)
      # so they actually live in the interaction subspace
      df$M[outlier_idx] <- 1L
    } else {
      shift_m <- max(mad(df$M, constant = 1.4826), 1e-4) * magnitude
      df$M[outlier_idx] <- df$M[outlier_idx] + shift_m
    }
  }
  if ("W" %in% names(df)) {
    shift_w <- max(mad(df$W, constant = 1.4826), 1e-4) * magnitude
    df$W[outlier_idx] <- df$W[outlier_idx] + shift_w
  }
  
  # 3b. For complex architecture, shift a few Z columns too so outliers are
  #     actually influential in the high-dimensional covariate space
  z_cols <- grep("^(Z_matrix\\.|X[0-9]+|Z[0-9]+)$", names(df), value = TRUE)
  if (length(z_cols) > 0) {
    for (zc in head(z_cols, 3)) {
      shift_z <- max(mad(df[[zc]], constant = 1.4826), 1e-4) * magnitude
      df[[zc]][outlier_idx] <- df[[zc]][outlier_idx] + shift_z
    }
  }
  
  # 3c. Recalculate Y to forcefully mask the target parameter
  clean_mod  <- lm(dgp_res$formula, data = df[-outlier_idx, ])
  expected_y <- predict(clean_mod, newdata = df[outlier_idx, ])
  scale_y    <- max(mad(df$y, constant = 1.4826), 1e-4)
  df$y[outlier_idx] <- expected_y - (sign(expected_y) * scale_y * magnitude * 2)
  
  dgp_res$data <- df
  dgp_res$true_outliers <- outlier_idx
  return(dgp_res)
}

# ------------------------------------------------------------------------------
# 4. Single-Iteration Worker (Called Inside future_map_dfr)
# ------------------------------------------------------------------------------
run_scaling_iteration <- function(iter_id, N, k, B, architecture,
                                  magnitude, iter_seed, rho = NA) {
  set.seed(iter_seed)
  
  tryCatch({
    
    # A. Generate Data
    dgp_clean <- generate_scaling_dgp(
      N = N, architecture = architecture,
      rho = if (!is.null(rho) && !is.na(rho)) rho else 0.8
    )
    
    # B. Inject Adversarial Outliers
    dgp_poisoned <- inject_safe_outliers(dgp_clean, k = k, magnitude = magnitude)
    
    # C. Fit base model on poisoned data
    mod  <- lm(dgp_poisoned$formula, data = dgp_poisoned$data)
    tpos <- dgp_poisoned$target_pos
    
    # D. Detect the influential set via influence::sens()
    detected_set <- fast_sens_topk(mod, 
      pos = tpos, sign = sign(coef(mod)[tpos]), k = k)
    
    # E. Extract FWL components: x = target column, Z = everything else
    #    Z must NEVER contain x (handoff rule #1)
    X_full   <- model.matrix(dgp_poisoned$formula, dgp_poisoned$data)
    x_target <- X_full[, tpos]
    Z_fwl    <- X_full[, -tpos, drop = FALSE]
    y_vec    <- dgp_poisoned$data$y
    
    # F. Run Greedy EVT
    t_greedy_start <- Sys.time()
    res_greedy <- evt_iter(
      y = y_vec, x = x_target, Z = Z_fwl,
      set = detected_set, block_count = B
    )
    t_greedy <- as.numeric(difftime(Sys.time(), t_greedy_start, units = "secs"))
    
    # G. Dual-sign detection (Oracle: pick the direction that finds more outliers)
    detected_pos <- fast_sens_topk(mod, pos = tpos, sign = 1, k = k)
    detected_neg <- fast_sens_topk(mod, pos = tpos, sign = -1, k = k)
    
    overlap_pos <- length(intersect(detected_pos, dgp_poisoned$true_outliers))
    overlap_neg <- length(intersect(detected_neg, dgp_poisoned$true_outliers))
    
    if (overlap_neg > overlap_pos) {
      detected_set_exact <- detected_neg
    } else {
      detected_set_exact <- detected_pos
    }
    
    # H. Run Exact Dinkelbach EVT on the best-direction set
    t_exact_start <- Sys.time()
    res_exact <- evt_iter_dm(
      y = y_vec, x = x_target, Z = Z_fwl,
      set = detected_set_exact, block_count = B
    )
    t_exact <- as.numeric(difftime(Sys.time(), t_exact_start, units = "secs"))
    
    # I. Detection metrics (greedy = single-sign, exact = oracle dual-sign)
    detection_rate       <- length(intersect(detected_set, dgp_poisoned$true_outliers)) / k
    detection_rate_exact <- length(intersect(detected_set_exact, dgp_poisoned$true_outliers)) / k
    
    data.frame(
      iter               = iter_id,
      detection_rate       = detection_rate,
      detection_rate_exact = detection_rate_exact,
      p_greedy         = res_greedy$p_value,
      converged_greedy = res_greedy$converged,
      p_exact          = res_exact$p_value,
      converged_exact  = res_exact$converged,
      cpu_greedy   = t_greedy,
      cpu_exact    = t_exact,
      error_msg    = NA_character_,
      stringsAsFactors = FALSE
    )
    
  }, error = function(e) {
    warning(sprintf("Iter %d failed: %s", iter_id, e$message))
    data.frame(
      iter               = iter_id,
      detection_rate       = NA_real_,
      detection_rate_exact = NA_real_,
      p_greedy         = NA_real_,
      converged_greedy = FALSE,
      p_exact          = NA_real_,
      converged_exact  = FALSE,
      cpu_greedy   = NA_real_,
      cpu_exact    = NA_real_,
      error_msg    = e$message,
      stringsAsFactors = FALSE
    )
  })
}

# ------------------------------------------------------------------------------
# 5. Define the Scenario Grid (iter is the INNER parallel dimension)
# ------------------------------------------------------------------------------
# --- Base grid ---
grid_base <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = c("simple", "complex", "interaction",
                   "triple_interaction", "nonlinear_nuisance",
                   "sparse_binary_interaction",
                   "polynomial_interaction"),
  stringsAsFactors = FALSE
)

# high_k_interaction with k up to 50 ---
grid_high_k <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20, 30, 50),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = "high_k_interaction",
  stringsAsFactors = FALSE
)

# --- Rho sweep for collinear_interaction ---
grid_collinear <- expand.grid(
  N            = c(500, 1000, 2000, 5000),
  k            = c(1, 3, 5, 10, 15, 20),
  B_type       = c("20", "50", "100", "sqrt"),
  architecture = "collinear_interaction",
  rho          = c(0.5, 0.7, 0.85, 0.95),
  stringsAsFactors = FALSE
)

# --- Combine and filter: drop rows where k/N > 0.05 ---
grid_base$rho    <- NA_real_
grid_high_k$rho  <- NA_real_
scenario_grid <- rbind(grid_base, grid_high_k, grid_collinear)
scenario_grid <- scenario_grid[scenario_grid$k / scenario_grid$N <= 0.05, ]
scenario_grid <- scenario_grid[order(scenario_grid$architecture,
                                     scenario_grid$N,
                                     scenario_grid$k), ]
rownames(scenario_grid) <- NULL

n_scenarios <- nrow(scenario_grid)
n_iters     <- sim_params$n_iters
total_rows  <- n_scenarios * n_iters

dir.create("../output/temp_03", recursive = TRUE, showWarnings = FALSE)

cat(sprintf(paste0(
  "Starting 03 Scaling Suite (Extended Boundary Tests).\n",
  "  N    = {500, 1000, 2000, 5000}\n",
  "  k    = {1, 3, 5, 10, 15, 20} + {30, 50} for high_k_interaction\n",
  "  B    = {20, 50, 100, sqrt(N)}\n",
  "  Arch = {simple, complex, interaction, triple_interaction, nonlinear_nuisance,\n",
  "          collinear_interaction (rho sweep), sparse_binary_interaction,\n",
  "          polynomial_interaction, high_k_interaction}\n",
  "  k/N filter: <= 0.05\n",
  "  Scenarios: %d | Iterations per scenario: %d | Total rows: %d\n\n"),
  n_scenarios, n_iters, total_rows))

# ------------------------------------------------------------------------------
# 6. Orchestrator Loop (Sequential scenarios, parallel iterations)
# ------------------------------------------------------------------------------
for (i in seq_len(n_scenarios)) {
  
  sc <- scenario_grid[i, ]
  chunk_file <- sprintf("../output/temp_03/03_chunk_%04d.rds", i)
  
  # Checkpoint: skip if already computed
  if (is_computed(chunk_file)) {
    if (i %% 20 == 0) cat(sprintf("[%04d/%04d] Skipping (cached): N=%d k=%d B=%s arch=%s\n",
                                  i, n_scenarios, sc$N, sc$k, sc$B_type, sc$architecture))
    next
  }
  
  cat(sprintf("[%04d/%04d] Computing: N=%d  k=%d  B=%s  arch=%s ... ",
              i, n_scenarios, sc$N, sc$k, sc$B_type, sc$architecture))
  
  # Resolve block count with adaptive cap (Fix C from diagnostic 04)
  B <- resolve_block_count(sc$B_type, sc$N, sc$k)
  B_raw <- if (sc$B_type == "sqrt") floor(sqrt(sc$N)) else as.numeric(sc$B_type)
  B_was_capped <- (B < B_raw)
  if (B_was_capped) {
    cat(sprintf("[B capped: %d->%d] ", B_raw, B))
  }
  
  # Parallel inner loop over iterations
  scenario_results <- furrr::future_map_dfr(
    seq_len(n_iters), function(iter_id) {
      
      iter_seed <- sim_params$seed + (i - 1) * n_iters + iter_id
      
      run_scaling_iteration(
        iter_id      = iter_id,
        N            = sc$N,
        k            = sc$k,
        B            = B,
        architecture = sc$architecture,
        magnitude    = sim_params$magnitude,
        iter_seed    = iter_seed,
        rho          = sc$rho
      )
      
    }, .options = furrr_options(seed = TRUE)
  )
  
  # Attach scenario-level columns
  scenario_results$N            <- sc$N
  scenario_results$k            <- sc$k
  scenario_results$B_type       <- sc$B_type
  scenario_results$B_actual     <- B
  scenario_results$B_raw        <- B_raw
  scenario_results$B_capped     <- B_was_capped
  scenario_results$architecture <- sc$architecture
  scenario_results$rho          <- sc$rho
  
  safe_save_rds(scenario_results, chunk_file)
  cat("Done.\n")
}

# ------------------------------------------------------------------------------
# 7. Compile Final Artifact
# ------------------------------------------------------------------------------
cat("\nAll scenarios completed. Assembling final dataset...\n")
final_data <- compile_checkpoints(
  temp_dir          = "../output/temp_03",
  pattern           = "^03_chunk_.*\\.rds$",
  final_output_path = "../output/03_scaling_results_master.rds",
  clear_temp        = FALSE
)
cat("Script 03 execution finished successfully.\n")

# ------------------------------------------------------------------------------
# 8. Sanity Check & Reporting Validation
# ------------------------------------------------------------------------------
res <- final_data

cat("\n============================\n")
cat("  03 SANITY CHECK\n")
cat("============================\n")

cat("\n=== 1. BASIC INTEGRITY ===\n")
cat(sprintf("Rows: %d / %d (Missing: %d)\n",
            nrow(res), total_rows, total_rows - nrow(res)))
cat("Missing CPU (greedy):", sum(is.na(res$cpu_greedy)), "\n")
cat("Missing CPU (exact):",  sum(is.na(res$cpu_exact)), "\n")
cat("Missing detection rates:", sum(is.na(res$detection_rate)), "\n")

# B-capping summary (Fix C diagnostic)
if ("B_capped" %in% names(res)) {
  n_capped <- sum(res$B_capped, na.rm = TRUE)
  cat(sprintf("\nB-cap applied: %d / %d rows (%.1f%%)\n",
              n_capped, nrow(res), n_capped / nrow(res) * 100))
  if (n_capped > 0) {
    cat("B-capped scenarios (B_type -> B_actual):\n")
    capped_summary <- unique(res[res$B_capped == TRUE,
                                 c("N","k","B_type","B_raw","B_actual")])
    capped_summary <- capped_summary[order(capped_summary$N, capped_summary$k), ]
    print(capped_summary, row.names = FALSE)
  }
}

# Grid coverage: every (N, k, B_type, architecture) should have n_iters rows
coverage <- res %>%
  group_by(N, k, B_type, architecture) %>%
  summarise(n = n(), .groups = "drop")
incomplete <- coverage %>% filter(n < n_iters)
if (nrow(incomplete) > 0) {
  cat(sprintf("WARNING: %d scenarios have fewer than %d iterations:\n",
              nrow(incomplete), n_iters))
  print(incomplete)
} else {
  cat(sprintf("All %d scenarios have full %d iterations.\n", n_scenarios, n_iters))
}

cat("\n=== 2. STABILITY & ERROR RATES ===\n")
total_errors <- sum(!is.na(res$error_msg))
cat("Total algorithm crashes caught:", total_errors, "\n")

if (total_errors > 0) {
  cat("\nErrors by architecture:\n")
  print(table(res$architecture, Has_Error = !is.na(res$error_msg)))
  
  cat("\nErrors by (architecture, k):\n")
  print(table(res$architecture, res$k, !is.na(res$error_msg)))
  
  cat("\nSample error messages:\n")
  print(head(unique(na.omit(res$error_msg)), 3))
} else {
  cat("Perfect stability. No crashes.\n")
}

cat("\n=== 3. CONVERGENCE ===\n")
cat(sprintf("Greedy EVT converged: %.1f%%\n",
            mean(res$converged_greedy, na.rm = TRUE) * 100))
cat(sprintf("Exact  EVT converged: %.1f%%\n",
            mean(res$converged_exact, na.rm = TRUE) * 100))

conv_by_arch <- res %>%
  group_by(architecture) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
print(conv_by_arch)

# Rho sweep: convergence by rho (collinear_interaction only)
conv_by_rho <- res %>%
  filter(architecture == "collinear_interaction", !is.na(rho)) %>%
  group_by(rho) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact,  na.rm = TRUE) * 100, 1),
    det_greedy  = round(mean(detection_rate,   na.rm = TRUE), 3),
    det_exact   = round(mean(detection_rate_exact, na.rm = TRUE), 3),
    .groups = "drop"
  )
if (nrow(conv_by_rho) > 0) {
  cat("\nCollinear interaction — by rho:\n")
  print(conv_by_rho)
}

conv_by_k <- res %>%
  group_by(k) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("\nConvergence by k:\n")
print(conv_by_k)

conv_by_N <- res %>%
  group_by(N) %>%
  summarise(
    conv_greedy = round(mean(converged_greedy, na.rm = TRUE) * 100, 1),
    conv_exact  = round(mean(converged_exact, na.rm = TRUE) * 100, 1),
    .groups = "drop"
  )
cat("\nConvergence by N:\n")
print(conv_by_N)

cat("\n=== 4. DETECTION PREVIEW ===\n")
det_table <- res %>%
  filter(!is.na(detection_rate)) %>%
  group_by(architecture, k) %>%
  summarise(
    det_greedy = round(mean(detection_rate), 3),
    det_exact  = round(mean(detection_rate_exact, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(architecture, k)
print(det_table, n = 30)

# Warn on trivial results
all_zero <- all(det_table$det_greedy == 0) && all(det_table$det_exact == 0)
all_one  <- all(det_table$det_greedy == 1) && all(det_table$det_exact == 1)
if (all_zero) cat("WARNING: All detection rates are 0 — algorithm may be blind.\n")
if (all_one)  cat("WARNING: All detection rates are 1 — problem may be too easy.\n")

# Detection by N to see if larger samples help
det_by_N <- res %>%
  filter(!is.na(detection_rate)) %>%
  group_by(N, k) %>%
  summarise(
    det_greedy = round(mean(detection_rate), 3),
    det_exact  = round(mean(detection_rate_exact, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  arrange(N, k)
cat("\nDetection by (N, k):\n")
print(det_by_N, n = 25)

cat("\n=== 5. REPORTING SANITY ===\n")
cat("Greedy detection rates outside [0,1]:",
    sum(res$detection_rate < 0 | res$detection_rate > 1, na.rm = TRUE), "\n")
cat("Exact  detection rates outside [0,1]:",
    sum(res$detection_rate_exact < 0 | res$detection_rate_exact > 1, na.rm = TRUE), "\n")
cat("Negative CPU (greedy):", sum(res$cpu_greedy < 0, na.rm = TRUE), "\n")
cat("Negative CPU (exact):",  sum(res$cpu_exact < 0, na.rm = TRUE), "\n")
cat("Implausibly large CPU (>600s, greedy):",
    sum(res$cpu_greedy > 600, na.rm = TRUE), "\n")
cat("Implausibly large CPU (>600s, exact):",
    sum(res$cpu_exact > 600, na.rm = TRUE), "\n")

# B_actual consistency (now accounts for adaptive capping)
res <- res %>%
  mutate(B_expected = mapply(resolve_block_count, B_type, N, k))
cat("B_actual mismatches:", sum(res$B_actual != res$B_expected, na.rm = TRUE), "\n")
res$B_expected <- NULL

# P-value sanity (converged rows should have p in [0,1])
p_bad_greedy <- res %>% filter(converged_greedy == TRUE & (p_greedy < 0 | p_greedy > 1))
p_bad_exact  <- res %>% filter(converged_exact  == TRUE & (p_exact  < 0 | p_exact  > 1))
cat("Greedy p-values outside [0,1] (converged):", nrow(p_bad_greedy), "\n")
cat("Exact  p-values outside [0,1] (converged):", nrow(p_bad_exact), "\n")

# NA in converged rows
na_greedy <- res %>% filter(converged_greedy == TRUE & is.na(p_greedy))
na_exact  <- res %>% filter(converged_exact  == TRUE & is.na(p_exact))
if (nrow(na_greedy) > 0) cat("WARNING:", nrow(na_greedy),
                             "converged greedy rows have NA p-values!\n")
if (nrow(na_exact) > 0)  cat("WARNING:", nrow(na_exact),
                             "converged exact rows have NA p-values!\n")

cat("\n=== 6. TIMING PROFILE ===\n")
timing_by_N <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(N) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("By sample size (N):\n")
print(timing_by_N)

timing_by_k <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(k) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("\nBy set size (k):\n")
print(timing_by_k)

timing_by_arch <- res %>%
  filter(!is.na(cpu_greedy) & !is.na(cpu_exact)) %>%
  group_by(architecture) %>%
  summarise(
    med_greedy = round(median(cpu_greedy), 3),
    med_exact  = round(median(cpu_exact), 3),
    ratio      = round(median(cpu_exact) / max(median(cpu_greedy), 1e-6), 1),
    .groups = "drop"
  )
cat("\nBy architecture:\n")
print(timing_by_arch)

# Cross-tabulation: N x k (the most informative scaling view)
timing_Nk <- res %>%
  filter(!is.na(cpu_exact)) %>%
  group_by(N, k) %>%
  summarise(
    med_exact = round(median(cpu_exact), 3),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = k, values_from = med_exact,
                     names_prefix = "k=")
cat("\nExact median CPU (sec) — N x k:\n")
print(timing_Nk)

cat("\nSanity check complete.\n")