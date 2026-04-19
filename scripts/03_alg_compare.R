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
library(ggplot2)

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
# 2. DGP
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
    exact_wins = max_exact > (max_greedy + 1e-8),
    percent_underestimation = ((max_exact - max_greedy) / max_exact) * 100
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

library(ggplot2)
p <- results %>%
  tidyr::pivot_longer(c(max_greedy, max_exact), names_to = "method", values_to = "value") %>%
  ggplot(aes(x = factor(k_size), y = value, fill = method)) +
  geom_boxplot(color = "#3300C0") +
  scale_fill_manual(values = c("mean_exact" = "#75FF6B", "mean_greedy" = "#FFB950"))+
  labs(title = "Exact vs Greedy Mean Influence by Set Size",
       x = "Set Size (k)", y = "Mean Block Maxima") +
  theme_minimal()

ggsave("../output/alg_compare_plot.png", plot = p, width = 7, height = 5, dpi = 300)
ggsave("../output/alg_compare_plot.pdf", plot = p, width = 7, height = 5)

saveRDS(results, "../output/algorithmic_bias_results.rds")



# ------------------------------------------------------------------------------
# 5. Size and Power Curves (Effect Size x-axis, Bootstrap Critical Value)
# ------------------------------------------------------------------------------
DELTA_GRID <- c(0, 2, 5, 10, 20, 50, 100)
N_ITER  <- 200
N_BOOT  <- 99
ALPHA   <- 0.05
k_fixed <- 20
set.seed(20260414)
cat("Starting Power Simulation...\n")

# Pre-allocate a list to store iteration results for memory efficiency
results_list <- vector("list", N_ITER)

for (i in seq_len(N_ITER)) {
  
  if (i %% 10 == 0) cat(sprintf("  Processing Power Iteration %d / %d...\n", i, N_ITER))
  
  # 5a. Generate BASELINE (Null) Data ONCE per iteration
  dat_null <- generate_sim_data(n_obs)
  df_null  <- data.frame(y = dat_null$y, x = dat_null$x,
                         Z1 = dat_null$Z[,1], Z2 = dat_null$Z[,2])
  
  # Identify the target set on the unshifted data
  base_mod   <- lm(y ~ x + Z1 + Z2, data = df_null)
  sens_obj   <- influence::sens(
    base_mod,
    lambda = influence::set_lambda("beta_i", pos = 2, sign = sign(coef(base_mod)[2]))
  )
  target_set <- sens_obj$influence$id[1:k_fixed]
  
  # 5b. FWL on NULL data -> Used exclusively for Bootstrap CV
  fwl_null <- testingMIS:::fwl(y = dat_null$y, X = dat_null$x, Z = dat_null$Z)
  X_null   <- fwl_null[, 2]
  R_null   <- residuals(lm(fwl_null[, 1] ~ X_null - 1))
  
  # 5c. Bootstrap CV from NULL data (Calculated ONCE per iteration)
  boot_greedy <- numeric(N_BOOT)
  boot_exact  <- numeric(N_BOOT)
  
  for (b in seq_len(N_BOOT)) {
    perm_set <- sample(n_obs, k_fixed)
    boot_greedy[b] <- max(abs(testingMIS::dfb_bmx(X_null, R_null, set = perm_set, block_count = block_count)))
    boot_exact[b]  <- max(abs(exact_dfb_bmx(X_null, R_null, set = perm_set, block_count = block_count)))
  }
  
  cv_greedy <- quantile(boot_greedy, 1 - ALPHA)
  cv_exact  <- quantile(boot_exact, 1 - ALPHA)
  
  # 5d. Inner loop over DELTA grid applying shift to the SAME baseline data
  iter_results <- data.frame(delta = DELTA_GRID, rej_greedy = 0, rej_exact = 0)
  
  for (j in seq_along(DELTA_GRID)) {
    delta <- DELTA_GRID[j]
    
    # Apply delta shift -> ALTERNATIVE data
    dat_alt <- dat_null
    dat_alt$y[target_set] <- dat_alt$y[target_set] + delta
    
    # FWL on SHIFTED data -> observed test statistic
    fwl_alt  <- testingMIS:::fwl(y = dat_alt$y, X = dat_alt$x, Z = dat_alt$Z)
    X_alt    <- fwl_alt[, 2]
    R_alt    <- residuals(lm(fwl_alt[, 1] ~ X_alt - 1))
    
    obs_greedy <- max(abs(testingMIS::dfb_bmx(X_alt, R_alt, set = target_set, block_count = block_count)))
    obs_exact  <- max(abs(exact_dfb_bmx(      X_alt, R_alt, set = target_set, block_count = block_count)))
    
    # Record binary rejection decisions
    iter_results$rej_greedy[j] <- as.integer(obs_greedy > cv_greedy)
    iter_results$rej_exact[j]  <- as.integer(obs_exact  > cv_exact)
  }
  
  results_list[[i]] <- iter_results
}

# 5e. Aggregate results across all iterations
power_results <- dplyr::bind_rows(results_list) %>%
  dplyr::group_by(delta) %>%
  dplyr::summarize(
    power_greedy = mean(rej_greedy),
    power_exact  = mean(rej_exact),
    .groups = "drop"
  )

# 5f. Plot
p_power <- power_results %>%
  tidyr::pivot_longer(c(power_greedy, power_exact),
                      names_to = "method", values_to = "rejection_rate") %>%
  dplyr::mutate(method = dplyr::recode(method,
                                       power_greedy = "Greedy",
                                       power_exact  = "Exact")) %>%
  ggplot(aes(x = delta, y = rejection_rate, color = method, linetype = method)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = ALPHA, linetype = "dashed", color = "gray60") +
  annotate("text", x = 0, y = ALPHA + 0.03,
           label = paste0("nominal α = ", ALPHA), hjust = 0, color = "gray50", size = 3) +
  scale_color_manual(values = c("Exact" = "#8FA8C8", "Greedy" = "#A8C8A8")) +
  # REPLACE WITH:
  scale_y_continuous(limits = c(0, NA), labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  labs(title    = "Size and Power Curves: Exact vs Greedy",
       subtitle = paste0("k = ", k_fixed, ", B = ", N_BOOT, " bootstrap draws, ",
                         N_ITER, " iterations per δ"),
       x = "Effect Size (δ)",  y = "Rejection Rate",
       color = "Method", linetype = "Method") +
  theme_minimal()

ggsave("../output/size_power_curves.png", plot = p_power, width = 7, height = 5, dpi = 300)
ggsave("../output/size_power_curves.pdf", plot = p_power, width = 7, height = 5)
cat("Size & Power curves saved.\n")