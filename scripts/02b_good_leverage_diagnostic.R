# ==============================================================================
# File: /scripts/02b_good_leverage_diagnostic.R
# Purpose: Diagnose WHY MIS falsely detects good-leverage points under
#          heavy-tailed errors. Uses FWL decomposition columns from sim_engine.
#
# Outputs:
#   ../output/02b_good_leverage_diagnostic.rds          (full iteration-level data)
#   ../output/tables/02b_tab_fwl_decomposition.tex      (mechanism table)
#   ../output/figures/02b_fig_overlap_by_error.pdf/.png (overlap distribution)
#   ../output/figures/02b_fig_dfb_ratio.pdf/.png        (DFBETA ratio comparison)
#   ../output/figures/02b_fig_lev_res_scatter.pdf/.png  (leverage vs residual driver)
#
# Inputs: ../R/ (all engine scripts)
# ==============================================================================

# 1. Dependencies
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(scales)
library(xtable)

source("../R/helpers_local.R")
source("../R/dgp_factory.R")
source("../R/influence_injector.R")
source("../R/utils_checkpoint.R")
source("../R/exact_dfb_bmx.R")
source("../R/evt_iter_dm.R")
source("../R/dinkelbach_topk.R")
source("../R/sim_engine.R")

# 2. Configuration
diag_params <- list(
  n_iters     = 200,
  n_obs       = 1000,
  k           = 10,
  magnitude   = 10,
  block_count = 50,
  alpha       = 0.05,
  seed        = 20260521
)

set.seed(diag_params$seed)

error_types <- c("normal", "mixed_normal", "skewed_t", "golm",
                 "beta_logistic", "gpd", "contaminated", "pareto")

# Output directories
dir.create("../output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("../output/figures", recursive = TRUE, showWarnings = FALSE)

# Shared theme
theme_thesis <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title       = element_text(face = "bold"),
      strip.background = element_rect(fill = "grey95", colour = "black", linewidth = 0.4),
      strip.text       = element_text(face = "bold"),
      panel.border     = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      panel.grid.minor = element_blank(),
      legend.position  = "bottom"
    )
}

# Error distribution ordering by theoretical tail weight
error_order <- c("normal", "beta_logistic", "mixed_normal", "skewed_t",
                 "contaminated", "golm", "pareto", "gpd")

# ==============================================================================
# 3. Run Simulation
# ==============================================================================
cat("=== 02b Good Leverage Diagnostic ===\n")
cat(sprintf("  n=%d, k=%d, mag=%d, iters=%d\n\n",
            diag_params$n_obs, diag_params$k,
            diag_params$magnitude, diag_params$n_iters))

all_results <- list()

for (err in error_types) {
  cat(sprintf("  %s ... ", err))
  
  iter_results <- purrr::map_dfr(seq_len(diag_params$n_iters), function(iter_id) {
    tryCatch({
      run_mis_iteration(
        iter           = iter_id,
        n              = diag_params$n_obs,
        p              = 1,
        x_type         = "normal",
        error_type     = err,
        outlier_method = "good_leverage",
        k              = diag_params$k,
        magnitude      = diag_params$magnitude,
        block_count    = diag_params$block_count,
        alpha          = diag_params$alpha
      )
    }, error = function(e) {
      warning(sprintf("Iter %d failed for %s: %s", iter_id, err, e$message))
      return(NULL)
    })
  })
  
  all_results[[err]] <- iter_results
  cat(sprintf("Done (%d rows)\n", nrow(iter_results)))
}

results <- bind_rows(all_results)
results$error_type <- factor(results$error_type, levels = error_order)

saveRDS(results, "../output/02b_good_leverage_diagnostic.rds")
cat(sprintf("\nSaved %d rows to ../output/02b_good_leverage_diagnostic.rds\n", nrow(results)))