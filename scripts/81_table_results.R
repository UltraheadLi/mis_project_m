# ==============================================================================
# File: script/81_table_results.R
# Purpose: Generates main-text LaTeX tables from master simulation datasets.
#          - Table 1: Diagnostic Breakdown (MIS vs Classical under Bad Leverage)
#          - Table 2: Algorithmic Performance (Exact vs Greedy by Architecture)
# Inputs:  ../output/02_final_distributions.rds
#          ../output/03_scaling_results_master.rds
# Outputs: ../output/tables/tab_01_diagnostics.tex
#          ../output/tables/tab_02_algorithmic_scaling.tex
# ==============================================================================

library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
library(tools)

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
input_02 <- "../output/02_final_distributions.rds"
input_03 <- "../output/03_scaling_results_master.rds"
out_dir <- "../output/tables"

if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Safely load data (using dummy data if missing for testing)
res_02 <- if (file.exists(input_02)) readRDS(input_02) else data.frame()
res_03 <- if (file.exists(input_03)) readRDS(input_03) else data.frame()

# ------------------------------------------------------------------------------
# DEFENSIVE COLUMN GUARD FOR res_03
# ------------------------------------------------------------------------------
if (nrow(res_03) > 0) {
  res_03_required <- c("converged_exact", "converged_greedy",
                       "detection_rate_exact", "detection_rate",
                       "p_exact", "p_greedy", "cpu_exact", "cpu_greedy")
  missing_03 <- setdiff(res_03_required, names(res_03))
  if (length(missing_03) > 0) {
    message(
      "WARNING [81_table_results.R]: res_03 is missing columns — likely stale chunks.\n",
      "  Missing: ", paste(missing_03, collapse = ", "), "\n",
      "  Fix: delete ../output/temp_03/ and re-run Script 03.\n",
      "  Patching with NA/FALSE for now so tables can still be generated."
    )
    for (col in missing_03) {
      if (col %in% c("converged_exact", "converged_greedy")) {
        res_03[[col]] <- FALSE
      } else {
        res_03[[col]] <- NA_real_
      }
    }
  }
}

# Helper: Clean names for publication
clean_names <- function(x) {
  x <- gsub("_", " ", x)
  x <- tools::toTitleCase(x)
  x <- gsub("Gpd", "GPD", x)
  x <- gsub("Golm", "GOLM", x)
  x <- gsub("High K", "High-k", x)
  return(x)
}

# Helper: Standard academic kable wrapper
thesis_table <- function(df, caption, label, align = NULL) {
  kable(df, format = "latex", booktabs = TRUE, caption = caption, label = label, align = align, linesep = "") %>%
    kable_styling(latex_options = c("hold_position"))
}

format_metrics <- function(df) {
  df %>%
    mutate(
      across(starts_with("Conv.") | starts_with("Det."), ~ sprintf("%.1f%%", .x * 100)),
      across(starts_with("CPU"), ~ sprintf("%.3f", .x))
    )
}

# ------------------------------------------------------------------------------
# Section 02 Tables
# ------------------------------------------------------------------------------

# Table 1: Detection Comparison (Focused on Bad Leverage)
tab_01_data <- res_02 %>%
  filter(outlier_method == "bad_leverage") %>%
  group_by(error_type) %>%
  summarise(
    `Cook's D` = mean(detect_cooks, na.rm = TRUE),
    Leverage   = mean(detect_lev, na.rm = TRUE),
    DFBETAS    = mean(detect_dfbetas, na.rm = TRUE),
    MIS        = mean(detection_success, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(error_type = clean_names(error_type)) %>%
  rename(`Error Distribution` = error_type) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.1f%%", . * 100)))

latex_tab_01 <- thesis_table(tab_01_data, 
                             caption = "Detection Accuracy: MIS vs. Classical Diagnostics (Masking Scenario)", 
                             label = "tab:detection_masking") %>%
  column_spec(5, bold = TRUE) # Highlight MIS

writeLines(latex_tab_01, file.path(out_dir, "tab_02a_diagnostics.tex"))


# Table 2: EVT Convergence Rate
tab_02_data <- res_02 %>%
  group_by(outlier_method, error_type) %>%
  summarise(Convergence = mean(converged, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    error_type = clean_names(error_type),
    outlier_method = clean_names(outlier_method)
  ) %>%
  pivot_wider(names_from = outlier_method, values_from = Convergence) %>%
  rename(`Error Distribution` = error_type) %>%
  mutate(across(where(is.numeric), ~ sprintf("%.1f%%", . * 100)))

latex_tab_02 <- thesis_table(tab_02_data, 
                             caption = "EVT Convergence Rates Across Outlier Scenarios", 
                             label = "tab:evt_convergence")

writeLines(latex_tab_02, file.path(out_dir, "tab_02b_convergence.tex"))


# ------------------------------------------------------------------------------
# Section 03 Tables
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Table 3a: Impact of Sample Size (N)
# Aggregated across k and B_type
# ------------------------------------------------------------------------------
tab_N <- res_03 %>%
  group_by(architecture, N) %>%
  summarise(
    `Conv. Exact`  = mean(converged_exact, na.rm = TRUE),
    `Conv. Greedy` = mean(converged_greedy, na.rm = TRUE),
    `Det. Exact`   = mean(detection_rate_exact, na.rm = TRUE),
    `Det. Greedy`  = mean(detection_rate, na.rm = TRUE),
    `CPU Exact`    = median(cpu_exact, na.rm = TRUE),
    `CPU Greedy`   = median(cpu_greedy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(architecture = clean_names(architecture)) %>%
  filter(N %in% c(500, 1000, 5000)) %>% # Display representative N values
  rename(Architecture = architecture, `Sample Size ($N$)` = N) %>%
  format_metrics()

latex_tab_N <- thesis_table(tab_N, 
                            caption = "Algorithmic Scaling by Sample Size ($N$)", 
                            label = "tab:scaling_N")

writeLines(latex_tab_N, file.path(out_dir, "tab_03a_scaling_N.tex"))

# ------------------------------------------------------------------------------
# Table 3b: Impact of Set Size (k)
# Aggregated across N and B_type
# ------------------------------------------------------------------------------
tab_k <- res_03 %>%
  group_by(architecture, k) %>%
  summarise(
    `Conv. Exact`  = mean(converged_exact, na.rm = TRUE),
    `Conv. Greedy` = mean(converged_greedy, na.rm = TRUE),
    `Det. Exact`   = mean(detection_rate_exact, na.rm = TRUE),
    `Det. Greedy`  = mean(detection_rate, na.rm = TRUE),
    `CPU Exact`    = median(cpu_exact, na.rm = TRUE),
    `CPU Greedy`   = median(cpu_greedy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(architecture = clean_names(architecture)) %>%
  filter(k %in% c(5, 10, 20)) %>% # Display representative k values
  rename(Architecture = architecture, `Set Size ($k$)` = k) %>%
  format_metrics()

latex_tab_k <- thesis_table(tab_k, 
                            caption = "Algorithmic Scaling by Influential Set Size ($k$)", 
                            label = "tab:scaling_k")

writeLines(latex_tab_k, file.path(out_dir, "tab_03b_scaling_k.tex"))

# ------------------------------------------------------------------------------
# Table 3c: Impact of Block Count (B_type)
# Aggregated across N and k
# ------------------------------------------------------------------------------
tab_B <- res_03 %>%
  group_by(architecture, B_type) %>%
  summarise(
    `Conv. Exact`  = mean(converged_exact, na.rm = TRUE),
    `Conv. Greedy` = mean(converged_greedy, na.rm = TRUE),
    `Det. Exact`   = mean(detection_rate_exact, na.rm = TRUE),
    `Det. Greedy`  = mean(detection_rate, na.rm = TRUE),
    `CPU Exact`    = median(cpu_exact, na.rm = TRUE),
    `CPU Greedy`   = median(cpu_greedy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(architecture = clean_names(architecture)) %>%
  # Optionally filter B_type if needed, or leave all (20, 50, 100, sqrt)
  rename(Architecture = architecture, `Blocks ($B$)` = B_type) %>%
  format_metrics()

latex_tab_B <- thesis_table(tab_B, 
                            caption = "EVT Sensitivity by Block Count ($B$)", 
                            label = "tab:scaling_B")

writeLines(latex_tab_B, file.path(out_dir, "tab_03c_scaling_B.tex"))

cat("LaTeX tables successfully generated in", out_dir, "\n")