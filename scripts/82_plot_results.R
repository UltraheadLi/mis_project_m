# ==============================================================================
# File: script/82_plot_results.R (Refactored)
# Purpose: Generates main-text and appendix figures with standard academic formatting.
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(tools)

# ------------------------------------------------------------------------------
# 0. Global Setup & Academic Theme
# ------------------------------------------------------------------------------
out_dir <- "../output/figures"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

input_02 <- "../output/02_final_distributions.rds"
input_03 <- "../output/03_scaling_results_master.rds"

# Okabe-Ito color-blind-safe palette
col_exact  <- "#0072B2"
col_greedy <- "#D55E00"

# Centralized Academic Theme
theme_thesis <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title       = element_text(face = "bold"),
      legend.position  = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#E8EAEB", color = "black", linewidth = 0.5),
      strip.text       = element_text(face = "bold", size = base_size - 2),
      panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# ------------------------------------------------------------------------------
# 1. Data Loading (with robust fallback)
# ------------------------------------------------------------------------------
load_or_simulate_02 <- function(path) {
  if (file.exists(path)) return(readRDS(path))
  warning("File 02 not found. Generating dummy data for testing.")
  expand.grid(
    error_type = c("normal", "mixed_normal", "skewed_t", "golm", "beta_logistic", "gpd", "contaminated", "pareto"),
    outlier_method = c("bad_leverage", "good_leverage", "vertical_outlier"),
    iter = 1:50
  ) %>% mutate(
    detection_success = runif(n(), 0.85, 1.0),
    detect_cooks      = runif(n(), 0.10, 0.40),
    detect_lev        = runif(n(), 0.20, 0.50),
    detect_dfbetas    = runif(n(), 0.30, 0.60),
    shape             = rnorm(n(), mean = 0.1, sd = 0.05),
    converged         = TRUE
  )
}

load_or_simulate_03 <- function(path) {
  if (file.exists(path)) return(readRDS(path))
  warning("File 03 not found. Generating dummy data for testing.")
  expand.grid(
    N = c(500, 1000, 2000, 5000), k = c(1, 3, 5, 10, 15, 20),
    B_actual = c(20, 50, 100), B_type = c("20", "50", "100", "sqrt"),
    architecture = c("simple", "complex", "interaction", "collinear_interaction", "triple_interaction", "nonlinear_nuisance"),
    iter = 1:5
  ) %>% mutate(
    detection_rate_exact = runif(n(), 0.90, 1.0),
    detection_rate       = pmax(0, runif(n(), 0.5, 0.9) - (k * 0.01)),
    converged_exact      = rbinom(n(), 1, 0.95),
    converged_greedy     = rbinom(n(), 1, 0.90),
    cpu_greedy           = runif(n(), 0.02, 0.08),
    cpu_exact            = (N / 1000) * (k / 5) * runif(n(), 0.8, 1.2)
  )
}

res_02 <- load_or_simulate_02(input_02)
res_03 <- load_or_simulate_03(input_03)

# ------------------------------------------------------------------------------
# DEFENSIVE COLUMN GUARD FOR res_03
# ------------------------------------------------------------------------------
if (nrow(res_03) > 0) {
  res_03_required <- c("converged_exact", "converged_greedy",
                       "detection_rate_exact", "detection_rate",
                       "cpu_exact", "cpu_greedy")
  missing_03 <- setdiff(res_03_required, names(res_03))
  if (length(missing_03) > 0) {
    message(
      "WARNING [82_plot_results.R]: res_03 is missing columns — likely stale chunks.\n",
      "  Missing: ", paste(missing_03, collapse = ", "), "\n",
      "  Fix: delete ../output/temp_03/ and re-run Script 03.\n",
      "  Patching with NA/FALSE so plots can still render (with NAs)."
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

# ------------------------------------------------------------------------------
# 2. Figure Generation
# ------------------------------------------------------------------------------

# --- Figure 1: Detection Heatmap ---
plot_heatmap <- res_02 %>%
  filter(outlier_method %in% c("bad_leverage")) %>%
  group_by(error_type, outlier_method) %>%
  summarise(across(c(detection_success, detect_cooks, detect_lev, detect_dfbetas), mean, na.rm = TRUE), .groups = "drop") %>%
  rename(MIS = detection_success, `Cook's D` = detect_cooks, Leverage = detect_lev, DFBETAS = detect_dfbetas) %>%
  pivot_longer(-c(error_type, outlier_method), names_to = "Method", values_to = "Rate") %>%
  mutate(Scenario = tools::toTitleCase(gsub("_", " ", outlier_method)))

p_heatmap <- ggplot(plot_heatmap, aes(x = Method, y = error_type, fill = Rate)) +
  geom_tile(color = "white") +
  geom_text(aes(label = percent(Rate, accuracy = 1)), size = 3.5, fontface = "bold") +
  facet_wrap(~ Scenario) +
  scale_fill_gradient2(low = col_greedy, mid = "#FFFFBF", high = col_exact, midpoint = 0.5, labels = percent) +
  labs(x = NULL, y = "Error Distribution") +
  theme_thesis() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- Figure 3: Greedy Trap (Bar) ---
p_bar <- res_03 %>%
  group_by(architecture) %>%
  summarise(Exact = mean(detection_rate_exact, na.rm=TRUE), Greedy = mean(detection_rate, na.rm=TRUE)) %>%
  pivot_longer(-architecture, names_to = "Algorithm", values_to = "Rate") %>%
  mutate(architecture = tools::toTitleCase(gsub("_", " ", architecture))) %>%
  ggplot(aes(x = reorder(architecture, -Rate), y = Rate, fill = Algorithm)) +
  geom_col(position = position_dodge(0.8), width = 0.7, color = "black") +
  scale_y_continuous(labels = percent_format()) +
  scale_fill_manual(values = c("Exact" = col_exact, "Greedy" = col_greedy)) +
  labs(x = NULL, y = "Mean Detection Rate") +
  theme_thesis() + theme(axis.text.x = element_text(angle = 20, hjust = 1))

# --- Figures 7 & 8: Convergence and Detection Grids ---
# Pre-process grid data to reduce redundancy
grid_data <- res_03 %>%
  group_by(N, B_type, k) %>%
  summarise(
    Conv_Exact = mean(converged_exact, na.rm = TRUE),
    Conv_Greedy = mean(converged_greedy, na.rm = TRUE),
    Det_Exact = mean(detection_rate_exact, na.rm = TRUE),
    Det_Greedy = mean(detection_rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    N_label = factor(paste("N =", N), levels = paste("N =", sort(unique(N)))),
    B_label = factor(paste("Blocks =", B_type), levels = paste("Blocks =", c("20", "50", "100", "sqrt")))
  )

plot_grid <- function(df, y_exact, y_greedy, title, y_label) {
  df %>%
    pivot_longer(cols = c({{y_exact}}, {{y_greedy}}), names_to = "Algorithm", values_to = "Value") %>%
    mutate(Algorithm = gsub(".*_", "", Algorithm)) %>%
    ggplot(aes(x = k, y = Value, color = Algorithm)) +
    geom_line(linewidth = 1) + geom_point(size = 2) +
    facet_grid(B_label ~ N_label) +
    scale_y_continuous(labels = percent_format(), limits = c(0, 1.05)) +
    scale_color_manual(values = c("Exact" = col_exact, "Greedy" = col_greedy)) +
    labs(x = "Influential Set Size (k)", y = y_label) +
    theme_thesis()
}

p_conv <- plot_grid(grid_data, Conv_Exact, Conv_Greedy, "EVT Convergence Rate", "Convergence Rate")
p_det  <- plot_grid(grid_data, Det_Exact, Det_Greedy, "MIS Detection Rate", "Detection Rate")

# ------------------------------------------------------------------------------
# 3. Output Generation
# ------------------------------------------------------------------------------
ggsave(file.path(out_dir, "fig_02_detection_heatmap.pdf"), p_heatmap, width = 12, height = 6)
ggsave(file.path(out_dir, "fig_03_greedy_bar.pdf"), p_bar, width = 10, height = 6)
ggsave(file.path(out_dir, "fig_03_convergence_grid.pdf"), p_conv, width = 10, height = 7)
ggsave(file.path(out_dir, "fig_03_detection_grid.pdf"), p_det, width = 10, height = 7)

ggsave(file.path(out_dir, "fig_02_detection_heatmap.png"), p_heatmap, width = 12, height = 6)
ggsave(file.path(out_dir, "fig_03_greedy_bar.png"), p_bar, width = 10, height = 6)
ggsave(file.path(out_dir, "fig_03_convergence_grid.png"), p_conv, width = 10, height = 7)
ggsave(file.path(out_dir, "fig_03_detection_grid.png"), p_det, width = 10, height = 7)

cat("Successfully generated and saved academic plots to", out_dir, "\n")