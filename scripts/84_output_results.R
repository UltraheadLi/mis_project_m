# ==============================================================================
# File: /script/04_output_results.R
# Purpose: Generate publication-quality tables and figures from the robust
#          comparison Monte Carlo simulation results.
#
# Outputs:
#   Tables (CSV + LaTeX):
#     04_table1_bias_main.csv / .tex   — Mean |Bias| by Estimator × Outlier × Error Dist
#     04_table2_coverage_and_k.csv / .tex — Coverage + k-selection accuracy
#
#   Figures (PDF):
#     04_fig1_coefficient_distributions.pdf — Coef estimates across DGPs
#     04_fig2_bias_coverage_tradeoff.pdf    — Bias vs Coverage scatter
#     04_fig3_contamination_illustration.pdf — Clean vs contaminated X
#
# Input: ../output/04_robust_comparison.rds  (compiled MC results)
#
# Usage: source("script/04_output_results.R")
# ==============================================================================

# ==============================================================================
# 0. Setup
# ==============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(xtable)

# --- Paths ---
input_path   <- "../output/04_robust_comparison_results.rds"
fig_dir      <- "../output/figures"
tab_dir      <- "../output/tables"

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---
sim <- readRDS(input_path)
cat(sprintf("Loaded %d simulation rows.\n", nrow(sim)))

# --- Source DGP + injector for Plot 3 ---
source("../R/dgp_factory.R")
source("../R/influence_injector.R")

# ==============================================================================
# Shared constants
# ==============================================================================

# Nice estimator labels (order matters for plotting)
estimator_levels <- c("full", "cd", "lev", "dfb",
                      "mis_alpha", "mis_peel", "mis_oracle",
                      "mm", "lts")

estimator_labels <- c(
  "full"       = "OLS (Full)",
  "cd"         = "Cook's D",
  "lev"        = "Leverage",
  "dfb"        = "DFBETAS",
  "mis_alpha"  = "MIS (alpha-k)",
  "mis_peel"   = "MIS (Peel)",
  "mis_oracle" = "MIS (Oracle)",
  "mm"         = "MM-Estimator",
  "lts"        = "LTS"
)

# Outlier type labels
outlier_labels <- c(
  "none"              = "No Contamination",
  "vertical_outlier"  = "Vertical Outlier",
  "good_leverage"     = "Good Leverage",
  "bad_leverage"      = "Bad Leverage"
)

# Error distribution ordering (increasing tail weight)
error_order <- c("normal", "mixed_normal", "beta_logistic", "skewed_t",
                 "contaminated", "golm", "pareto", "gpd")

error_labels <- c(
  "normal"         = "Normal",
  "mixed_normal"   = "Mixed Normal",
  "beta_logistic"  = "Beta-Logistic",
  "skewed_t"       = "Skewed-t",
  "contaminated"   = "Contaminated",
  "golm"           = "GOLM",
  "pareto"         = "Pareto",
  "gpd"            = "GPD"
)

# Colour palette for 9 estimators
estimator_colours <- c(
  "OLS (Full)"      = "#999999",
  "Cook's D"        = "#E69F00",
  "Leverage"        = "#F0E442",
  "DFBETAS"         = "#D55E00",
  "MIS (alpha-k)"   = "#56B4E9",
  "MIS (Peel)"      = "#0072B2",
  "MIS (Oracle)"    = "#009E73",
  "MM-Estimator"    = "#CC79A7",
  "LTS"             = "#882255"
)

# ==============================================================================
# Helper: Reshape simulation data to long format
# ==============================================================================
reshape_bias_long <- function(df) {
  df %>%
    select(iter, x_type, error_type, outlier_method, set_size,
           starts_with("bias_")) %>%
    pivot_longer(cols = starts_with("bias_"),
                 names_to = "estimator",
                 names_prefix = "bias_",
                 values_to = "abs_bias") %>%
    mutate(
      estimator = factor(estimator, levels = estimator_levels,
                         labels = estimator_labels[estimator_levels]),
      outlier_method = factor(outlier_method,
                              levels = names(outlier_labels),
                              labels = outlier_labels),
      error_type = factor(error_type,
                          levels = error_order,
                          labels = error_labels[error_order])
    )
}

reshape_coef_long <- function(df) {
  df %>%
    select(iter, x_type, error_type, outlier_method, set_size,
           starts_with("coef_")) %>%
    pivot_longer(cols = starts_with("coef_"),
                 names_to = "estimator",
                 names_prefix = "coef_",
                 values_to = "coefficient") %>%
    mutate(
      estimator = factor(estimator, levels = estimator_levels,
                         labels = estimator_labels[estimator_levels]),
      outlier_method = factor(outlier_method,
                              levels = names(outlier_labels),
                              labels = outlier_labels),
      error_type = factor(error_type,
                          levels = error_order,
                          labels = error_labels[error_order])
    )
}

reshape_coverage_long <- function(df) {
  df %>%
    select(iter, x_type, error_type, outlier_method, set_size,
           starts_with("cov_")) %>%
    pivot_longer(cols = starts_with("cov_"),
                 names_to = "estimator",
                 names_prefix = "cov_",
                 values_to = "coverage") %>%
    mutate(
      estimator = factor(estimator, levels = estimator_levels,
                         labels = estimator_labels[estimator_levels]),
      outlier_method = factor(outlier_method,
                              levels = names(outlier_labels),
                              labels = outlier_labels),
      error_type = factor(error_type,
                          levels = error_order,
                          labels = error_labels[error_order])
    )
}

# ==============================================================================
# Shared ggplot theme for publication
# ==============================================================================
theme_thesis <- function(base_size = 10) {
  theme_minimal(base_size = base_size) +
    theme(
      text              = element_text(family = ""),
      plot.title         = element_text(face = "bold", size = base_size + 2,
                                        hjust = 0.5),
      plot.subtitle      = element_text(size = base_size, hjust = 0.5,
                                        color = "grey40"),
      strip.text         = element_text(face = "bold", size = base_size),
      strip.background   = element_rect(fill = "grey95", colour = NA),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x        = element_text(angle = 35, hjust = 1, size = base_size - 1),
      axis.title         = element_text(size = base_size),
      legend.position    = "bottom",
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 1),
      legend.key.size    = unit(0.4, "cm"),
      plot.margin        = margin(10, 10, 10, 10)
    )
}

# ==============================================================================
# TABLE 1: Mean Absolute Bias — Estimator × Outlier Type × Error Distribution
# ==============================================================================
cat("Generating Table 1: Bias...\n")

bias_long <- reshape_bias_long(sim)

tab1 <- bias_long %>%
  group_by(estimator, outlier_method, error_type) %>%
  summarise(mean_bias = mean(abs_bias, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = c(outlier_method, error_type),
              values_from = mean_bias,
              names_sep = " | ")

# --- Save CSV ---
write.csv(tab1, file.path(tab_dir, "04_table1_bias_main.csv"),
          row.names = FALSE)

# --- Save LaTeX ---
# For thesis: a condensed version aggregated over error distributions
tab1_condensed <- bias_long %>%
  group_by(estimator, outlier_method) %>%
  summarise(
    mean_bias = sprintf("%.4f", mean(abs_bias, na.rm = TRUE)),
    sd_bias   = sprintf("%.4f", sd(abs_bias, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(cell = paste0(mean_bias, " (", sd_bias, ")")) %>%
  select(estimator, outlier_method, cell) %>%
  pivot_wider(names_from = outlier_method, values_from = cell)

# Also generate the full breakdown (estimator × outlier × error) in wide form
tab1_full <- bias_long %>%
  group_by(estimator, outlier_method, error_type) %>%
  summarise(mean_bias = sprintf("%.4f", mean(abs_bias, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_wider(names_from = error_type, values_from = mean_bias)

# Save both
write.csv(tab1_condensed, file.path(tab_dir, "04_table1_bias_condensed.csv"),
          row.names = FALSE)

latex1 <- xtable(tab1_condensed,
                 caption = "Mean Absolute Bias (SD) by Estimator and Outlier Type",
                 label = "tab:bias_main",
                 align = c("l", "l", rep("r", ncol(tab1_condensed) - 1)))
print(latex1,
      file = file.path(tab_dir, "04_table1_bias_main.tex"),
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      caption.placement = "top")

cat("  Table 1 saved.\n")

# ==============================================================================
# TABLE 2: Coverage Rates + k-Selection Summary
# ==============================================================================
cat("Generating Table 2: Coverage & k-selection...\n")

# --- Left half: Coverage rates ---
cov_long <- reshape_coverage_long(sim)

cov_summary <- cov_long %>%
  group_by(estimator, outlier_method) %>%
  summarise(coverage = sprintf("%.3f", mean(coverage, na.rm = TRUE)),
            .groups = "drop") %>%
  pivot_wider(names_from = outlier_method, values_from = coverage)

# --- Right half: k-selection accuracy (MIS variants only) ---
k_summary <- sim %>%
  filter(outlier_method != "none") %>%
  group_by(outlier_method) %>%
  summarise(
    true_k      = sprintf("%.1f", mean(set_size, na.rm = TRUE)),
    k_alpha_mean = sprintf("%.1f", mean(k_alpha, na.rm = TRUE)),
    k_peel_mean  = sprintf("%.1f", mean(k_peel, na.rm = TRUE)),
    k_oracle_mean = sprintf("%.1f", mean(k_oracle, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    outlier_method = factor(outlier_method,
                            levels = names(outlier_labels),
                            labels = outlier_labels)
  )

# Save coverage
write.csv(cov_summary, file.path(tab_dir, "04_table2_coverage.csv"),
          row.names = FALSE)

latex_cov <- xtable(cov_summary,
                    caption = "95\\% CI Coverage Rates by Estimator and Outlier Type",
                    label = "tab:coverage",
                    align = c("l", "l", rep("r", ncol(cov_summary) - 1)))
print(latex_cov,
      file = file.path(tab_dir, "04_table2_coverage.tex"),
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      caption.placement = "top")

# Save k-selection
write.csv(k_summary, file.path(tab_dir, "04_table2_k_selection.csv"),
          row.names = FALSE)

latex_k <- xtable(k_summary,
                  caption = "Mean Estimated $k$ by MIS Variant and Outlier Type",
                  label = "tab:k_selection",
                  align = c("l", "l", rep("r", ncol(k_summary) - 1)))
print(latex_k,
      file = file.path(tab_dir, "04_table2_k_selection.tex"),
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      caption.placement = "top")

cat("  Table 2 saved.\n")

# ==============================================================================
# FIGURE 1: Coefficient Estimates Across Error Distributions
#
# X-axis:  Error distribution (ordered by tail weight)
# Y-axis:  Estimated coefficient (true beta = 1)
# Colour:  Estimator
# Facets:  Outlier type
# Geom:    Median line + Q5/Q95 ribbon with visible boundary lines
#
# Key fix: use 5th/95th percentiles (resistant to extreme MC draws)
#          + coord_cartesian to clip y-axis to informative region
# ==============================================================================
cat("Generating Figure 1: Coefficient distributions...\n")

coef_long <- reshape_coef_long(sim)

# --- Trimmed summaries: median + 5th/95th percentiles ---
coef_ribbon <- coef_long %>%
  group_by(estimator, outlier_method, error_type) %>%
  summarise(
    median_coef = median(coefficient, na.rm = TRUE),
    q05 = quantile(coefficient, 0.05, na.rm = TRUE),
    q95 = quantile(coefficient, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

# --- Determine sensible y-axis limits per facet ---
# Use the 1st/99th percentile of ALL medians + ribbon edges within each
# outlier type, with a small pad, so the plot zooms to the action region.
y_limits <- coef_ribbon %>%
  group_by(outlier_method) %>%
  summarise(
    ylo = quantile(q05, 0.01, na.rm = TRUE),
    yhi = quantile(q95, 0.99, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    pad = (yhi - ylo) * 0.10,
    ylo = ylo - pad,
    yhi = yhi + pad
  )

# Build per-facet blank data for scale limits (ggplot trick)
dummy_limits <- y_limits %>%
  select(outlier_method, ylo, yhi) %>%
  pivot_longer(cols = c(ylo, yhi), values_to = "coefficient") %>%
  mutate(error_type = levels(coef_ribbon$error_type)[1]) %>%
  select(outlier_method, error_type, coefficient)

fig1 <- ggplot(coef_ribbon,
               aes(x = error_type, group = estimator)) +
  # Shaded ribbon (90% interval)
  geom_ribbon(aes(ymin = q05, ymax = q95, fill = estimator),
              alpha = 0.10) +
  # Boundary lines for the ribbon
  geom_line(aes(y = q05, colour = estimator),
            linewidth = 0.3, linetype = "dotted") +
  geom_line(aes(y = q95, colour = estimator),
            linewidth = 0.3, linetype = "dotted") +
  # Median line + points
  geom_line(aes(y = median_coef, colour = estimator),
            linewidth = 0.7) +
  geom_point(aes(y = median_coef, colour = estimator),
             size = 1.5) +
  # True beta reference
  geom_hline(yintercept = 1, linetype = "dashed", colour = "black",
             linewidth = 0.5) +
  # Invisible points to set per-facet y limits
  geom_blank(data = dummy_limits,
             aes(x = error_type, y = coefficient),
             inherit.aes = FALSE) +
  facet_wrap(~ outlier_method, ncol = 2, scales = "free_y") +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  scale_fill_manual(values = estimator_colours, name = "Estimator") +
  labs(
    title = "Coefficient Estimates Across Error Distributions",
    subtitle = expression("Median (solid) with 90% interval (dotted); dashed line = true " * beta * " = 1"),
    x = "Error Distribution",
    y = expression(hat(beta))
  ) +
  theme_thesis(base_size = 10) +
  guides(colour = guide_legend(nrow = 2,
                               override.aes = list(linewidth = 1.2,
                                                   linetype = "solid")),
         fill   = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "04_fig1_coefficient_distributions.pdf"),
       plot = fig1, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "04_fig1_coefficient_distributions.png"),
       plot = fig1, width = 10, height = 8, dpi = 300)

cat("  Figure 1 saved.\n")

# ==============================================================================
# FIGURE 2: Bias–Coverage Trade-off Scatter
#
# X-axis: Median absolute bias (log10 scale — handles huge range)
# Y-axis: 95% CI coverage rate
# Points: One per estimator, labelled
# Facets: Outlier type
#
# Key fix: use MEDIAN bias (not mean) to resist extreme MC draws,
#          log10 x-axis so all estimators are spatially separated,
#          cleaner label placement.
# ==============================================================================
cat("Generating Figure 2: Bias-coverage tradeoff...\n")

cov_long <- reshape_coverage_long(sim)

bias_summary <- bias_long %>%
  group_by(estimator, outlier_method) %>%
  summarise(median_bias = median(abs_bias, na.rm = TRUE), .groups = "drop")

cov_summary_num <- cov_long %>%
  group_by(estimator, outlier_method) %>%
  summarise(mean_cov = mean(coverage, na.rm = TRUE), .groups = "drop")

tradeoff <- inner_join(bias_summary, cov_summary_num,
                       by = c("estimator", "outlier_method"))

# Assign short labels for plot annotation (avoids clutter)
est_short <- c(
  "OLS (Full)"      = "OLS",
  "Cook's D"        = "Cook",
  "Leverage"        = "Lev",
  "DFBETAS"         = "DFB",
  "MIS (alpha-k)"   = "MIS-\u03b1",
  "MIS (Peel)"      = "MIS-P",
  "MIS (Oracle)"    = "MIS-O",
  "MM-Estimator"    = "MM",
  "LTS"             = "LTS"
)
tradeoff$short_label <- est_short[as.character(tradeoff$estimator)]

# Estimator grouping for visual distinction
est_group <- c(
  "OLS (Full)"      = "Baseline",
  "Cook's D"        = "Classical",
  "Leverage"        = "Classical",
  "DFBETAS"         = "Classical",
  "MIS (alpha-k)"   = "MIS",
  "MIS (Peel)"      = "MIS",
  "MIS (Oracle)"    = "MIS",
  "MM-Estimator"    = "Robust",
  "LTS"             = "Robust"
)
tradeoff$group <- est_group[as.character(tradeoff$estimator)]

# Shape palette by group
group_shapes <- c(
  "Baseline"  = 16,
  "Classical" = 17,
  "MIS"       = 15,
  "Robust"    = 18
)
tradeoff$group <- factor(tradeoff$group,
                         levels = c("Baseline", "Classical", "MIS", "Robust"))

fig2 <- ggplot(tradeoff,
               aes(x = median_bias, y = mean_cov,
                   colour = estimator, shape = group)) +
  # Ideal zone shading
  annotate("rect",
           xmin = -Inf, xmax = Inf, ymin = 0.94, ymax = 0.96,
           fill = "grey90", alpha = 0.5) +
  geom_hline(yintercept = 0.95, linetype = "dashed", colour = "grey40",
             linewidth = 0.4) +
  geom_point(size = 3.5, stroke = 0.7) +
  ggrepel::geom_text_repel(
    aes(label = short_label),
    size = 2.8, fontface = "bold",
    max.overlaps = 20, seed = 42,
    segment.colour = "grey60", segment.size = 0.3,
    min.segment.length = 0.2,
    box.padding = 0.4, point.padding = 0.3,
    show.legend = FALSE
  ) +
  facet_wrap(~ outlier_method, ncol = 2, scales = "free_x") +
  scale_x_log10(labels = scales::label_number(accuracy = 0.001)) +
  scale_colour_manual(values = estimator_colours, name = "Estimator") +
  scale_shape_manual(values = group_shapes, name = "Method Class") +
  labs(
    title = "Bias\u2013Coverage Trade-off by Outlier Type",
    subtitle = "Ideal region: low bias (left) + coverage near 0.95 (grey band); x-axis on log scale",
    x = "Median Absolute Bias (log scale)",
    y = "95% CI Coverage Rate"
  ) +
  theme_thesis(base_size = 10) +
  theme(panel.grid.major.x = element_line(colour = "grey92")) +
  guides(colour = guide_legend(nrow = 2, order = 1),
         shape  = guide_legend(nrow = 1, order = 2))

ggsave(file.path(fig_dir, "04_fig2_bias_coverage_tradeoff.pdf"),
       plot = fig2, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "04_fig2_bias_coverage_tradeoff.png"),
       plot = fig2, width = 10, height = 8, dpi = 300)
cat("  Figure 2 saved.\n")

# ==============================================================================
# FIGURE 3: Contamination Illustration (Clean vs Injected X)
#
# Shows one simulated dataset per outlier type to visualise what each
# contamination topology looks like in (x, y) space.
# ==============================================================================
cat("Generating Figure 3: Contamination illustration...\n")

set.seed(2024)

# Generate one clean dataset
n_demo <- 200
dat_clean <- generate_complex_data(n = n_demo, p = 1,
                                   x_type = "normal",
                                   error_type = "normal")

# Create panels for each outlier type
outlier_methods <- c("vertical_outlier", "good_leverage", "bad_leverage")
k_demo <- 10
mag_demo <- 5

panels <- list()

# Panel 1: Clean
df_clean <- data.frame(x = dat_clean$X[, 1], y = dat_clean$y,
                       status = "Clean", panel = "No Contamination")
panels[[1]] <- df_clean

# Panels 2-4: Each injection type
for (om in outlier_methods) {
  set.seed(2024)  # Reset so clean data is identical
  dat_c <- generate_complex_data(n = n_demo, p = 1,
                                 x_type = "normal",
                                 error_type = "normal")
  dat_inj <- apply_influence_shift(dat_c, method = om,
                                   k = k_demo, magnitude = mag_demo)
  
  status <- rep("Clean", n_demo)
  status[dat_inj$outlier_indices] <- "Contaminated"
  
  df_panel <- data.frame(x = dat_inj$X[, 1], y = dat_inj$y,
                         status = status,
                         panel = outlier_labels[om])
  panels[[length(panels) + 1]] <- df_panel
}

df_illust <- bind_rows(panels)
df_illust$panel <- factor(df_illust$panel,
                          levels = c("No Contamination",
                                     outlier_labels[outlier_methods]))

fig3 <- ggplot(df_illust, aes(x = x, y = y)) +
  # Clean OLS fit (on clean obs only)
  # Clean OLS fit (on clean obs only)
  geom_smooth(data = df_illust %>% filter(status == "Clean"),
              aes(linetype = "Clean OLS"),
              method = "lm", se = FALSE, colour = "grey40", linewidth = 1) +
  # Contaminated OLS fit (on ALL obs) — shows how the line shifts
  geom_smooth(data = df_illust %>% filter(panel != "No Contamination"),
              aes(linetype = "Contaminated OLS"),
              method = "lm", se = FALSE, colour = "#D55E00", linewidth = 1) +
  scale_linetype_manual(name = "Models", 
                        values = c("Clean OLS" = "solid", "Contaminated OLS" = "dashed")) +
  # Points
  geom_point(aes(colour = status, size = status, alpha = status)) +
  facet_wrap(~ panel, ncol = 2, scales = "free") +
  scale_colour_manual(values = c("Clean" = "#0072B2",
                                 "Contaminated" = "#D55E00"),
                      name = "Observation") +
  scale_size_manual(values = c("Clean" = 1, "Contaminated" = 2.5),
                    guide = "none") +
  scale_alpha_manual(values = c("Clean" = 0.35, "Contaminated" = 0.9),
                     guide = "none") +
  labs(
    title = "Contamination Topologies: Clean vs. Injected Observations",
    x = "X (Predictor)",
    y = "Y (Response)"
  ) +
  theme_thesis(base_size = 10) +
  guides(
    colour = guide_legend(override.aes = list(size = 4, alpha = 1)),
    linetype = guide_legend(override.aes = list(linewidth = 2))
  ) +
  theme(
    legend.position = "bottom",
    legend.key.width = unit(1.5, "cm"),
    plot.title = element_text(size = 18, face = "bold"),    # Makes title larger
    legend.title = element_text(size = 14, face = "bold"),  # Makes legend titles larger
    legend.text = element_text(size = 14),                 # Makes legend text larger
    strip.text = element_text(size = 14, face = "bold")
  )

ggsave(file.path(fig_dir, "04_fig3_contamination_illustration.pdf"),
       plot = fig3, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "04_fig3_contamination_illustration.png"),
       plot = fig3, width = 10, height = 8, dpi = 300)

cat("  Figure 3 saved.\n")

# ==============================================================================
# Summary
# ==============================================================================
cat("\n========================================\n")
cat("All outputs generated successfully.\n")
cat("========================================\n")
cat(sprintf("Tables:  %s\n", tab_dir))
cat(sprintf("Figures: %s\n", fig_dir))
cat("Files:\n")
cat("  Tables:\n")
cat("    04_table1_bias_main.csv / .tex\n")
cat("    04_table1_bias_condensed.csv\n")
cat("    04_table2_coverage.csv / .tex\n")
cat("    04_table2_k_selection.csv / .tex\n")
cat("  Figures:\n")
cat("    04_fig1_coefficient_distributions.pdf\n")
cat("    04_fig2_bias_coverage_tradeoff.pdf\n")
cat("    04_fig3_contamination_illustration.pdf\n")
cat("========================================\n")