# ==============================================================================
# Script: /script/01_demo_masking.R
# Purpose: Demonstrate classical LOO failure vs MIS under masking.
# Outputs: ../output/fig1_masking_demo.pdf&png
# ==============================================================================

# 1. Setup & Dependencies
library(influence)
source("../R/diagnostics_classical.R") 
dir.create("../output", showWarnings = FALSE)

set.seed(123)

# 2. DGP
n_base <- 50
x_base <- rnorm(n_base)
y_base <- 1.5 * x_base + rnorm(n_base, sd = 0.5)

# Inject outlier pairs:
# A/B/C — masking demo (High, Medium, Low intensity)
# D     — sign-flip group
x_out <- c(
  5.0, 5.1,          # Pair A (Extreme)
  3.0, 3.1,          # Pair B (Moderate)
  1.5, 1.6,          # Pair C (Mild)
  4.5, 4.7, 4.9, 5.2 # Group D (Sign-flip)
)
y_out <- c(
  -5.0, -5.1,        # Pair A
  -2.0, -2.1,        # Pair B
  -0.5, -0.6,        # Pair C
  -9.0, -9.4, -9.8, -10.4  # Group D — y deeply negative at high x
)

df <- data.frame(
  id = 1:60,
  x  = c(x_base, x_out),
  y  = c(y_base, y_out)
)

mdl_base <- lm(y ~ x, data = df)

# 3. Detection (Targeting k=10)
ids_cooks <- get_classical_set(mdl_base, target_var = "x",k=10,  metric = "cooks_d")
ids_lev   <- get_classical_set(mdl_base, target_var = "x",k=10,  metric = "leverage")
ids_dfb   <- get_classical_set(mdl_base, target_var = "x",k=10,  metric = "dfbetas_target")

mis_run <- sens(mdl_base, lambda = set_lambda("beta_i", pos = 2, sign = -1), options = set_options("all"))
ids_mis <- mis_run$influence$id[1:10]

# 4. Generate Output 
# Build agreement grid data for Plot 3
all_ids    <- sort(unique(unlist(flagged)))
tool_names <- names(flagged)

agreement_mat <- matrix(
  FALSE,
  nrow = length(all_ids),
  ncol = length(tool_names),
  dimnames = list(as.character(all_ids), tool_names)
)
for (tool in tool_names) {
  agreement_mat[as.character(flagged[[tool]]), tool] <- TRUE
}

# Draw 2×2 panel set (for both pdf and png)
draw_panels <- function() {
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  
  # --- Plot 1: Classical flags ---
  plot(df$x, df$y,
       main = "Classical Tools: Flagged Points",
       xlab = "X", ylab = "Y",
       pch = 16, col = "gray80", cex = 0.9)
  classical_tools <- c("cooks", "dfb", "lev")
  for (tool in classical_tools) {
    ids <- flagged[[tool]]
    points(df$x[ids], df$y[ids],
           col = tool_cols[tool],
           pch = tool_pchs[tool],
           cex = 2, lwd = 2)
  }
  legend("topleft",
         legend = tool_labels[classical_tools],
         col    = tool_cols[classical_tools],
         pch    = tool_pchs[classical_tools],
         pt.cex = 1.6, pt.lwd = 2, cex = 0.8, bg = "white")
  
  # --- Plot 2: MIS flags ---
  plot(df$x, df$y,
       main = "MIS: Flagged Points",
       xlab = "X", ylab = "Y",
       pch = 16, col = "gray80", cex = 0.9)
  points(df$x[ids_mis], df$y[ids_mis],
         col = tool_cols["mis"],
         pch = tool_pchs["mis"],
         cex = 2, lwd = 2)
  legend("topleft",
         legend = tool_labels["mis"],
         col    = tool_cols["mis"],
         pch    = tool_pchs["mis"],
         pt.cex = 1.6, pt.lwd = 2, cex = 0.8, bg = "white")
  
  # --- Plot 3: Agreement grid ---
  n_obs  <- nrow(agreement_mat)
  n_tool <- ncol(agreement_mat)
  col_order <- c("cooks", "dfb", "lev", "mis")
  
  # Bottom margin for angled labels
  par(mar = c(6, 4, 3, 1))
  
  plot.new()
  plot.window(xlim = c(0.5, n_tool + 0.5), ylim = c(0.5, n_obs + 0.5))
  title(main = "Flag Agreement: Observations × Tool", cex.main = 1)
  
  # Y-axis: observation IDs
  axis(2, at = seq_len(n_obs), labels = all_ids, las = 1, cex.axis = 0.65, tick = FALSE)
  
  # X-axis: tool labels at 45 degrees using text() instead of axis()
  text(x      = seq_along(col_order),
       y      = 0.5 - (n_obs * 0.06),
       labels = tool_labels[col_order],
       srt    = 45,
       adj    = c(1, 1),
       xpd    = TRUE,
       cex    = 0.80)
  
  # Grid background
  abline(v = seq(0.5, n_tool + 0.5, 1), col = "gray90", lwd = 0.8)
  abline(h = seq(0.5, n_obs + 0.5, 1), col = "gray90", lwd = 0.8)
  
  # Cells
  for (j in seq_along(col_order)) {
    tool <- col_order[j]
    for (i in seq_len(n_obs)) {
      obs_id <- all_ids[i]
      if (agreement_mat[as.character(obs_id), tool]) {
        points(j, i, pch = 15, cex = 2.2, col = tool_cols[tool])
      } else {
        points(j, i, pch = 15, cex = 2.2, col = "gray93")
      }
    }
  }
  
  # Legend
  legend("bottomright",
         legend = c("Flagged", "Not flagged"),
         pch    = c(15, 15),
         col    = c("gray40", "gray93"),
         pt.cex = 1.4, cex = 0.72, bg = "white")

  
  # --- Plot 4: Fitted lines after removing each flagged set ---
  # Margin for Plot 4
  par(mar = c(4, 4, 3, 1))
  
  plot(df$x, df$y,
       main = "Fitted Lines After Removing Flagged Sets",
       xlab = "X", ylab = "Y",
       pch = 16, col = "gray80", cex = 0.9)
  abline(mdl_base, col = "black", lty = 2, lwd = 1.5)
  for (tool in names(flagged)) {
    ids     <- flagged[[tool]]
    mdl_new <- lm(y ~ x, data = df[-ids, ])
    abline(mdl_new, col = tool_cols[tool], lty = 1, lwd = 2)
  }
  legend("topleft",
         legend = c("Original Fit", tool_labels),
         col    = c("black", tool_cols),
         lty    = c(2, rep(1, 4)),
         lwd    = c(1.5, rep(2, 4)),
         cex    = 0.8, bg = "white")
}

# PDF
pdf("../output/fig1_masking_demo.pdf", width = 14, height = 12)
draw_panels()
dev.off()

# PNG
png("../output/fig1_masking_demo.png", width = 14, height = 12,
    units = "in", res = 300)
draw_panels()
dev.off()

cat("Simulation complete.")