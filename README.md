# Classical Diagnostics Baseline (`/R/diagnostics_classical.R`)

## Overview
This module provides a standardized interface for computing classical, leave-one-out (LOO) regression diagnostics. It is designed to serve as the baseline comparison for evaluating the performance of exact Most Influential Sets (MIS) detection algorithms, particularly in the presence of masking and joint influence.

## Functions Included
The script is modular. You can extract individual metrics or aggregate them:
* `get_leverage(model)`: Extracts diagonal hat matrix values.
* `get_cooks_d(model)`: Computes Cook's Distance for all observations.
* `get_dfbetas(model, target_var)`: Computes DFBETAS for a specific covariate.
* `get_all_classical(model, target_var)`: Aggregates the above into a clean, dependency-free `data.frame`.
* `get_classical_set(model, target_var, k, metric)`: Extracts an "influential set" using classical heuristics, either by taking the top `k` absolute values or by applying standard statistical thresholds ($4/n$, $2p/n$, $2/\sqrt{n}$).

## Quick Usage

```R
# Source the functions
source("R/diagnostics_classical.R")

# Fit a standard linear model
fit <- lm(mpg ~ wt + hp, data = mtcars)

# 1. Get all raw diagnostics in a neat table
metrics_df <- get_all_classical(fit, target_var = "wt")

# 2. Extract the IDs of the top 3 most influential points based on Cook's D
top_3_ids <- get_classical_set(fit, target_var = "wt", k = 3, metric = "cooks_d")

# 3. Extract IDs that exceed theoretical statistical thresholds for DFBETAS
flagged_ids <- get_classical_set(fit, target_var = "wt", k = NULL, metric = "dfbetas_target")