# Phase 1 — Classical vs. MIS Detection

| Path | Role |
|---|---|
| `R/diagnostics_classical.R` | LOO baseline toolkit |
| `script/01_classical_vs_mis_detection.R` | Simulation script |
| `output/fig1_classical_vs_mis_detection.png` | Output figure |

## R Module: `diagnostics_classical.R`

Standardised interface for computing LOO diagnostics, designed for direct comparison against MIS output.

- `get_leverage(model)`: Hat matrix diagonal.
- `get_cooks_d(model)`: Cook's Distance for all observations.
- `get_dfbetas(model, target_var)`: DFBETAS for a specific covariate.
- `get_all_classical(model, target_var)`: All three metrics combined in a single `data.frame`.
- `get_classical_set(model, target_var, k, metric)`: Top-$k$ flagged set, or threshold-based flagging using standard cutoffs ($4/n$, $2p/n$, $2/\sqrt{n}$) when `k = NULL`.

```r
source("R/diagnostics_classical.R")
fit <- lm(mpg ~ wt + hp, data = mtcars)

metrics_df <- get_all_classical(fit, target_var = "wt")
top3       <- get_classical_set(fit, target_var = "wt", k = 3, metric = "cooks_d")
flagged    <- get_classical_set(fit, target_var = "wt", k = NULL, metric = "dfbetas_target")
```

## Script: `01_classical_vs_mis_detection.R`

Constructs a synthetic DGP with a true positive slope ($Y = 1.5X + \varepsilon$), then injects two pathological structures:

- **Masked pairs** — tightly clustered outliers whose marginal LOO influence is near zero, blinding Cook's D, DFBETAS, and Leverage.
- **Sign-flippers** — a high-leverage cluster that inverts the estimated slope from positive to negative when classical tools fail to flag it.

Runs all four detection methods (top $k = 10$) and writes a 2×2 figure comparing flagged points, MIS-only flags, tool agreement, and post-removal fitted lines.

---

# Phase 2 — Exact EVT Inference & Algorithmic Comparison

Rough Structure Plan.
testingmis/
│
├── R/                              # THE "FAT ENGINE" (All reusable math & logic)
│   ├── dgp_factory.R               # Houses `generate_complex_data()` for all 9 distributions
│   ├── influence_injector.R        # Houses `apply_influence_shift()` (+/- mu, / phi logic)
│   ├── sim_engine.R                # Wrappers: `run_distribution_sim()`, `run_scaling_sim()`
│   └── utils_checkpoint.R          # Helpers to safely save and resume local .rds files
│
├── script/                         # THE "SLIM CONTROLLERS" (Executes the pipelines)
│   ├── 02_run_sim_dist.R           # Defines the DGP grid and loops through distributions
│   ├── 03_run_scaling_grid.R       # Defines the n, blocks, delta grid and runs the loop
│   ├── 81_table_results.R          # Reads .rds files -> Generates LaTeX/CSV tables
│   └── 82_plot_results.R           # Reads .rds files -> Generates ggplot2 PDFs
│
└── output/                         # THE ARTIFACTS (Never manually edit anything here)
    ├── temp/                       # Crucial for laptops! Stores chunked .rds files during loops
    ├── figures/                    # Final plots (e.g., fig_evd_robustness.pdf)
    ├── tables/                     # Final tables (e.g., tab_scaling_matrix.tex)
    ├── 02_final_distributions.rds  # The fully merged dataset from script 02
    └── 03_final_scaling.rds        # The fully merged dataset from script 03

THE 5 DEGREES OF FREEDOM / SHAPE VARIANTS
├── 1.5 = Infinite Variance (Extreme breakdown)
├── 2.1 = Barely finite variance
├── 3.0 = Standard heavy tail
├── 5.0 = Moderate tail
├── 10.0 = Almost normal

Using 14 threads (7 cores) of CPU.