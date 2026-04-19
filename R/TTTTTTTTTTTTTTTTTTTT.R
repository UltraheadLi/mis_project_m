# ==============================================================================
# File: demo_one_call.R
# Purpose: Strip away furrr, tryCatch, and everything else. Call evt_iter_dm
#          ONCE in the most naked possible way and show what it does.
#
# Run line-by-line from the R console (not via source()) so you can see
# every message.
# ==============================================================================

library(testingMIS)
library(influence)

source("dgp_factory.R")

# ---- 1. Does evt_iter_dm even exist? Where does it live? --------------------
cat("\n== Function existence check ==\n")
cat("exists('evt_iter_dm'): ", exists("evt_iter_dm"), "\n")
if (exists("evt_iter_dm")) {
  cat("class(evt_iter_dm) : ", class(evt_iter_dm), "\n")
  cat("environment        : "); print(environment(evt_iter_dm))
  cat("formal arguments   :\n"); print(formals(evt_iter_dm))
  # where is it attached from?
  cat("find('evt_iter_dm'):\n"); print(find("evt_iter_dm"))
}

# ---- 2. Build one clean normal dataset --------------------------------------
set.seed(20260415)
dat <- generate_complex_data(n = 1000, p = 1,
                             x_type = "normal", error_type = "normal")
set_idx <- sample(seq_along(dat$y), 3)
Z_intercept <- matrix(1, nrow = 1000, ncol = 1)

# ---- 3. Call it NAKEDLY. No tryCatch. -------------------------------------
# If it errors, R will print the error and stop here - that is informative.
# If it warns, the warning will print.
# If it succeeds, we print the return object in full.
cat("\n== Calling evt_iter_dm naked ==\n")
options(warn = 1)  # print warnings immediately, do not buffer

res <- evt_iter_dm(
  y           = dat$y,
  x           = dat$X[, 1],
  Z           = Z_intercept,
  set         = set_idx,
  block_count = 25
)

cat("\n== Return object introspection ==\n")
cat("class(res):\n");   print(class(res))
cat("\nnames(res):\n"); print(names(res))
cat("\nstr(res):\n");   str(res)
cat("\nPrint res:\n");  print(res)

# ---- 4. If res has a 'converged' field, flag it explicitly ------------------
if ("converged" %in% names(res)) {
  cat("\n>>> res$converged =", res$converged, "  (class:", class(res$converged), ")\n")
}

# ---- 5. Try varying block_count - maybe 25 is too few for GEV fit ----------
cat("\n== Sweep block_count ==\n")
for (bc in c(10, 25, 50, 100)) {
  r <- tryCatch(
    evt_iter_dm(y = dat$y, x = dat$X[,1], Z = Z_intercept,
                set = set_idx, block_count = bc),
    error = function(e) list(ERROR = conditionMessage(e))
  )
  cat(sprintf("  block_count=%3d  ->  ", bc))
  if (!is.null(r$ERROR)) {
    cat("ERROR:", r$ERROR, "\n")
  } else if ("converged" %in% names(r)) {
    cat("converged =", r$converged,
        "  shape =", if(is.null(r$shape)) NA else round(r$shape, 3), "\n")
  } else {
    cat("returned names: ", paste(names(r), collapse=", "), "\n")
  }
}