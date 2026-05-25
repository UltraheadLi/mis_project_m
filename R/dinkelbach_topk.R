# ==============================================================================
# File: /R/dinkelbach_topk.R
# Purpose: General-purpose exact MIS detection via Dinkelbach's method for
#          linear-fractional programming. Finds the k observations whose joint
#          removal maximally shifts a target coefficient in a given direction.
#
# Theory:
#   The set-DFBETA for removing set S from OLS(y ~ X) on coefficient j is:
#
#       DFBETA_j(S) = e_j' (X'X)^{-1} X_S' (I - H_SS)^{-1} r_S
#
#   For the univariate FWL-orthogonalised case (regress out all but x_j),
#   this simplifies to a linear-fractional program:
#
#       max_S  sum_{i in S} x_i r_i  /  (sum_{all} x_i^2 - sum_{i in S} x_i^2)
#
#   Dinkelbach's method solves this exactly by iterating:
#     1. Form weights w_i = n_i - lambda * d_i
#     2. Select the top-k by w_i
#     3. Update lambda = sum(n[top_k]) / (C + sum(d[top_k]))
#   until convergence. This is guaranteed to converge to the global optimum
#   in O(k * n) per iteration, typically 3-8 iterations.
#
# Interface:
#   Two entry points are provided:
#     1. dinkelbach_topk()    — low-level, works on raw vectors (X, R)
#     2. dinkelbach_topk_lm() — high-level, works on an lm object
#                               (drop-in replacement for fast_sens_topk)
#
# Performance:
#   - O(n log n) per Dinkelbach iteration (dominated by partial sort)
#   - Typically 3-8 iterations to converge (guaranteed monotone)
#   - Comparable wall-clock to fast_sens_topk for n < 5000
#   - Exact: no approximation, no refinement loop needed
#
# Called by: /R/sim_robust_engine.R, /script/04_compare_robust.R
# ==============================================================================


#' Exact MIS Detection via Dinkelbach's Method (Low-Level)
#'
#' Finds the k indices from candidate vectors (x, r) that maximise the
#' linear-fractional objective:
#'
#'   sgn * sum(x[S] * r[S]) / (sum_x2_all - sum(x[S]^2))
#'
#' where sum_x2_all is the total sum of squares of the full predictor
#' (including observations not in the candidate set).
#'
#' This is the core Dinkelbach solver. It operates on pre-orthogonalised
#' (FWL) vectors and assumes a simple regression structure. For the
#' multivariate case, the caller must perform FWL projection first.
#'
#' @param x       Numeric vector; FWL-orthogonalised predictor values for
#'                the candidate observations (length n_candidates).
#' @param r       Numeric vector; OLS residuals for the candidate
#'                observations (same length as x).
#' @param k       Integer; size of the influential set to find.
#' @param sgn     Integer; +1 or -1, the direction of influence to maximise.
#' @param sum_x2  Numeric; total sum of x^2 over ALL observations (not just
#'                candidates). This is the denominator anchor. If NULL
#'                (default), computed as sum(x^2) — correct when candidates
#'                are the full sample.
#' @param max_iter Integer; maximum Dinkelbach iterations (default = 50).
#'                Convergence is guaranteed and typically occurs in 3-8.
#' @param tol     Numeric; convergence tolerance for lambda (default = 1e-9).
#'
#' @return A list with components:
#'   \item{indices}{Integer vector of length k — positions in x/r of the
#'                  most influential observations.}
#'   \item{dfbeta}{Numeric; the exact DFBETA value of the selected set
#'                 (signed, in the direction of sgn).}
#'   \item{lambda}{Numeric; the converged Dinkelbach parameter.}
#'   \item{iterations}{Integer; number of Dinkelbach iterations used.}
#' @export
dinkelbach_topk <- function(x, r, k, sgn = 1L,
                            sum_x2 = NULL,
                            max_iter = 50L, tol = 1e-9) {

  n <- length(x)

  # --- Input validation ---
  if (length(r) != n) {
    stop("x and r must have the same length.")
  }
  if (k < 1L || k > n) {
    stop(sprintf("k must be between 1 and %d (got %d).", n, k))
  }
  if (!sgn %in% c(1L, -1L)) {
    stop("sgn must be +1 or -1.")
  }

  # --- Precompute numerator and denominator components ---
  # Objective: max_S  sgn * sum(x[S]*r[S]) / (sum_x2_all - sum(x[S]^2))
  #
  # Dinkelbach form:  max_S  sum(n_val[S]) / (C + sum(d_val[S]))
  #   where  n_val_i = sgn * x_i * r_i
  #          d_val_i = -x_i^2
  #          C       = sum_x2_all

  n_val <- sgn * (x * r)
  d_val <- -(x^2)

  if (is.null(sum_x2)) {
    sum_x2 <- sum(x^2)
  }

  # --- Dinkelbach iteration ---
  lambda <- 0
  idx <- integer(k)
  n_iter <- 0L

  for (iter in seq_len(max_iter)) {
    n_iter <- iter

    # Step 1: Form parametric weights
    w <- n_val - lambda * d_val

    # Step 2: Select top-k by weight (partial sort for efficiency)
    # order() is O(n log n); for very large n, a partial sort could be
    # used, but order() is well-optimised in R and sufficient for n < 50k
    idx <- order(w, decreasing = TRUE)[seq_len(k)]

    # Step 3: Update lambda (the Dinkelbach ratio)
    num <- sum(n_val[idx])
    den <- sum_x2 + sum(d_val[idx])

    # Guard against degenerate denominator (all x are in S)
    if (abs(den) < 1e-15) {
      warning("Dinkelbach denominator near zero — degenerate design.")
      break
    }

    new_lambda <- num / den

    # Step 4: Convergence check
    if (abs(new_lambda - lambda) < tol) {
      lambda <- new_lambda
      break
    }
    lambda <- new_lambda
  }

  # --- Compute the signed DFBETA ---
  dfbeta_val <- sgn * lambda

  return(list(
    indices    = idx,
    dfbeta     = dfbeta_val,
    lambda     = lambda,
    iterations = n_iter
  ))
}


#' Exact MIS Detection via Dinkelbach's Method (lm Interface)
#'
#' Drop-in replacement for \code{fast_sens_topk}. Accepts a fitted lm object
#' and returns the k row indices of the most influential observations for
#' a target coefficient, using the exact Dinkelbach solver.
#'
#' Internally performs FWL (Frisch-Waugh-Lovell) projection to reduce the
#' multivariate problem to a univariate linear-fractional program, then
#' calls \code{dinkelbach_topk} on the projected vectors.
#'
#' @param mod   A fitted lm object.
#' @param pos   Integer; column position of the target coefficient in the
#'              design matrix (e.g., 2 for the first slope when an intercept
#'              is present). Default = 2.
#' @param sign  Integer; +1 or -1 direction for influence maximisation.
#'              +1 finds the set whose removal most increases beta[pos];
#'              -1 finds the set whose removal most decreases it.
#' @param k     Integer; number of most influential observations to return.
#'
#' @return Integer vector of length k — the original row indices of the k
#'         most influential observations. Same return type as
#'         \code{fast_sens_topk} for drop-in compatibility.
#' @export
dinkelbach_topk_lm <- function(mod, pos = 2L, sign = 1L, k = 1L) {

  X <- stats::model.matrix(mod)
  y <- stats::model.response(stats::model.frame(mod))
  N <- nrow(X)
  p <- ncol(X)

  if (pos < 1L || pos > p) {
    stop(sprintf("pos must be between 1 and %d (got %d).", p, pos))
  }
  if (k < 1L || k > N) {
    stop(sprintf("k must be between 1 and %d (got %d).", N, k))
  }

  # ------------------------------------------------------------------
  # FWL Projection: partial out all columns except pos
  #
  # Let X = [Z | x_j] where x_j is the target column.
  # FWL gives:
  #   x_fwl = M_Z x_j       (residuals of x_j on Z)
  #   y_fwl = M_Z y          (residuals of y on Z)
  #   r_fwl = y_fwl - x_fwl * beta_j   (= full OLS residuals)
  #
  # The DFBETA for beta_j from removing set S depends only on
  # (x_fwl, r_fwl), reducing the problem to the univariate case.
  # ------------------------------------------------------------------

  if (p == 1L) {
    # No nuisance regressors — trivial case
    x_fwl <- X[, 1]
    r_fwl <- stats::residuals(mod)
  } else {
    Z_cols <- setdiff(seq_len(p), pos)
    Z <- X[, Z_cols, drop = FALSE]

    # Project out Z from both x_j and y
    qr_Z <- qr(Z)
    x_fwl <- qr.resid(qr_Z, X[, pos])
    y_fwl <- qr.resid(qr_Z, y)

    # Residuals of the FWL regression (= full model residuals)
    # beta_j_fwl = sum(x_fwl * y_fwl) / sum(x_fwl^2)
    r_fwl <- y_fwl - x_fwl * (sum(x_fwl * y_fwl) / sum(x_fwl^2))
  }

  # Total sum of squares of x_fwl (full sample — the denominator anchor)
  sum_x2_full <- sum(x_fwl^2)

  # --- Call the core Dinkelbach solver ---
  result <- dinkelbach_topk(
    x      = x_fwl,
    r      = r_fwl,
    k      = k,
    sgn    = as.integer(sign),
    sum_x2 = sum_x2_full
  )

  return(result$indices)
}