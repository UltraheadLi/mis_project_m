# ==============================================================================
# File: /R/helpers_local.R
# Purpose: Local implementations of testingMIS internal functions.
#          This file MUST be sourced FIRST, before any other R/ files.
#          Here are: testingMIS:::fwl, testingMIS:::make_blocks,
#                    testingMIS::dfbeta_numeric, testingMIS::dfb_bmx,
#                    testingMIS::estimate_dfb_evd
# ==============================================================================

# ---------- Core helpers (were unexported in testingMIS) ----------

#' Frisch-Waugh-Lovell projection: residualizes y and X on Z
#' @noRd
fwl <- function(y, X, Z) {
  Q <- qr(as.matrix(Z))
  cbind(
    Yp = y - qr.fitted(Q, y),
    Xp = as.matrix(X) - qr.fitted(Q, as.matrix(X))
  )
}

#' Splits a vector into equal-sized blocks (drops remainder)
#' @noRd
make_blocks <- function(X, block_size) {
  n <- length(X) %/% block_size
  matrix(X[seq_len(n * block_size)], nrow = block_size, ncol = n)
}

# ---------- Exported functions from testingMIS ----------

#' Set DFBETA for numeric inputs
dfbeta_numeric <- function(y, X, set, col_X = 1) {
  X <- as.matrix(X)[, col_X]
  sum_xy_all <- sum(X * y)
  sum_xsq_all <- sum(X^2)
  c(
    dfbeta = -((sum_xy_all - sum(X[set] * y[set])) /
      (sum_xsq_all - sum(X[set]^2)) -
      sum_xy_all / sum_xsq_all)
  )
}

#' Set Influence Block Maxima (Greedy)
dfb_bmx <- function(X, R, set, block_count) {
  sgn <- sign(sum(X[set] * R[set]))
  if (sgn == 0) {
    stop("dfbeta of set is exactly zero")
  }
  which.extr <- if (sgn > 0) which.max else which.min

  X <- X[-set]
  R <- R[-set]

  sumX2 <- sum(X^2)
  nS <- length(set)
  block_size <- length(X) %/% block_count
  Xbl <- make_blocks(X, block_size)
  Rbl <- make_blocks(R, block_size)

  dfb_mat <- Xbl * Rbl / (sumX2 - Xbl^2)
  m_list <- list(apply(dfb_mat, 2, which.extr))

  if (nS > 1) {
    for (s in 2:nS) {
      dfb_mat <- vapply(
        seq_len(block_count),
        function(i) {
          selected_idx <- vapply(m_list, function(m) m[i], integer(1))
          remaining <- setdiff(seq_len(block_size), selected_idx)
          result <- rep(-sgn * Inf, block_size)
          result[remaining] <- Xbl[remaining, i] *
            Rbl[remaining, i] /
            (sumX2 - sum(Xbl[selected_idx, i]^2) - Xbl[remaining, i]^2)
          result
        },
        numeric(block_size)
      )
      m_list[[s]] <- apply(dfb_mat, 2, which.extr)
    }
  }

  all_selected <- do.call(rbind, m_list)

  vapply(
    seq_len(block_count),
    function(i) {
      idx <- all_selected[, i]
      sum(Xbl[idx, i] * Rbl[idx, i]) / (sumX2 - sum(Xbl[idx, i]^2))
    },
    numeric(1)
  )
}

#' Estimate EVD of Set Influence (Greedy version)
estimate_dfb_evd <- function(y, x, Z, set, block_count = 20, verbose = TRUE) {
  fwl_vars <- fwl(y = y, X = x, Z = Z)

  Y <- fwl_vars[, 1]
  X <- fwl_vars[, 2]

  fwl_lm <- lm(Y ~ X - 1)
  R <- fwl_lm |> residuals()

  set_dfb <- dfbeta_numeric(Y, X, set)

  bm_X <- apply(make_blocks(X[-set], length(X[-set]) %/% block_count), 2, max)
  bm_R <- apply(make_blocks(R[-set], length(R[-set]) %/% block_count), 2, max)

  x_evd <- evd::fgev(bm_X)
  r_evd <- evd::fgev(bm_R)

  is_x_frechet <- x_evd$estimate["shape"] - 1.96 * x_evd$std.err["shape"] > 0
  is_r_frechet <- r_evd$estimate["shape"] - 1.96 * r_evd$std.err["shape"] > 0

  tail_coef <- max(
    is_x_frechet * x_evd$estimate["shape"],
    is_r_frechet * r_evd$estimate["shape"]
  )

  if (isTRUE(verbose)) {
    cat(sprintf(
      "\nX: Shape = %.4f \t| Pr(>x) = %.4f\nR: Shape = %.4f \t| Pr(>x) = %.4f\n",
      x_evd$estimate["shape"],
      pnorm(-x_evd$estimate["shape"], 0, x_evd$std.err["shape"]),
      r_evd$estimate["shape"],
      pnorm(-r_evd$estimate["shape"], 0, r_evd$std.err["shape"])
    ))
  }

  Delta_bmx <- abs(dfb_bmx(X, R, set = set, block_count))
  fit_evd_bm <- evd::fgev(Delta_bmx, shape = tail_coef)
  fit_evd_bm$estimate["shape"] <- tail_coef

  list(
    params = fit_evd_bm$estimate,
    set_dfb = set_dfb,
    block_maxima = Delta_bmx
  )
}