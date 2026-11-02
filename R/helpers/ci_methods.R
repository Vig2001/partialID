# Confidence intervals for a partially identified band fused from two
# independent studies. Port of src/helpers/ci_methods.py and
# optimal_weight_grid in src/helpers/weight_optimisers.py, with separate
# variances at the lower and upper ends.
#
# Inputs are R x 2 matrices: row r is replicate r, column k is study k.
# Each end uses critical value z; z = qnorm(1 - alpha / 2) gives a
# Bonferroni joint CI for the band. Each method returns a list with
#   lo, hi        CI endpoints;
#   L_hat, U_hat  fused bound estimates;
#   w_L, w_U      weight on study 1 at the lower (omega_0) and upper (omega_1)
#                 end.

W_GRID <- seq(0, 1, length.out = 1001)  # Candidate weights, spacing 1e-3.


#' Stops unless the bound matrices share shape R x 2, variances are
#' non-negative and z > 0.
validate_bounds <- function(L, U, s2_L, s2_U, z) {
  if (!is.matrix(L) || ncol(L) != 2) {
    stop(sprintf("L must be an R x 2 matrix, got dim %s",
                 paste(dim(L), collapse = " x ")))
  }
  for (x in list(U, s2_L, s2_U)) {
    if (!identical(dim(x), dim(L))) {
      stop("L, U, s2_L and s2_U must share a shape")
    }
  }
  if (any(s2_L < 0) || any(s2_U < 0)) stop("variances must be non-negative")
  if (z <= 0) stop(sprintf("z must be positive, got z = %s", z))
}


#' Stops unless w is a scalar or a length-R vector in [0, 1].
validate_weight <- function(w, n_reps, name) {
  if (!(length(w) %in% c(1, n_reps))) {
    stop(sprintf("%s must have length 1 or %d, got %d", name, n_reps,
                 length(w)))
  }
  if (any(w < 0 | w > 1)) {
    stop(sprintf("%s must lie in [0, 1], got range [%s, %s]", name,
                 min(w), max(w)))
  }
}


#' Naive intersection CI: largest lower bound, smallest upper bound.
#'
#' Each end uses the standard error of the study attaining it, so the CI
#' ignores the selection step and can undercover when the two bounds at an
#' end are close.
intersection_ci <- function(L, U, s2_L, s2_U, z) {
  validate_bounds(L, U, s2_L, s2_U, z)
  rows <- seq_len(nrow(L))
  k_L <- max.col(L, ties.method = "first")
  k_U <- max.col(-U, ties.method = "first")
  L_hat <- L[cbind(rows, k_L)]
  U_hat <- U[cbind(rows, k_U)]
  list(
    lo = L_hat - z * sqrt(s2_L[cbind(rows, k_L)]),
    hi = U_hat + z * sqrt(s2_U[cbind(rows, k_U)]),
    L_hat = L_hat,
    U_hat = U_hat,
    w_L = as.numeric(k_L == 1),
    w_U = as.numeric(k_U == 1)
  )
}


#' Convex-combination CI with weight w on study 1 at each end.
#'
#' The studies are independent, so w L_1 + (1 - w) L_2 has variance
#' w^2 s2_L1 + (1 - w)^2 s2_L2; likewise at the upper end.
#'
#' @param w_L,w_U Weights on study 1; scalars or length-R vectors in [0, 1].
convex_ci <- function(L, U, s2_L, s2_U, w_L, w_U, z) {
  validate_bounds(L, U, s2_L, s2_U, z)
  validate_weight(w_L, nrow(L), "w_L")
  validate_weight(w_U, nrow(L), "w_U")
  L_hat <- w_L * L[, 1] + (1 - w_L) * L[, 2]
  U_hat <- w_U * U[, 1] + (1 - w_U) * U[, 2]
  se_L <- sqrt(w_L^2 * s2_L[, 1] + (1 - w_L)^2 * s2_L[, 2])
  se_U <- sqrt(w_U^2 * s2_U[, 1] + (1 - w_U)^2 * s2_U[, 2])
  list(
    lo = L_hat - z * se_L,
    hi = U_hat + z * se_U,
    L_hat = L_hat,
    U_hat = U_hat,
    w_L = rep_len(w_L, nrow(L)),
    w_U = rep_len(w_U, nrow(L))
  )
}


#' Weight on bound 1 optimised over w_grid, vectorised over replicates.
#'
#' lower = TRUE maximises the lower CI limit
#'   w m1 + (1 - w) m2 - z sqrt(w^2 v1 + (1 - w)^2 v2);
#' lower = FALSE minimises the upper limit, with + z in place of - z.
#' The objective is concave (convex) in w, so the grid maximum (minimum) is
#' the global optimum up to the grid spacing.
#'
#' @param m1,m2 Bound estimates, length R.
#' @param v1,v2 Variances used in the penalty, length R.
#' @return Weights, length R.
optimal_weight_grid <- function(m1, m2, v1, v2, lower, z, w_grid = W_GRID) {
  mean_part <- outer(m1, w_grid) + outer(m2, 1 - w_grid)
  penalty <- z * sqrt(outer(v1, w_grid^2) + outer(v2, (1 - w_grid)^2))
  idx <- if (lower) {
    max.col(mean_part - penalty, ties.method = "first")
  } else {
    max.col(-(mean_part + penalty), ties.method = "first")
  }
  w_grid[idx]
}
