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


#' Cross-fitted convex CI: fold-wise convex bounds, averaged over K folds.
#'
#' Port of cross_fit_convex_ci in src/helpers/ci_methods.py. Fold f
#' contributes w_f bound_1f + (1 - w_f) bound_2f, with weights chosen on the
#' other folds. The folds are disjoint, so the average of the K fold bounds
#' has variance (1 / K^2) sum_f [w_f^2 s2_1f + (1 - w_f)^2 s2_2f] at each
#' end. Each weight is treated as fixed given its training folds; the
#' dependence between a fold's estimate and the weights the other folds draw
#' from it is ignored.
#'
#' @param L,U,s2_L,s2_U Lists of K per-fold R x 2 matrices.
#' @param w_L,w_U Lists of K per-fold weight vectors on study 1, length R.
#' @return CI list; w_L and w_U are the weights averaged over folds.
cross_fit_convex_ci <- function(L, U, s2_L, s2_U, w_L, w_U, z) {
  n_folds <- length(L)
  lengths <- vapply(list(U, s2_L, s2_U, w_L, w_U), length, integer(1))
  if (n_folds < 2 || any(lengths != n_folds)) {
    stop(sprintf("need >= 2 folds and one entry per fold, got %s",
                 paste(c(n_folds, lengths), collapse = ", ")))
  }
  folds <- lapply(seq_len(n_folds), function(f) {
    convex_ci(L[[f]], U[[f]], s2_L[[f]], s2_U[[f]], w_L[[f]], w_U[[f]], z)
  })
  fold_mean <- function(name) {
    Reduce(`+`, lapply(folds, `[[`, name)) / n_folds
  }
  L_hat <- fold_mean("L_hat")
  U_hat <- fold_mean("U_hat")
  # Each fold's convex_ci half-width is z times its standard error.
  var_L <- Reduce(`+`, lapply(folds, function(f) ((f$L_hat - f$lo) / z)^2))
  var_U <- Reduce(`+`, lapply(folds, function(f) ((f$hi - f$U_hat) / z)^2))
  list(
    lo = L_hat - z * sqrt(var_L) / n_folds,
    hi = U_hat + z * sqrt(var_U) / n_folds,
    L_hat = L_hat,
    U_hat = U_hat,
    w_L = fold_mean("w_L"),
    w_U = fold_mean("w_U")
  )
}


#' Row-wise min_j (theta_j + k s_j) and the column attaining it.
#'
#' @param theta,s R x J matrices of estimates and standard errors.
#' @return List with value (length R) and index (length R).
smallest_corrected <- function(theta, s, k) {
  corrected <- theta + k * s
  index <- max.col(-corrected, ties.method = "first")
  list(value = corrected[cbind(seq_len(nrow(theta)), index)], index = index)
}


#' CLR upper end for min_j theta_j with every bound kept (no selection).
#'
#' @return List with ci_end (critical value k_p), estimate (k_half, the
#'   half-median-unbiased estimate) and index (column attaining ci_end).
clr_upper <- function(theta, s, k_p, k_half) {
  end <- smallest_corrected(theta, s, k_p)
  list(ci_end = end$value,
       estimate = smallest_corrected(theta, s, k_half)$value,
       index = end$index)
}


#' The level-quantile of max(Z_1, Z_2), Z_1 and Z_2 independent N(0, 1).
#'
#' P(max(Z_1, Z_2) <= k) = Phi(k)^2, so k = Phi^{-1}(sqrt(level)).
max_of_two_quantile <- function(level) qnorm(sqrt(level))


#' Stops unless p lies in (lower, 1).
validate_level <- function(p, lower = 0) {
  if (!(p > lower && p < 1)) {
    stop(sprintf("p must lie in (%s, 1), got p = %s", lower, p))
  }
}


#' Intersection CI of Chernozhukov, Lee and Rosen (2013) for two studies.
#'
#' Port of clr_intersection_ci in src/helpers/ci_methods.py, with separate
#' variances at each end. At the upper end each U_k is corrected upward,
#' U_k + k_p se_k, and the CI end is the smaller corrected estimate. There
#' is no adaptive inequality selection step: both studies are kept, so k_p
#' is the p-quantile of the larger of two independent N(0, 1),
#' Phi^{-1}(sqrt(p)) (2.239 at p = 0.975). The lower end is the mirror
#' image.
#'
#' @param p Coverage level of each one-sided end, e.g. 1 - alpha / 2.
#' @return CI list; L_hat and U_hat are the half-median-unbiased estimates
#'   (level 1/2); w_L and w_U are 1 where study 1 attains that CI end.
clr_intersection_ci <- function(L, U, s2_L, s2_U, p) {
  validate_bounds(L, U, s2_L, s2_U, 1)
  validate_level(p)
  k_p <- max_of_two_quantile(p)
  k_half <- max_of_two_quantile(0.5)
  upper <- clr_upper(U, sqrt(s2_U), k_p, k_half)
  # max_k L_k = -min_k (-L_k): the lower end is the upper end of -L.
  lower <- clr_upper(-L, sqrt(s2_L), k_p, k_half)
  list(
    lo = -lower$ci_end,
    hi = upper$ci_end,
    L_hat = -lower$estimate,
    U_hat = upper$estimate,
    w_L = as.numeric(lower$index == 1),
    w_U = as.numeric(upper$index == 1)
  )
}


#' P(sup_w a(w)'Z <= k) for k >= 0, with Z ~ N(0, I_2).
#'
#' a(w) is the standardised-error direction of weight w; over w in [0, 1] it
#' sweeps the positive quadrant, so the sup is ||Z|| when both Z_k > 0 and
#' max(Z_1, Z_2) otherwise (see _quadrant_sup_cdf in the Python module).
quadrant_sup_cdf <- function(k) pnorm(k) - exp(-k^2 / 2) / 4

QUADRANT_MIN_LEVEL <- 0.25  # cdf of the sup at k = 0
QUADRANT_K_MAX <- 8  # upper bracket for its quantile


#' The level-quantile of the sup in quadrant_sup_cdf(), for level > 1/4:
#' 2.353 at 0.975 and 0.563 at 1/2.
quadrant_sup_quantile <- function(level) {
  validate_level(level, QUADRANT_MIN_LEVEL)
  uniroot(function(k) quadrant_sup_cdf(k) - level,
          c(0, QUADRANT_K_MAX), tol = 1e-12)$root
}


#' w theta_1 + (1 - w) theta_2 and its standard error over the weight grid.
#'
#' @return List with theta and s, each R x length(w_grid).
convex_family <- function(theta, s2, w_grid) {
  list(theta = outer(theta[, 1], w_grid) + outer(theta[, 2], 1 - w_grid),
       s = sqrt(outer(s2[, 1], w_grid^2) + outer(s2[, 2], (1 - w_grid)^2)))
}


#' Convex CI without a sample split, calibrated by CLR over all weights.
#'
#' Port of clr_convex_ci in src/helpers/ci_methods.py, with separate
#' variances at each end. Every convex combination w U_1 + (1 - w) U_2
#' bounds min(U_1, U_2) from above, so the upper end is an intersection bound
#' over w in [0, 1]. Each is corrected by k_p s(w), and the CI end is the
#' smallest corrected bound. With no adaptive inequality selection step,
#' k_p is the p-quantile of the sup of the standardised errors over all
#' weights (quadrant_sup_quantile(), 2.353 at p = 0.975). The lower end is
#' the mirror image.
#'
#' @param p Coverage level of each one-sided end, in (1/4, 1).
#' @param w_grid Equally spaced weights from 0 to 1.
#' @return CI list; L_hat and U_hat are the half-median-unbiased estimates
#'   (level 1/2); w_L and w_U are the weights attaining the CI ends.
clr_convex_ci <- function(L, U, s2_L, s2_U, p, w_grid = W_GRID) {
  validate_bounds(L, U, s2_L, s2_U, 1)
  if (any(s2_L <= 0) || any(s2_U <= 0)) stop("variances must be positive")
  steps <- diff(w_grid)
  if (w_grid[1] != 0 || w_grid[length(w_grid)] != 1 || any(steps <= 0) ||
      max(abs(steps - steps[1])) > 1e-12) {
    stop("w_grid must be equally spaced from 0 to 1")
  }
  k_p <- quadrant_sup_quantile(p)
  k_half <- quadrant_sup_quantile(0.5)
  fam_U <- convex_family(U, s2_U, w_grid)
  upper <- clr_upper(fam_U$theta, fam_U$s, k_p, k_half)
  # max_w L(w) = -min_w (-L(w)), with the same weight on study 1.
  fam_L <- convex_family(-L, s2_L, w_grid)
  lower <- clr_upper(fam_L$theta, fam_L$s, k_p, k_half)
  list(
    lo = -lower$ci_end,
    hi = upper$ci_end,
    L_hat = -lower$estimate,
    U_hat = upper$estimate,
    w_L = w_grid[lower$index],
    w_U = w_grid[upper$index]
  )
}
