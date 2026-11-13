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
#   w_L, w_U      weight on study 1 at the lower (omega_L) and upper (omega_U)
#                 end.

W_GRID <- seq(0, 1, length.out = 1001)  # Candidate weights, spacing 1e-3.
# Weights of a single study: 0 (study 2) or 1 (study 1). Optimising the
# convex objective over these two alone selects, at each end, the study with
# the tighter CI limit, i.e. intersects the two CIs (CORNER_WEIGHTS in
# src/toy_sim.py).
CORNER_WEIGHTS <- c(0, 1)


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


# ------------- calibrated CIs and CLR with selection (2026-11-13) -----------
# Ports of calibrated_ci_intersection, clr_intersection_selection_ci and
# clr_convex_selection_ci in src/helpers/ci_methods.py (same algorithms, so
# the two agree to rounding error), with separate variances at each end.
# The selection step is CLR (2013) eqs. (3.4)-(3.5), p. 687, with gamma_n =
# 1 - 0.1 / log n (l_n := log n, p. 686); the CI end is eq. (3.7), p. 688.
# Critical values come from Hotelling's tube formula for d = 1, exact here
# (arc_exceedance()), instead of CLR's simulation (Definition 3, p. 690).

SELECTION_GAMMA_CONSTANT <- 0.1  # gamma_n = 1 - 0.1 / log n (CLR)
SIMPSON_PANELS <- 64  # composite Simpson panels for the dropped-arc integral
BISECTION_STEPS <- 60  # bisection steps for the kept-set critical value


#' Intersection of the two studies' CIs, each with critical value z.
#'
#' lo = max_k (L_k - z se_k), hi = min_k (U_k + z se_k). With the Sidak value
#' this equals clr_intersection_ci() at p = 1 - alpha / 2; it is written out
#' separately so that the tables show the equality in the results.
#'
#' @return CI list; L_hat and U_hat are the estimates of the studies attaining
#'   the CI ends; w_L and w_U are 1 where study 1 attains them.
calibrated_ci_intersection <- function(L, U, s2_L, s2_U, z) {
  validate_bounds(L, U, s2_L, s2_U, z)
  rows <- seq_len(nrow(L))
  lower_limits <- L - z * sqrt(s2_L)
  upper_limits <- U + z * sqrt(s2_U)
  k_L <- max.col(lower_limits, ties.method = "first")
  k_U <- max.col(-upper_limits, ties.method = "first")
  list(
    lo = lower_limits[cbind(rows, k_L)],
    hi = upper_limits[cbind(rows, k_U)],
    L_hat = L[cbind(rows, k_L)],
    U_hat = U[cbind(rows, k_U)],
    w_L = as.numeric(k_L == 1),
    w_U = as.numeric(k_U == 1)
  )
}


#' CLR's selection level gamma_n = 1 - 0.1 / log(n); n = n_1 + n_2 here.
clr_selection_level <- function(n) {
  if (n < 2) stop(sprintf("n must be >= 2, got %s", n))
  1 - SELECTION_GAMMA_CONSTANT / log(n)
}


#' CLR's adaptive inequality selection at an upper end: bound j is kept if
#' theta_j <= min_i (theta_i + k_gamma s_i) + 2 k_gamma s_j.
#'
#' @return Logical R x J matrix.
kept_by_selection <- function(theta, s, k_gamma) {
  floor <- apply(theta + k_gamma * s, 1, min)
  theta <= floor + 2 * k_gamma * s
}


#' CLR upper end for min(theta_1, theta_2) with selection: with m kept
#' bounds the critical value is Phi^{-1}(p^{1/m}).
#'
#' @return List with ci_end, estimate (level 1/2) and index.
clr_upper_selected_intersection <- function(theta, s, k_gamma, p) {
  m <- rowSums(kept_by_selection(theta, s, k_gamma))
  end <- smallest_corrected(theta, s, qnorm(p^(1 / m)))
  list(ci_end = end$value,
       estimate = smallest_corrected(theta, s, qnorm(0.5^(1 / m)))$value,
       index = end$index)
}


#' Stops unless gamma lies in (1/2, 1).
validate_gamma <- function(gamma) {
  if (!(gamma > 0.5 && gamma < 1)) {
    stop(sprintf("gamma must lie in (0.5, 1), got gamma = %s", gamma))
  }
}


#' CLR intersection CI with adaptive inequality selection, two studies.
#'
#' Selection uses k_gamma = Phi^{-1}(sqrt(gamma)); the critical value is
#' Phi^{-1}(p^{1/m}) for m kept bounds (2.239 or 1.960 at p = 0.975); the CI
#' end is the smallest corrected bound over both studies.
#'
#' @param gamma Selection level, e.g. clr_selection_level(n_1 + n_2).
#' @return CI list; L_hat and U_hat are the half-median-unbiased estimates.
clr_intersection_selection_ci <- function(L, U, s2_L, s2_U, p, gamma) {
  validate_bounds(L, U, s2_L, s2_U, 1)
  validate_level(p, 0.5)
  validate_gamma(gamma)
  k_gamma <- max_of_two_quantile(gamma)
  upper <- clr_upper_selected_intersection(U, sqrt(s2_U), k_gamma, p)
  lower <- clr_upper_selected_intersection(-L, sqrt(s2_L), k_gamma, p)
  list(
    lo = -lower$ci_end,
    hi = upper$ci_end,
    L_hat = -lower$estimate,
    U_hat = upper$estimate,
    w_L = as.numeric(lower$index == 1),
    w_U = as.numeric(upper$index == 1)
  )
}


#' (1/pi) int_0^{delta/2} exp(-k^2 / (2 cos^2 t)) dt per dropped arc delta,
#' by composite Simpson; an R x G matrix, zero for a zero angle. k has
#' length R, dropped_arcs is R x G.
dropped_arc_integral <- function(k, dropped_arcs) {
  nodes <- seq(0, 1, length.out = 2 * SIMPSON_PANELS + 1)
  weights <- rep(1, length(nodes))
  weights[seq(2, length(nodes) - 1, by = 2)] <- 4
  weights[seq(3, length(nodes) - 2, by = 2)] <- 2
  total <- matrix(0, nrow(dropped_arcs), ncol(dropped_arcs))
  for (i in seq_along(nodes)) {
    t <- dropped_arcs / 2 * nodes[i]
    total <- total + weights[i] * exp(-k^2 / (2 * cos(t)^2))
  }
  step <- dropped_arcs / 2 / (2 * SIMPSON_PANELS)
  step / 3 * total / pi
}


#' P(sup over a union of arcs S of u'Z > k), Z ~ N(0, I_2), k >= 0, S
#' inside the positive quadrant: Q(k) + |S| / (2 pi) exp(-k^2 / 2) plus the
#' dropped-arc integrals (see _arc_exceedance in the Python module).
#'
#' Hotelling's tube formula for a one-dimensional index set (d = 1), exact
#' here since Z is two-dimensional and S is an arc no longer than pi / 2. It
#' is the "exact version of the tube formula" that CLR's simulated critical
#' values approximate (CLR 2013, p. 698, Remark 2: "A sharper analytical
#' approach can be based on Hotelling's tube method ... That approach is
#' tractable for the case of d = 1, but does not immediately extend to
#' d > 1. Note that the simulation based approach is effectively a numeric
#' version of the exact version of the tube formula and is less conservative
#' than using simplified tube formulas."). Derivation in
#' partialID_Overleaf/Notes/Simulation_Notes/toy_sim_explained.tex.
arc_exceedance <- function(k, arc_length, dropped_arcs) {
  pnorm(k, lower.tail = FALSE) + arc_length / (2 * pi) * exp(-k^2 / 2) +
    rowSums(dropped_arc_integral(k, dropped_arcs))
}


#' The level-quantile of that sup per row, by bisection between the
#' quantiles for one direction and for the whole quadrant.
arc_quantile <- function(level, arc_length, dropped_arcs) {
  lo <- rep(qnorm(level), length(arc_length))
  hi <- rep(quadrant_sup_quantile(level), length(arc_length))
  for (step in seq_len(BISECTION_STEPS)) {
    mid <- (lo + hi) / 2
    above <- arc_exceedance(mid, arc_length, dropped_arcs) > 1 - level
    lo <- ifelse(above, mid, lo)
    hi <- ifelse(above, hi, mid)
  }
  (lo + hi) / 2
}


#' Total angle and dropped_arcs of the standardised-error directions of the kept
#' weights; weight w has angle atan2((1 - w) se_2, w se_1).
#'
#' @param kept Logical R x W matrix over w_grid.
#' @param se R x 2 standard errors of the two bounds.
#' @return List with arc_length (length R) and dropped_arcs (R x G,
#'   zero-padded).
kept_arcs <- function(kept, se, w_grid) {
  n_w <- length(w_grid)
  phi <- atan2(outer(se[, 2], 1 - w_grid), outer(se[, 1], w_grid))
  both_kept <- kept[, -n_w, drop = FALSE] & kept[, -1, drop = FALSE]
  arc_length <- rowSums((phi[, -n_w, drop = FALSE] - phi[, -1, drop = FALSE]) *
                          both_kept)
  # A dropped arc ends at kept column j when the previous kept column is
  # < j - 1.
  previous <- matrix(0L, nrow(kept), n_w)
  last <- integer(nrow(kept))
  for (j in seq_len(n_w)) {
    previous[, j] <- last
    last <- ifelse(kept[, j], j, last)
  }
  ends_dropped <- kept & previous > 0 & previous < col(kept) - 1
  where <- which(ends_dropped, arr.ind = TRUE)
  if (nrow(where) == 0) {
    return(list(arc_length = arc_length,
                dropped_arcs = matrix(0, nrow(kept), 1)))
  }
  sizes <- phi[cbind(where[, 1], previous[where])] - phi[where]
  rank <- ave(seq_len(nrow(where)), where[, 1], FUN = seq_along)
  dropped_arcs <- matrix(0, nrow(kept), max(rank))
  dropped_arcs[cbind(where[, 1], rank)] <- sizes
  list(arc_length = arc_length, dropped_arcs = dropped_arcs)
}


#' CLR upper end over all convex weights with selection.
#'
#' @return List with ci_end, estimate (level 1/2) and index.
clr_upper_selected_convex <- function(theta, s, se, w_grid, k_gamma, p) {
  arcs <- kept_arcs(kept_by_selection(theta, s, k_gamma), se, w_grid)
  k_p <- arc_quantile(p, arcs$arc_length, arcs$dropped_arcs)
  k_half <- arc_quantile(0.5, arcs$arc_length, arcs$dropped_arcs)
  end <- smallest_corrected(theta, s, k_p)
  list(ci_end = end$value,
       estimate = smallest_corrected(theta, s, k_half)$value,
       index = end$index)
}


#' Convex CLR CI with adaptive inequality selection over the weights.
#'
#' Selection uses k_gamma = quadrant_sup_quantile(gamma); the critical value
#' is the p-quantile of the sup over the kept weights (arc_quantile()); the
#' CI end is the smallest corrected bound over all weights.
#'
#' @param gamma Selection level, e.g. clr_selection_level(n_1 + n_2).
#' @return CI list; L_hat and U_hat are the half-median-unbiased estimates.
clr_convex_selection_ci <- function(L, U, s2_L, s2_U, p, gamma,
                                    w_grid = W_GRID) {
  validate_bounds(L, U, s2_L, s2_U, 1)
  if (any(s2_L <= 0) || any(s2_U <= 0)) stop("variances must be positive")
  validate_level(p, 0.5)
  validate_gamma(gamma)
  k_gamma <- quadrant_sup_quantile(gamma)
  fam_U <- convex_family(U, s2_U, w_grid)
  upper <- clr_upper_selected_convex(fam_U$theta, fam_U$s, sqrt(s2_U), w_grid,
                                     k_gamma, p)
  fam_L <- convex_family(-L, s2_L, w_grid)
  lower <- clr_upper_selected_convex(fam_L$theta, fam_L$s, sqrt(s2_L), w_grid,
                                     k_gamma, p)
  list(
    lo = -lower$ci_end,
    hi = upper$ci_end,
    L_hat = -lower$estimate,
    U_hat = upper$estimate,
    w_L = w_grid[lower$index],
    w_U = w_grid[upper$index]
  )
}
