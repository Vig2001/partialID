# Closed-form sharp bounds of the DGK continuous DGP (Dorn, Guo and Kallus;
# 2025, Sec. 7.1) under the marginal sensitivity model (MSM), and the sample
# splits of the two-study design. Data generation itself is DGK's code
# (R/helpers/dgk.R). Derivation of the bounds:
# partialID_Overleaf/Results/Causal Simulations/causal_sim_explained.tex.

MEAN_SIGMA <- 4 / 3  # E[sigma(X)] = E[1 + X4^2] for X4 ~ Unif[-1, 1].


#' Stops unless every Lambda is finite and at least 1.
validate_lambdas <- function(lambdas) {
  if (length(lambdas) == 0 || any(!is.finite(lambdas)) || any(lambdas < 1)) {
    stop(sprintf("lambdas must be finite and >= 1, got %s",
                 paste(lambdas, collapse = ", ")))
  }
}


#' Restricts a study to the units indexed by idx.
subset_study <- function(study, idx) {
  list(X = study$X[idx, , drop = FALSE], Z = study$Z[idx], Y = study$Y[idx])
}


#' Random partition of 1, ..., n into fold A and fold B.
#'
#' Fold A (floor(frac_a n) units) selects the convex weights; fold B gives
#' the reported CI.
#'
#' @return List with index vectors a and b.
split_folds <- function(n, frac_a) {
  if (!(frac_a > 0 && frac_a < 1)) {
    stop(sprintf("frac_a must lie in (0, 1), got frac_a = %s", frac_a))
  }
  idx_a <- sort(sample.int(n, floor(frac_a * n)))
  list(a = idx_a, b = setdiff(seq_len(n), idx_a))
}


#' K_Lambda = (Lambda^2 - 1) / Lambda * phi(Phi^{-1}(Lambda / (Lambda + 1))).
k_lambda <- function(lambdas) {
  validate_lambdas(lambdas)
  (lambdas^2 - 1) / lambdas * dnorm(qnorm(lambdas / (lambdas + 1)))
}


#' Sharp ATE bounds of one study, [-(4/3) K_Lambda, (4/3) K_Lambda].
#'
#' Same formula as DGK's dvds_plotting.r ("true_upper", "true_lower").
#'
#' @return data.frame with columns lambda, lower and upper.
true_band <- function(lambdas) {
  half_width <- MEAN_SIGMA * k_lambda(lambdas)
  data.frame(lambda = lambdas, lower = -half_width, upper = half_width)
}
