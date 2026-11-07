# The two-study (data fusion) design built around DGK's code: the fold
# split used by the convex method, and the true sharp bounds of the DGK DGPs
# (Dorn, Guo and Kallus; 2025, Sec. 7.1) under the marginal sensitivity
# model (MSM), closed form for the continuous outcome and Monte Carlo for the
# binary one. Everything here is ours; DGK's own code is in R/helpers/dgk.R.
# Derivations:
# partialID_Overleaf/Results/Causal Simulations/causal_sim_explained.tex.

MEAN_SIGMA <- 4 / 3  # E[sigma(X)] = E[1 + X4^2] for X4 ~ Unif[-1, 1].
# Binary truth by Monte Carlo over X, as in DGK's dvds_plotting.r.
TRUTH_MC_SIZE <- 500000
TRUTH_MC_SEED <- 1


#' Stops unless every Lambda is finite and at least 1.
validate_lambdas <- function(lambdas) {
  if (length(lambdas) == 0 || any(!is.finite(lambdas)) || any(lambdas < 1)) {
    stop(sprintf("lambdas must be finite and >= 1, got %s",
                 paste(lambdas, collapse = ", ")))
  }
}


#' Conditional mean of DGK's continuous outcome, m(x) (dvds_sim.R line 26).
continuous_mean <- function(X) {
  2 * sign(X[, 1]) + X[, 2] + X[, 2] * X[, 3]
}


#' DGK's linear index of the binary outcome, eta(x) (dvds_sim.R line 25).
binary_index <- function(X) {
  0.5 * X[, 1] + X[, 2] + 0.25 * X[, 2] * X[, 3]
}


#' Stops unless the binary outcome's denominator k is finite and >= 1.
validate_binary_denominator <- function(denominator) {
  if (!(length(denominator) == 1 && is.finite(denominator) &&
        denominator >= 1)) {
    stop(sprintf("binary denominator must be finite and >= 1, got %s",
                 paste(denominator, collapse = ", ")))
  }
}


#' P(Y = 1 | x) of the binary outcome, p_k(x) = 1 / (k + exp(eta(x))).
#'
#' k = 1 is DGK's outcome (prevalence 0.5); larger k makes it rarer (about
#' 0.062 at k = 15).
binary_mean <- function(X, denominator = 1) {
  validate_binary_denominator(denominator)
  1 / (denominator + exp(binary_index(X)))
}


#' Draws one study: DGK's DGP, with the continuous outcome noise scaled and
#' the binary outcome optionally made rarer.
#'
#' Calls DGK's verbatim simulate_dgk() and then rescales the noise around
#' m(X), so Y = m(X) + noise_scale (1 + X4^2) eps. With binary_denominator
#' k > 1, DGK's binary outcome is thinned: each Y = 1 is kept with
#' probability p_k(X) / p_1(X) <= 1, so P(Y = 1 | X) = p_k(X) exactly, Y
#' still does not depend on Z given X, and the new outcome never exceeds
#' DGK's. Thinning draws n uniforms after DGK's draws; with k = 1 there are
#' no extra draws, and noise_scale = 1, k = 1 return DGK's study unchanged.
#'
#' @param n Sample size.
#' @param noise_scale Multiplier c > 0 of the continuous outcome's noise.
#' @param binary_denominator k >= 1 in p_k(x) (binary_mean()).
#' @return List with X, Z, Ybinary and Ycontin, as simulate_dgk().
simulate_study <- function(n, noise_scale = 1, binary_denominator = 1) {
  if (!(is.finite(noise_scale) && noise_scale > 0)) {
    stop(sprintf("noise_scale must be positive, got %s", noise_scale))
  }
  validate_binary_denominator(binary_denominator)
  study <- simulate_dgk(n)
  if (noise_scale != 1) {
    m <- continuous_mean(study$X)
    study$Ycontin <- m + noise_scale * (study$Ycontin - m)
  }
  if (binary_denominator != 1) {
    keep_prob <- binary_mean(study$X, binary_denominator) /
      binary_mean(study$X, 1)
    keep <- runif(n) < keep_prob
    study$Ybinary <- as.integer(study$Ybinary == 1 & keep)
  }
  study
}


#' Restricts a study to the units indexed by idx.
subset_study <- function(study, idx) {
  list(X = study$X[idx, , drop = FALSE], Z = study$Z[idx],
       Ybinary = study$Ybinary[idx], Ycontin = study$Ycontin[idx])
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


#' Sharp ATE bounds, continuous outcome: [-c (4/3) K_Lambda, c (4/3) K_Lambda].
#'
#' Same formula as DGK's dvds_plotting.r ("true_upper", "true_lower") when
#' c = 1. Scaling the noise by c scales sigma(x), and hence the CVaR term and
#' the band, by c (causal_sim.tex, continuous derivation).
#'
#' @param noise_scale Multiplier c of the outcome noise (simulate_study()).
#' @return data.frame with columns lambda, lower and upper.
continuous_true_band <- function(lambdas, noise_scale = 1) {
  half_width <- noise_scale * MEAN_SIGMA * k_lambda(lambdas)
  data.frame(lambda = lambdas, lower = -half_width, upper = half_width)
}


#' Sharp ATE bounds, binary outcome, by Monte Carlo over X.
#'
#' Computed as in DGK's dvds_plotting.r, with their seed and generator, so
#' the truth is theirs. With mu(x) = P(Y = 1 | x) and l = Lambda,
#'   psi1+ = E[e mu + (1 - e) min(1 - 1/l + mu/l, l mu)],
#'   psi0- = E[(1 - e) mu + e max(1 - l + l mu, mu/l)],
#' and ATE+ = psi1+ - psi0-; ATE- likewise with min and max exchanged. The
#' caller's random-number state is restored afterwards.
#'
#' @param denominator k in mu(x) = p_k(x) (binary_mean()); 1 is DGK's.
#' @return data.frame with lambda, lower, upper and their Monte Carlo s.e.
binary_true_band <- function(lambdas, denominator = 1) {
  validate_lambdas(lambdas)
  validate_binary_denominator(denominator)
  old_kind <- RNGkind()
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (is.null(old_seed)) {
      rm(".Random.seed", envir = globalenv())
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  set.seed(TRUTH_MC_SEED, kind = "Mersenne-Twister", normal.kind = "Inversion",
           sample.kind = "Rejection")
  n <- TRUTH_MC_SIZE
  X <- matrix(runif(n * 5, -1, 1), nrow = n)
  e <- 1 / (1 + exp(X[, 1] + 0.5 * (X[, 2] > 0) + 0.5 * X[, 2] * X[, 3]))
  mu <- binary_mean(X, denominator)
  # Per-draw contributions to ATE+ and ATE-; their means are the bounds.
  ate_plus <- sapply(lambdas, function(l) {
    (e * mu + (1 - e) * pmin(1 - 1 / l + mu / l, mu * l)) -
      ((1 - e) * mu + e * pmax(1 - l + mu * l, mu / l))
  })
  ate_minus <- sapply(lambdas, function(l) {
    (e * mu + (1 - e) * pmax(1 - l + mu * l, mu / l)) -
      ((1 - e) * mu + e * pmin(1 - 1 / l + mu / l, mu * l))
  })
  data.frame(lambda = lambdas,
             lower = colMeans(ate_minus), upper = colMeans(ate_plus),
             lower_mcse = apply(ate_minus, 2, sd) / sqrt(n),
             upper_mcse = apply(ate_plus, 2, sd) / sqrt(n))
}


#' Sharp ATE bounds for both outcomes, one row per (outcome, Lambda).
#'
#' @param lambdas Sensitivity parameters; lambda_index is their position.
#' @param noise_scale Multiplier of the continuous outcome noise; the binary
#'   band does not depend on it.
#' @param binary_denominator k of the binary outcome (binary_mean()); the
#'   continuous band does not depend on it.
#' @return data.frame with outcome, lambda_index, lambda, lower, upper,
#'   lower_mcse and upper_mcse (0 for the closed-form continuous band).
true_bands <- function(lambdas, noise_scale = 1, binary_denominator = 1) {
  continuous <- continuous_true_band(lambdas, noise_scale)
  continuous$lower_mcse <- 0
  continuous$upper_mcse <- 0
  rbind(
    data.frame(outcome = "binary", lambda_index = seq_along(lambdas),
               binary_true_band(lambdas, binary_denominator)),
    data.frame(outcome = "continuous", lambda_index = seq_along(lambdas),
               continuous))
}
