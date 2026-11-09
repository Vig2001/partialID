# DVDS bounds of Dorn, Guo and Kallus (DGK) with parametric nuisances, for
# the lean causal simulation (R/causal_sim_quick.R).
#
# The nuisances are fitted by logistic and linear regression on the true
# functional forms of DGK's DGP (R/helpers/dgk.R, simulate_dgk()), in place
# of their forests:
#   e(x)          logistic regression of Z on X1, 1(X2 > 0) and X2 X3: the
#                 true model;
#   mu_z(x)       binary outcome: logistic regression of Y on X1, X2 and
#                 X2 X3 within arm z, exact for k = 1 and nearly so for
#                 k > 1 (binary_mean() in data_fusion_design.R);
#   m_z, s_z      continuous outcome: Y = m(x) + s(x) eps with eps ~ N(0, 1)
#                 within arm z; m by least squares on sign(X1), X2 and X2 X3,
#                 s by least squares of |residual| on 1 and X4^2, times
#                 sqrt(pi / 2) (E|eps| = sqrt(2 / pi)). Both are the true
#                 forms, also with study 2's noise scaled by c.
# Every unit's nuisances are predicted out of fold (cross_fit_nuisances(),
# 5 folds within each fit by default), as DGK's forests predict out of bag:
# a unit's own outcome never enters the fit that predicts for it. In-sample
# fits shrink the residuals and bias both bounds inwards by O(p / n). For
# the cross-fitted CI, a held-out fold's nuisances come from models fitted
# on the other folds (dvds_bounds_lean_train_test(), standard cross-fitting).
# The binary bounds come from DGK's dvds.binary(), sourced unmodified from
# their clone (it takes mu0, mu1 and e only). The continuous bounds use
# DGK's influence function (simplified_dvds.R, phi.continuous()), written
# out here with the normal quantile and CVaR of the fitted location-scale
# model in place of their quantile and CVaR forests. Notes:
# partialID_Overleaf/Notes/Simulation_Notes/causal_sim_quick_explained.tex.

SCALE_FLOOR <- 1e-8  # Smallest fitted s(x); a guard, not reached in practice.
NUISANCE_FOLDS <- 5  # Default folds for out-of-fold nuisance predictions.


#' Features of the true models: propensity, binary mean, continuous mean.
#'
#' @param X n x 5 covariate matrix.
#' @return data.frame with columns x1, x2, x23 (X2 X3), x4sq (X4^2), x2pos
#'   (1(X2 > 0)) and sx1 (sign(X1)).
lean_features <- function(X) {
  data.frame(x1 = X[, 1], x2 = X[, 2], x23 = X[, 2] * X[, 3],
             x4sq = X[, 4]^2, x2pos = as.numeric(X[, 2] > 0),
             sx1 = sign(X[, 1]))
}


#' Propensity score e(x) = P(Z = 1 | x) by logistic regression.
#'
#' @param features,Z Training data.
#' @param newdata Features of the units to predict for.
fit_propensity <- function(features, Z, newdata) {
  fit <- stats::glm(Z ~ x1 + x2pos + x23, family = stats::binomial(),
                    data = cbind(features, Z = Z))
  unname(stats::predict(fit, newdata = newdata, type = "response"))
}


#' Arm-wise binary outcome means mu_z(x), by logistic regression.
#'
#' @return List with mu0 and mu1, predicted for the rows of newdata.
fit_binary_means <- function(features, Z, Y, newdata) {
  data <- cbind(features, Y = Y)
  arm_mean <- function(z) {
    fit <- stats::glm(Y ~ x1 + x2 + x23, family = stats::binomial(),
                      data = data[Z == z, ])
    unname(stats::predict(fit, newdata = newdata, type = "response"))
  }
  list(mu0 = arm_mean(0), mu1 = arm_mean(1))
}


#' Arm-wise normal location-scale fit of the continuous outcome.
#'
#' m_z(x) by least squares on sign(X1), X2 and X2 X3; s_z(x) by least
#' squares of the absolute residual on 1 and X4^2, rescaled by
#' sqrt(pi / 2) since E|eps| = sqrt(2 / pi) for eps ~ N(0, 1).
#'
#' @return List with m0, m1, s0 and s1, predicted for the rows of newdata.
fit_location_scale <- function(features, Z, Y, newdata) {
  data <- cbind(features, Y = Y)
  arm_fit <- function(z) {
    arm <- data[Z == z, ]
    mean_fit <- stats::lm(Y ~ sx1 + x2 + x23, data = arm)
    arm$abs_resid <- abs(stats::residuals(mean_fit))
    scale_fit <- stats::lm(abs_resid ~ x4sq, data = arm)
    s <- sqrt(pi / 2) * stats::predict(scale_fit, newdata = newdata)
    list(m = unname(stats::predict(mean_fit, newdata = newdata)),
         s = pmax(unname(s), SCALE_FLOOR))
  }
  fit0 <- arm_fit(0)
  fit1 <- arm_fit(1)
  list(m0 = fit0$m, m1 = fit1$m, s0 = fit0$s, s1 = fit1$s)
}


#' Nuisance models fitted on one sample and predicted at other units.
#'
#' @param train List from simulate_study() (or a subset of one): the units
#'   the models are fitted on.
#' @param test_features lean_features() of the units to predict for.
#' @param outcomes Subset of "binary", "continuous".
#' @return List with e, and mu0, mu1 (binary) and m0, m1, s0, s1
#'   (continuous), each of length nrow(test_features).
predict_nuisances <- function(train, test_features, outcomes) {
  train_features <- lean_features(train$X)
  pred <- list(e = fit_propensity(train_features, train$Z, test_features))
  if ("binary" %in% outcomes) {
    pred <- c(pred, fit_binary_means(train_features, train$Z, train$Ybinary,
                                     test_features))
  }
  if ("continuous" %in% outcomes) {
    pred <- c(pred, fit_location_scale(train_features, train$Z,
                                       train$Ycontin, test_features))
  }
  pred
}


#' Nuisance predictions for every unit of a study, out of fold.
#'
#' The units are split at random into n_folds folds; each fold's predictions
#' come from models fitted on the other folds. With n_folds = 1 the models
#' are fitted on, and predict for, every unit (in-sample).
#'
#' @param study List from simulate_study() (or a fold of it).
#' @param outcomes Subset of "binary", "continuous".
#' @param n_folds Number of folds, >= 1.
#' @return List with e, and mu0, mu1 (binary) and m0, m1, s0, s1
#'   (continuous), each of length n.
cross_fit_nuisances <- function(study, outcomes, n_folds) {
  features <- lean_features(study$X)
  n <- nrow(features)
  folds <- sample(rep_len(seq_len(n_folds), n))
  nuisance_names <- c("e", if ("binary" %in% outcomes) c("mu0", "mu1"),
                      if ("continuous" %in% outcomes) c("m0", "m1", "s0", "s1"))
  nu <- lapply(setNames(nuisance_names, nuisance_names),
               function(name) numeric(n))
  for (f in seq_len(n_folds)) {
    test <- folds == f
    train <- if (n_folds == 1) test else !test
    pred <- predict_nuisances(subset_study(study, which(train)),
                              features[test, , drop = FALSE], outcomes)
    for (name in nuisance_names) nu[[name]][test] <- pred[[name]]
  }
  nu
}


#' DGK's influence function for the upper bound on E[Y(1)], one column per
#' Lambda, with the normal quantile and CVaR of a location-scale fit.
#'
#' Mirrors phi.continuous() in DGK's simplified_dvds.R with tau = Lambda /
#' (Lambda + 1), Q_tau(x) = m(x) + s(x) Phi^{-1}(tau) in place of their
#' quantile forest and CVaR_tau(x) = m(x) + s(x) phi(Phi^{-1}(tau)) /
#' (1 - tau) in place of their CVaR forest. A is the arm indicator and p its
#' probability, so other bounds follow by the same substitutions as in their
#' dvds.continuous(): -Y, -mu, -m for a lower bound and 1 - Z, 1 - e for the
#' control arm.
#'
#' @param Y Outcomes, length n.
#' @param A Arm indicator (0/1), length n.
#' @param mu Arm mean predictions, length n.
#' @param m,s Location and scale predictions, length n.
#' @param p P(A = 1 | x), length n, in (0, 1).
#' @param lambdas Sensitivity parameters, each >= 1.
#' @return n x length(lambdas) matrix of influence function values.
phi_continuous_lean <- function(Y, A, mu, m, s, p, lambdas) {
  vapply(lambdas, function(lambda) {
    tau <- lambda / (lambda + 1)
    q <- qnorm(tau)
    quantile <- m + s * q
    cvar <- m + s * dnorm(q) / (1 - tau)
    check <- quantile + pmax(0, Y - quantile) / (1 - tau)
    kappa <- mu / lambda + (1 - 1 / lambda) * cvar
    pseudo <- Y / lambda + (1 - 1 / lambda) * check
    A * Y + (1 - A) * kappa + (1 - p) / p * A * (pseudo - kappa)
  }, numeric(length(Y)))
}


#' Continuous-outcome DVDS bounds on the ATE with lean nuisances.
#'
#' Combines four calls of phi_continuous_lean() as DGK's dvds.continuous()
#' combines phi.continuous(): ATE^+ = phi_1^+ - phi_0^- and ATE^- =
#' phi_1^- - phi_0^+, each bound the mean of its influence values, with
#' variance var / n.
#'
#' @param nu List with e, m0, m1, s0 and s1 (cross_fit_nuisances()).
#' @return data.frame as bounds_from_dgk(): outcome, lambda_index, lambda,
#'   lower, upper, s2_lower, s2_upper, n.
dvds_continuous_lean <- function(study, lambdas, nu) {
  Y <- study$Ycontin
  Z <- study$Z
  phi_1_plus <- phi_continuous_lean(Y, Z, nu$m1, nu$m1, nu$s1, nu$e, lambdas)
  phi_1_minus <- -phi_continuous_lean(-Y, Z, -nu$m1, -nu$m1, nu$s1, nu$e,
                                      lambdas)
  phi_0_plus <- phi_continuous_lean(Y, 1 - Z, nu$m0, nu$m0, nu$s0, 1 - nu$e,
                                    lambdas)
  phi_0_minus <- -phi_continuous_lean(-Y, 1 - Z, -nu$m0, -nu$m0, nu$s0,
                                      1 - nu$e, lambdas)
  ate_plus <- phi_1_plus - phi_0_minus
  ate_minus <- phi_1_minus - phi_0_plus
  n <- length(Y)
  data.frame(outcome = "continuous", lambda_index = seq_along(lambdas),
             lambda = lambdas, lower = colMeans(ate_minus),
             upper = colMeans(ate_plus),
             s2_lower = apply(ate_minus, 2, stats::var) / n,
             s2_upper = apply(ate_plus, 2, stats::var) / n, n = n)
}


#' DVDS estimates of the sharp ATE bounds and their variances, per Lambda,
#' with lean nuisances; same output as dvds_bounds() in dgk.R.
#'
#' @param study List from simulate_study() (or a fold of it).
#' @param lambdas Sensitivity parameters, each >= 1.
#' @param outcomes Non-empty subset of "binary", "continuous".
#' @param dgk Environment from load_dgk(), which supplies dvds.binary().
#' @param nuisance_folds Folds for out-of-fold nuisance predictions; 1 fits
#'   in-sample.
#' @return data.frame, one row per (outcome, Lambda); see bounds_from_dgk().
dvds_bounds_lean <- function(study, lambdas, outcomes, dgk,
                             nuisance_folds = NUISANCE_FOLDS) {
  validate_lambdas(lambdas)
  if (min(sum(study$Z), sum(1 - study$Z)) < MIN_ARM_SIZE) {
    stop(sprintf("each arm needs >= %d units, got %d treated, %d control",
                 MIN_ARM_SIZE, sum(study$Z), sum(1 - study$Z)))
  }
  nu <- cross_fit_nuisances(study, outcomes, nuisance_folds)
  bounds_from_nuisances(study, lambdas, outcomes, nu, dgk)
}


#' DVDS bounds of a training sample and a held-out sample, with lean
#' nuisances fitted on the training sample (standard cross-fitting).
#'
#' The training sample's bounds are those of dvds_bounds_lean(), with
#' out-of-fold nuisances within it; they only set the convex weights. The
#' held-out sample's nuisances are predicted by models fitted on the whole
#' training sample, so no unit's nuisances use its own outcome, and its
#' bounds are means over the held-out units given the training sample.
#'
#' @param train,test Lists from simulate_study() (or subsets of one).
#' @param lambdas Sensitivity parameters, each >= 1.
#' @param outcomes Non-empty subset of "binary", "continuous".
#' @param dgk Environment from load_dgk(), which supplies dvds.binary().
#' @param nuisance_folds Folds for the training sample's own out-of-fold
#'   nuisance predictions.
#' @return List with data.frames train and test; see bounds_from_dgk().
dvds_bounds_lean_train_test <- function(train, test, lambdas, outcomes, dgk,
                                        nuisance_folds = NUISANCE_FOLDS) {
  if (length(test$Z) < 2) {
    stop(sprintf("the held-out sample needs >= 2 units, got %d",
                 length(test$Z)))
  }
  train_bounds <- dvds_bounds_lean(train, lambdas, outcomes, dgk,
                                   nuisance_folds)
  nu_test <- predict_nuisances(train, lean_features(test$X), outcomes)
  list(train = train_bounds,
       test = bounds_from_nuisances(test, lambdas, outcomes, nu_test, dgk))
}


#' DVDS bounds of one sample from nuisance predictions at its units.
#'
#' @param study List from simulate_study() (or a subset of one).
#' @param nu List from cross_fit_nuisances() or predict_nuisances(), for the
#'   units of study.
#' @return data.frame, one row per (outcome, Lambda); see bounds_from_dgk().
bounds_from_nuisances <- function(study, lambdas, outcomes, nu, dgk) {
  stopifnot(all(nu$e > 0 & nu$e < 1))  # DGK's code divides by e and 1 - e.
  n <- length(study$Z)
  rows <- list()
  if ("binary" %in% outcomes) {
    binary <- do.call(rbind, lapply(lambdas, function(l) {
      dgk$dvds.binary(Y = study$Ybinary, Z = study$Z, mu0 = nu$mu0,
                      mu1 = nu$mu1, e = nu$e, Lambda = l)
    }))
    rows[[length(rows) + 1]] <- bounds_from_dgk(binary, "binary", n)
  }
  if ("continuous" %in% outcomes) {
    rows[[length(rows) + 1]] <- dvds_continuous_lean(study, lambdas, nu)
  }
  do.call(rbind, rows)
}
