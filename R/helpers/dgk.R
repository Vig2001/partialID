# DGK's simulation code (Dorn, Guo and Kallus; 2025, Section 7.1), used as
# is. Their repository is cloned at a pinned commit in external/ (README).
# Data generation and nuisance fits are copied verbatim from their
# dvds_sim.R (line numbers cited); the DVDS estimator is their
# dvds.continuous(), sourced unmodified from simplified_dvds.R.

DGK_DIR <- file.path("external", "DoublyValidDoublySharpSensitivity")
DGK_COMMIT <- "26cbbf6e90dcd0648c837c0d58568375291dd522"  # 2024-03-22
DGK_SCRIPTS <- c("dvds_sim.R", "simplified_dvds.R", "dvds_plotting.r")
DGK_ALPHA <- 0.05  # Default alpha of dvds.continuous(); sets its CI width.
GRF_FORESTS <- c("probability_forest", "regression_forest", "quantile_forest")
MIN_ARM_SIZE <- 20  # Fewest treated or control units for which forests are fit.


#' Stops unless the packages DGK's scripts load are installed.
#'
#' Named literally so that renv records them in renv.lock.
check_dgk_packages <- function() {
  ok <- c(grf = requireNamespace("grf", quietly = TRUE),
          dplyr = requireNamespace("dplyr", quietly = TRUE),
          ggplot2 = requireNamespace("ggplot2", quietly = TRUE),
          progress = requireNamespace("progress", quietly = TRUE))
  if (!all(ok)) {
    stop(sprintf("missing packages %s; run renv::restore()",
                 paste(names(ok)[!ok], collapse = ", ")))
  }
}


#' Runs git inside the clone; returns its output lines.
git_dgk <- function(...) {
  system2("git", c("-C", DGK_DIR, ...), stdout = TRUE, stderr = FALSE)
}


#' Stops unless the clone is at the pinned commit with unedited scripts, and
#' the packages they load are installed.
check_clone <- function() {
  check_dgk_packages()
  if (!dir.exists(DGK_DIR)) {
    stop(sprintf("clone not found; run: git clone %s %s",
                 "https://github.com/CausalML/DoublyValidDoublySharpSensitivity",
                 DGK_DIR))
  }
  head <- git_dgk("rev-parse", "HEAD")
  if (!identical(head, DGK_COMMIT)) {
    stop(sprintf("clone is at %s, expected %s", head, DGK_COMMIT))
  }
  edited <- git_dgk("status", "--porcelain", "--", DGK_SCRIPTS)
  if (length(edited) > 0) {
    stop(sprintf("DGK scripts are edited: %s", paste(edited, collapse = "; ")))
  }
}


#' Environment holding DGK's simplified_dvds.R functions.
#'
#' Their code never sets grf's thread count. The forest constructors are
#' masked here by versions with num.threads fixed, so parallel workers do
#' not oversubscribe the cores; their file is untouched, and grf output does
#' not depend on the thread count (checked, grf 2.6.1).
#'
#' @param num_threads grf threads per forest; NULL uses every core.
load_dgk <- function(num_threads = NULL) {
  dgk <- new.env(parent = globalenv())
  for (name in GRF_FORESTS) {
    dgk[[name]] <- local({
      forest <- getExportedValue("grf", name)
      function(...) forest(..., num.threads = num_threads)
    })
  }
  suppressPackageStartupMessages(
    sys.source(file.path(DGK_DIR, "simplified_dvds.R"), envir = dgk))
  dgk
}


#' Draws one study from the DGK DGP: dvds_sim.R lines 23-26, verbatim.
#'
#' Ybinary is drawn and discarded so that the use of random numbers matches
#' theirs.
#'
#' @param n Sample size.
#' @return List with X (n x 5), Z (0/1) and Y (continuous outcome).
simulate_dgk <- function(n) {
  if (n < 2) stop(sprintf("n must be >= 2, got n = %s", n))
  d = 5
  X = matrix(runif(n*d, min=-1, max=1), nrow=n)
  Z = rbinom(n, 1, 1/(1 + exp(X[,1] + 0.5*(X[,2]>0) + 0.5*X[,2]*X[,3])))
  Ybinary = rbinom(n, 1, 1/(1 + exp(0.5*X[,1] + X[,2] + 0.25*X[,2]*X[,3])))
  Ycontin = 2*sign(X[,1]) + X[,2] + X[,2]*X[,3] + (1 + X[,4]^2)*rnorm(n)
  list(X = X, Z = Z, Y = Ycontin)
}


#' Propensity and arm-wise outcome means, all out-of-bag where in sample.
#'
#' dvds_sim.R lines 29, 30 (continuous part), 33-34 and 39-42, verbatim.
#' Their binary-outcome forests (lines 31-32) are not fitted.
#'
#' @param dgk Environment from load_dgk(), which supplies the forests.
#' @return List with e, mu0 and mu1, each of length n.
dgk_nuisances <- function(X, Z, Ycontin, dgk) {
  scope <- new.env(parent = dgk)
  scope$X <- X
  scope$Z <- Z
  scope$Ycontin <- Ycontin
  scope$n <- length(Z)
  local({
    e = predict(probability_forest(X=X, Y=factor(Z)))$prediction[,"1"]
    mu0contin = mu1contin = numeric(n)
    fit0contin = regression_forest(X = X[Z==0,], Y = Ycontin[Z==0])
    fit1contin = regression_forest(X = X[Z==1,], Y = Ycontin[Z==1])
    mu0contin[Z==0] = predict(fit0contin)$predictions
    mu0contin[Z==1] = predict(fit0contin, X[Z==1,])$predictions
    mu1contin[Z==0] = predict(fit1contin, X[Z==0,])$predictions
    mu1contin[Z==1] = predict(fit1contin)$predictions
  }, envir = scope)
  list(e = scope$e, mu0 = scope$mu0contin, mu1 = scope$mu1contin)
}


#' DVDS estimates of the sharp ATE bounds and their variances, per Lambda.
#'
#' Calls DGK's dvds.continuous(). It returns bounds and CIs, not variances;
#' its CI is mean -/+ z sd / sqrt(n) with z = qnorm(1 - DGK_ALPHA / 2), so
#' the variance of each bound is recovered exactly as ((CI - mean) / z)^2.
#'
#' @param study List with X, Z and Y.
#' @param lambdas Sensitivity parameters, each >= 1.
#' @param num_threads grf threads per forest; NULL uses every core.
#' @return data.frame, one row per Lambda: lambda_index, lambda, lower,
#'   upper, s2_lower, s2_upper and n.
dvds_bounds <- function(study, lambdas, num_threads = NULL) {
  validate_lambdas(lambdas)
  if (min(sum(study$Z), sum(1 - study$Z)) < MIN_ARM_SIZE) {
    stop(sprintf("each arm needs >= %d units, got %d treated, %d control",
                 MIN_ARM_SIZE, sum(study$Z), sum(1 - study$Z)))
  }
  dgk <- load_dgk(num_threads)
  nu <- dgk_nuisances(study$X, study$Z, study$Y, dgk)
  # Their code divides by e and 1 - e without a check.
  stopifnot(all(nu$e > 0 & nu$e < 1))
  out <- dgk$dvds.continuous(X = study$X, Y = study$Y, Z = study$Z,
                             mu0 = nu$mu0, mu1 = nu$mu1, e = nu$e,
                             Lambdas = lambdas)
  z <- qnorm(1 - DGK_ALPHA / 2)
  data.frame(
    lambda_index = seq_along(lambdas),
    lambda = out$Lambda,
    lower = out$lower,
    upper = out$upper,
    s2_lower = ((out$lower - out$lower.CI) / z)^2,
    s2_upper = ((out$upper.CI - out$upper) / z)^2,
    n = length(study$Y)
  )
}
