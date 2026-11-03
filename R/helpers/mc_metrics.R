# Monte Carlo summaries of CI performance against a known true band.
# Port of src/helpers/mc_metrics.py; column names match.


#' Mean of x and its Monte Carlo standard error.
#'
#' @param x One value per replicate, length >= 2.
#' @return Named vector c(mean, mcse).
mean_se <- function(x) {
  x <- as.numeric(x)
  if (length(x) < 2) {
    stop(sprintf("x needs at least two entries, got %d", length(x)))
  }
  c(mean = mean(x), mcse = sd(x) / sqrt(length(x)))
}


#' Coverage, CI width, estimated band width, bias and weight summaries.
#'
#' One method in one cell. The estimated band width is U_hat - L_hat, the
#' width of the fused bound estimates without the sampling margin. Each
#' per-replicate quantity q gives columns q (mean) and q_mcse.
#'
#' @param ci List with lo, hi, L_hat, U_hat, w_L, w_U (see ci_methods.R).
#' @param true_lo,true_hi Ends of the true band.
#' @param theta True estimand, inside [true_lo, true_hi].
#' @return One-row data.frame.
summarise_ci <- function(ci, true_lo, true_hi, theta) {
  width <- ci$hi - ci$lo
  per_replicate <- list(
    cov_joint = ci$lo <= true_lo & ci$hi >= true_hi,
    cov_lo = ci$lo <= true_lo,
    cov_hi = ci$hi >= true_hi,
    cov_theta = ci$lo <= theta & ci$hi >= theta,
    width = width,
    band_width = ci$U_hat - ci$L_hat,
    bias_L = ci$L_hat - true_lo,
    bias_U = ci$U_hat - true_hi,
    wL = ci$w_L,
    wU = ci$w_U
  )
  row <- list()
  for (name in names(per_replicate)) {
    estimate <- mean_se(per_replicate[[name]])
    row[[name]] <- estimate[["mean"]]
    row[[paste0(name, "_mcse")]] <- estimate[["mcse"]]
  }
  row$width_median <- median(width)
  row$wL_sd <- sd(ci$w_L)
  row$wU_sd <- sd(ci$w_U)
  row$wL_interior <- mean(ci$w_L > 0 & ci$w_L < 1)
  row$wU_interior <- mean(ci$w_U > 0 & ci$w_U < 1)
  as.data.frame(row)
}
