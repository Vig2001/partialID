# Scratch: Gaussian prediction of M3 / M6 widths for candidate designs.
# Bound estimates are drawn as N(true bound, sd^2 / m), with sd measured in
# the pilot (pilot_se.csv) and m the size of the sample that is used, so
# DVDS bias and forest quality in small folds are ignored. M3 and M6 are the
# causal sim's own functions (R/helpers/ci_methods.R).
# Run from the repo root: Rscript <this file>
source(file.path("R", "helpers", "data_fusion_design.R"))
source(file.path("R", "helpers", "ci_methods.R"))
SP <- Sys.getenv("PILOT_OUT")
LAMBDAS <- seq(1, 3, length.out = 20)
R <- 1000
Z <- qnorm(0.975); Z_CAL <- qnorm(1 - 0.05 / 4)

# --- sd per (outcome, gamma, lambda, end), pooled over reps and n ---
pilot <- read.csv(file.path(SP, "pilot_se.csv"))
pilot$sd_lower <- sqrt(pilot$s2_lower * pilot$n)
pilot$sd_upper <- sqrt(pilot$s2_upper * pilot$n)
# The sd grows with n when gamma > 1 (forests smooth e at small n), so each
# study uses the pilot n nearest its own size on the log scale.
sds <- aggregate(cbind(sd_lower, sd_upper) ~ outcome + gamma + n +
                   lambda_index, data = pilot, FUN = median)
pilot_n <- function(n) {
  ns <- sort(unique(pilot$n)); ns[which.min(abs(log(ns) - log(n)))]
}
by_n <- aggregate(cbind(sd_lower, sd_upper) ~ outcome + gamma + n,
                  data = pilot, FUN = median)
cat("Median sd by n (checks 1/sqrt(n) scaling):\n"); print(by_n)

# --- truth: binary band depends on gamma through e(X) ---
binary_band_gamma <- function(gamma) {
  set.seed(1); n <- 5e5
  X <- matrix(runif(n * 5, -1, 1), nrow = n)
  e <- 1 / (1 + exp(gamma * (X[, 1] + 0.5 * (X[, 2] > 0) +
                               0.5 * X[, 2] * X[, 3])))
  mu <- 1 / (1 + exp(0.5 * X[, 1] + X[, 2] + 0.25 * X[, 2] * X[, 3]))
  up <- sapply(LAMBDAS, function(l) mean(
    (e * mu + (1 - e) * pmin(1 - 1 / l + mu / l, mu * l)) -
      ((1 - e) * mu + e * pmax(1 - l + mu * l, mu / l))))
  lo <- sapply(LAMBDAS, function(l) mean(
    (e * mu + (1 - e) * pmax(1 - l + mu * l, mu / l)) -
      ((1 - e) * mu + e * pmin(1 - 1 / l + mu / l, mu * l))))
  list(lower = lo, upper = up)
}
truth <- function(outcome, gamma) {
  if (outcome == "continuous") {
    b <- continuous_true_band(LAMBDAS); list(lower = b$lower, upper = b$upper)
  } else binary_band_gamma(gamma)
}
TRUTH <- list()
for (o in c("binary", "continuous")) for (g in c(1, 2, 3)) {
  TRUTH[[paste(o, g)]] <- truth(o, g)
}

sd_of <- function(o, g, end, n) {
  s <- sds[sds$outcome == o & sds$gamma == g & sds$n == pilot_n(n), ]
  s[order(s$lambda_index), paste0("sd_", end)]
}

# --- one design: study 1 gamma = 1, study 2 gamma = g2 ---
predict_design <- function(n1, n2, g2, frac_a, outcome) {
  set.seed(2)
  n <- c(n1, n2); n_a <- floor(frac_a * n); n_b <- n - n_a
  t1 <- TRUTH[[paste(outcome, 1)]]; t2 <- TRUTH[[paste(outcome, g2)]]
  sd1 <- list(lower = sd_of(outcome, 1, "lower", n1),
              upper = sd_of(outcome, 1, "upper", n1))
  sd2 <- list(lower = sd_of(outcome, g2, "lower", n2),
              upper = sd_of(outcome, g2, "upper", n2))
  rows <- list()
  for (i in seq_along(LAMBDAS)) for (j in seq_along(LAMBDAS)) {
    draw <- function(m) {
      out <- list()
      for (end in c("lower", "upper")) {
        mean <- c(t1[[end]][i], t2[[end]][j])
        sd <- c(sd1[[end]][i], sd2[[end]][j])
        out[[end]] <- cbind(rnorm(R, mean[1], sd[1] / sqrt(m[1])),
                            rnorm(R, mean[2], sd[2] / sqrt(m[2])))
        out[[paste0("s2_", end)]] <- matrix(sd^2 / m, R, 2, byrow = TRUE)
      }
      out
    }
    full <- draw(n); a <- draw(n_a); b <- draw(n_b)
    # Fold-A variances rescaled to the fold-B size, as in pair_cis().
    sc <- matrix(n_a / n_b, R, 2, byrow = TRUE)
    w_L <- optimal_weight_grid(a$lower[, 1], a$lower[, 2],
                               a$s2_lower[, 1] * sc[, 1],
                               a$s2_lower[, 2] * sc[, 2], TRUE, Z)
    w_U <- optimal_weight_grid(a$upper[, 1], a$upper[, 2],
                               a$s2_upper[, 1] * sc[, 1],
                               a$s2_upper[, 2] * sc[, 2], FALSE, Z)
    m3 <- convex_ci(b$lower, b$upper, b$s2_lower, b$s2_upper, w_L, w_U, Z)
    m6 <- intersection_ci(full$lower, full$upper, full$s2_lower,
                          full$s2_upper, Z_CAL)
    tl <- max(t1$lower[i], t2$lower[j]); tu <- min(t1$upper[i], t2$upper[j])
    # Single-study CI (z = 1.96, full sample); valid only if that study
    # attains both ends of the intersection band.
    w5 <- c(t1$upper[i] - t1$lower[i] +
              Z * (sd1$lower[i] + sd1$upper[i]) / sqrt(n1),
            t2$upper[j] - t2$lower[j] +
              Z * (sd2$lower[j] + sd2$upper[j]) / sqrt(n2))
    valid5 <- c(t1$lower[i] >= t2$lower[j] & t1$upper[i] <= t2$upper[j],
                t2$lower[j] >= t1$lower[i] & t2$upper[j] <= t1$upper[i])
    rows[[length(rows) + 1]] <- data.frame(
      lambda_1 = LAMBDAS[i], lambda_2 = LAMBDAS[j],
      width_3 = mean(m3$hi - m3$lo), width_6 = mean(m6$hi - m6$lo),
      width_5_valid = if (any(valid5)) min(w5[valid5]) else Inf,
      wL = mean(w_L), wU = mean(w_U),
      cov_3 = mean(m3$lo <= tl & m3$hi >= tu),
      cov_6 = mean(m6$lo <= tl & m6$hi >= tu))
  }
  d <- do.call(rbind, rows); d$ratio <- d$width_3 / d$width_6; d
}

designs <- expand.grid(frac_a = c(0.1, 0.15, 0.2, 0.3), g2 = c(1, 2, 3),
                       n = c("2000,500", "5000,1000"),
                       outcome = c("binary", "continuous"),
                       stringsAsFactors = FALSE)
fits <- parallel::mclapply(seq_len(nrow(designs)), function(k) {
  ds <- designs[k, ]; nn <- as.integer(strsplit(ds$n, ",")[[1]])
  predict_design(nn[1], nn[2], ds$g2, ds$frac_a, ds$outcome)
}, mc.cores = 8)
out <- list(); detail <- list()
for (k in seq_len(nrow(designs))) {
  ds <- designs[k, ]; nn <- as.integer(strsplit(ds$n, ",")[[1]])
  d <- fits[[k]]
  best <- d[which.min(d$ratio), ]
  out[[k]] <- data.frame(ds, fold_a_2 = floor(ds$frac_a * nn[2]),
                         m3_narrower = sum(d$ratio < 0.99),
                         fusion_wins = sum(d$ratio < 0.99 &
                                             d$width_3 < d$width_5_valid),
                         best_fusion = round(min(c(Inf, d$ratio[
                           d$width_3 < d$width_5_valid])), 3),
                         min_ratio = round(min(d$ratio), 3),
                         median_ratio = round(median(d$ratio), 3),
                         best_l1 = round(best$lambda_1, 2),
                         best_l2 = round(best$lambda_2, 2),
                         min_cov_3 = round(min(d$cov_3), 3),
                         min_cov_6 = round(min(d$cov_6), 3))
  detail[[k]] <- cbind(ds, d)
}
out <- do.call(rbind, out)
print(out, row.names = FALSE)
write.csv(out, file.path(SP, "predict_designs.csv"), row.names = FALSE)
write.csv(do.call(rbind, detail), file.path(SP, "predict_pairs.csv"),
          row.names = FALSE)
