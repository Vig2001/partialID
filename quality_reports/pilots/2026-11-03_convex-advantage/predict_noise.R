# Scratch: Gaussian prediction of M3 / M6 / M5 widths, continuous outcome,
# study 2 noise scaled by c (band +/- c (4/3) K_Lambda). Bound estimates are
# N(true bound, sd^2 / m), sd from pilot_noise.csv at the study's own n,
# interpolated in Lambda; DVDS bias and small-fold forest quality ignored.
# Run from the repo root: PILOT_OUT=<dir> Rscript <this file>
source(file.path("R", "helpers", "data_fusion_design.R"))
source(file.path("R", "helpers", "ci_methods.R"))
SP <- Sys.getenv("PILOT_OUT")
LAMBDAS <- seq(1, 3, by = 0.05)
R <- 1000
Z <- qnorm(0.975); Z_CAL <- qnorm(1 - 0.05 / 4)

pilot <- read.csv(file.path(SP, "pilot_noise.csv"))
pilot$sd_lower <- sqrt(pilot$s2_lower * pilot$n)
pilot$sd_upper <- sqrt(pilot$s2_upper * pilot$n)
sds <- aggregate(cbind(sd_lower, sd_upper) ~ c_noise + n + lambda,
                 data = pilot, FUN = median)
check <- aggregate(cbind(sd_lower, sd_upper) ~ c_noise + n, data = pilot,
                   FUN = median)
check$ratio_to_c1 <- NA
for (k in seq_len(nrow(check))) {
  base <- check[check$c_noise == 1 & check$n == check$n[k], ]
  check$ratio_to_c1[k] <- (check$sd_lower[k] + check$sd_upper[k]) /
    (base$sd_lower + base$sd_upper)
}
cat("Median sd by c and n (ratio to c = 1):\n"); print(check, digits = 3)

sd_curve <- function(c_noise, n, end) {
  s <- sds[sds$c_noise == c_noise & sds$n == n, ]
  approx(s$lambda, s[[paste0("sd_", end)]], xout = LAMBDAS, rule = 2)$y
}
band <- function(c_noise) {
  h <- c_noise * MEAN_SIGMA * k_lambda(LAMBDAS)
  list(lower = -h, upper = h)
}

predict_design <- function(n1, n2, c2, frac_a) {
  set.seed(2)
  n <- c(n1, n2); n_a <- floor(frac_a * n); n_b <- n - n_a
  t1 <- band(1); t2 <- band(c2)
  sd1 <- list(lower = sd_curve(1, n1, "lower"), upper = sd_curve(1, n1, "upper"))
  sd2 <- list(lower = sd_curve(c2, n2, "lower"),
              upper = sd_curve(c2, n2, "upper"))
  rows <- vector("list", length(LAMBDAS)^2); r <- 0
  for (i in seq_along(LAMBDAS)) for (j in seq_along(LAMBDAS)) {
    draw <- function(m) {
      out <- list()
      for (end in c("lower", "upper")) {
        mu <- c(t1[[end]][i], t2[[end]][j]); sd <- c(sd1[[end]][i], sd2[[end]][j])
        out[[end]] <- cbind(rnorm(R, mu[1], sd[1] / sqrt(m[1])),
                            rnorm(R, mu[2], sd[2] / sqrt(m[2])))
        out[[paste0("s2_", end)]] <- matrix(sd^2 / m, R, 2, byrow = TRUE)
      }
      out
    }
    full <- draw(n); a <- draw(n_a); b <- draw(n_b)
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
    # Bands are nested (both symmetric about 0): the tighter study is the
    # only valid single-study CI.
    k <- if (t1$upper[i] <= t2$upper[j]) 1 else 2
    sd_k <- if (k == 1) c(sd1$lower[i], sd1$upper[i]) else
      c(sd2$lower[j], sd2$upper[j])
    tu <- min(t1$upper[i], t2$upper[j])
    r <- r + 1
    rows[[r]] <- data.frame(
      lambda_1 = LAMBDAS[i], lambda_2 = LAMBDAS[j], tighter = k,
      rel_gap = abs(t1$upper[i] - t2$upper[j]) / sqrt(sd1$upper[i]^2 / n1 +
                                                       sd2$upper[j]^2 / n2),
      width_3 = mean(m3$hi - m3$lo), width_6 = mean(m6$hi - m6$lo),
      width_5 = 2 * tu + Z * sum(sd_k) / sqrt(n[k]),
      wL = mean(w_L), wU = mean(w_U))
  }
  d <- do.call(rbind, rows)
  d$r36 <- d$width_3 / d$width_6; d$r35 <- d$width_3 / d$width_5; d
}

designs <- expand.grid(frac_a = c(0.2, 0.3), c2 = c(1, 1.5, 2, 2.5),
                       n = c("2000,500", "1000,500"), stringsAsFactors = FALSE)
fits <- parallel::mclapply(seq_len(nrow(designs)), function(k) {
  nn <- as.integer(strsplit(designs$n[k], ",")[[1]])
  predict_design(nn[1], nn[2], designs$c2[k], designs$frac_a[k])
}, mc.cores = 8)
out <- do.call(rbind, lapply(seq_len(nrow(designs)), function(k) {
  d <- fits[[k]]; win <- d$r36 < 0.99 & d$r35 < 0.99
  best <- d[win, ][which.min(pmax(d$r36, d$r35)[win]), ]
  data.frame(designs[k, ], fusion_wins = sum(win),
             best_vs_M6 = if (any(win)) round(best$r36, 3) else NA,
             best_vs_M5 = if (any(win)) round(best$r35, 3) else NA,
             best_l1 = if (any(win)) best$lambda_1 else NA,
             best_l2 = if (any(win)) best$lambda_2 else NA,
             wins_l1_ge_2 = sum(win & d$lambda_1 >= 2 & d$lambda_2 >= 1.5),
             median_r36 = round(median(d$r36), 3))
}))
print(out, row.names = FALSE)
detail <- do.call(rbind, lapply(seq_len(nrow(designs)), function(k) {
  cbind(designs[k, ], fits[[k]])
}))
write.csv(detail, file.path(SP, "predict_noise_pairs.csv"), row.names = FALSE)
