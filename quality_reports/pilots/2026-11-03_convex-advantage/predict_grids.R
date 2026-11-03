# Scratch: as predict_noise.R, but separate grids: study 1 (c = 1, n1) on
# Lambda in [3, 5], study 2 (noise c, n2) on [1, 3]. sd from the pilots.
source(file.path("R", "helpers", "data_fusion_design.R"))
source(file.path("R", "helpers", "ci_methods.R"))
SP <- Sys.getenv("PILOT_OUT")
L1 <- seq(3, 5, by = 0.05); L2 <- seq(1, 3, by = 0.05)
R <- 1000; Z <- qnorm(0.975); Z_CAL <- qnorm(1 - 0.05 / 4)
n1 <- 2000; n2 <- 500
p1 <- read.csv(file.path(SP, "pilot_study1_high.csv"))
p2 <- read.csv(file.path(SP, "pilot_noise.csv")); p2 <- p2[p2$n == n2, ]
sd_curve <- function(p, grid, end, c_noise = NULL) {
  if (!is.null(c_noise)) p <- p[p$c_noise == c_noise, ]
  p$sd <- sqrt(p[[paste0("s2_", end)]] * p$n)
  s <- aggregate(sd ~ lambda, data = p, FUN = median)
  approx(s$lambda, s$sd, xout = grid, rule = 2)$y
}
# Study 1 at large Lambda: SD and bias of the estimates (10 reps).
t1 <- MEAN_SIGMA * k_lambda(p1$lambda)
b <- aggregate(cbind(sd_up = sqrt(s2_upper * n), bias_up = upper - t1,
                     bias_lo = lower + t1) ~ lambda, data = p1, FUN = mean)
b$bias_up_in_se <- b$bias_up / (b$sd_up / sqrt(n1))
b$bias_lo_in_se <- b$bias_lo / (b$sd_up / sqrt(n1))
cat("Study 1, n = 2000, Lambda 3-5 (means over 10 reps):\n")
print(b[b$lambda %in% c(3, 3.5, 4, 4.5, 5), ], digits = 3, row.names = FALSE)

predict <- function(c2, frac_a) {
  set.seed(2)
  n <- c(n1, n2); n_a <- floor(frac_a * n); n_b <- n - n_a
  h1 <- MEAN_SIGMA * k_lambda(L1); h2 <- c2 * MEAN_SIGMA * k_lambda(L2)
  s1 <- list(lower = sd_curve(p1, L1, "lower"), upper = sd_curve(p1, L1, "upper"))
  s2 <- list(lower = sd_curve(p2, L2, "lower", c2),
             upper = sd_curve(p2, L2, "upper", c2))
  rows <- list()
  for (i in seq_along(L1)) for (j in seq_along(L2)) {
    # Only pairs near the matching curve can show interior weights.
    if (abs(h1[i] - h2[j]) > 0.25) next
    tt <- list(lower = c(-h1[i], -h2[j]), upper = c(h1[i], h2[j]))
    sd <- list(lower = c(s1$lower[i], s2$lower[j]),
               upper = c(s1$upper[i], s2$upper[j]))
    draw <- function(m) {
      out <- list()
      for (end in c("lower", "upper")) {
        out[[end]] <- cbind(rnorm(R, tt[[end]][1], sd[[end]][1] / sqrt(m[1])),
                            rnorm(R, tt[[end]][2], sd[[end]][2] / sqrt(m[2])))
        out[[paste0("s2_", end)]] <- matrix(sd[[end]]^2 / m, R, 2, byrow = TRUE)
      }
      out
    }
    full <- draw(n); a <- draw(n_a); bb <- draw(n_b)
    sc <- matrix(n_a / n_b, R, 2, byrow = TRUE)
    w_L <- optimal_weight_grid(a$lower[, 1], a$lower[, 2], a$s2_lower[, 1] * sc[, 1],
                               a$s2_lower[, 2] * sc[, 2], TRUE, Z)
    w_U <- optimal_weight_grid(a$upper[, 1], a$upper[, 2], a$s2_upper[, 1] * sc[, 1],
                               a$s2_upper[, 2] * sc[, 2], FALSE, Z)
    m3 <- convex_ci(bb$lower, bb$upper, bb$s2_lower, bb$s2_upper, w_L, w_U, Z)
    m6 <- intersection_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, Z_CAL)
    k <- if (h1[i] <= h2[j]) 1 else 2
    w5 <- 2 * min(h1[i], h2[j]) + Z * (sd$lower[k] + sd$upper[k]) / sqrt(n[k])
    rows[[length(rows) + 1]] <- data.frame(
      c2 = c2, frac_a = frac_a, lambda_1 = L1[i], lambda_2 = L2[j],
      w3 = mean(m3$hi - m3$lo), w6 = mean(m6$hi - m6$lo), w5 = w5,
      wL = mean(w_L), wU = mean(w_U))
  }
  d <- do.call(rbind, rows); d$r36 <- d$w3 / d$w6; d$r35 <- d$w3 / d$w5; d
}
designs <- expand.grid(c2 = c(2, 2.5), frac_a = c(0.2, 0.3))
fits <- parallel::mclapply(seq_len(nrow(designs)), function(k)
  predict(designs$c2[k], designs$frac_a[k]), mc.cores = 4)
d <- do.call(rbind, fits)
d$win <- d$r36 < 0.99 & d$r35 < 0.99
write.csv(d, file.path(SP, "predict_grids_pairs.csv"), row.names = FALSE)
s <- do.call(rbind, lapply(split(d, list(d$c2, d$frac_a)), function(x) {
  w <- x[x$win, ]
  data.frame(c2 = x$c2[1], frac_a = x$frac_a[1], wins = nrow(w),
             best_r36 = round(min(w$r36), 3), med_r36 = round(median(w$r36), 3),
             med_r35 = round(median(w$r35), 3),
             lam1 = paste(range(w$lambda_1), collapse = "-"),
             lam2 = paste(range(w$lambda_2), collapse = "-"),
             med_w = round(median(w$wL), 2))
}))
print(s, row.names = FALSE)
x <- d[d$c2 == 2.5 & d$frac_a == 0.3 & d$lambda_1 %in% c(3, 4, 5) & d$win, ]
print(x[, c("lambda_1", "lambda_2", "r36", "r35", "wL", "wU")], digits = 3,
      row.names = FALSE)
