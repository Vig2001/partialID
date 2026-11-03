# Pilot (scratch): SD of DGK's DVDS bound estimators when the continuous
# outcome noise is scaled by c: Y = m(X) + c (1 + X4^2) eps. Propensity and
# binary outcome as in DGK. sd = sqrt(n * s2).
# Run from the repo root: PILOT_OUT=<dir> Rscript <this file> <reps> <cores>
source(file.path("R", "helpers", "data_fusion_design.R"))
source(file.path("R", "helpers", "dgk.R"))

args <- commandArgs(trailingOnly = TRUE)
reps <- as.integer(args[1]); cores <- as.integer(args[2])
LAMBDAS <- seq(3, 5, by = 0.1)

simulate_noise <- function(n, c_noise) {
  d <- 5
  X <- matrix(runif(n * d, min = -1, max = 1), nrow = n)
  Z <- rbinom(n, 1, 1 / (1 + exp(X[, 1] + 0.5 * (X[, 2] > 0) +
                                   0.5 * X[, 2] * X[, 3])))
  Ybinary <- rbinom(n, 1, 1 / (1 + exp(0.5 * X[, 1] + X[, 2] +
                                         0.25 * X[, 2] * X[, 3])))
  Ycontin <- 2 * sign(X[, 1]) + X[, 2] + X[, 2] * X[, 3] +
    c_noise * (1 + X[, 4]^2) * rnorm(n)
  list(X = X, Z = Z, Ybinary = Ybinary, Ycontin = Ycontin)
}

grid <- expand.grid(rep = seq_len(reps), c_noise = 1,
                    n = 2000)
RNGkind("L'Ecuyer-CMRG"); set.seed(12)
start <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  study <- simulate_noise(g$n, g$c_noise)
  out <- tryCatch(dvds_bounds(study, LAMBDAS, num_threads = 1),
                  error = function(e) NULL)
  if (is.null(out)) return(NULL)
  out <- out[out$outcome == "continuous", ]
  cbind(out, c_noise = g$c_noise, rep = g$rep)
}, mc.cores = cores, mc.set.seed = TRUE)
cat(sprintf("failed fits: %d\n", sum(vapply(res, is.null, logical(1)))))
res <- do.call(rbind, res)
cat(sprintf("time: %.1f min\n", as.numeric(Sys.time() - start, units = "mins")))
write.csv(res, file.path(Sys.getenv("PILOT_OUT"), "pilot_study1_high.csv"),
          row.names = FALSE)
