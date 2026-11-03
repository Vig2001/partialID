# Pilot (scratch): SD of DGK's DVDS bound estimators when the propensity
# logit is scaled by gamma. sd = sqrt(n * s2), so SE at any n is sd / sqrt(n).
# Run from the repo root: Rscript <this file> <reps> <cores>
source(file.path("R", "helpers", "data_fusion_design.R"))
source(file.path("R", "helpers", "dgk.R"))

args <- commandArgs(trailingOnly = TRUE)
reps <- as.integer(args[1]); cores <- as.integer(args[2])
LAMBDAS <- seq(1, 3, length.out = 20)

simulate_gamma <- function(n, gamma) {
  d <- 5
  X <- matrix(runif(n * d, min = -1, max = 1), nrow = n)
  lin <- X[, 1] + 0.5 * (X[, 2] > 0) + 0.5 * X[, 2] * X[, 3]
  Z <- rbinom(n, 1, 1 / (1 + exp(gamma * lin)))
  Ybinary <- rbinom(n, 1, 1 / (1 + exp(0.5 * X[, 1] + X[, 2] +
                                         0.25 * X[, 2] * X[, 3])))
  Ycontin <- 2 * sign(X[, 1]) + X[, 2] + X[, 2] * X[, 3] +
    (1 + X[, 4]^2) * rnorm(n)
  list(X = X, Z = Z, Ybinary = Ybinary, Ycontin = Ycontin)
}

grid <- expand.grid(rep = seq_len(reps), gamma = c(1, 2, 3),
                    n = c(1000, 5000))
RNGkind("L'Ecuyer-CMRG"); set.seed(11)
start <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, ]
  study <- simulate_gamma(g$n, g$gamma)
  out <- tryCatch(dvds_bounds(study, LAMBDAS, num_threads = 1),
                  error = function(e) NULL)
  if (is.null(out)) return(NULL)
  cbind(out, gamma = g$gamma, rep = g$rep, treated = mean(study$Z))
}, mc.cores = cores, mc.set.seed = TRUE)
cat(sprintf("failed fits: %d\n", sum(vapply(res, is.null, logical(1)))))
res <- do.call(rbind, res)
cat(sprintf("time: %.1f min\n", as.numeric(Sys.time() - start, units = "mins")))
write.csv(res, file.path(Sys.getenv("PILOT_OUT"), "pilot_se.csv"),
          row.names = FALSE)
