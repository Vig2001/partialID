# Two observational studies (OS) from the DGK DGP, differing only in sample
# size n_k and sensitivity parameter Lambda_k. Each study has a binary and a
# continuous outcome (DGK's two designs, drawn from the same X and Z). For
# each outcome, each study gives DVDS estimates of its sharp ATE bounds
# [L_k, U_k], computed with DGK's own code (R/helpers/dgk.R; their repository
# cloned in external/). CIs for the fused band [max_k L_k, min_k U_k] are
# compared over all pairs (Lambda_1, Lambda_2), separately per outcome:
#   M6  calibrated intersection       (full-sample fits, z = Phi^{-1}(1 - a/4))
#   M5  each study alone              (full-sample fits; DGK's CI)
#   M3  convex split                  (weights from fold A, CI from fold B;
#                                      DVDS fitted separately in each fold)
# The naive intersection (M1) is omitted: it undercovers at small gaps, so
# width comparisons with it are not like-for-like.
#
# Usage, from the repo root:
#   Rscript R/causal_sim.R [--seed 1] [--reps 1000] [--n1 2000] [--n2 500]
#     [--lambdas 1,1.5,...] [--frac-a 0.3] [--alpha 0.05] [--workers W]
#     [--backend multisession]
#     [--out-dir src/outputs/causal_sim] [--summarise-only]
#
# Plans:   quality_reports/plans/2026-11-02_dvds-r-two-os.md,
#          2026-11-02_dgk-exact-replication.md, 2026-11-02_causal-sim-binary.md
# Notes:   partialID_Overleaf/Notes/Project_Notes/causal_sim_code_explained.tex
# Outputs: src/outputs/causal_sim/

# Paths are relative to the repo root, where .Rprofile also activates renv.
if (!file.exists(file.path("R", "causal_sim.R"))) {
  stop("run from the repo root: Rscript R/causal_sim.R")
}
for (helper in c("data_fusion_design.R", "dgk.R", "ci_methods.R",
                 "mc_metrics.R", "run_time.R")) {
  source(file.path("R", "helpers", helper))
}

# Default Lambda grid for both studies: DGK's (dvds_sim.R line 10). The
# number of values does not affect run time (forests are fitted once per fit).
LAMBDAS <- seq(1, 3, length.out = 20)
FITS <- c("full", "fold_a", "fold_b")
OUTCOMES <- c("binary", "continuous")
# Progress-bar labels; short, so the bar fits an 80-column console.
FIT_LABELS <- c(full = "intersection", fold_a = "convex split A",
                fold_b = "convex split B")
M3 <- "M3 convex (split)"
M5_D1 <- "M5 dataset 1 only"
M5_D2 <- "M5 dataset 2 only"
M6 <- "M6 intersection (calibrated)"
MIN_FOLD_SIZE <- 100  # Smallest fold; DVDS fits forests within each arm.
THETA <- 0  # True ATE: neither outcome depends on Z in the DGK DGP.


# ----------------------------- configuration ------------------------------
parse_args <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opt <- function(...) optparse::make_option(...)
  options <- list(
    opt("--seed", type = "integer", default = 1,
        help = "Seed for the per-replicate streams [default %default]."),
    opt("--reps", type = "integer", default = 1000,
        help = "Monte Carlo replicates [default %default]."),
    opt("--n1", type = "integer", default = 2000,
        help = "Size of study 1 [default %default]."),
    opt("--n2", type = "integer", default = 500,
        help = "Size of study 2 [default %default]."),
    opt("--lambdas", type = "character",
        default = NA_character_,
        help = paste("Comma-separated Lambda grid for both studies; must",
                     "include 1 [default: DGK's seq(1, 3, length.out = 20)].")),
    opt("--frac-a", dest = "frac_a", type = "double", default = 0.3,
        help = "Fraction of each study in fold A [default %default]."),
    opt("--alpha", type = "double", default = 0.05,
        help = paste("Nominal level; each end uses z = qnorm(1 - alpha / 2)",
                     "[default %default].")),
    opt("--workers", type = "integer",
        default = parallelly::availableCores(omit = 1),
        help = "Parallel workers [default: available cores minus 1]."),
    opt("--backend", type = "character", default = "multisession",
        help = paste("future backend: multisession, multicore or sequential",
                     "[default %default].")),
    opt("--out-dir", dest = "out_dir", type = "character",
        default = file.path("src", "outputs", "causal_sim"),
        help = "Output directory [default %default]."),
    opt("--summarise-only", dest = "summarise_only", action = "store_true",
        default = FALSE,
        help = "Recompute CIs and summaries from estimates.csv; no fitting.")
  )
  optparse::parse_args(optparse::OptionParser(option_list = options),
                       args = argv)
}


#' Validated settings, with the critical values implied by alpha.
make_config <- function(args) {
  n <- c(args$n1, args$n2)
  n_a <- floor(args$frac_a * n)
  lambdas <- if (is.na(args$lambdas)) LAMBDAS else
    as.numeric(strsplit(args$lambdas, ",")[[1]])
  validate_lambdas(lambdas)
  # Lambda = 1 anchors the collapse check (check_collapse_at_one).
  if (!(1 %in% lambdas) || is.unsorted(lambdas, strictly = TRUE)) {
    stop(sprintf("lambdas must be strictly increasing and include 1, got %s",
                 args$lambdas))
  }
  if (args$reps < 2) stop(sprintf("reps must be >= 2, got %d", args$reps))
  if (!(args$alpha > 0 && args$alpha < 1)) {
    stop(sprintf("alpha must lie in (0, 1), got %s", args$alpha))
  }
  if (!(args$frac_a > 0 && args$frac_a < 1)) {
    stop(sprintf("frac_a must lie in (0, 1), got %s", args$frac_a))
  }
  if (min(n_a, n - n_a) < MIN_FOLD_SIZE) {
    stop(sprintf("every fold needs >= %d units, got fold sizes %s",
                 MIN_FOLD_SIZE, paste(c(n_a, n - n_a), collapse = ", ")))
  }
  if (!(args$backend %in% c("multisession", "multicore", "sequential"))) {
    stop(sprintf("unknown backend %s", args$backend))
  }
  if (args$workers < 1) stop(sprintf("workers must be >= 1, got %d",
                                     args$workers))
  list(seed = args$seed, reps = args$reps, n = n, frac_a = args$frac_a,
       alpha = args$alpha, z = qnorm(1 - args$alpha / 2),
       z_calibrated = qnorm(1 - args$alpha / 4), lambdas = lambdas,
       workers = args$workers, backend = args$backend,
       out_dir = args$out_dir)
}


# ----------------------------- simulation ------------------------------
#' One L'Ecuyer-CMRG stream per replicate, derived from seed.
#'
#' Replicate r always uses stream r, so results do not depend on the number
#' of workers or the order in which replicates finish.
replicate_seeds <- function(seed, reps) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  seeds <- vector("list", reps)
  current <- .Random.seed
  for (r in seq_len(reps)) {
    seeds[[r]] <- current
    current <- parallel::nextRNGStream(current)
  }
  seeds
}


#' DVDS estimates for one replicate: 2 studies x 3 fits (full, fold A, B),
#' each fit giving both outcomes.
#'
#' @param rep Replicate number.
#' @param progress progressr progressor; receives a label before each fit
#'   and one step after it.
#' @return data.frame of dvds_bounds() rows with columns rep, study and fit.
simulate_replicate <- function(rep, cfg, progress, num_threads) {
  studies <- lapply(cfg$n, simulate_dgk)
  rows <- list()
  for (k in seq_along(studies)) {
    folds <- split_folds(cfg$n[k], cfg$frac_a)
    samples <- list(full = studies[[k]],
                    fold_a = subset_study(studies[[k]], folds$a),
                    fold_b = subset_study(studies[[k]], folds$b))
    for (fit in FITS) {
      label <- sprintf("rep %d/%d: %s, study %d", rep, cfg$reps,
                       FIT_LABELS[[fit]], k)
      progress(message = label, amount = 0)
      est <- dvds_bounds(samples[[fit]], cfg$lambdas, num_threads)
      rows[[length(rows) + 1]] <- data.frame(rep = rep, study = k,
                                             fit = fit, est)
      progress(message = label, amount = 1)
    }
  }
  do.call(rbind, rows)
}


#' Sets the future backend. Sequential runs give grf every thread; parallel
#' runs give each worker one. grf output does not depend on its thread count
#' (checked for grf 2.6.1), so results are identical across backends.
#'
#' @return Number of grf threads per fit (NULL = all).
set_backend <- function(cfg) {
  if (cfg$workers == 1 || cfg$backend == "sequential") {
    future::plan(future::sequential)
    return(NULL)
  }
  strategy <- switch(cfg$backend, multisession = future::multisession,
                     multicore = future::multicore)
  future::plan(strategy, workers = cfg$workers)
  1
}


#' Runs every replicate in parallel behind one progress bar.
#'
#' @return data.frame with one row per (rep, study, fit, outcome, Lambda).
run_simulation <- function(cfg) {
  num_threads <- set_backend(cfg)
  on.exit(future::plan(future::sequential), add = TRUE)
  seeds <- replicate_seeds(cfg$seed, cfg$reps)
  # The bar is drawn by the progress package, which progressr only suggests;
  # naming it here also records it in renv.lock.
  if (!requireNamespace("progress", quietly = TRUE)) {
    stop("package 'progress' is required; run renv::restore()")
  }
  # progressr disables handlers outside interactive sessions unless enabled.
  bar <- progressr::handler_progress(
    format = "[:bar] :percent ETA :eta | :message",
    clear = FALSE, enable = TRUE)
  results <- progressr::with_progress(handlers = bar, enable = TRUE, {
    progress <- progressr::progressor(
      steps = cfg$reps * length(cfg$n) * length(FITS))
    future.apply::future_lapply(
      seq_len(cfg$reps), simulate_replicate, cfg = cfg,
      progress = progress, num_threads = num_threads,
      future.seed = seeds, future.chunk.size = 1, future.packages = "grf")
  })
  do.call(rbind, results)
}


# ----------------------------- confidence intervals ------------------------
#' Rearranges the estimates of one outcome into arrays
#' [replicate, lambda_index, study].
#'
#' @param estimates Rows of estimates.csv for a single outcome.
#' @return Nested list: arrays[[fit]][[column]].
bound_arrays <- function(estimates) {
  reps <- sort(unique(estimates$rep))
  n_lambdas <- max(estimates$lambda_index)
  columns <- c("lower", "upper", "s2_lower", "s2_upper", "n")
  arrays <- list()
  for (fit in FITS) {
    sub <- estimates[estimates$fit == fit, ]
    idx <- cbind(match(sub$rep, reps), sub$lambda_index, sub$study)
    arrays[[fit]] <- lapply(setNames(columns, columns), function(column) {
      out <- array(NA_real_, c(length(reps), n_lambdas, 2))
      out[idx] <- sub[[column]]
      stopifnot(!anyNA(out))
      out
    })
  }
  arrays
}


#' CIs of every method for study 1 at Lambda index i, study 2 at index j.
#'
#' @return Named list of CI lists (see ci_methods.R), one per method.
pair_cis <- function(arrays, i, j, z, z_calibrated) {
  # R x 2 matrices: column 1 is study 1 at Lambda_i, column 2 study 2 at
  # Lambda_j.
  pick <- function(fit) {
    lapply(arrays[[fit]], function(a) cbind(a[, i, 1], a[, j, 2]))
  }
  full <- pick("full")
  a <- pick("fold_a")
  b <- pick("fold_b")
  # Weights on fold A: fold-A means, variances rescaled to fold-B size, i.e.
  # var(phi_A) / n_B (derivations.tex, eq. est-oracleweight).
  scale <- a$n / b$n
  w_L <- optimal_weight_grid(a$lower[, 1], a$lower[, 2],
                             a$s2_lower[, 1] * scale[, 1],
                             a$s2_lower[, 2] * scale[, 2], lower = TRUE, z)
  w_U <- optimal_weight_grid(a$upper[, 1], a$upper[, 2],
                             a$s2_upper[, 1] * scale[, 1],
                             a$s2_upper[, 2] * scale[, 2], lower = FALSE, z)
  cis <- list(
    convex_ci(b$lower, b$upper, b$s2_lower, b$s2_upper, w_L, w_U, z),
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 1, 1, z),
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 0, 0, z),
    intersection_ci(full$lower, full$upper, full$s2_lower, full$s2_upper,
                    z_calibrated)
  )
  setNames(cis, c(M3, M5_D1, M5_D2, M6))
}


#' Summary rows for every (outcome, Lambda_1, Lambda_2, method).
#'
#' The target is the fused population band [max_k L_k, min_k U_k]. The
#' bands are nested, so it equals the band of min(Lambda_1, Lambda_2). The
#' gaps g_L = |L(Lambda_1) - L(Lambda_2)| and g_U = |U(Lambda_1) -
#' U(Lambda_2)| are equal for the continuous outcome (symmetric band) but
#' not for the binary one.
#'
#' @param truth Output of true_bands() for the Lambda grid of the estimates.
summarise_estimates <- function(estimates, truth, cfg) {
  rows <- list()
  for (outcome in unique(estimates$outcome)) {
    arrays <- bound_arrays(estimates[estimates$outcome == outcome, ])
    band <- truth[truth$outcome == outcome, ]
    for (i in band$lambda_index) {
      for (j in band$lambda_index) {
        true_lo <- max(band$lower[i], band$lower[j])
        true_hi <- min(band$upper[i], band$upper[j])
        cis <- pair_cis(arrays, i, j, cfg$z, cfg$z_calibrated)
        for (method in names(cis)) {
          rows[[length(rows) + 1]] <- data.frame(
            outcome = outcome,
            lambda_1 = band$lambda[i], lambda_2 = band$lambda[j],
            gap_L = abs(band$lower[i] - band$lower[j]),
            gap_U = abs(band$upper[i] - band$upper[j]), method = method,
            true_lo = true_lo, true_hi = true_hi,
            summarise_ci(cis[[method]], true_lo, true_hi, THETA))
        }
      }
    }
  }
  do.call(rbind, rows)
}


#' Coverage of each study's own band by its own CI, per outcome and Lambda
#' (DGK's criterion). With --n1 1000, study 1 replicates DGK's design.
single_study_coverage <- function(estimates, truth, z) {
  full <- merge(estimates[estimates$fit == "full", ],
                truth[c("outcome", "lambda_index", "lower", "upper")],
                by = c("outcome", "lambda_index"),
                suffixes = c("", "_true"))
  covered <- (full$lower - z * sqrt(full$s2_lower) <= full$lower_true) &
    (full$upper + z * sqrt(full$s2_upper) >= full$upper_true)
  out <- aggregate(list(coverage = covered),
                   list(outcome = full$outcome, study = full$study,
                        n = full$n, lambda = full$lambda),
                   mean)
  out$coverage_mcse <- sqrt(out$coverage * (1 - out$coverage) /
                              length(unique(full$rep)))
  out[order(out$outcome, out$study, out$lambda), ]
}


# ----------------------------- sanity checks ------------------------------
#' At Lambda = 1 both DVDS estimators reduce to AIPW, so ATE^+ = ATE^- in
#' every fit (continuous: kappa = mu, pseudo = Y; binary: tau = 1/2,
#' kappa = mu).
check_collapse_at_one <- function(estimates) {
  at_one <- estimates[estimates$lambda == 1, ]
  worst <- max(abs(at_one$upper - at_one$lower))
  list(name = "Band collapses to a point at Lambda = 1 (every fit)",
       ok = worst <= 1e-8, detail = sprintf("max |U - L| = %.2e", worst))
}


#' Mean bounds widen with Lambda, for every (outcome, study, fit).
check_monotone <- function(estimates) {
  groups <- split(estimates,
                  list(estimates$outcome, estimates$study, estimates$fit))
  ok <- vapply(groups, function(g) {
    lo <- tapply(g$lower, g$lambda_index, mean)
    hi <- tapply(g$upper, g$lambda_index, mean)
    all(diff(lo) <= 0) && all(diff(hi) >= 0)
  }, logical(1))
  list(name = "Mean L decreasing and mean U increasing in Lambda",
       ok = all(ok),
       detail = sprintf("%d of %d (outcome, study, fit) groups monotone",
                        sum(ok),
                        length(ok)))
}


#' M6 is valid whichever study attains each end: joint coverage
#' >= 1 - alpha up to 3 MC s.e. in every pair.
check_calibrated_valid <- function(summary, alpha) {
  m6 <- summary[summary$method == M6, ]
  # Cells with coverage 1 have MC s.e. 0 and give z = +Inf, which passes.
  zs <- (m6$cov_joint - (1 - alpha)) / m6$cov_joint_mcse
  worst <- tapply(m6$cov_joint, m6$outcome, min)
  list(name = sprintf("M6 joint coverage >= %g in every pair", 1 - alpha),
       ok = min(zs) >= -3,
       detail = sprintf("min coverage %s; min z = %.2f over %d pairs",
                        paste(sprintf("%s %.4f", names(worst), worst),
                              collapse = ", "),
                        min(zs), nrow(m6)))
}


#' Prints PASS/FAIL for each check; returns TRUE if all pass.
run_sanity_checks <- function(estimates, summary, alpha) {
  checks <- list(check_collapse_at_one(estimates), check_monotone(estimates),
                 check_calibrated_valid(summary, alpha))
  cat("\nSanity checks\n")
  for (check in checks) {
    cat(sprintf("[%s] %s: %s\n", if (check$ok) "PASS" else "FAIL",
                check$name, check$detail))
  }
  n_passed <- sum(vapply(checks, `[[`, logical(1), "ok"))
  cat(sprintf("%d/%d checks passed\n", n_passed, length(checks)))
  n_passed == length(checks)
}


# ----------------------------- output ------------------------------
#' Prints M3 and M6 for a few pairs on the diagonal and at Lambda_1 = max.
print_highlights <- function(summary) {
  lambdas <- sort(unique(summary$lambda_1))
  show <- lambdas[round(seq(1, length(lambdas), length.out = 5))]
  keep <- summary$lambda_1 %in% show & summary$lambda_2 %in% show &
    (summary$lambda_1 == summary$lambda_2 |
       summary$lambda_1 == max(lambdas)) &
    summary$method %in% c(M3, M6)
  columns <- c("outcome", "lambda_1", "lambda_2", "method", "cov_joint",
               "cov_lo", "cov_hi", "width", "wL", "wU")
  old <- options(width = 150)
  on.exit(options(old), add = TRUE)
  cat("\nDiagonal and Lambda_1 = max (full results in summary.csv)\n")
  print(format(summary[keep, columns], digits = 3), row.names = FALSE)
}


#' Short git hash of HEAD, with "-dirty" if R/ or renv.lock has changes.
git_hash <- function() {
  hash <- tryCatch(
    system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE,
            stderr = FALSE),
    error = function(e) "unknown")
  dirty <- tryCatch(
    system2("git", c("status", "--porcelain", "--", "R", "renv.lock"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) character(0))
  paste0(hash, if (length(dirty) > 0) "-dirty" else "")
}


#' Records the command, settings and software versions of a run.
write_run_info <- function(cfg, start, path) {
  packages <- c("grf", "dplyr", "future", "future.apply", "progressr",
                "optparse")
  versions <- vapply(packages, function(p) {
    as.character(utils::packageVersion(p))
  }, character(1))
  lines <- c(
    sprintf("date: %s", format(start, "%Y-%m-%d %H:%M")),
    sprintf("git: %s", git_hash()),
    sprintf("DGK commit: %s (estimator and DGP code, unmodified)", DGK_COMMIT),
    sprintf("command: Rscript R/causal_sim.R %s",
            paste(commandArgs(trailingOnly = TRUE), collapse = " ")),
    sprintf("seed: %d | reps: %d | n: %s | frac_a: %g | alpha: %g",
            cfg$seed, cfg$reps, paste(cfg$n, collapse = ", "), cfg$frac_a,
            cfg$alpha),
    sprintf("outcomes: %s", paste(OUTCOMES, collapse = ", ")),
    sprintf("lambdas: %s", paste(cfg$lambdas, collapse = ", ")),
    sprintf("backend: %s | workers: %d", cfg$backend, cfg$workers),
    sprintf("simulation run time: %s", format_elapsed(start)),
    R.version.string,
    paste(sprintf("%s %s", packages, versions), collapse = " | ")
  )
  writeLines(lines, path)
}


# ----------------------------- main ------------------------------
main <- function() {
  start <- Sys.time()
  args <- parse_args()
  cfg <- make_config(args)
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  estimates_path <- file.path(cfg$out_dir, "estimates.csv")
  if (args$summarise_only) {
    estimates <- read.csv(estimates_path)
    # Runs before 2026-11-02 estimated the continuous outcome only.
    if (is.null(estimates$outcome)) estimates$outcome <- "continuous"
  } else {
    check_clone()
    estimates <- run_simulation(cfg)
    write.csv(estimates, estimates_path, row.names = FALSE)
    write_run_info(cfg, start, file.path(cfg$out_dir, "run_info.txt"))
  }
  # True bands, once per (outcome, Lambda); binary by Monte Carlo.
  lambdas <- as.numeric(tapply(estimates$lambda, estimates$lambda_index,
                               `[`, 1))
  truth <- true_bands(lambdas)
  truth <- truth[truth$outcome %in% estimates$outcome, ]
  write.csv(truth, file.path(cfg$out_dir, "true_bands.csv"), row.names = FALSE)
  summary <- summarise_estimates(estimates, truth, cfg)
  write.csv(summary, file.path(cfg$out_dir, "summary.csv"), row.names = FALSE)
  write.csv(single_study_coverage(estimates, truth, cfg$z),
            file.path(cfg$out_dir, "single_study_coverage.csv"),
            row.names = FALSE)
  run_sanity_checks(estimates, summary, cfg$alpha)
  print_highlights(summary)
  cat(sprintf("\nWritten to %s\n", cfg$out_dir))
  cat(sprintf("Total run time: %s\n", format_elapsed(start)))
}


if (sys.nframe() == 0L) main()
