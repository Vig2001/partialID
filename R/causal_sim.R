# Two observational studies (OS) from the DGK DGP, differing in sample size
# n_k, sensitivity parameter grid Lambda_k and (optionally) the scale c_2 of
# study 2's continuous outcome noise and the rarity k_1 of study 1's binary
# outcome, P(Y = 1 | x) = 1 / (k_1 + exp(eta(x))) (DGK's is k = 1). Each
# study has a binary and a continuous outcome (DGK's two designs, drawn from
# the same X and Z; --outcomes can run one of them only). For
# each outcome, each study gives DVDS estimates of its sharp ATE bounds
# [L_k, U_k], computed with DGK's own code (R/helpers/dgk.R; their repository
# cloned in external/). CIs for the fused band [max_k L_k, min_k U_k] are
# compared over all pairs (Lambda_1, Lambda_2), separately per outcome:
#   dataset 1, dataset 2        each study alone (full-sample fits; DGK's CI)
#   intersection (calibrated)   full-sample fits, z = Phi^{-1}(1 - alpha/4)
#   convex (split)              weights from fold A, CI from fold B (DVDS
#                               fitted separately in each fold)
#
# Usage, from the repo root:
#   Rscript R/causal_sim.R [--seed 1] [--reps 1000] [--n1 2000] [--n2 500]
#     [--lambdas-1 1,1.5,... | 3:5:0.05] [--lambdas-2 ...] [--noise-scale-2 1]
#     [--binary-denominator-1 1] [--outcomes binary,continuous]
#     [--frac-a 0.3] [--alpha 0.05] [--workers W] [--backend multisession]
#     [--out-dir src/outputs/causal_sim] [--summarise-only]
# The defaults reproduce the baseline run (both studies on DGK's grid,
# c = 1, k = 1, both outcomes).
#
# Plans:   quality_reports/plans/done/2026-11-02_dvds-r-two-os.md,
#          2026-11-02_dgk-exact-replication.md, 2026-11-02_causal-sim-binary.md,
#          2026-11-03_convex-advantage-noise-design.md,
#          2026-11-05_binary-rare-outcome-k15.md
# Notes:   partialID_Overleaf/Notes/Simulation_Notes/causal_sim_code_explained.tex
# Outputs: src/outputs/causal_sim/ (summary.csv: coverage and widths only;
#          summary_full.csv: every summary column)

# Paths are relative to the repo root, where .Rprofile also activates renv.
if (!file.exists(file.path("R", "causal_sim.R"))) {
  stop("run from the repo root: Rscript R/causal_sim.R")
}
for (helper in c("data_fusion_design.R", "dgk.R", "ci_methods.R",
                 "mc_metrics.R", "run_time.R")) {
  source(file.path("R", "helpers", helper))
}

# Default Lambda grid for both studies: DGK's (dvds_sim.R line 10). The
# number of values hardly affects the fitting time (forests are fitted once
# per fit), but the summary stage grows with the number of pairs.
LAMBDAS <- seq(1, 3, length.out = 20)
FITS <- c("full", "fold_a", "fold_b")
OUTCOMES <- c("binary", "continuous")
# Progress-bar labels; short, so the bar fits an 80-column console.
FIT_LABELS <- c(full = "intersection", fold_a = "convex split A",
                fold_b = "convex split B")
# Methods, in the order they are reported.
M5_D1 <- "dataset 1"
M5_D2 <- "dataset 2"
M6 <- "intersection (calibrated)"
M3 <- "convex (split)"
METHODS <- c(M5_D1, M5_D2, M6, M3)
# Columns of the short summary.csv; every column is in summary_full.csv.
SHORT_COLUMNS <- c("outcome", "lambda_1", "lambda_2", "method", "cov_joint",
                   "ci_width", "est_band_width", "true_band_width")
SHORT_DIGITS <- 3
MIN_FOLD_SIZE <- 100  # Smallest fold; DVDS fits forests within each arm.
THETA <- 0  # True ATE: neither outcome depends on Z in the DGK DGP.


# ----------------------------- configuration ------------------------------
#' @param extra_options optparse options appended by scripts that reuse this
#'   one (R/causal_sim_quick.R).
#' @param default_out_dir Default of --out-dir.
parse_args <- function(argv = commandArgs(trailingOnly = TRUE),
                       extra_options = list(),
                       default_out_dir = file.path("src", "outputs",
                                                   "causal_sim")) {
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
    opt("--lambdas-1", dest = "lambdas_1", type = "character",
        default = NA_character_,
        help = paste("Lambda grid of study 1: comma-separated values or",
                     "from:to:by [default: DGK's seq(1, 3, length.out = 20)]."
                     )),
    opt("--lambdas-2", dest = "lambdas_2", type = "character",
        default = NA_character_,
        help = "Lambda grid of study 2, as --lambdas-1 [default: DGK's]."),
    opt("--noise-scale-2", dest = "noise_scale_2", type = "double",
        default = 1,
        help = paste("Scale c of study 2's continuous outcome noise; its true",
                     "band scales by c [default %default].")),
    opt("--binary-denominator-1", dest = "binary_denominator_1",
        type = "double", default = 1,
        help = paste("k >= 1 for study 1's binary outcome, P(Y = 1 | x) =",
                     "1 / (k + exp(eta(x))); study 2 keeps DGK's k = 1",
                     "[default %default].")),
    opt("--outcomes", type = "character", default = "binary,continuous",
        help = paste("Comma-separated outcomes to estimate; 'binary' alone",
                     "skips the slow continuous DVDS [default %default].")),
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
        default = default_out_dir,
        help = "Output directory [default %default]."),
    opt("--summarise-only", dest = "summarise_only", action = "store_true",
        default = FALSE,
        help = "Recompute CIs and summaries from estimates.csv; no fitting.")
  )
  parser <- optparse::OptionParser(option_list = c(options, extra_options))
  optparse::parse_args(parser, args = argv)
}


#' A Lambda grid from a flag: NA (DGK's grid), "a,b,c" or "from:to:by".
#'
#' from:to:by values are rounded to 10 decimals, so 3:5:0.05 gives exactly
#' 3, 3.05, ..., 5.
parse_lambdas <- function(text) {
  if (is.na(text)) return(LAMBDAS)
  if (grepl(":", text, fixed = TRUE)) {
    parts <- as.numeric(strsplit(text, ":", fixed = TRUE)[[1]])
    if (length(parts) != 3 || anyNA(parts) || parts[3] <= 0) {
      stop(sprintf("a grid from:to:by needs three numbers, by > 0; got %s",
                   text))
    }
    lambdas <- round(seq(parts[1], parts[2] + 1e-9, by = parts[3]), 10)
  } else {
    lambdas <- as.numeric(strsplit(text, ",", fixed = TRUE)[[1]])
  }
  validate_lambdas(lambdas)
  if (is.unsorted(lambdas, strictly = TRUE)) {
    stop(sprintf("lambdas must be strictly increasing, got %s", text))
  }
  lambdas
}


#' Validated settings, with the critical values implied by alpha.
make_config <- function(args) {
  n <- c(args$n1, args$n2)
  n_a <- floor(args$frac_a * n)
  lambdas <- list(parse_lambdas(args$lambdas_1), parse_lambdas(args$lambdas_2))
  if (!(is.finite(args$noise_scale_2) && args$noise_scale_2 > 0)) {
    stop(sprintf("noise-scale-2 must be positive, got %s",
                 args$noise_scale_2))
  }
  validate_binary_denominator(args$binary_denominator_1)
  outcomes <- strsplit(args$outcomes, ",", fixed = TRUE)[[1]]
  if (length(outcomes) == 0 || !all(outcomes %in% OUTCOMES) ||
      anyDuplicated(outcomes)) {
    stop(sprintf("outcomes must be distinct values from %s, got %s",
                 paste(OUTCOMES, collapse = ", "), args$outcomes))
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
       noise_scale = c(1, args$noise_scale_2),
       binary_denominator = c(args$binary_denominator_1, 1),
       outcomes = OUTCOMES[OUTCOMES %in% outcomes],
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
#' each fit giving both outcomes on that study's Lambda grid.
#'
#' @param rep Replicate number.
#' @param progress progressr progressor; receives a label before each fit
#'   and one step after it.
#' @return data.frame of dvds_bounds() rows with columns rep, study and fit.
simulate_replicate <- function(rep, cfg, progress, num_threads) {
  studies <- Map(simulate_study, cfg$n, cfg$noise_scale,
                 cfg$binary_denominator)
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
      est <- dvds_bounds(samples[[fit]], cfg$lambdas[[k]], num_threads,
                         cfg$outcomes)
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
#' Rearranges the estimates of one outcome into one matrix
#' [replicate, lambda_index] per study (the studies' grids may differ).
#'
#' @param estimates Rows of estimates.csv for a single outcome.
#' @param fits Names of the fits to arrange.
#' @return Nested list: arrays[[fit]][[column]][[study]].
bound_arrays <- function(estimates, fits = FITS) {
  reps <- sort(unique(estimates$rep))
  columns <- c("lower", "upper", "s2_lower", "s2_upper", "n")
  arrays <- list()
  for (fit in fits) {
    arrays[[fit]] <- lapply(setNames(columns, columns), function(column) {
      lapply(1:2, function(k) {
        sub <- estimates[estimates$fit == fit & estimates$study == k, ]
        out <- matrix(NA_real_, length(reps), max(sub$lambda_index))
        out[cbind(match(sub$rep, reps), sub$lambda_index)] <- sub[[column]]
        stopifnot(!anyNA(out))
        out
      })
    })
  }
  arrays
}


#' CIs of every method for study 1 at Lambda index i, study 2 at index j.
#'
#' @return Named list of CI lists (see ci_methods.R), one per method, in the
#'   order of METHODS.
pair_cis <- function(arrays, i, j, z, z_calibrated) {
  # R x 2 matrices: column 1 is study 1 at Lambda_i, column 2 study 2 at
  # Lambda_j.
  pick <- function(fit) {
    lapply(arrays[[fit]], function(a) cbind(a[[1]][, i], a[[2]][, j]))
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
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 1, 1, z),
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 0, 0, z),
    intersection_ci(full$lower, full$upper, full$s2_lower, full$s2_upper,
                    z_calibrated),
    convex_ci(b$lower, b$upper, b$s2_lower, b$s2_upper, w_L, w_U, z)
  )
  setNames(cis, METHODS)
}


#' True bands of both studies, one row per (study, outcome, Lambda).
#'
#' Each study has its own Lambda grid; study 2's continuous band is scaled
#' by its noise scale c (continuous_true_band()), and each study's binary
#' band uses its own denominator k (binary_true_band()).
study_true_bands <- function(cfg) {
  do.call(rbind, lapply(1:2, function(k) {
    data.frame(study = k,
               true_bands(cfg$lambdas[[k]], cfg$noise_scale[k],
                          cfg$binary_denominator[k]))
  }))
}


#' Summary rows for every (outcome, Lambda_1, Lambda_2, method).
#'
#' The target is the fused population band [max_k L_k, min_k U_k] of study
#' 1 at Lambda_1 and study 2 at Lambda_2. The gaps g_L = |L_1 - L_2| and
#' g_U = |U_1 - U_2| are equal for the continuous outcome (symmetric bands)
#' but not for the binary one.
#'
#' @param truth Output of study_true_bands() for the grids of the estimates.
#' @return data.frame with every summary column (summary_full.csv).
summarise_estimates <- function(estimates, truth, cfg) {
  rows <- list()
  for (outcome in unique(estimates$outcome)) {
    arrays <- bound_arrays(estimates[estimates$outcome == outcome, ])
    band_1 <- truth[truth$outcome == outcome & truth$study == 1, ]
    band_2 <- truth[truth$outcome == outcome & truth$study == 2, ]
    for (i in band_1$lambda_index) {
      for (j in band_2$lambda_index) {
        true_lo <- max(band_1$lower[i], band_2$lower[j])
        true_hi <- min(band_1$upper[i], band_2$upper[j])
        cis <- pair_cis(arrays, i, j, cfg$z, cfg$z_calibrated)
        for (method in names(cis)) {
          rows[[length(rows) + 1]] <- data.frame(
            outcome = outcome,
            lambda_1 = band_1$lambda[i], lambda_2 = band_2$lambda[j],
            gap_L = abs(band_1$lower[i] - band_2$lower[j]),
            gap_U = abs(band_1$upper[i] - band_2$upper[j]), method = method,
            true_lo = true_lo, true_hi = true_hi,
            true_band_width = true_hi - true_lo,
            summarise_ci(cis[[method]], true_lo, true_hi, THETA))
        }
      }
    }
  }
  do.call(rbind, rows)
}


#' The short summary: joint coverage, CI width, and estimated against true
#' band width, rounded for reading.
short_summary <- function(summary) {
  out <- summary
  out$ci_width <- out$width
  out$est_band_width <- out$band_width
  out <- out[SHORT_COLUMNS]
  numeric <- vapply(out, is.numeric, logical(1))
  out[numeric] <- lapply(out[numeric], round, SHORT_DIGITS)
  out
}


#' Coverage of each study's own band by its own CI, per outcome and Lambda
#' (DGK's criterion). With --n1 1000, study 1 replicates DGK's design.
single_study_coverage <- function(estimates, truth, z) {
  full <- merge(estimates[estimates$fit == "full", ],
                truth[c("study", "outcome", "lambda_index", "lower", "upper")],
                by = c("study", "outcome", "lambda_index"),
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


#' Stops unless each study's estimated Lambdas match its true-band grid.
#'
#' Guards --summarise-only runs whose flags differ from the original run.
check_grids_match <- function(estimates, truth) {
  for (k in 1:2) {
    est <- unique(estimates[estimates$study == k, c("lambda_index", "lambda")])
    tru <- unique(truth[truth$study == k, c("lambda_index", "lambda")])
    both <- merge(est, tru, by = "lambda_index", suffixes = c("_est", "_true"))
    if (nrow(both) != nrow(est) || nrow(est) != nrow(tru) ||
        max(abs(both$lambda_est - both$lambda_true)) > 1e-8) {
      stop(sprintf(paste("study %d: the Lambda grid of the estimates does",
                         "not match the flags; pass the run's --lambdas-%d"),
                   k, k))
    }
  }
}


# ----------------------------- sanity checks ------------------------------
#' At Lambda = 1 both DVDS estimators reduce to AIPW, so ATE^+ = ATE^- in
#' every fit (continuous: kappa = mu, pseudo = Y; binary: tau = 1/2,
#' kappa = mu).
check_collapse_at_one <- function(estimates) {
  at_one <- estimates[estimates$lambda == 1, ]
  if (nrow(at_one) == 0) {
    return(list(name = "Band collapses to a point at Lambda = 1",
                ok = TRUE, detail = "skipped: no grid contains Lambda = 1"))
  }
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


#' The calibrated intersection is valid whichever study attains each end:
#' joint coverage >= 1 - alpha up to 3 MC s.e. in every pair.
check_calibrated_valid <- function(summary, alpha) {
  m6 <- summary[summary$method == M6, ]
  # Cells with coverage 1 have MC s.e. 0 and give z = +Inf, which passes.
  zs <- (m6$cov_joint - (1 - alpha)) / m6$cov_joint_mcse
  worst <- tapply(m6$cov_joint, m6$outcome, min)
  list(name = sprintf("Calibrated intersection joint coverage >= %g, all pairs",
                      1 - alpha),
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
#' Prints the short summary for five evenly spaced Lambda_1 values, each
#' paired with the Lambda_2 at which the two true bands are closest.
#'
#' For fixed Lambda_1 the fused true band widens with Lambda_2 until study
#' 2's band contains study 1's, then stays flat; the first Lambda_2 on the
#' flat part is where the two bands are closest.
print_highlights <- function(short) {
  lambdas <- sort(unique(short$lambda_1))
  show <- lambdas[round(seq(1, length(lambdas), length.out = 5))]
  rows <- list()
  for (outcome in unique(short$outcome)) {
    for (l1 in show) {
      sub <- short[short$outcome == outcome & short$lambda_1 == l1, ]
      one_method <- sub[sub$method == M5_D1, ]
      l2 <- one_method$lambda_2[which.max(one_method$true_band_width)]
      rows[[length(rows) + 1]] <- sub[sub$lambda_2 == l2, ]
    }
  }
  old <- options(width = 150)
  on.exit(options(old), add = TRUE)
  cat("\nSelected pairs (all pairs in summary.csv)\n")
  print(do.call(rbind, rows), row.names = FALSE)
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
    sprintf("outcomes: %s", paste(cfg$outcomes, collapse = ", ")),
    sprintf("lambdas study 1: %s", paste(cfg$lambdas[[1]], collapse = ", ")),
    sprintf("lambdas study 2: %s", paste(cfg$lambdas[[2]], collapse = ", ")),
    sprintf("noise scale (continuous outcome): study 1 %g, study 2 %g",
            cfg$noise_scale[1], cfg$noise_scale[2]),
    sprintf("binary denominator k: study 1 %g, study 2 %g",
            cfg$binary_denominator[1], cfg$binary_denominator[2]),
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
  # True bands, once per (study, outcome, Lambda); binary by Monte Carlo.
  # With --summarise-only the flags must match the run (see run_info.txt).
  truth <- study_true_bands(cfg)
  truth <- truth[truth$outcome %in% estimates$outcome, ]
  check_grids_match(estimates, truth)
  write.csv(truth, file.path(cfg$out_dir, "true_bands.csv"), row.names = FALSE)
  summary <- summarise_estimates(estimates, truth, cfg)
  short <- short_summary(summary)
  write.csv(summary, file.path(cfg$out_dir, "summary_full.csv"),
            row.names = FALSE)
  write.csv(short, file.path(cfg$out_dir, "summary.csv"), row.names = FALSE)
  write.csv(single_study_coverage(estimates, truth, cfg$z),
            file.path(cfg$out_dir, "single_study_coverage.csv"),
            row.names = FALSE)
  run_sanity_checks(estimates, summary, cfg$alpha)
  print_highlights(short)
  cat(sprintf("\nWritten to %s\n", cfg$out_dir))
  cat(sprintf("Total run time: %s\n", format_elapsed(start)))
}


if (sys.nframe() == 0L) main()
