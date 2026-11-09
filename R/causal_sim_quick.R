# Lean version of R/causal_sim.R: the same two-study DGK design and CI
# comparison, with parametric nuisances (R/helpers/dvds_lean.R) in place of
# DGK's forests, so that every valid method of the toy, and cross-fitting,
# runs in minutes. Per outcome and pair (Lambda_1, Lambda_2):
#   dataset 1, dataset 2        each study alone (full-sample fits)
#   intersection (calibrated)   full-sample fits, z = Phi^{-1}(1 - alpha/4)
#   convex (split)              weights from fold A, CI from fold B
#   convex (cross-fit)          K folds per study, standard cross-fitting:
#                               for each fold f, nuisance models fitted on
#                               the other K - 1 folds give the weights (from
#                               the complement's own bounds, variances
#                               scaled to the full n) and predict fold f's
#                               nuisances, for fold f's bounds; the fold
#                               bounds are averaged
#   intersection (CLR)          CLR intersection bounds, every bound kept
#   convex (CLR)                CLR over all convex combinations, no split
# M8 and M9 use level p = 1 - alpha/2 at each end, as in src/toy_sim.py.
# Fits per study and replicate: full, fold A and fold B of both studies
# first, then for each fold f of the K one fit on its complement, which also
# predicts for fold f (3 + K in all). Drawing the cross-fitting folds last
# keeps the other methods' fits unchanged when K changes.
#
# --nuisances forest (binary outcome only) uses DGK's forests instead, as
# R/causal_sim.R does. Full, fold A and fold B are fitted in R/causal_sim.R's
# order, so they equal that script's fits for the same seed and settings.
# For M7 the forests fitted on fold f's complement give the complement's
# bounds (out of bag) and predict fold f's nuisances.
#
# Usage, from the repo root (flags as R/causal_sim.R, plus those below):
#   Rscript R/causal_sim_quick.R [--folds 5] [--reps 1000] [--n1 3000]
#     [--nuisance-folds 5] [--n2 500] [--noise-scale-2 2.5]
#     [--binary-denominator-1 15] [--nuisances lean|forest]
#     [--outcomes binary|continuous] [--lambdas-1 ...] [--lambdas-2 ...]
#     [--workers W] [--backend multicore] [--out-dir ...] [--summarise-only]
#
# Plans:   quality_reports/plans/done/2026-11-07_causal-sim-quick.md,
#          2026-11-08_causal-sim-forest-binary-crossfit.md (forest mode)
# Notes:   partialID_Overleaf/Notes/Simulation_Notes/causal_sim_quick_explained.tex
# Outputs: src/outputs/causal_sim_quick/ by default (lean), or
#          src/outputs/causal_sim_forest/ (forest)
# TODO: after the CLeaR deadline, move the parts shared with causal_sim.R
#   (summary loop, run info) into R/helpers/ instead of overriding them here.

if (!file.exists(file.path("R", "causal_sim_quick.R"))) {
  stop("run from the repo root: Rscript R/causal_sim_quick.R")
}
source(file.path("R", "causal_sim.R"))  # Its main() runs only when called.
source(file.path("R", "helpers", "dvds_lean.R"))

M7 <- "convex (cross-fit)"
M8 <- "intersection (CLR)"
M9 <- "convex (CLR)"
METHODS_QUICK <- c(M5_D1, M5_D2, M6, M3, M7, M8, M9)
DEFAULT_FOLDS <- 5
DEFAULT_OUT_DIR <- file.path("src", "outputs", "causal_sim_quick")
DEFAULT_OUT_DIR_FOREST <- file.path("src", "outputs", "causal_sim_forest")
NUISANCES <- c("lean", "forest")
SINGLE_STUDY_SES <- 3  # Single-study coverage check: MC s.e. below nominal.


# ----------------------------- configuration ------------------------------
#' Names of the fits per study: full, the split folds, then each
#' cross-fitting fold held out and its complement.
quick_fits <- function(n_folds) {
  c("full", "fold_a", "fold_b", sprintf("cf_held_%d", seq_len(n_folds)),
    sprintf("cf_train_%d", seq_len(n_folds)))
}


#' Validated --nuisances, and the output directory.
#'
#' A forest run left at the lean default directory is moved to
#' DEFAULT_OUT_DIR_FOREST, so it never overwrites a lean run.
#'
#' @return List with nuisances and out_dir.
nuisance_settings <- function(args) {
  if (!(args$nuisances %in% NUISANCES)) {
    stop(sprintf("nuisances must be one of %s, got %s",
                 paste(NUISANCES, collapse = ", "), args$nuisances))
  }
  out_dir <- args$out_dir
  if (args$nuisances == "forest" && identical(out_dir, DEFAULT_OUT_DIR)) {
    out_dir <- DEFAULT_OUT_DIR_FOREST
  }
  list(nuisances = args$nuisances, out_dir = out_dir)
}


#' Settings of R/causal_sim.R plus the number of cross-fitting folds and the
#' nuisance models.
quick_config <- function(argv = commandArgs(trailingOnly = TRUE)) {
  extra <- list(
    optparse::make_option(
      "--folds", type = "integer", default = DEFAULT_FOLDS,
      help = "Cross-fitting folds K per study (M7) [default %default]."),
    optparse::make_option(
      "--nuisance-folds", dest = "nuisance_folds", type = "integer",
      default = NUISANCE_FOLDS,
      help = paste("Lean nuisances only: folds for out-of-fold nuisance",
                   "predictions within each fit; 1 fits in-sample",
                   "[default %default].")),
    optparse::make_option(
      "--nuisances", type = "character", default = "lean",
      help = paste("lean (parametric, R/helpers/dvds_lean.R) or forest",
                   "(DGK's forests; needs --outcomes binary)",
                   "[default %default].")))
  args <- parse_args(argv, extra_options = extra,
                     default_out_dir = DEFAULT_OUT_DIR)
  nuisance <- nuisance_settings(args)
  args$out_dir <- nuisance$out_dir
  cfg <- make_config(args)
  if (nuisance$nuisances == "forest" && !identical(cfg$outcomes, "binary")) {
    stop(sprintf("forest nuisances need --outcomes binary, got %s",
                 paste(cfg$outcomes, collapse = ", ")))
  }
  cfg$nuisances <- nuisance$nuisances
  if (args$folds < 2) stop(sprintf("folds must be >= 2, got %d", args$folds))
  smallest <- min(floor(cfg$n / args$folds))
  if (smallest < MIN_FOLD_SIZE) {
    stop(sprintf("every cross-fitting fold needs >= %d units, got %d",
                 MIN_FOLD_SIZE, smallest))
  }
  if (args$nuisance_folds < 1) {
    stop(sprintf("nuisance-folds must be >= 1, got %d", args$nuisance_folds))
  }
  cfg$folds <- args$folds
  cfg$nuisance_folds <- args$nuisance_folds
  cfg$fits <- quick_fits(args$folds)
  cfg$p <- 1 - cfg$alpha / 2  # Level of each one-sided CLR end.
  cfg$summarise_only <- args$summarise_only
  cfg
}


# ----------------------------- simulation ------------------------------
#' Random partition of 1, ..., n into K folds of sizes differing by at most 1.
#'
#' @return Fold label in 1, ..., K for each unit.
cross_fit_labels <- function(n, n_folds) {
  sample(rep_len(seq_len(n_folds), n))
}


#' Lean fits of study k on its full sample, fold A and fold B.
split_fits_lean <- function(study, k, rep, cfg, dgk) {
  folds <- split_folds(cfg$n[k], cfg$frac_a)
  samples <- list(full = study, fold_a = subset_study(study, folds$a),
                  fold_b = subset_study(study, folds$b))
  rows <- lapply(FITS, function(fit) {
    est <- dvds_bounds_lean(samples[[fit]], cfg$lambdas[[k]], cfg$outcomes,
                            dgk, cfg$nuisance_folds)
    data.frame(rep = rep, study = k, fit = fit, est)
  })
  do.call(rbind, rows)
}


#' Lean fits of study k for M7, by standard cross-fitting: for each fold f,
#' models fitted on its complement give the complement's bounds (cf_train_f,
#' for the weights) and fold f's bounds (cf_held_f).
cross_fits_lean <- function(study, k, rep, cfg, dgk) {
  labels <- cross_fit_labels(cfg$n[k], cfg$folds)
  rows <- lapply(seq_len(cfg$folds), function(f) {
    est <- dvds_bounds_lean_train_test(
      subset_study(study, which(labels != f)),
      subset_study(study, which(labels == f)), cfg$lambdas[[k]],
      cfg$outcomes, dgk, cfg$nuisance_folds)
    rbind(data.frame(rep = rep, study = k, fit = sprintf("cf_train_%d", f),
                     est$train),
          data.frame(rep = rep, study = k, fit = sprintf("cf_held_%d", f),
                     est$test))
  })
  do.call(rbind, rows)
}


#' Lean DVDS estimates for one replicate: full, fold A and fold B of both
#' studies, then the cross-fitting fits of both.
#'
#' @return data.frame of dvds_bounds_lean() rows with rep, study and fit.
simulate_replicate_quick <- function(rep, cfg, progress) {
  dgk <- load_dgk(1)
  studies <- Map(simulate_study, cfg$n, cfg$noise_scale,
                 cfg$binary_denominator)
  split_rows <- lapply(seq_along(studies), function(k) {
    split_fits_lean(studies[[k]], k, rep, cfg, dgk)
  })
  cross_rows <- lapply(seq_along(studies), function(k) {
    cross_fits_lean(studies[[k]], k, rep, cfg, dgk)
  })
  progress(message = sprintf("rep %d/%d", rep, cfg$reps))
  do.call(rbind, c(split_rows, cross_rows))
}


#' Forest fits of study k on its full sample, fold A and fold B, drawn and
#' fitted exactly as in simulate_replicate() of R/causal_sim.R.
split_fits_forest <- function(study, k, rep, cfg) {
  folds <- split_folds(cfg$n[k], cfg$frac_a)
  samples <- list(full = study, fold_a = subset_study(study, folds$a),
                  fold_b = subset_study(study, folds$b))
  rows <- lapply(FITS, function(fit) {
    est <- dvds_bounds(samples[[fit]], cfg$lambdas[[k]], cfg$num_threads,
                       cfg$outcomes)
    data.frame(rep = rep, study = k, fit = fit, est)
  })
  do.call(rbind, rows)
}


#' Forest fits of study k for M7, by standard cross-fitting: for each fold
#' f, forests fitted on its complement give the complement's bounds
#' (cf_train_f, for the weights) and fold f's bounds (cf_held_f).
cross_fits_forest <- function(study, k, rep, cfg) {
  labels <- cross_fit_labels(cfg$n[k], cfg$folds)
  rows <- lapply(seq_len(cfg$folds), function(f) {
    held <- subset_study(study, which(labels == f))
    train <- subset_study(study, which(labels != f))
    est <- dvds_bounds_train_test(train, held, cfg$lambdas[[k]],
                                  cfg$num_threads, cfg$outcomes)
    rbind(data.frame(rep = rep, study = k, fit = sprintf("cf_train_%d", f),
                     est$train),
          data.frame(rep = rep, study = k, fit = sprintf("cf_held_%d", f),
                     est$test))
  })
  do.call(rbind, rows)
}


#' Forest DVDS estimates for one replicate.
#'
#' Full, fold A and fold B of both studies are fitted first, in the order of
#' R/causal_sim.R: grf draws its seeds from R's random-number stream, so
#' these fits equal that script's for the same seed and settings. Only then
#' are the cross-fitting folds drawn and fitted.
#'
#' @return data.frame of dvds_bounds() rows with rep, study and fit.
simulate_replicate_forest <- function(rep, cfg, progress) {
  studies <- Map(simulate_study, cfg$n, cfg$noise_scale,
                 cfg$binary_denominator)
  split_rows <- lapply(seq_along(studies), function(k) {
    split_fits_forest(studies[[k]], k, rep, cfg)
  })
  cross_rows <- lapply(seq_along(studies), function(k) {
    cross_fits_forest(studies[[k]], k, rep, cfg)
  })
  progress(message = sprintf("rep %d/%d", rep, cfg$reps))
  do.call(rbind, c(split_rows, cross_rows))
}


#' Runs every replicate in parallel behind one progress bar.
#'
#' @return data.frame with one row per (rep, study, fit, outcome, Lambda).
run_simulation_quick <- function(cfg) {
  # grf threads per forest fit (forest nuisances only).
  cfg$num_threads <- set_backend(cfg)
  on.exit(future::plan(future::sequential), add = TRUE)
  forest <- cfg$nuisances == "forest"
  replicate_fit <- if (forest) {
    simulate_replicate_forest
  } else {
    simulate_replicate_quick
  }
  seeds <- replicate_seeds(cfg$seed, cfg$reps)
  if (!requireNamespace("progress", quietly = TRUE)) {
    stop("package 'progress' is required; run renv::restore()")
  }
  bar <- progressr::handler_progress(
    format = "[:bar] :percent ETA :eta | :message", clear = FALSE,
    enable = TRUE)
  results <- progressr::with_progress(handlers = bar, enable = TRUE, {
    progress <- progressr::progressor(steps = cfg$reps)
    future.apply::future_lapply(
      seq_len(cfg$reps), replicate_fit, cfg = cfg, progress = progress,
      future.seed = seeds, future.chunk.size = 1,
      future.packages = if (forest) "grf" else NULL)
  })
  do.call(rbind, results)
}


# ----------------------------- confidence intervals ------------------------
#' The R x 2 matrices (study 1 at Lambda index i, study 2 at index j) of
#' one fit, for each column of bound_arrays().
pick_fit <- function(arrays, fit, i, j) {
  lapply(arrays[[fit]], function(a) cbind(a[[1]][, i], a[[2]][, j]))
}


#' Plug-in weights (omega_0, omega_1) from one fit, with its variances
#' rescaled by `scale` (an R x 2 matrix) to the size the final CI uses.
fit_weights <- function(fit, scale, z) {
  list(w_L = optimal_weight_grid(fit$lower[, 1], fit$lower[, 2],
                                 fit$s2_lower[, 1] * scale[, 1],
                                 fit$s2_lower[, 2] * scale[, 2],
                                 lower = TRUE, z),
       w_U = optimal_weight_grid(fit$upper[, 1], fit$upper[, 2],
                                 fit$s2_upper[, 1] * scale[, 1],
                                 fit$s2_upper[, 2] * scale[, 2],
                                 lower = FALSE, z))
}


#' Cross-fitted convex CI at one pair: fold f's weights from its complement
#' (variances scaled to the full n, the size of the final CI's variance),
#' its bounds from fold f itself.
cross_fit_pair_ci <- function(arrays, i, j, cfg) {
  full <- pick_fit(arrays, "full", i, j)
  held <- lapply(seq_len(cfg$folds), function(f) {
    pick_fit(arrays, sprintf("cf_held_%d", f), i, j)
  })
  weights <- lapply(seq_len(cfg$folds), function(f) {
    train <- pick_fit(arrays, sprintf("cf_train_%d", f), i, j)
    fit_weights(train, train$n / full$n, cfg$z)
  })
  column <- function(x, name) lapply(x, `[[`, name)
  cross_fit_convex_ci(column(held, "lower"), column(held, "upper"),
                      column(held, "s2_lower"), column(held, "s2_upper"),
                      column(weights, "w_L"), column(weights, "w_U"), cfg$z)
}


#' CIs of every method in METHODS_QUICK at study 1's Lambda index i and
#' study 2's index j.
pair_cis_quick <- function(arrays, i, j, cfg) {
  full <- pick_fit(arrays, "full", i, j)
  a <- pick_fit(arrays, "fold_a", i, j)
  b <- pick_fit(arrays, "fold_b", i, j)
  # Fold-A variances rescaled to fold-B size (derivations.tex,
  # eq. est-oracleweight).
  split <- fit_weights(a, a$n / b$n, cfg$z)
  cis <- list(
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 1, 1,
              cfg$z),
    convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper, 0, 0,
              cfg$z),
    intersection_ci(full$lower, full$upper, full$s2_lower, full$s2_upper,
                    cfg$z_calibrated),
    convex_ci(b$lower, b$upper, b$s2_lower, b$s2_upper, split$w_L, split$w_U,
              cfg$z),
    cross_fit_pair_ci(arrays, i, j, cfg),
    clr_intersection_ci(full$lower, full$upper, full$s2_lower, full$s2_upper,
                        cfg$p),
    clr_convex_ci(full$lower, full$upper, full$s2_lower, full$s2_upper,
                  cfg$p)
  )
  setNames(cis, METHODS_QUICK)
}


#' Summary rows for one outcome and study 1's Lambda index i, every j.
summarise_lambda_1 <- function(i, arrays, band_1, band_2, outcome, cfg) {
  rows <- list()
  for (j in band_2$lambda_index) {
    true_lo <- max(band_1$lower[i], band_2$lower[j])
    true_hi <- min(band_1$upper[i], band_2$upper[j])
    cis <- pair_cis_quick(arrays, i, j, cfg)
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
  do.call(rbind, rows)
}


#' Summary rows for every (outcome, Lambda_1, Lambda_2, method), as
#' summarise_estimates() in R/causal_sim.R, in parallel over Lambda_1.
summarise_estimates_quick <- function(estimates, truth, cfg) {
  set_backend(cfg)
  on.exit(future::plan(future::sequential), add = TRUE)
  rows <- list()
  for (outcome in unique(estimates$outcome)) {
    arrays <- bound_arrays(estimates[estimates$outcome == outcome, ],
                           cfg$fits)
    band_1 <- truth[truth$outcome == outcome & truth$study == 1, ]
    band_2 <- truth[truth$outcome == outcome & truth$study == 2, ]
    rows <- c(rows, future.apply::future_lapply(
      band_1$lambda_index, summarise_lambda_1, arrays = arrays,
      band_1 = band_1, band_2 = band_2, outcome = outcome, cfg = cfg))
  }
  do.call(rbind, rows)
}


# ----------------------------- sanity checks ------------------------------
#' Each study's own CI covers its own band, at every Lambda, up to
#' SINGLE_STUDY_SES MC s.e. (DGK's criterion; lean nuisances are correctly
#' specified, so this should hold).
check_single_study <- function(single, alpha) {
  # Coverage 1 (0) has MC s.e. 0 and gives z = +Inf (-Inf).
  zs <- (single$coverage - (1 - alpha)) / single$coverage_mcse
  worst <- single[which.min(zs), ]
  list(name = sprintf("Single-study coverage >= %g, every study and Lambda",
                      1 - alpha),
       ok = min(zs) >= -SINGLE_STUDY_SES,
       detail = sprintf("lowest %.4f (%s, study %d, Lambda %g); %d of %d below",
                        worst$coverage, worst$outcome, worst$study,
                        worst$lambda, sum(zs < -SINGLE_STUDY_SES),
                        nrow(single)))
}


#' A method's joint coverage >= 1 - alpha up to 3 MC s.e. in every pair.
check_method_valid <- function(summary, alpha, method) {
  rows <- summary[summary$method == method, ]
  zs <- (rows$cov_joint - (1 - alpha)) / rows$cov_joint_mcse
  list(name = sprintf("%s joint coverage >= %g, all pairs", method,
                      1 - alpha),
       ok = min(zs) >= -3,
       detail = sprintf("min coverage %.4f; %d of %d pairs below by > 3 MC s.e.",
                        min(rows$cov_joint), sum(zs < -3), nrow(rows)))
}


#' Prints PASS/FAIL for each check; returns TRUE if all pass.
run_sanity_checks_quick <- function(estimates, summary, single, alpha) {
  checks <- list(check_collapse_at_one(estimates), check_monotone(estimates),
                 check_single_study(single, alpha),
                 check_calibrated_valid(summary, alpha),
                 check_method_valid(summary, alpha, M3),
                 check_method_valid(summary, alpha, M8),
                 check_method_valid(summary, alpha, M9),
                 check_method_valid(summary, alpha, M7))
  cat(paste("\nSanity checks (the cross-fit check may fail at ties, as in",
            "the toy; with forests at k = 15, study 1 is known to undercover",
            "its own band, as in the forest run behind Table 5)\n"))
  for (check in checks) {
    cat(sprintf("[%s] %s: %s\n", if (check$ok) "PASS" else "FAIL",
                check$name, check$detail))
  }
  n_passed <- sum(vapply(checks, `[[`, logical(1), "ok"))
  cat(sprintf("%d/%d checks passed\n", n_passed, length(checks)))
  n_passed == length(checks)
}


# ----------------------------- output ------------------------------
#' One line describing the nuisance models and M7's fold bounds.
nuisance_description <- function(cfg) {
  cross_fit <- paste("cross-fit: standard (fold f's nuisances from models",
                     "fitted on the other folds)")
  if (cfg$nuisances == "lean") {
    return(sprintf(paste("nuisances: parametric, true functional forms",
                         "(R/helpers/dvds_lean.R), out of fold with %d",
                         "folds | %s"), cfg$nuisance_folds, cross_fit))
  }
  sprintf(paste("nuisances: DGK's forests (dvds_sim.R lines 29-42),",
                "out of bag | %s"), cross_fit)
}


#' Records the command, settings, run times and software versions.
write_run_info_quick <- function(cfg, times, path) {
  packages <- c("future", "future.apply", "progressr", "optparse")
  if (cfg$nuisances == "forest") packages <- c("grf", packages)
  versions <- vapply(packages, function(p) {
    as.character(utils::packageVersion(p))
  }, character(1))
  estimator <- if (cfg$nuisances == "forest") {
    "estimator and DGP code, unmodified"
  } else {
    "dvds.binary(), unmodified"
  }
  lines <- c(
    sprintf("date: %s", format(times$start, "%Y-%m-%d %H:%M")),
    sprintf("git: %s", git_hash()),
    sprintf("DGK commit: %s (%s)", DGK_COMMIT, estimator),
    nuisance_description(cfg),
    sprintf("command: Rscript R/causal_sim_quick.R %s",
            paste(commandArgs(trailingOnly = TRUE), collapse = " ")),
    sprintf("seed: %d | reps: %d | n: %s | frac_a: %g | folds: %d | alpha: %g",
            cfg$seed, cfg$reps, paste(cfg$n, collapse = ", "), cfg$frac_a,
            cfg$folds, cfg$alpha),
    sprintf("outcomes: %s", paste(cfg$outcomes, collapse = ", ")),
    sprintf("lambdas study 1: %s", paste(cfg$lambdas[[1]], collapse = ", ")),
    sprintf("lambdas study 2: %s", paste(cfg$lambdas[[2]], collapse = ", ")),
    sprintf("noise scale (continuous outcome): study 1 %g, study 2 %g",
            cfg$noise_scale[1], cfg$noise_scale[2]),
    sprintf("binary denominator k: study 1 %g, study 2 %g",
            cfg$binary_denominator[1], cfg$binary_denominator[2]),
    sprintf("backend: %s | workers: %d", cfg$backend, cfg$workers),
    sprintf("fitting run time: %s", times$fitting),
    sprintf("summary run time: %s", times$summary),
    R.version.string,
    paste(sprintf("%s %s", packages, versions), collapse = " | ")
  )
  writeLines(lines, path)
}


# ----------------------------- main ------------------------------
main_quick <- function() {
  start <- Sys.time()
  cfg <- quick_config()
  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  estimates_path <- file.path(cfg$out_dir, "estimates.csv")
  if (cfg$summarise_only) {
    estimates <- read.csv(estimates_path)
    fitting <- "not rerun (--summarise-only)"
  } else {
    check_clone()
    estimates <- run_simulation_quick(cfg)
    write.csv(estimates, estimates_path, row.names = FALSE)
    fitting <- format_elapsed(start)
  }
  summary_start <- Sys.time()
  truth <- study_true_bands(cfg)
  truth <- truth[truth$outcome %in% estimates$outcome, ]
  check_grids_match(estimates, truth)
  write.csv(truth, file.path(cfg$out_dir, "true_bands.csv"), row.names = FALSE)
  summary <- summarise_estimates_quick(estimates, truth, cfg)
  write.csv(summary, file.path(cfg$out_dir, "summary_full.csv"),
            row.names = FALSE)
  write.csv(short_summary(summary), file.path(cfg$out_dir, "summary.csv"),
            row.names = FALSE)
  single <- single_study_coverage(estimates, truth, cfg$z)
  write.csv(single, file.path(cfg$out_dir, "single_study_coverage.csv"),
            row.names = FALSE)
  times <- list(start = start, fitting = fitting,
                summary = format_elapsed(summary_start))
  if (!cfg$summarise_only) {
    write_run_info_quick(cfg, times, file.path(cfg$out_dir, "run_info.txt"))
  }
  run_sanity_checks_quick(estimates, summary, single, cfg$alpha)
  print_highlights(short_summary(summary))
  cat(sprintf("\nWritten to %s\n", cfg$out_dir))
  cat(sprintf("Fitting: %s | summary: %s | total: %s\n", times$fitting,
              times$summary, format_elapsed(start)))
}


if (sys.nframe() == 0L) main_quick()
