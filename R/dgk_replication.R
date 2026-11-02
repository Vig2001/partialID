# Replication of the DGK simulation (Dorn, Guo and Kallus, Section 7.1) with
# their code, unmodified. Their repository is cloned at a pinned commit in
# external/ (see README); this script
#   1. checks the commit and that their scripts are unedited;
#   2. saves their committed outputs (the results behind the paper);
#   3. runs their dvds_sim.R and dvds_plotting.r, each in a fresh R process
#      with the project's renv library;
#   4. copies our outputs out, then restores theirs in the clone;
#   5. compares the two runs per (outcome, method, Lambda).
# Exact agreement is not expected: the data draws match (same set.seed(1)),
# but forest fits differ across grf versions. Agreement is judged within
# Monte Carlo error.
#
# Usage, from the repo root:
#   Rscript R/dgk_replication.R [--check-only]
#
# Plan:    quality_reports/plans/2026-11-02_dgk-exact-replication.md
# Outputs: src/outputs/dgk_replication/

if (!file.exists(file.path("R", "dgk_replication.R"))) {
  stop("run from the repo root: Rscript R/dgk_replication.R")
}
source(file.path("R", "helpers", "dgk.R"))  # DGK_DIR, DGK_COMMIT, check_clone

DGK_OUTPUT <- file.path("Simulations", "dvds_output")
RESULT_FILES <- c("binary_dvds.csv", "binary_zsb.csv", "continuous_dvds.csv",
                  "continuous_zsb.csv", "binary_coverageRates.tex",
                  "continuous_coverageRates.tex",
                  file.path("art", "binary_comparison.pdf"),
                  file.path("art", "continuous_comparison.pdf"))
OUT_DIR <- file.path("src", "outputs", "dgk_replication")
TRUTH_MC_SIZE <- 500000  # As in dvds_plotting.r (binary truth by simulation).
Z_FLAG <- 3  # |z| above which a run-to-run difference is reported.


# ----------------------------- clone ------------------------------
#' Copies the listed result files from one folder to another.
copy_results <- function(from, to) {
  for (file in RESULT_FILES) {
    dir.create(dirname(file.path(to, file)), recursive = TRUE,
               showWarnings = FALSE)
    stopifnot(file.copy(file.path(from, file), file.path(to, file),
                        overwrite = TRUE))
  }
}


#' Runs one DGK script, unmodified, in a fresh R process.
#'
#' The child starts in the repo root, so .Rprofile activates renv; it then
#' moves into the clone, where the script's relative paths resolve.
run_dgk_script <- function(script) {
  cat(sprintf("\nRunning %s (unmodified)\n", script))
  code <- sprintf("setwd('%s'); source('%s')", DGK_DIR, script)
  status <- system2("Rscript", c("-e", shQuote(code)))
  if (status != 0) stop(sprintf("%s failed with status %d", script, status))
}


# ----------------------------- comparison ------------------------------
#' True sharp bounds per (outcome, Lambda), computed as in dvds_plotting.r.
#'
#' Continuous: closed form. Binary: Monte Carlo over X, with set.seed(1) as
#' in their script, so both runs are judged against the same truth.
true_bounds <- function(lambdas) {
  k <- (1 + 1 / 3) * (lambdas^2 - 1) / lambdas *
    dnorm(qnorm(lambdas / (lambdas + 1)))
  continuous <- data.frame(outcome = "continuous", Lambda = lambdas,
                           true_lower = -k, true_upper = k)
  set.seed(1)
  X <- matrix(runif(TRUTH_MC_SIZE * 5, -1, 1), nrow = TRUTH_MC_SIZE)
  e <- 1 / (1 + exp(X[, 1] + 0.5 * (X[, 2] > 0) + 0.5 * X[, 2] * X[, 3]))
  mu <- 1 / (1 + exp(0.5 * X[, 1] + X[, 2] + 0.25 * X[, 2] * X[, 3]))
  upper <- sapply(lambdas, function(l) {
    mean(e * mu + (1 - e) * pmin(1 - 1 / l + mu / l, mu * l)) -
      mean((1 - e) * mu + e * pmax(1 - l + mu * l, mu / l))
  })
  lower <- sapply(lambdas, function(l) {
    mean(e * mu + (1 - e) * pmax(1 - l + mu * l, mu / l)) -
      mean((1 - e) * mu + e * pmin(1 - 1 / l + mu / l, mu * l))
  })
  binary <- data.frame(outcome = "binary", Lambda = lambdas,
                       true_lower = lower, true_upper = upper)
  rbind(continuous, binary)
}


#' All four result tables of one run, stacked, with outcome and method.
read_run <- function(folder) {
  parts <- list()
  for (outcome in c("binary", "continuous")) {
    for (method in c("dvds", "zsb")) {
      path <- file.path(folder, sprintf("%s_%s.csv", outcome, method))
      table <- read.table(path, header = TRUE)
      parts[[length(parts) + 1]] <- data.frame(
        outcome = outcome, method = method,
        table[c("Lambda", "lower", "upper", "lower.CI", "upper.CI")])
    }
  }
  do.call(rbind, parts)
}


#' Per (outcome, method, Lambda): mean bounds and DVDS coverage for one run.
#'
#' Coverage uses DGK's strict rule: lower.CI < true_lower and
#' true_upper < upper.CI. ZSB has no CI, so its coverage is NA.
summarise_run <- function(run, truth) {
  run <- merge(run, truth, by = c("outcome", "Lambda"))
  run$covered <- run$lower.CI < run$true_lower &
    run$true_upper < run$upper.CI
  groups <- split(run, list(run$outcome, run$method, run$Lambda), drop = TRUE)
  rows <- lapply(groups, function(g) {
    data.frame(outcome = g$outcome[1], method = g$method[1],
               Lambda = g$Lambda[1], reps = nrow(g),
               mean_lower = mean(g$lower), sd_lower = sd(g$lower),
               mean_upper = mean(g$upper), sd_upper = sd(g$upper),
               coverage = mean(g$covered))
  })
  do.call(rbind, rows)
}


#' Joins the two runs and adds two-sample z statistics for each difference.
#'
#' The runs use the same data draws, so these z statistics are conservative.
compare_runs <- function(original, replication) {
  both <- merge(original, replication,
                by = c("outcome", "method", "Lambda"),
                suffixes = c("_original", "_replication"))
  z_mean <- function(end) {
    diff <- both[[paste0("mean_", end, "_replication")]] -
      both[[paste0("mean_", end, "_original")]]
    se <- sqrt(both[[paste0("sd_", end, "_original")]]^2 /
                 both$reps_original +
                 both[[paste0("sd_", end, "_replication")]]^2 /
                 both$reps_replication)
    diff / se
  }
  both$z_lower <- z_mean("lower")
  both$z_upper <- z_mean("upper")
  p1 <- both$coverage_original
  p2 <- both$coverage_replication
  se <- sqrt(p1 * (1 - p1) / both$reps_original +
               p2 * (1 - p2) / both$reps_replication)
  both$z_coverage <- ifelse(se > 0, (p2 - p1) / se, 0)
  both[order(both$outcome, both$method, both$Lambda), ]
}


#' Prints the largest |z| per outcome and method, and both coverage summaries.
print_comparison <- function(comparison, folders) {
  cat("\nReplication vs DGK's committed results (|z| > 3 flagged)\n")
  for (key in unique(paste(comparison$outcome, comparison$method))) {
    rows <- comparison[paste(comparison$outcome, comparison$method) == key, ]
    zs <- abs(c(rows$z_lower, rows$z_upper))
    z_cov <- abs(rows$z_coverage)
    cat(sprintf("%-16s max |z| bounds %.2f, coverage %s%s\n", key, max(zs),
                if (all(is.na(z_cov))) "n/a" else
                  sprintf("%.2f", max(z_cov, na.rm = TRUE)),
                if (max(c(zs, z_cov), na.rm = TRUE) > Z_FLAG) "  <-- check"
                else ""))
  }
  for (outcome in c("continuous", "binary")) {
    file <- sprintf("%s_coverageRates.tex", outcome)
    cat(sprintf("\n%s DVDS coverage\n  DGK:  %s\n  ours: %s\n", outcome,
                readLines(file.path(folders$original, file)),
                readLines(file.path(folders$replication, file))))
  }
}


#' Records the commits, settings and software versions of the run.
write_run_info <- function(path, start, ran) {
  ours <- tryCatch(system2("git", c("rev-parse", "--short", "HEAD"),
                           stdout = TRUE, stderr = FALSE),
                   error = function(e) "unknown")
  packages <- c("grf", "dplyr", "ggplot2", "progress")
  versions <- vapply(packages, function(p) {
    as.character(utils::packageVersion(p))
  }, character(1))
  lines <- c(
    sprintf("date: %s", format(start, "%Y-%m-%d %H:%M")),
    sprintf("projectone git: %s", ours),
    sprintf("DGK commit: %s (unmodified dvds_sim.R, dvds_plotting.r)",
            DGK_COMMIT),
    sprintf("settings: DGK's own (set.seed(1), 1000 reps, n = 1000, 20 Lambdas)"),
    if (ran) sprintf("run time: %.1f min",
                     as.numeric(difftime(Sys.time(), start, units = "mins")))
    else "run time: comparison only (--check-only)",
    R.version.string,
    paste(sprintf("%s %s", packages, versions), collapse = " | "))
  writeLines(lines, path)
}


# ----------------------------- main ------------------------------
main <- function() {
  parser <- optparse::OptionParser(option_list = list(
    optparse::make_option("--check-only", dest = "check_only",
                          action = "store_true", default = FALSE,
                          help = "Compare saved outputs; do not rerun DGK.")))
  args <- optparse::parse_args(parser)
  check_clone()
  folders <- list(original = file.path(OUT_DIR, "original"),
                  replication = file.path(OUT_DIR, "replication"))
  start <- Sys.time()
  if (!args$check_only) {
    # Their committed results, taken before their scripts overwrite them.
    copy_results(file.path(DGK_DIR, DGK_OUTPUT), folders$original)
    run_dgk_script("dvds_sim.R")
    run_dgk_script("dvds_plotting.r")
    copy_results(file.path(DGK_DIR, DGK_OUTPUT), folders$replication)
    copy_results(folders$original, file.path(DGK_DIR, DGK_OUTPUT))
    write_run_info(file.path(OUT_DIR, "run_info.txt"), start, ran = TRUE)
  }
  original <- read_run(folders$original)
  truth <- true_bounds(sort(unique(original$Lambda)))
  comparison <- compare_runs(summarise_run(original, truth),
                             summarise_run(read_run(folders$replication),
                                           truth))
  write.csv(comparison, file.path(OUT_DIR, "comparison.csv"),
            row.names = FALSE)
  print_comparison(comparison, folders)
  cat(sprintf("\nWritten to %s\n", OUT_DIR))
}


if (sys.nframe() == 0L) main()
