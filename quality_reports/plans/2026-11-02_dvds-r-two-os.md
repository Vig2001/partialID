# Plan: two-OS DVDS simulation in R (naive intersection vs convex split)

Date: 2026-11-02. Status: **approved 2026-11-02**, with one addition: a note
on study 2 and an explanation of the code in a `.tex` file in
`partialID_Overleaf/Notes/Project_Notes/`.
**Renamed 2026-11-02 (after the full run):** `R/dvds_sim.R` → `R/causal_sim.R`,
outputs `src/outputs/causal_sim/`, notes `causal_sim_code_explained.tex`; M1
(naive intersection) dropped from the comparators. Names below are as planned.
Builds on `2026-10-28_dvds-simulation.md`; the decisions below settle its open
questions for this first simulation.

## Decisions (Q&A, 2026-11-02)
- **Language and layout:** R, in a top-level `R/` folder as in README "Planned
  structure"; `renv` at the repo root; outputs to `src/outputs/dvds_sim/`.
  No `R/stata/` yet (CLR is out of scope).
- **DGP:** DGK's continuous DGP, unchanged, for both studies: all five
  covariates observed, ATE = 0, no unobserved confounding. The studies differ
  only in $n_k$ and $\Lambda_k$.
- **Sizes:** study 1 $n_1 = 2000$, study 2 $n_2 = 500$ (SE ratio about 2).
- **$\Lambda$:** each study estimated on the grid $1, 1.25, \ldots, 3$ (9
  values; flag `--lambdas`); CIs formed for all 81 pairs
  $(\Lambda_1, \Lambda_2)$. Changed on 2026-11-02 from DGK's 20-value grid
  (400 pairs) at the user's request; the grid size does not affect run time
  (10.6 s per n = 2000 fit for 6, 9 or 20 values).
- **Methods:** naive intersection (M1), convex split (M3), calibrated
  intersection (M6), study 1 only, study 2 only. Labels as in `toy_sim.py`.
- **Split:** each study split at random into fold A (30%) and fold B (70%);
  DVDS run separately within each fold.
- **Run:** 1000 replicates, parallel over replicates on the laptop (one grf
  thread per worker).
- **Progress bar:** one bar over all DVDS fits (replicates x 6), label naming
  the replicate, method, study and fold most recently started, with ETA.
- **Packages (renv.lock):** renv, grf, progress, optparse, future,
  future.apply, progressr. Everything else in base R.
- **Comments:** academic, minimal words, each step explained.

## What the design implies
With identical DGPs, study $k$'s population band is
$[L_k, U_k] = [-\tfrac43 K_{\Lambda_k}, +\tfrac43 K_{\Lambda_k}]$, with
$K_\Lambda = \tfrac{\Lambda^2-1}{\Lambda}\phi(\Phi^{-1}(\tfrac{\Lambda}{\Lambda+1}))$
(`causal_sim_explained.tex`, formerly `dvds_sim_explained.tex`). The bands are nested, so:
- the true fused band $[\max_k L_k, \min_k U_k]$ is the band of
  $\min(\Lambda_1, \Lambda_2)$;
- the gap is the same at both ends, $g = \tfrac43(K_{\Lambda_{\max}} - K_{\Lambda_{\min}})$,
  i.e. the diagonal $g_L = g_U$ of the toy;
- $\Lambda_1 = \Lambda_2$ is the zero-gap case (naive intersection
  undercovers); $\Lambda_2 < \Lambda_1$ (small study narrower but noisier) is
  where the convex weights face a real trade-off.

## Estimator (port of DGK `simplified_dvds.R`, continuous case)
- Nuisances with grf, out-of-bag: propensity (probability forest), arm-wise
  means (regression forests), arm-wise $\tau$-quantiles (quantile forest), and
  CVaR from the weights of one regression forest fitted at the median $\tau$
  (DGK's shortcut, kept).
- For each $\Lambda$: influence-function values for $\mathrm{ATE}^\pm$;
  $\hat L = \overline{\varphi^-}$, $\hat U = \overline{\varphi^+}$, variances
  $\widehat{\mathrm{var}}(\varphi^\pm)/n$.
- Per-arm bounds ($\psi_1^\pm, \psi_0^\pm$) also returned, for sanity checks
  only. Note: DGK's own per-arm columns (`upper.1 = mean(phi1plus)`) average
  over all $\Lambda$; this is not ported.

## CI methods (port of `src/helpers/ci_methods.py`, `optimal_weight_grid`)
- $z = \Phi^{-1}(1-\alpha/2)$ per end, $\alpha = 0.05$.
- **M1:** $\max_k \hat L_k - z\,\widehat{se}$ of the attaining study;
  $\min_k \hat U_k + z\,\widehat{se}$ likewise.
- **M6:** as M1 with $z = \Phi^{-1}(1-\alpha/4)$.
- **Single study:** $\hat L_k - z\,\widehat{se}_k$, $\hat U_k + z\,\widehat{se}_k$ (DGK's CI).
- **M3 convex split:** $\omega_0$ maximises
  $\omega \hat L_{1A} + (1-\omega)\hat L_{2A} - z\sqrt{\omega^2 \tilde s_1 + (1-\omega)^2 \tilde s_2}$
  on a 1001-point grid, with $\tilde s_k = \widehat{\mathrm{var}}(\varphi_{kA})/n_{kB}$
  (fold-A variance at fold-B size, as in `toy_sim.py` and
  `derivations.tex` eq. `est-oracleweight`). The lower limit is
  $\omega_0 \hat L_{1B} + (1-\omega_0)\hat L_{2B} - z\sqrt{\omega_0^2 s_{1B} + (1-\omega_0)^2 s_{2B}}$.
  Upper end analogous with $\omega_1$ (minimised).

## Truth and metrics
- Truth: $L_k, U_k$ in closed form (above); target band
  $[\max_k L_k, \min_k U_k]$; $\theta = 0$.
- Per $(\Lambda_1, \Lambda_2, \text{method})$, as `mc_metrics.summarise_ci`:
  `cov_joint, cov_lo, cov_hi, cov_theta, width, bias_L, bias_U, wL, wU`, each
  with `_mcse`, plus `width_median`.

## Files
```
R/
├── dvds_sim.R            # entry: flags, parallel loop, progress bar, save
└── helpers/
    ├── dgp.R             # simulate_dgk(), split_folds(), true_band()
    ├── dvds.R            # fit_nuisances(), phi_continuous(), dvds_bounds()
    ├── ci_methods.R      # intersection_ci(), convex_ci(), optimal_weight_grid()
    └── mc_metrics.R      # mean_se(), summarise_ci()
renv.lock, renv/ (library gitignored by renv), .Rprofile, .renvignore
```
- Run from the repo root: `Rscript R/dvds_sim.R [flags]`.
- Flags: `--seed 1`, `--reps 1000`, `--n1 2000`, `--n2 500`, `--frac-a 0.3`,
  `--alpha 0.05`, `--workers` (default: available cores minus 1),
  `--out-dir src/outputs/dvds_sim`, `--summarise-only` (recompute CIs and
  metrics from a saved `estimates.csv`, no refitting).
- `.renvignore` excludes `scripts_misc/`, `partialID_Overleaf/`,
  `project1venv/` and `src/` from renv's dependency scan.

## Outputs (`src/outputs/dvds_sim/`)
- `estimates.csv`: one row per (replicate, study, fit in {full, fold_a,
  fold_b}, $\Lambda$): bounds, variances, $n$, per-arm bounds. About 120k rows
  (~22 MB measured); gitignored.
- `summary.csv`: one row per $(\Lambda_1, \Lambda_2, \text{method})$, 405 rows.
- `run_info.txt`: git hash, command, seed, workers, R and package versions,
  run time.

## Reproducibility
- R has no generator objects (deviation from AGENTS.md §5): each replicate
  gets its own L'Ecuyer-CMRG stream derived from `--seed`, so results do not
  depend on the number of workers. grf seeds are drawn from that stream.
- AGENTS.md is written for Python; R code follows its principles (small pure
  functions, no hard-coded paths, CLI flags, validation with `stop()`), with
  roxygen-style `#'` docs in place of Google docstrings.

## Sanity checks (printed PASS/FAIL, as in `toy_sim.py`)
1. At $\Lambda = 1$, $\hat L = \hat U$ exactly in every fit.
2. Mean $\hat L$ non-increasing and mean $\hat U$ non-decreasing in $\Lambda$.
3. Per-arm means (study 1, full sample) against per-arm truth
   $\psi_1^+ = K_\Lambda E[(1-e)\sigma]$, $\psi_0^- = -K_\Lambda E[e\sigma]$
   (one-off large Monte Carlo). These are asymmetric, so they catch sign errors
   that the symmetric ATE band hides. Reported, not PASS/FAIL (DVDS has
   finite-sample bias).
4. M6 joint coverage $\geq 0.95$ up to 3 MC s.e. in every pair.
5. Optional replication: `--n1 1000` makes "study 1 only" directly comparable
   with DGK's reported 90.5–98.2%.

## Caveats
- Study 2's folds are small (fold A 150, about 66 treated; fold B 350). DVDS
  already undercovers at $n = 1000$ near $\Lambda = 1$ (DGK: 90.5%), so M3 may
  undercover there because of the per-study estimator, not the weighting. The
  single-study CIs help separate the two.
- Width comparisons only between methods with valid coverage.
- Variances are estimated: outside the known-variance Gaussian section.

## Timing (measured 2026-11-02, grf 2.6.1, this laptop)
One replicate (6 fits) takes 8.6 s with all threads and 24.7 s with one
thread. With W single-thread workers, about $24.7/W$ s per replicate, so
roughly 45–60 min for 1000 replicates. Homebrew R compiles all packages
from source (grf about 80 s).

**Measured after implementation:** 20 replicates on 9 multisession workers
(10 cores) took 2.1 min: ~25 s worker start-up, then ~5 s per replicate.
1000 replicates: about 85 min. The summary step takes seconds.

## Steps
1. renv: bare init at the repo root; install the seven packages; snapshot.
   renv's global cache lies outside the sandbox, so the cache is disabled and
   packages go into `renv/library/`.
2. Write `R/helpers/*.R` and `R/dvds_sim.R`.
3. Smoke test: `--reps 4 --workers 2`; check the bar, outputs and checks.
4. Time 20 replicates on all workers; report the projected full run. You run
   the full simulation yourself (e.g. `caffeinate -is Rscript R/dvds_sim.R`).
5. Update README (R setup, run command, timing), `.gitignore`, and the CLAUDE.md
   "Current state" section.
