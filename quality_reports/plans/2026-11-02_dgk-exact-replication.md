# Plan: run DGK's code exactly; replicate their simulation

Date: 2026-11-02. Status: **approved 2026-11-02** (dplyr + ggplot2 into renv;
gitignored clone in `external/` at a pinned commit; RHC real-data analysis out
of scope).

## Their repository (github.com/CausalML/DoublyValidDoublySharpSensitivity)
- `dvds_sim.R` (paper Section 7.1): `set.seed(1)`, 1000 replicates, n = 1000,
  `Lambdas = seq(1, 3, length.out = 20)`, binary and continuous outcomes, DVDS
  and ZSB; one process, grf on all cores; writes `Simulations/dvds_output/`.
- `simplified_dvds.R`: `phi.continuous`, `dvds.continuous`, `dvds.binary`,
  `extrema.aipw` (ZSB).
- `dvds_plotting.r`: closed-form truth, box plots, coverage sentences
  (`continuous_coverageRates.tex`).
- Their output CSVs are committed, so the replication is compared with their
  actual numbers.
- Packages: grf, progress, dplyr, ggplot2. No LICENSE file: use privately, do
  not redistribute; a public release fetches it at the pinned commit.
- Exact numeric agreement is unlikely (R 4.2.1 and an older grf); the target
  is distributional agreement.

## Steps
1. **Vendor:** `git clone` into `external/DoublyValidDoublySharpSensitivity/`;
   pin the commit hash (README, run_info); `external/` gitignored. Their files
   are never edited.
2. **Replicate unmodified:** from inside the clone, `Rscript dvds_sim.R` then
   `Rscript dvds_plotting.r` (user runs it; ~1 h estimated). Their committed
   outputs stay recoverable from git in the clone.
   `R/dgk_replication_check.R` compares ours with theirs per Lambda (mean and
   sd of bounds, coverage); results to `src/outputs/dgk_replication/`.
3. **causal_sim uses their code:** `R/helpers/dvds.R` (port) replaced by an
   adapter that `source()`s their `simplified_dvds.R` and calls
   `dvds.continuous()` unchanged; data generation and nuisance lines copied
   verbatim from `dvds_sim.R` (including the `Ybinary` draw, so random-number
   use matches), marked with source line numbers. Variances recovered exactly
   as s2 = ((upper.CI - upper) / z)^2. Forest functions masked inside each
   worker with 1-thread versions (their file untouched; grf output does not
   depend on threads). Their arm-wise columns (averaged over all Lambda) and
   `report_arm_bounds` dropped.
4. **Defaults:** Lambda grid = their 20 values (flag kept). The two-study
   design (folds, CIs, summaries) is unchanged.
5. **Verify, rerun, document:** port vs adapter on 2 replicates; rerun the
   causal simulation (user, ~85-90 min); update notes .tex, README, CLAUDE.md.

The parametric-nuisance plan becomes a later `--nuisance` option, not a
replacement.
