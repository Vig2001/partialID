# Plan: binary and continuous outcomes in the causal simulation

Date: 2026-11-02. Status: **approved 2026-11-02**, with one addition: both R
scripts print how long they took.

## Why it is cheap (measured, n = 2000, one thread)
- Continuous fit as now (11 forests): 11.1 s.
- Binary adds two probability forests for mu0, mu1 (e is shared): ~2.6 s.
- `dvds.binary()` over 20 Lambdas: 0.04 s.
- Binary truth: DGK's Monte Carlo (500,000 draws), once per Lambda, outside
  the loop: seconds.
So about +20-25% run time (~105 min for the default run).

Fitting the binary forests also makes our nuisance code a full verbatim copy
of DGK's `dvds_sim.R` lines 29-42, so random-number use matches theirs.

## Steps
1. `R/helpers/dgk.R`: `simulate_dgk()` keeps `Ybinary`; `dgk_nuisances()`
   runs lines 29-42 verbatim (binary forests included); `dvds_bounds()`
   returns rows for both outcomes (`outcome` column). Binary uses their
   `dvds.binary()` per Lambda (their lines 45-49); variances recovered from
   the CIs as for continuous.
2. Truth: binary bands from their `dvds_plotting.r` Monte Carlo (500,000
   draws, `set.seed(1)`, RNG state restored), with its Monte Carlo s.e.;
   computed once in the summary stage, saved to `true_bands.csv`. Shared with
   `R/dgk_replication.R` (moved to `R/helpers/dgp.R`, since renamed `data_fusion_design.R`).
3. Summary by outcome. Binary bands are asymmetric, so `gap` becomes
   `gap_L`, `gap_U` (toy notation). `summary.csv`: 3200 rows.
4. Checks per outcome: collapse at Lambda = 1 (binary: tau = 0.5, kappa = mu,
   so DVDS reduces to AIPW), monotone in Lambda, M6 validity, single-study
   coverage (compare DGK binary 94.6-97.5%).
5. Run time printed at the end of `causal_sim.R` and `dgk_replication.R`
   (and kept in `run_info.txt`).
6. Smoke test; update notes .tex, README, CLAUDE.md. Full run by the user.

Deviation from the draft: no `--outcomes` flag. All nuisance lines are fitted
anyway (to stay verbatim), so a flag would save almost nothing.
