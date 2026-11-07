# Causal sim: binary outcome with a rare outcome in study 1 (k = 15)

Date: 2026-11-05. Status: **code changes and checks done 2026-11-07
(Claude, at the user's request); the full run, table and paper steps are
the user's.**

## Outcome (2026-11-07)

Code changes 1-3 implemented as below; `simulate_study()` also gained
`binary_index()`, `binary_mean()` and `validate_binary_denominator()`.
Checks (scratch, sequential backend because the sandbox blocks the sockets
of multisession workers):
1. Defaults unchanged: 2 replicates with the n = 3000 noise settings
   reproduce reps 1-2 of `causal_sim_noise_n3000_reps1000/estimates.csv`
   exactly (984 rows, max difference 0).
2. Thinning (n = 1e6, k = 15): prevalence 0.0616 (model 0.0618); decile
   bins of eta within 1.4 binomial s.e. of 1 / (15 + exp(eta)); arm means
   0.0611 / 0.0623 against model 0.0612 / 0.0624; Y_15 <= Y_1 always, and
   X, Z, Ycontin unchanged.
3. Truth: k = 15 widths 0.051, 0.093, 0.130, 0.165; `binary_mean(X, 1)` is
   bit-identical to DGK's formula (the saved true_bands.csv differs by
   5e-16, CSV rounding only).
4. Smoke run (`--reps 2 --binary-denominator-1 15 --outcomes binary`, grids
   as in the full run): 3/3 sanity checks pass; about 10 s per replicate
   on one process with every grf thread (39 CPU-s), so the full run is
   about 1.5 h on 9 workers.

Docs updated: header of `R/causal_sim.R`, README, CLAUDE.md, and
`causal_sim_code_explained.tex` (compiled, no errors or overfull boxes).

## Why

With the DGK binary outcome unchanged, the SE ratio of the bound estimates
is just sqrt(n1 / n2) (about 2.45 at 3000/500), and the convex split never
beats the calibrated intersection (1000-rep run, best ratio 0.997). The
Bernoulli variance is already near its maximum (p ranges over 0.15-0.85,
mean 0.5), so study 2 cannot be made noisier; instead make study 1's outcome
rare, p_1(x) = 1 / (k + exp(eta(x))), eta(x) = 0.5 x1 + x2 + 0.25 x2 x3.

A quick scratch check (2026-11-05, 200 reps, n = 3000/500, not saved;
scripts `binary_k_check.R`, `binary_k_analyse.R` in that session's
scratchpad) gave:

| k | prevalence | SE ratio, best pairs | split narrower (> 2 paired s.e.) | best ratio | where |
|---|---|---|---|---|---|
| 1 | 0.50 | 2.4-2.55 | 25 (0) | 0.994 | Lambda_1 2.25-3.25 |
| 9 | 0.10 | 3.4-3.85 | 19 (9) | 0.908 | Lambda_1 1.5-2, Lambda_2 1.0-1.2 |
| 20 | 0.05 | 4.6-5.4 | 38 (20) | 0.825 | Lambda_1 1.5-2.25, Lambda_2 1.0-1.1 |

k = 15 sits between (prevalence 0.062). Its true band width for study 1 is
0.051, 0.093, 0.130, 0.165 at Lambda = 1.5, 2, 2.5, 3; study 2's (k = 1) is
0.164 at Lambda = 1.25 and 0.293 at 1.5, so the near-ties are at Lambda_2
close to 1.

Caveats to carry into the paper:
- Study 1's DVDS SE grows with Lambda when the outcome is rare (k = 9:
  0.0128 at Lambda 2, 0.0231 at Lambda 5), so the SE ratio is large only at
  small Lambda_1.
- Ties need Lambda_2 near 1, i.e. study 2 close to assuming no unmeasured
  confounding.
- The two studies share the ATE (0) only because DGK's Y does not depend on
  Z; with outcome prevalence 6% vs 50% they are quite different populations.

## Code changes (projectone)

1. `R/helpers/data_fusion_design.R`
   - `binary_index(X)`: DGK's linear index eta(X) (dvds_sim.R line 24),
     beside `continuous_mean()`.
   - `simulate_study(n, noise_scale = 1, binary_denominator = 1)`: after the
     verbatim `simulate_dgk(n)`, if `binary_denominator != 1`, thin the binary
     outcome: keep each Y = 1 with probability p_k(X) / p_1(X), where
     p_k = 1 / (k + exp(eta)). Then P(Y = 1 | X) = p_k(X) exactly, and the
     new outcome is coupled to DGK's (Y_k <= Y_1). Validate k >= 1. No extra
     random draws when k = 1, so the defaults reproduce the baseline exactly.
   - `binary_true_band(lambdas, denominator = 1)`: mu = 1 / (denominator +
     exp(eta)); `true_bands(lambdas, noise_scale = 1, binary_denominator = 1)`
     passes it on.
2. `R/causal_sim.R`
   - Flag `--binary-denominator-1` (double, default 1): k for study 1's
     binary outcome; study 2 stays at 1. Validated in `make_config`;
     `cfg$binary_denominator = c(k, 1)`.
   - `simulate_replicate`: `Map(simulate_study, cfg$n, cfg$noise_scale,
     cfg$binary_denominator)`.
   - `study_true_bands`: pass `cfg$binary_denominator[k]`.
   - Flag `--outcomes` (default "binary,continuous"): with "binary",
     `dvds_bounds()` skips `dvds.continuous()` (the slow part; DGK's verbatim
     nuisance lines stay as they are). `write_run_info` prints the outcomes
     run and the binary denominators. Check that `summarise_estimates`,
     `single_study_coverage` and the run checks only loop over outcomes
     present in the estimates.
3. Docs: header comment of `R/causal_sim.R`, README (flags),
   `partialID_Overleaf/Notes/Simulation_Notes/causal_sim_code_explained.tex`
   (thinning, new flags), CLAUDE.md current state.

## Checks before the full run

1. Defaults unchanged: `Rscript R/causal_sim.R --reps 2 --n1 3000 --n2 500
   --noise-scale-2 2.5 --lambdas-1 3:5:0.05 --lambdas-2 1:3:0.05
   --out-dir <scratch>`; estimates.csv must equal reps 1-2 of
   `src/outputs/causal_sim_noise_n3000_reps1000/estimates.csv` (replicate r
   always uses stream r).
2. Thinning: for n = 1e6 and k = 15, mean(Ybinary) is about 0.062, and the
   mean within bins of eta matches 1 / (15 + exp(eta)).
3. Truth: `binary_true_band(c(1.5, 2, 2.5, 3), 15)` widths are about 0.051,
   0.093, 0.130, 0.165 (scratch values above).
4. Smoke run: `--reps 2 --binary-denominator-1 15 --outcomes binary`.

## Full run (from the repo root, 9 cores)

```
Rscript R/causal_sim.R --reps 1000 --n1 3000 --n2 500 \
  --binary-denominator-1 15 --outcomes binary \
  --lambdas-1 1:4:0.05 --lambdas-2 1:2:0.05 \
  --workers 9 --out-dir src/outputs/causal_sim_binary_k15
```

61 x 21 = 1281 pairs. Rough time: 1.5-2 h on 9 workers (binary only; the
scratch check took 30 min for 200 reps and 3 values of k). Keep the machine
awake (e.g. `caffeinate -i` before the command on macOS).

Expected: split narrower than the calibrated intersection at Lambda_1 about
1.5-2.5, Lambda_2 about 1.0-1.2, best ratio between the k = 9 and k = 20
values (0.83-0.91), both methods valid. If not, stop and compare with the
scratch check before writing anything up.

## Table and paper

1. `python -m src.latex.causal_sim_tables --summary
   src/outputs/causal_sim_binary_k15/summary_full.csv --outcome binary
   --lambdas-1 <values with wins, e.g. 1.5 1.75 2 2.25 2.5> --top 3
   --macro-decimals 4` (no code change; the Lambda_1 rows are chosen after
   seeing the summary).
2. Via /update-paper-figs: `_<hash>` copies in `Tables/causal_sim/`,
   PROVENANCE.md section.
3. Decide (user): replace the current binary table, or add the k = 15 table
   beside it. The caption must say study 1's outcome is rare (prevalence
   about 6%, via k = 15), that the ATE is shared because Y does not depend
   on Z, and that the near-ties sit at Lambda_2 close to 1.
4. Compile; no commits by Claude.

## Open questions

- Grids: is Lambda_1 up to 4 worth the extra pairs, given the SE ratio falls
  with Lambda_1?
- Should k also be offered for study 2 (`--binary-denominator-2`) for
  symmetry with `--noise-scale-2`? Not needed for this run.
