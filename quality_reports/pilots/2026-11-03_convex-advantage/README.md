# Pilots for the convex-advantage design (2026-11-03)

These are the evidence behind
`quality_reports/plans/2026-11-03_convex-advantage-noise-design.md`. They
are exploratory scripts, not part of `R/`. They source the repo's helpers
(`R/helpers/dgk.R`, `data_fusion_design.R`, `ci_methods.R`) and need DGK's
clone in `external/`.

Run everything from the repo root, with `PILOT_OUT` set to this folder:

```
D=quality_reports/pilots/2026-11-03_convex-advantage
PILOT_OUT=$D Rscript $D/pilot_noise.R 10 8
```

The arguments to the pilot scripts are reps and cores. Forking via
`parallel::mclapply` means results depend on the core count, so the pilots
should reproduce only with 8 cores (only the final predictor has been
rerun from this folder and checked). The predictors are seeded and deterministic
given the pilot CSVs.

| Script | What it does | Output | Time |
|---|---|---|---|
| `pilot_se.R` | Forest pilot: SD of DGK's DVDS bounds with study 2's propensity logit scaled by gamma = 1, 2, 3; n = 1000, 5000; 10 reps | `pilot_se.csv` | 6 min |
| `predict_widths.R` | Gaussian predictor for the gamma designs: M3/M6/M5 widths over all 400 pairs, for n in {(2000, 500), (5000, 1000)} and frac_a in {0.1, 0.15, 0.2, 0.3} | `predict_designs.csv`, `predict_pairs.csv` | 2 min |
| `pilot_noise.R` | Forest pilot: SD with the continuous outcome noise scaled by c = 1, 1.5, 2, 2.5; n = 500, 1000, 2000; 10 reps | `pilot_noise.csv` | 3.5 min |
| `predict_noise.R` | Gaussian predictor for the noise designs on one 0.05-step grid on [1, 3] | `predict_noise_pairs.csv` | 3 min |
| `pilot_study1_high.R` | Forest pilot: study 1 (c = 1, n = 2000) at Lambda in [3, 5]: SD and bias | `pilot_study1_high.csv` | 1 min |
| `predict_grids.R` | Gaussian predictor for the final design: study 1 on [3, 5], study 2 on [1, 3], c in {2, 2.5}, frac_a in {0.2, 0.3}; pairs within 0.25 of the matching curve | `predict_grids_pairs.csv` | 1 min |

`predict_widths.R` needs `pilot_se.csv`; `predict_noise.R` needs
`pilot_noise.csv`; `predict_grids.R` needs `pilot_noise.csv` and
`pilot_study1_high.csv`. The large pair-level CSVs (3-4 MB) are not kept;
the predictors regenerate them.

## Gaussian predictor: what it assumes
- Bound estimates are drawn as N(true bound, sd^2 / m), with sd taken from
  the forest pilots and m the size of the sample used (full, fold A or
  fold B).
- M3 and M6 use the causal sim's own functions.
- It ignores DVDS bias, forest quality in small folds, and the correlation
  between the lower and upper estimates.
- At the baseline design it reproduces the real run's median M3/M6 ratio:
  1.029 vs 1.028 (binary) and 1.033 vs 1.038 (continuous).
- In `predict_widths.R`, `min_cov_3` and `min_cov_6` are minima over 400
  pairs of 1000 draws each. Values near 0.93 are MC noise, not
  undercoverage.
