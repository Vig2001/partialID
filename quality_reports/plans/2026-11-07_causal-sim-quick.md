# Plan: `R/causal_sim_quick.R`, the causal simulation with lean nuisances and every method

Date: 2026-11-07. Status: awaiting approval. No earlier plan was saved for this. The idea appears only as the "parametric-nuisance option" in `done/2026-11-02_dgk-exact-replication.md` and the session logs.

## Why

The forest version (`R/causal_sim.R`) takes 2 h at n = (2000, 500) and took 9 h 45 min at (3000, 500). That leaves no room for CLR (M8, M9), cross-fitting (M7) or trying several designs before 16 Nov.

The lean script uses the same two-study DGK design and the same CI comparison, with parametric nuisances fitted by logistic and linear regression. With those, fitting costs milliseconds, so every method in the toy can run in minutes per configuration.

## Design

1. **Data:** unchanged. Use `simulate_study(n, noise_scale, binary_denominator)` from `R/helpers/data_fusion_design.R`, so `--noise-scale-2` and `--binary-denominator-1` work as in the main script. The true bands are the same functions.

2. **Nuisances** go in a new `R/helpers/dvds_lean.R`, fitted per fit (full sample, fold or cross-fit fold).
   - **e(x):** logistic regression on X1, 1(X2 > 0) and X2 X3. This is the true functional form, so e is correctly specified.
   - **Binary outcome, mu_z(x):** logistic regression per arm on X1, X2 and X2 X3. This is exact for k = 1. For k > 1, the true logit is -log(k - 1 + exp(eta)), so logistic regression is slightly misspecified. Because e is correct, DGK's double validity still makes the bounds valid, possibly conservative. Then call DGK's `dvds.binary()` verbatim: it only takes mu0, mu1 and e.
   - **Continuous outcome:** a Gaussian location-scale fit per arm.
     - m_z(x) by `lm` on sign(X1), X2 and X2 X3.
     - sigma(x) by `lm` of |residual| on 1 and X4^2, times sqrt(pi / 2).
     - Then Q_tau = m + sigma Phi^{-1}(tau) and CVaR_tau = m + sigma phi(Phi^{-1}(tau)) / (1 - tau).
     - These go into DGK's influence function (the `phi.continuous` formula, reimplemented with parametric Q and CVaR in place of their forests). This is exact for the DGK DGP.
   - **Why the true functional form:** the simulation tests the fusion layer, not nuisance learning. It also removes the single-study undercoverage of the forest DVDS at Lambda = 1 (about 0.91), which made M6's validity check fail there.
   - **Alternative:** `--nuisance-features linear` would use main effects X1..X5 only, as a misspecified variant.

3. **Fits per study and replicate:**
   - full sample;
   - fold A and fold B (M3);
   - for M7 (K = 5, as in `src/toy_sim.py::cross_fit_ci`): each fold k is a held-out fit, and its weights come from a fit on the other K - 1 folds, with variances scaled to the full n.

   That is 3 + 2K = 13 fits per study, 26 per replicate. With n_2 = 500, a held-out fold has 100 units (about 44 treated). This is at `MIN_FOLD_SIZE`, so I'll check it.

4. **Methods:** all the toy's, except the oracle M2, which needs the true standard errors of the bounds.
   - M1 naive intersection;
   - M3 split;
   - M4 no split (plug-in weights; expected to undercover);
   - M5 dataset 1 and dataset 2;
   - M6 calibrated intersection;
   - M7 cross-fit;
   - M8 CLR intersection;
   - M9 convex CLR.

   M8 and M9 use p = 1 - alpha / 2 per end, with no inequality selection, as in the toy.

5. **R ports in `R/helpers/ci_methods.R`:** `cross_fit_convex_ci`, `clr_intersection_ci` and `clr_convex_ci`, from the Python, with separate lower and upper variances.
   - The quadrant-sup quantile is found with `uniroot`.
   - Each port is checked against the Python function on the same random inputs: a scratch script writes the inputs, and the outputs must agree to about 1e-12.

6. **Structure:** `R/causal_sim_quick.R` sources `R/causal_sim.R`, whose `main()` only runs when it is called directly. It reuses the flags, seeds, parallel runner, summaries, checks and output, and replaces only:
   - the fitting step (`simulate_replicate`, `FITS`);
   - `METHODS`;
   - `pair_cis`.

   This leaves the main script untouched before the deadline. TODO, after submission: move the shared parts into `R/helpers/`. The default output directory is `src/outputs/causal_sim_quick/`.
   - New flags: `--folds` (K, default 5) and `--nuisance-features` (`true` or `linear`).

7. **Faster summaries:** the summary stage runs in parallel over Lambda_1, with `future_lapply` and the same backend. With 9 methods, it would otherwise dominate the run time (see below).

## Checks

1. The R ports match the Python (step 5).
2. At Lambda = 1 the lean continuous and binary bounds collapse to AIPW, so L = U.
3. Lean and forest DVDS on the same 2 replicates give close bounds: differences within about one s.e.
4. The existing sanity checks pass, plus three new ones:
   - single-study coverage is at least 0.95 minus 3 MC s.e. at every Lambda;
   - M8 and M9 are valid in every pair;
   - M1's coverage near ties is printed.
5. Smoke run with 2 replicates; then 100 replicates to measure the run time before any 1000-rep run.

## Estimated run time (measured pieces, n = (3000, 500))

- **Lean nuisance fit** (e plus both arms, both outcomes): about 5 ms.
- **DGK's `dvds.binary()`** over 61 Lambdas: about 54 ms per fit. The lean continuous step should take a similar time (vectorised over Lambda).
- **So about 0.1 s per fit**, or about 2.6 s per replicate (26 fits) on one core.
  - **Fitting, 1000 replicates:** about 45 min on one core; about 5-6 min on 9 workers.
- **Summary stage:** about 18 weight-grid evaluations per pair (M3 2, M4 2, M7 10, M9 about 4), at 21 ms each (R = 1000 replicates, 1001 weights). That is about 0.4 s per pair:
  - about 11 min per outcome on the 41 x 41 grid single-threaded;
  - about 1.5 min per outcome in parallel on 9 workers.
- **Total:** about 10 min for 1000 replicates of both outcomes on 9 workers, about 30 min without the parallel summary. The forest version takes 2-10 h for the same design.

## Out of scope

- Changing `R/causal_sim.R`'s results.
- The oracle M2.
- Paper tables; those come later through `causal_sim_tables`, after you choose a design.
