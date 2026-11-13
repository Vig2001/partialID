# Confidence Intervals for Data-Fusion-Based Partial Identification

Point identification (ID) of a causal estimand relies on strong, untestable assumptions, such as no unobserved confounding in an observational study (OS) or study exchangeability of potential outcomes (POs). Relaxing such an assumption gives partial ID: a band that contains the causal quantity instead of a point estimate. Most such bands relax one assumption in one dataset (OS or RCT). \cite{lanners_data_2025} narrow the identification band by intersecting the bands from an OS and an RCT. This project extends that work by minimising the width of the confidence interval for (a) the fused ID band and (b) the causal estimand. We motivate the method with simulations, then apply it to sensitivity regimes.

## Core principles

- **Plan first:** for non-trivial tasks, present a plan and wait for my approval. Save approved plans to `quality_reports/plans/YYYY-MM-DD_short-name.md`; once enacted, move them to `quality_reports/plans/done/` and update references to them.
- **Verify after:** after changing a `.tex` file, compile it and check the log for errors, undefined references and missing citations. After changing code, run it and check the output before saying it is done.
- **The paper is a separate repo:** `partialID_Overleaf/` syncs with Overleaf. Commit and push paper changes from inside that folder; never add it to the projectone repo. The paper is `partialID_Overleaf/compiled.tex`, which inputs the sections from their folders.
- **Everything the paper uses lives inside `partialID_Overleaf/`** (figures, tables, `.bib`), or Overleaf can't see it.
- **Numbers come from code:** never type a result into the paper by hand; it must trace to code in `src/` and the file that code wrote.
- **Learn from corrections:** when I correct you, add a one-line entry under Corrections.
- **Log before leaving:** when I say let's wrap up or summarise the day, run /session-log first.
- **Coding:** follow @AGENTS.md.

## Folder structure

- `src/`: the main codebase (e.g. `toy_sim.py`).
    - `helpers/`: functions; `plotting/`: plotting functions; `latex/`: scripts that write paper and slide figures and tables
    - `outputs/`: results written by scripts, one folder per script
- `R/` (helpers in `R/helpers/`): the causal simulations; `tests/`: pytest tests; `checks/`: one-off numerical checks
- `scripts_misc/`: old rough working. Not the source of truth: don't edit, import from or review it unless I ask
- `project1venv/`: virtual environment (never edit)
- `partialID_Overleaf/`: the paper (`compiled.tex`)
- `partialID_Overleaf/Slides/`: Beamer talks, one folder per talk (`Slides/YYYY-MM_<name>/YYYY-MM_<name>.tex`), using the shared UCL theme (`beamerthemeUCL.sty`) and `\input{../preamble.tex}` unless told otherwise
- `quality_reports/`: plans, session logs, reviews

settings.json blocks reading `private_data/` (for future data-use agreements); open-source data goes in a separate `public_data/` folder.

## Commands

- Compile paper: `cd partialID_Overleaf && latexmk -pdf compiled.tex` (`-xelatex` if Overleaf uses XeLaTeX)
- Python: `source project1venv/bin/activate`; dependencies in `requirements.txt`
- Run analysis: `python -m src.[script]` from the repo root (`src` is installed editable)
- Tests: `pytest` from the repo root
- R: `Rscript R/causal_sim.R [--reps ...]`, `Rscript R/causal_sim_quick.R ...` or `Rscript R/dgk_replication.R` from the repo root. DGK's repo must be cloned in `external/` (see README). `.Rprofile` activates renv; restore packages with `Rscript -e 'renv::restore()'`
- Bibliography: `partialID_Overleaf/references.bib`

## Paper conventions

- Citations: `natbib` with `plainnat` (author–year); `\citet{}` in text, `\citep{}` in parentheses, never bare `\cite{}`
- Notation: `partialID_Overleaf/Methodology/notation_setup.tex` (Section `sec:notation`, kept lean): estimand $\tau$ = ATE, two independent datasets, bands $[L_k, U_k]$, true band, gaps and ties, weights $\omega_L$, $\omega_U$. The CVaR quantile level in the DGK appendix is $\eta$.
- Spelling: British English

## Current state

### Method and framing
- Convex combination method: stress-testing it and proving theorems. Brieuc's shorter proof of Prop 4.2 (`prop42_proof.pdf`) is under review; notes in `quality_reports/reviews/2026-11-02_prop42-proof-notes.md` (open: keep Remark 2 or not).
- Abstract draft (2026-11-06, `Introduction/abstract.tex`): two independent datasets with different ID assumptions; the estimated intersection is biased inwards, so naive inference undercovers and the calibrated version is conservative. Contributions:
    1. convex fusion that converges to the intersection of population bands and is asymptotically valid under sample splitting;
    2. the problem as a special case of CLR (2013), extended to pool variance across estimators (convex CLR);
    3. toy and causal simulations;
    4. K-fold cross-fitting as an extension of the split.

### Method codes
The code labels methods M1–M14. Use the names in chat and in the paper.

| Code | Name | Notes |
| --- | --- | --- |
| M1 | naive intersection | |
| M2 | convex, oracle | |
| M3 | convex, split | |
| M4 | convex, no split | called "plug-in" before 2026-11-02 |
| M5 | dataset 1 / dataset 2 alone | "large OS" / "small OS" in the causal sims |
| M6 | calibrated intersection | Šidák z = 2.239 (below) |
| M7 | convex, cross-fit | K = 5 |
| M8 | CLR intersection | no selection; k = 2.239 |
| M9 | convex CLR | no selection; k = 2.353 |
| M10 | split intersection | paper: "intersection CIs, split" |
| M11 | split band intersection | |
| M12 | CLR intersection with selection | adaptive inequality selection |
| M13 | convex CLR with selection | tube formula for d = 1 |
| M14 | intersection of calibrated CIs | equals M8 |

### Toy simulation
`python -m src.toy_sim`, about 1 min 40 s; outputs in `src/outputs/toy_sim/`. Known true band [1, 3]; r in {1, ..., 5}, gaps g in {0, 0.05, 0.1, 0.2, 0.5, 1}, 180 cells.

- **Cross-fit** (2026-11-03, plan `done/2026-11-03_toy-cross-fit.md`; K = 2 tried and dropped): M7 is narrower than M3 everywhere, but undercovers at exact ties (min joint coverage 0.925).
- **Early findings:**
    - Dataset 1 alone is always valid, and beats M3 wherever M3 beats M6.
    - The oracle M2 is at most 2–10% narrower than the best valid comparator.
    - The naive intersection undercovers at small gaps, so width comparisons with it are not like-for-like.
- **Estimated band ends** (`--metric ends`; bold `\inward{}` = end inside [1, 3] by more than 2 MC s.e.):
    - M6 is pulled in at both ends by similar amounts (up to [1.200, 2.793] at r = 5, g = 0).
    - M8 and M9 are pulled in at both ends near ties for r ≥ 2.
    - M3 and M7 are never pulled in; at large r they are conservative at the upper end, where the noisy dataset 2 sets U.
- **CLR intersection, M8** (plans `done/2026-11-05_toy-clr.md`, `done/2026-11-06_clr-no-selection.md`):
    - Setup: p = 1 − α/2 per end, estimates at p = 1/2.
    - No selection step since 2026-11-06: I found the selection-dependent tube formulas hard to follow, and per-replicate simulation was too slow. So k = Φ⁻¹(√p) = 2.239.
    - It equals M6 at r = 1 and is narrower at r > 1, because it corrects each bound before taking the max.
    - Valid in all 180 cells (min 0.947). M3 is wider in every cell (by 0.8–11%).
    - M7 is narrower than M8, by at most 3.1%, in the cells where both cover: 25/25 at r = 1, 19/28 at r = 2, 9/33 at r = 3, 1/36 at r = 4. Elsewhere M7 is up to 5% wider.
    - Its half-median-unbiased band is narrower than the truth near ties at large r (1.883 at r = 5, g = 0).
- **Convex CLR, M9** (plan `done/2026-11-05_toy-convex-clr.md`; `ci_methods.clr_convex_ci`): CLR over all weights w in [0, 1], i.e. a convex CI with no split.
    - Critical value: k is the quantile of the sup over weights. The sup is ‖Z‖ if both errors are positive, else max(Z₁, Z₂), so P(sup ≤ k) = Φ(k) − exp(−k²/2)/4. This gives k = 2.353 at p = 0.975 and 0.563 at 1/2.
    - Plugging M4's weight into M8's k undercovers at ties (0.937).
    - Valid in all 180 cells (min 0.944).
    - Against M8: narrower in 138/180 cells by up to 1.5%, wider by up to 1.3%. In excess width it is up to 11% narrower at r = 1 and 2.5% narrower at r = 5, and up to 5% wider elsewhere.
    - Against M3: width ratio 0.895–1.001.
    - Width over the per-cell narrowest method with validity theory, median / max: M9 1.003 / 1.022, M8 1.010 / 1.019, dataset 1 1.009 / 1.423, M6 1.064 / 1.228, M3 1.072 / 1.130.
- **Split intersection** (2026-11-11, plan `done/2026-11-11_split-selection-comparator.md`):
    - M10 is M3 with weights restricted to {0, 1}. It intersects the fold-A CIs (variances at fold-B size) and reports the winning dataset's fold-B CI at each end. It answers the referee question "why not select one dataset on the training fold?"
    - M11 (2026-11-12) selects by estimated bound alone (z = 0, `BAND_INTERSECTION_Z`).
    - Both are valid in all 180 cells (min 0.946 and 0.944).
    - M10 / M3 width: median by r 1.013–1.018 (range 0.996–1.029). At equal gaps M10 is narrower only at (r, g) = (1, 0.5) and (2, 1).
    - M11 equals M10 at r = 1. M11 / M3 median by r: 1.013, 1.049, 1.108, 1.170, 1.229 (max 1.267). It is narrower than both only at g = 1, r = 2, 3.
    - Oracle weights never lose to selection, but estimated weights can: training noise pulls the weight off the corner.
- **Šidák** (2026-11-12, supervisor): M6 uses z = Φ⁻¹(√(1 − α/2)) = 2.2390 at each end, which is Šidák over the two independent datasets and Bonferroni over the two ends. It was Φ⁻¹(1 − α/4) = 2.2414 before.
    - Code: `ToyConfig.z_calibrated` and `make_config()` in `R/causal_sim.R`.
    - Only M6 changed (at most 0.16% narrower); no table row changed.
- **Convergence** (2026-11-04, plan `done/2026-11-04_toy-bound-convergence.md`): `python -m src.toy_convergence` plots estimated bands against n₁ (100 to 1e5, n₂ = n₁/4, r = 2) at fixed population gaps g.
    - Split: outward bias (1 − E[w])g, which vanishes as the weights reach the corners (more slowly for small g), and none at a tie. Its upper-end bias / SD peaks at about 0.65 for every g, so convergence is not uniform in g.
    - Intersection: inward bias of about +0.68 SD at a tie, at every n.
    - `uniform_convergence.pdf` is kept for a possible extension on uniform convergence.

### Causal simulation (R)
- **Scripts:**
    - `R/causal_sim.R` runs DGK's own code (`R/helpers/dgk.R`): their DGP and nuisance lines 23–42 verbatim, with `dvds.binary()` and `dvds.continuous()` sourced from the clone in `external/` (commit `26cbbf6`; gitignored; it has no licence, so never redistribute it). It reproduces their replicate 1 to 5e-15.
    - `R/causal_sim.R` flags:
        - `--noise-scale-2`: scales study 2's continuous noise, and so its true band, by c.
        - `--lambdas-1` / `--lambdas-2`: a grid per study, `from:to:by` allowed.
        - `--binary-denominator-1`: k for study 1's binary outcome, P(Y = 1 | x) = 1/(k + exp(η)), drawn by thinning; the default 1 is DGK's outcome.
        - `--outcomes`: which outcomes to run.
        - The defaults reproduce DGK's design exactly.
    - `R/causal_sim_quick.R` + `R/helpers/dvds_lean.R` (2026-11-08, plan `done/2026-11-07_causal-sim-quick.md`, notes `causal_sim_quick_explained.tex`) run the paper's designs.
        - Nuisances are parametric models of the true form, predicted out of fold (`--nuisance-folds 5`). In-sample fits biased the bounds inwards and broke cross-fitting; those runs are kept in `*_insample/`.
        - `--nuisances forest` (binary only) uses DGK's forests.
        - Run time: 5–10 min per design; the forest run took 3 h 21 min on 9 workers.
    - Cross-fitting (M7) is standard since 2026-11-09: models fitted on the other folds predict the held-out fold. The fold-alone runs are kept in `*_foldalone/`.
    - `--summarise-only` resummarises without refitting; in the sandbox add `--workers 1 --backend sequential`.
    - `summary.csv` is short; `summary_full.csv` has everything.
    - The R ports match the Python methods (M7–M9 to 2e-13, M12–M14 to 1e-14).
- **Paper designs:** n = (3000, 500), Λ step 0.05, 1000 reps.
    - Continuous: c = 2.5, Λ₁ in [3, 5], Λ₂ in [1, 3].
    - Binary: study 1's rare outcome k = 15 (prevalence 0.0616), Λ₁ in [1, 4], Λ₂ in [1, 2].
- **Runs** (in `src/outputs/`):
    - Main text, parametric nuisances: `causal_sim_quick_continuous`, `causal_sim_quick_binary_k15`.
    - Appendix, forest nuisances: `causal_sim_forest_binary_k15`, and `causal_sim_noise_n3000_reps1000` (continuous).
- **Findings:**
    - With DGK's design unchanged, the convex split is never narrower than the calibrated intersection (SE ratio ≈ 2 is too small to beat the fold-B cost). DVDS undercovers at Λ = 1 (≈ 0.91). Hence the scaled noise and the rare outcome.
    - Continuous, c = 2.5: the split is narrower than the calibrated intersection in 304/1681 cells (best ratio 0.963). The bands are always nested, and dataset 1 alone is narrower than the split in every such cell.
    - Binary, k = 15: study 1's DVDS undercovers its own band from Λ₁ ≈ 1.4 (0.92 at Λ₁ = 4).
        - Cause (checked 2026-11-09): the estimates are unbiased and the SEs right on average, but estimate and SE are correlated (−0.84 lower, +0.86 upper at Λ₁ = 3.5), so each end misses 3.5–4.4%.
        - Convex CIs inherit this where their weight is on study 1 (Λ₂ > 1.5). M6, M8 and M9 survive through their 2.24 critical value.
        - About 40% of binary cells sit at the validity cut-off, so valid counts move by tens of cells under a fold redraw.
    - Valid cells:
        - Continuous parametric: all methods in all cells (M7 min 0.943).
        - Binary parametric, of 1281: dataset 1 426, split 725, cross-fit 605, split intersection 366.
        - Binary forest, of 1281: dataset 1 537, split 559, cross-fit 588, split intersection 429, M8 1279; M6 and M9 all.
    - Width over the best valid CI at near ties (2026-11-08), median:
        - Parametric: CLR methods ≤ 1.016; M7 1.028 (continuous) and 1.066 (binary); split 1.067 and 1.129; calibrated 1.118–1.168.
        - Forest: M9 1.000, M8 1.018, M7 1.057, M3 1.142, M6 1.157.
    - The convex split is never narrower than the CLR intersection (0/1681, 0/1281).
    - Excess-width medians over valid cells (2026-11-13; one-off check, not yet written by code):
        - The convex CLR is 4.5% wider than CLR in both designs. It is narrower in only 21% (continuous) and 27% (binary) of cells, all within 2 pooled SE of a tie.
        - CLR is 18% narrower than the convex split, in every cell.
    - Selection (M12, M13) is valid everywhere, except the CLR intersection with selection in 3 parametric binary, 2 forest binary and 123 continuous forest cells. The 123 have Λ₂ in [1, 1.1]: study 1 is dropped, and study 2's own forest DVDS CI undercovers at n = 500.
- **Other:**
    - `R/dgk_replication.R` reruns DGK's simulation unmodified (~40 min); not run yet.
    - DGK's DGP and closed-form band: `Results/Causal Simulations/causal_sim_explained.tex`.
    - Next, with no plan saved yet: larger n₁ and a finer grid.
    - renv at the repo root (cache disabled, library in `renv/library/`).
    - `R/stata/` (CLR via RStata) not created yet; a draft Python CLR is in `quality_reports/plans/2026-10-30_clr-implementation-draft.md`.
    - Old ZSB/NIW code is in `scripts_misc/outdated_help/`.

### Comparator sets and selection (2026-11-13)
Plan: `quality_reports/plans/2026-11-13_comparator-sets.md`. Stage 1 is done. Stage 2, CLR with adaptive inequality selection, is implemented and in the paper.

- **Three sets** in every width table:
    1. Sample splitting: split intersection, convex split, convex cross-fit.
    2. Simultaneous inference: intersection of calibrated CIs and CLR, as two identical columns; convex CLR.
    3. Adaptive inequality selection (header on two lines): CLR, convex CLR.
- **Selection methods:**
    - M13's critical value is Hotelling's tube formula for d = 1, which is exact here (CLR 2013, p. 698, Remark 2).
    - `checks/clr_selection_check.py` compares the analytical critical values with CLR's simulated ones in one causal cell: they agree.
    - In the toy, selection matters only at gaps ≥ 0.5.
- **Excess width:** tables report mean CI width minus the true band width. The causal tables add a true-band-width column.
- **Causal table rows:** `causal_sim_tables --row-rule convex-wins`. For each Λ₁, they are the (up to) three Λ₂ where the convex split beats the split intersection and the convex CLR beats CLR, by the largest margin, with all four valid. The text states this rule.
- **Removed:** the calibrated intersection is out of every table and figure (still computed and checked in code). The naive intersection is only in `fig_undercoverage_mt.pdf`.
- **Open:** whether to keep the causal Results opening remark about DGK's original design.

### Where things are in the paper
- **Toy, main text:**
    - `tab:toy-width-coverage` (three sets, excess width).
    - `fig:toy-coverage-by-gap` (`fig_coverage_by_gap_mt.pdf`, from `src.latex.toy_sim_paper_figs`).
- **Toy, appendix (`app:sim`):**
    - Oracle paragraph (`app:toy-oracle`).
    - `tab:toy-split-intersection`: split band intersection, convex split, oracle convex (`toy_sim_slide_figs --band-oracle`).
- **Appendix "Estimated Bands" (`app:estimated-bands`):**
    - Band-ends table `tab:toy-band-width`, with the calibrated column labelled "estimated bands intersected".
    - Convergence figure `fig_convergence_app.pdf`.
- **Causal, main text:**
    - `tab:causal-continuous-param` and `tab:causal-binary-param` (eleven columns, scaled with `\resizebox`).
    - `fig:causal-width-ratio`: one figure with continuous and binary panels side by side (`fig_width_mt.pdf`, from `src.latex.causal_sim_box_panels`, which reads the `width_ratio_cells_*.csv` files written by `causal_sim_paper_figs`).
- **Causal, appendix (`app:sim-forest`):** binary forest figure and table, continuous forest table.
- **Files:**
    - Tables: `Tables/{toy_sim,causal_sim}/{total_width,excess_width}/`, with the same names in both; the paper inputs `excess_width`.
    - Box plots: `Figures/causal_sim/total_width/`.
    - Commands and history: `Figures/{toy_sim,causal_sim}/PROVENANCE.md`.
    - PDFs are gitignored in the paper repo, so upload them to Overleaf by hand.
    - Slide versions of the toy figures and tables come from `python -m src.latex.toy_sim_slide_figs`.

### Code
- Code follows AGENTS.md.
    - CIs: `src/helpers/ci_methods.py`.
    - MC summaries: `src/helpers/mc_metrics.py`.
    - Plots: `src/plotting/`.
    - Paper and slide outputs: `src/latex/`.
- `tests/` has 49 tests, all passing. pytest 9.2.0 is installed but not yet in `requirements.txt`.

## Skills

| Command | What it does |
| --- | --- |
| `/compile-paper [file]` | Compile the paper (or a deck) and report errors and warnings |
| `/proofread [file]` | Grammar, typos, consistency report |
| `/review-paper` | Proofreader + methods reviewer + citation checker + referee, verified and combined (only runs when typed) |
| `/sync-overleaf` | Guided pull/commit/push for the paper |
| `/session-log` | Save what we did and what is next |
| `/simulation-study [script]` | Run a simulation and check coverage, bias, width |
| `/next-steps` | Prioritised next steps from logs, To_Do and pending plans |
| `/challenge-idea [idea or plan]` | Stress-test an idea, or review a plan in `quality_reports/plans/`, before investing time |
| `/make-slides [topic]` | Plan, build, compile and audit a Beamer talk |
| `/write-tests [module]` | Propose tests for a `src/` module, then write them after approval |
| `/update-paper-figs [script]` | Regenerate a figure or table, copy it into `partialID_Overleaf/`, record provenance |

Agents (in `.claude/agents/`): code-reviewer, proofreader, methods-reviewer, citation-checker, counterexample-hunter, referee, slide-auditor.

## Corrections

<!-- Claude adds entries here, e.g. "Use \citet not \cite for in-text citations" -->
- Compare widths only between methods with valid coverage; never call a valid CI "wider" than one that undercovers.
- Each new talk gets its own folder inside `Slides/`.
- Respect each section's assumptions: in the known-variance Gaussian section, never suggest estimating $\sigma_k^2$.
- Write derivations in the paper as connected prose between displays, never with bold "Step 1", "Result." headings.
- When variances are estimated, don't contrast with the known-variance section; say validity relies on the variance estimators being consistent.
- Table captions above, figure captions below (`compiled.tex` sets `\captionsetup[table]{position=top}`). Keep captions lean: shared conventions (nominal level, replications, coverage and bold, folds, selection, row rule) go in the simulation's setup text; MC s.e. as one "at most" clause from the `...Max` macros.
- Don't rank methods by how often each is narrowest per cell (that uses hindsight). Compare fixed rules by median and worst-case width relative to the per-cell best.
- Cross-fitting with nuisances means standard cross-fitting: fit on the other folds, predict on the held-out fold. Don't port the toy's fold-alone fits.
- Related work: the closest papers are Lanners et al. (2025), Chernozhukov, Lee & Rosen (2013) and Yu et al. (2026; JRSS-B vol. 88, published version of arXiv 2412.03731). Han & McCloskey (2024) targets a parameter selected from data, so it is not close. Keep econometric references minimal; the contrast is that they don't combine data sources for partial ID.
- Since 2026-11-11 (supervisor): no single-dataset comparators in any paper table or figure. Causal box plots show CI width / true band width on a log scale; the binary panel uses Λ₁ in [2, 4], Λ₂ in (1, 2], stated in the panel title. Plan `done/2026-11-11_drop-single-dataset-truth-ratio.md`.
- Paper input files are named `fig_|tab_<content>_<mt|app>` (main text or appendix), with no hash. The hash, command and old name go in PROVENANCE.md. This replaces the hash-in-name rule of /update-paper-figs.
- In the paper use the Methods notation, not code names: weights $\omega_L$, $\omega_U$; oracle $\omega^{\ast}$ (`\ast`; estimate $\hat\omega^{\ast}$); population bounds $L_k, U_k$.
