# Plan: CausaliTea reading-group talk (2026-10-12)

## Goal
30-minute Beamer talk for the CausaliTea reading group, inspired by the JICI 2026
lightning talk and `partialID_Overleaf/Notes/causalitea_plan.tex`.
Main message: when fusing two partial-ID bands, optimise the confidence interval,
not the estimated band.

## Decisions (agreed with user)
- Setting: two OS datasets with different unobserved confounding; say the method
  generalises to one OS + one RCT. Subscripts 1, 2 (matches band_toy).
- Sensitivity model: ZSB / MSM only. No NIW. Dorn, Guo & Kallus later.
- Numbers and figures: from `src/outputs/band_toy/results.csv` (current), not
  the outdated toy table in `Methods/Simulations/toy_example.tex`.
- No single-weight argument. No toy numbers before the toy section.
- Claim on width: an interior oracle weight gives a narrower expected CI than
  any fixed per-end choice of dataset; the split costs efficiency. Not claimed
  against the naive intersection (which undercovers at small gaps).
- Preamble: opacity fading (true PDF opacity), 2pt progress-bar footline,
  drop miniframes, `\band` TikZ macro, remove duplicate `\usetheme`,
  `no footline` -> `nofootline`.
- Paper fixes done: `derivations.tex` omega_1 arg max -> arg min (and removed
  the incorrect `:=`), `niw.tex` Lambda -> Gamma in the constraint.

## Outline (30 frames + title)
1. Why fuse partially identified bands? (1-3)
2. Sensitivity analysis in one dataset: ZSB (4-8)
3. Intersecting the two bands (9-12, forgetful statistician)
4. Convex combination targeting the CI (13-17)
5. Theory under Gaussian, known variance (18-22)
6. Toy simulation (23-26; figures from results.csv)
7. Discussion: why split + alternatives (winner's curse lit. incl. Inference on
   Winners, CLR; cross-fitting aggregation); next steps (formalise width gain,
   asymptotic validity for IF-based estimates e.g. Dorn-Guo-Kallus, compare
   with winner's curse lit.); future work (non-Gaussian, >2 datasets,
   Imbens-Manski); take-home.

## Pending
- Approval of `src/band_toy_slide_figs.py` (reads results.csv, writes
  single-panel PDFs with git hash to `partialID_Overleaf/Figures/`).
- Bib entries for CLR and Inference on Winners (not in references.bib).

## Files
- `partialID_Overleaf/Slides/2026-10_causalitea/2026-10_causalitea.tex`
- `partialID_Overleaf/Slides/preamble.tex`
