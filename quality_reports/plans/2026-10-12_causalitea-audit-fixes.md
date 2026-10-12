# Plan: slide-audit fixes for the CausaliTea deck (2026-10-12)

Deck: `partialID_Overleaf/Slides/2026-10_causalitea/2026-10_causalitea.tex`.
From the slide-auditor run on 2026-10-12, checked against `results.csv`.
Status: **awaiting the user's choice of which items to apply.** Edit nothing until then.

Audit found no critical errors. Proposition, Theorem, weight conditions, precision weight
and all running-example numbers match `Methods/derivations.tex` and `results.csv`.
Rejected audit findings, already checked on rendered pages: heatmap legibility (readable
at 1080p) and the table slide overflowing (it fits above the progress bar).

## Major
1. **Table slide, plug-in bullet.** "Sometimes more biased than intersection" doesn't fit
   the cells shown, where the plug-in is less biased (r=5, g=0: 0.066 vs 0.200). It is more
   biased in only 24/108 cells, but undercovers more than the intersection in 93/108.
   Fix: "Plug-in weights (no split) undercover even more than the intersection".
2. **Concavity slide.** Add "The comparison with the intersection is by simulation (next
   section)". The theory compares only with fixed per-end choices of one dataset.
3. **Bootstrap slide.** Say "asymptotically covers the whole band with probability at
   least 1−α". ZSB's percentile bootstrap guarantee is asymptotic; check the exact
   statement in ZSB.
4. **Undefined symbols.** Define B_eval(ω1) = Û_eval(ω1) + z√v_U(ω1) and v_U ("defined
   analogously"). Link [L_true, U_true] to the band [L, U]. On the separate-weights slide,
   say L_k are population ends and L̂_k their estimates. Mention A2's upper-end analogue.
5. **Heatmap bullet.** "r = 1 shown; worst over all r: intersection 89.8, split 94.4
   (MC s.e. 0.3–0.4)". Split minimum is at r=5; intersection MC s.e. is 0.4.
6. **Word count over ~40.**
   - Forgetful statistician (~75 words). E.g. "Too wide: the winner may be the noisier
     estimate"; the current "values in the CI which would not be chosen" is unclear.
   - Toy set-up (~80): move the "Compared: …" line into a speaker note.
   - Proposition (~75): merge the two proof bullets.

## Minor
7. Frame 9 title: "better noisy estimate" → "more extreme noisy estimate".
8. Separate-weights slide: add "a single weight collapses to picking one dataset"
   (`Methods/convcomb.tex`).
9. Critical-region slide: label z_{1-α/2}, not a bare z.
10. Precision-weight slide: "the gain shrinks" → "the gain can shrink" (not in the paper).
11. Worst-case tilt title: add "(upper end; the lower end flips them)".
12. Toy picture: "data 1/2" → "dataset 1/2"; add v_k = SE_k^2.
13. MSM and ZSB slides: footnote citations appear from click 1. Use `\only<n->{\footnote{\fullcite{..}}}`.
14. Merge Next steps and Future work into one slide (29 content frames for 30 min).

## Other open items
- Bib entries to add by hand: CLR (2013), Andrews–Kitagawa–McCloskey (2024), cross-fitting.
  These are only needed if the Discussion section returns; it was removed.
- After committing projectone: rerun `python -m src.band_toy_slide_figs` and rename the
  `_wip` files to the short hash. Paths in the deck and PROVENANCE.md need updating.
