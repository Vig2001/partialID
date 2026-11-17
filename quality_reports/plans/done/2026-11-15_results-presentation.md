# Results presentation changes (2026-11-15)

Goal: change tables and figures so the CLeaR Overleaf project needs only the
`Tables/` and `Figures/` folders uploaded, plus a one-time replacement of each
affected `tabular` block by an `\input` line and one new figure environment.

Decisions (user, 2026-11-15):
- Tables become self-contained: each file in `Tables/` holds the whole
  `tabular` (column spec, header, rules, body).
- Legend names for the selection CIs: "CLR (AIS)" and "convex, CLR (AIS)";
  the simultaneous CLR is "CLR (Šidák)" in every table and figure.
- Excess-width figure: hollow markers where coverage is more than 2 MC s.e.
  below nominal.
- The band-ends appendix table keeps its "estimated bands intersected" column.
- Main prose is not changed; captions only if needed, kept as lean.

Steps:
1. `src.toy_sim --gaps`: rerun with g = 2 added into `src/outputs/toy_sim_g2/`;
   check the 180 existing cells equal the current `results.csv`.
2. Toy Table 1 (excess and total width): drop the calibrated column, rename
   CLR, add g = 2 rows; self-contained.
3. Figure 2 (`fig_coverage_by_gap_mt.pdf`): add CLR (AIS), convex CLR (AIS)
   and g = 2; same file name.
4. New `fig_excess_width_by_gap_mt.pdf`: same layout for excess width.
5. Causal tables (two main text, two forest appendix): drop the calibrated
   column, rename CLR; self-contained; rows unchanged.
6. Causal box plots (main and forest): add the AIS CIs, rename CLR; same file
   names; report if a caption validity entry is needed.
7. Copy into `partialID_Overleaf/` under the same names, make the matching
   Results/Appendix edits there, update both PROVENANCE.md files, run tests,
   compile.

Flagged, not changed: toy_example.tex line 9 says the tables report both the
calibrated intersection and CLR.
