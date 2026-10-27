# Plan: replace \hyperref with typed \ref (2026-10-26)

Scope: files input by `partialID_Overleaf/compiled.tex` only. `Outdated/` and
`Notes/Sensitivity_Analysis/` left unchanged (user decision).

1. `MSM/equivalence.tex:10` — Proof of Theorem~\ref{thm:msm_equi}
2. `MSM/equivalence.tex:32` — (see Section~\ref{sec:zsb})
3. `PartialIDConditional/dg_ipw.tex:5` — the ZSB bounds (Section~\ref{sec:zsb})
4. `Methods/derivations.tex:3` — Section~\ref{sec:int_and_hyp}
5. `Methods/derivations.tex:169` — Theorem~\ref{thm:band_validity}
6. `Methods/derivations.tex:173` — Theorem~\ref{thm:band_validity}
7. `Methods/convcomb.tex:45` (comment) — old (2) -> Equation~\ref{eq:oracle-objective},
   old (3) -> Equation~\ref{eq:oracle-weight}; new label added to the
   $\omega_0^\star$ argmax equation in `derivations.tex`.

Status: done; compiled cleanly, no undefined references.
