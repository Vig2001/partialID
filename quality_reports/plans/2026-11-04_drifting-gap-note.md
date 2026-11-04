# Note: fixed versus drifting gaps in the convergence toy

Date: 2026-11-04. Status: **note, no action needed.** Background for the
convergence figure (plan `done/2026-11-04_toy-bound-convergence.md`) and for a
possible extension on uniform convergence.

## The gap is fixed in any real problem
The gap g is the distance between the two datasets' population bounds
(g_L = L_1 - L_2, g_U = U_1 - U_2; in the equal-gap toy also
mu_1 - mu_2 = g). For a given data-generating process it is fixed. Only the
estimates change with n: the estimated gap g_hat = L_hat_1 - L_hat_2 equals
g plus noise with SD sqrt(SE_1^2 + SE_2^2), so it concentrates on g as n
grows. g is fixed but unknown to the analyst; if it were known, the weights
would not need estimating.

## What a drifting gap is for
A gap that shrinks with n, g_n = h / sqrt(n) (equivalently g_n = h SE(n)),
is a mathematical device, not a model of data. It asks how a method behaves
when, at the sample size at hand, the gap is comparable to the standard
error.

- With g fixed and n -> infinity, g / SE -> infinity, so fixed-g
  (pointwise) asymptotics always describe the easy, separated case. If at
  the actual n the gap is about one SE, that limit misdescribes the
  situation.
- Holding h = g / SE fixed as n grows keeps the difficulty of the actual
  sample, so the limit is a better approximation of finite-sample
  behaviour.
- Analogy: local (Pitman) alternatives in power analysis. At a fixed
  alternative the power of a consistent test tends to 1, which is
  uninformative; local alternatives theta_0 + h / sqrt(n) give useful
  approximations. Nobody believes the parameter moves with n.
- In partial identification the same device is used to study whether CIs
  cover uniformly over the parameter space, including near ties, which is
  where intersection-type methods fail. Pointers to check against the
  library before citing: Imbens and Manski (2004), Stoye (2009), Andrews
  and Soares (2010).

## Two different questions
1. **Consistency**: do the estimated bounds converge to the true bounds?
   A statement about one data-generating process at a time, so fixed g is
   the right design. This is what the convergence figure in the toy
   section shows: for g > 0 the weights reach the corners and the outward
   bias (1 - E[w_L]) g vanishes; for g = 0 the estimates are unbiased at
   every n.
2. **Uniform behaviour over g**: at a given n, how badly can the method do
   over all gaps the analyst cannot rule out? This needs either a drifting
   gap or many fixed gaps compared at each n. It matters for the validity
   of inference (coverage), not for consistency.

## What was tried and what was kept (2026-11-04)
- First version: a "near tie" panel with g_n = 2 SE_1(n). Dropped: a
  population gap does not change with n, and the panel answered question 2
  in a figure about question 1.
- Second version: a separate figure (`uniform_convergence.pdf`) of the
  standardised outward bias, bias / SD, at fixed small gaps g in
  {0, 0.01, 0.02, 0.05, 0.1, 0.2}. It shows non-uniformity with fixed gaps
  only: for the convex split's upper end, each fixed g rises to a peak of
  about 0.65 SD and then falls, with the peak at larger n for smaller g;
  the intersection is biased inwards, about -0.68 SD at a tie at every n.
  Taken out of the notes as confusing for a consistency point; still
  written by `python -m src.toy_convergence`, with the copy
  `Figures/toy_sim/uniform_convergence_wip.pdf` marked unused in
  `Figures/toy_sim/PROVENANCE.md`.

## Possible extension
A section on uniform validity: show coverage (not only bias) uniformly over
g, link the intersection's inward bias near ties to its undercoverage and
the split's outward bias to its conservativeness, and compare with an
adaptive critical value such as Chernozhukov, Lee and Rosen (2013).
`uniform_convergence.pdf` and the coverage heatmaps are the starting
material.

## Lesson
Match the asymptotic framing to the claim: fixed g for consistency;
drifting or many fixed gaps only for uniformity. Raise the second as a
separate, optional point rather than building it into a figure meant for
the first.
