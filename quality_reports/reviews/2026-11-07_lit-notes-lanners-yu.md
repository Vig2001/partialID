# Literature notes: Lanners et al. (2025) and Yu et al. (2025)

Written 2026-11-07 from the copies in Zotero, for the literature review
(`partialID_Overleaf/Literature Review/df_sensitivity.tex`). Page numbers are
PDF pages of those copies. Nothing here has been put into the paper.

---

## 1. Lanners et al. (2025), *Data Fusion for Partial Identification of Causal Effects*

- Bib key `lanners_data_2025`; Zotero `NQLKB42F` (a duplicate, `QKIABS55`, also exists).

### What they do

- **Setting.** An RCT and an OS, each with its own sensitivity parameter. The
  fused bound combines two components:
  - `v`, estimated from the experimental samples, with parameter γ;
  - `w`, estimated from the observational samples, with parameter ρ.

  The assignment of γ and ρ follows the influence-function terms on p. 7.
- **The problem.** The sharp bounds contain a max and a min, which are not
  pathwise differentiable, so efficient influence functions (EIFs) cannot be
  derived directly (Section 5, p. 6).
- **The fix: the Boltzmann operator** (Section 5.1, p. 6). For two values, it is
  λ₁x₁ + λ₂x₂ with λ₁ = e^{αx₁}/(e^{αx₁} + e^{αx₂}) and λ₂ = 1 − λ₁. It tends to
  max(x₁, x₂) as α → ∞ and to min(x₁, x₂) as α → −∞.
- **Applied pointwise in X, inside the expectation** (Lemma 2, p. 6). The bound
  for each covariate value is g₁(X)μ(X,1,t) + g₀(X){λ₁v + λ₂w}. For any α > 0
  (under a stated condition), the smoothed bounds are valid bounds on τ, and they
  converge to the sharp bounds as α → ∞.
- **Inference** (p. 7). They give an EIF for the smoothed functional, a
  cross-fitted bias-corrected estimator, and asymptotic normality for fixed α
  (Theorem 2). The data dependence of the weight appears in the EIF through the
  factor {λ₁ + αλ₁λ₂(v − w)}.
- **Choosing α** (Appendix E.1, p. 52). α is user-chosen.
  - Larger α approximates the max and min more closely.
  - "Especially in small sample sizes", large α can make the bias-correction
    term unstable where the max or min switches between its arguments, because
    the operator has steep gradients there.
  - They use α = 10 as "a moderate value" and say it "may be increased in larger
    datasets". The simulation settings on p. 56 also use α = 10.
  - This is a remark: I did not see an experiment that illustrates the
    instability, so say they *note* it rather than *illustrate* it.

### How it relates to our work: it is a choice of ω

λ is a data-dependent convex weight, so the Boltzmann operator is a particular
rule for choosing our weights. It differs from ours in three ways:

1. **Purpose.** λ approximates the max and min so that an influence function
   exists. It is not chosen to narrow the confidence interval.
2. **Level.** λ is applied pointwise in X, inside the expectation, so the weight
   can change with the covariates. Our ω₀ and ω₁ act on the ATE-level bounds.
3. **How the data-dependent weight is handled.** Their route is smoothness and
   the influence function: fix α, so the smoothed bound is a regular functional
   whose estimator is asymptotically normal. This is **a third route to valid
   inference with a data-dependent weight**, alongside our two:
   - **sample splitting**: choose the weight on one fold and evaluate on the
     other (convex split, M3, with the uniform-in-ω proof);
   - **a simultaneous critical value**: bound the sup over all weights, as in
     CLR without inequality selection (convex CLR, M9).

   This fits the "menu of methods" framing: the core difficulty is the
   data-dependent weight, and each route pays for it differently.
   - Splitting pays in sample size.
   - CLR pays in the critical value (2.353 instead of 1.96 at level 0.975).
   - Smoothing pays in an outward bias (the smoothed bound is wider than the
     sharp one) and becomes unstable as α grows.

**Not implementing it is defensible**, because it uses its own sensitivity model
and works at the covariate level, provided the paper says why. A referee may
still ask for the ATE-level Boltzmann weighting rule as a comparator. It would
be cheap in the toy:
- **Lower end:** B_α(L̂₁, L̂₂) − z·se.
- **Upper end:** B_{−α}(Û₁, Û₂) + z·se.
- **Standard error, by the delta method:** the gradient is
  ∂B/∂x₁ = λ₁ + αλ₁λ₂(x₁ − x₂), matching their EIF factor.

**Wording to check.** `Methods/convcomb.tex` currently says their weight "is
deterministically chosen". Only α is fixed; the weight moves with the estimated
bounds. A referee who knows the paper may object.

### Draft paragraph for `Methods/convcomb.tex` (replacing the Lanners paragraph; not applied)

> \citet{lanners_data_2025} face the same difficulty: the intersection of two
> bounds involves a maximum and a minimum, which are not pathwise
> differentiable. They replace each with the Boltzmann operator, a convex
> combination of the two bounds whose weights are a softmax of the bounds
> themselves, scaled by a user-chosen smoothing parameter. Their approach can
> therefore be read as a particular data-dependent choice of weights, which
> differs from ours in three respects. The weights are chosen to approximate the
> maximum and minimum, so that an efficient influence function exists, rather
> than to narrow the confidence interval. They are applied pointwise in the
> covariates, inside the expectation, whereas our weights $\omega_0$ and
> $\omega_1$ act on the ATE-level bounds. Finally, the dependence of the weights
> on the data is handled by smoothness rather than by sample splitting: for a
> fixed smoothing parameter the smoothed bounds are valid but wider than the
> sharp bounds, and admit standard influence-function-based inference, while
> larger values, which approximate the sharp bounds more closely, make the bias
> correction unstable where the two bounds cross
> \citep[Appendix~E.1]{lanners_data_2025}. Choosing covariate-dependent weights
> $\omega_0(x)$ and $\omega_1(x)$ with the inference tools developed here is a
> natural extension, which we leave to future work.

Their smoothing parameter is called α, which clashes with our CI level, so the
draft does not name it.

### Conservativeness: ATE level versus covariate level

There are two separate sources of conservativeness, and the paper should keep
them apart.

- **The target band itself.** E[max(l₁(X), l₂(X))] ≥ max(E l₁, E l₂). So
  whenever the dataset that sets the bound changes with X, the ATE-level
  intersection is wider than the covariate-level one, and no CI method recovers
  that.
  - In the toy there are no covariates, so this gap is zero.
  - Whether the bound functions cross in X in the causal sim has not been
    checked.
- **The inference cost:** sample splitting and critical values.

Going to the covariate level has two caveats:
- **Extra assumptions.** The datasets must share covariates and bound the same
  conditional estimand given X. Lanners' γ is what handles this for the RCT.
- **The same problem one level down.** The max inside the expectation is
  non-regular, which is exactly why Lanners smooth.

### The conditional method (clarified 2026-11-07)

Meaning: keep the ATE as the estimand, but let the weight depend on the
covariates. For the lower end, L(ω₀) = E[ω₀(X) l₁(X) + (1 − ω₀(X)) l₂(X)];
the upper end is analogous with ω₁(x).

- **It is Lanners' level.** Their λ(X) is one choice of ω(x). The conditional
  method is their structure with a variance-aware weight rule and our inference.
- **What it buys:**
  - a sharper target. The best weight is 1{l₁(x) ≥ l₂(x)}, giving
    E[max(l₁, l₂)(X)] ≥ max(E l₁, E l₂);
  - variance pooling where l₁(x) ≈ l₂(x).
- **What it needs** (strictly stronger than A2): both conditional bounds contain
  the same τ(x) at every x, which requires shared covariates, conditional
  exchangeability across datasets, and a common target distribution of X.
- **Inference:**
  - **Split:** extends naturally. For fixed ω(·), L(ω) is linear: weight each
    unit's doubly robust signal by ω(Xᵢ) (check DVDS's signal has the conditional
    bound as its conditional mean). Choose ω̂(·) on fold A and evaluate on
    fold B. The proof structure is the same, plus a CLT uniform over weight
    functions, to check.
  - **CLR:** does not carry over directly; the sup is over a function class.
  - **Smoothing:** Lanners' route.
  - **Plug-in indicator:** non-regular where l₁ ≈ l₂ on a set of positive
    probability (margin conditions; see the covariate-assisted bounds
    literature, e.g. Levis et al., from memory, so check).
- **Our simulations:**
  - The toy has no covariates, so there is no gain.
  - In the causal sim, the conditional bounds may be nested at every x, since
    the OS differ only in Λ. Then the only gain is pooling. Check whether they
    cross.
- **Status:** future work, as in the draft Lanners paragraph.

### PPI as an alternative to the convex combination (open)

Prediction-powered inference (Angelopoulos et al., 2023; this is from memory, so
check it) gains from the correlation between predictions and labels *on the same
units*, for a *single* target. Our two datasets are independent and bound
*different* quantities (L₁ ≠ L₂), so there is no correlation to exploit and no
shared target. PPI could sharpen each dataset's own bound estimates, for example
with extra unlabelled covariate data, but that is separate from the fusion step.
It only becomes an alternative with a concrete construction.

---

## 2. Yu et al. (2025), *Using a Two-Parameter Sensitivity Analysis Framework to Efficiently Combine Randomized and Non-randomized Studies*

- Ruoqi Yu, Bikram Karmakar, Jessica Vandeleest and Eleanor Bimla Schwarz;
  arXiv:2412.03731. The Zotero copy `Z7YPKMRC` is dated 2025-11-29; `B8ZAN38G`
  is another copy.
- **Fix the bib entry before citing.** `references.bib` has `yu_using_nodate`,
  which lists only Ruoqi Yu and has no date. Zotero's citation key for the
  preprint is `yu_using_2025`.

### Setting

- **Data:** one RCT and one OS.
- **Target:** the average treatment effect on the OS treated population (ATOT),
  β⋆ = E{Y(1) − Y(0) | Z = 1} in the OS (eq. 1, p. 8).
- **Selection:** selection into the RCT depends only on covariates
  (Assumption 1, p. 6), and the RCT's covariate support may be smaller than the
  OS's (Assumption 2, p. 6).
- **Design:** a triplet matching algorithm aligns OS treated, OS control and RCT
  units (Section 3, p. 8).

### Two sensitivity models, one per study

- **OS (unmeasured confounding):** a Rosenbaum-type model with parameter Γ
  (model (2), Section 4.1). The OS p-value (eq. 7, p. 20) is a sensitivity
  analysis p-value: it is computed at the worst case Γ allows, so it is not
  uniform and errs on the conservative side.
- **RCT (generalisability bias):** Δ ≥ 0 bounds the gap between the ATOT on the
  common support and the full ATOT (eq. 6, p. 18). Δ is on the outcome scale,
  and they suggest standardising it by √(S_t² + S_c²) (p. 19).
- **RCT inference** is by randomisation and assumes a constant additive
  treatment effect for RCT units (p. 20). The p-value is eq. 8 (p. 21), with a
  large-sample alternative (Theorem 3).

### The combination (the "p-value correction")

Fix (Γ, Δ). For each candidate value β₀, compute the one-sided sensitivity
p-values p_o(β₀) from the OS and p_r(β₀) from the RCT. Then combine them with
**Fisher's product** (eq. 9, p. 21):

- **Upper limit:** β̂_U = sup{β₀ : p_o(β₀) × p_r(β₀) ≥ κ_α}, with
  κ_α = exp(−χ²_{4;1−α}/2), using the χ² distribution with 4 degrees of freedom
  (details in Supplement S2).
- **Lower limit:** the same with the alternative reversed.
- **Two-sided interval:** both limits at α/2 (p. 23).
- **Monotonicity:** a supremum in each p-value's definition makes it monotone in
  β₀, which the efficiency proof needs (pp. 20–21).
- **The correction:** the sensitivity p-values are not uniform, so the standard
  Fisher calibration is not exact. Theorem 4 (p. 21) shows the combined interval
  is still asymptotically valid, and conservative when Γ > 1 or Δ > 0.

### Theory

- **Theorem 4:** validity, as above.
- **Theorem 5** (p. 22, under Assumption 4): for large enough samples, the
  combined upper limit is strictly below both single-study upper limits.
  - The asymptotics follow Bahadur efficiency and design sensitivity: the sample
    sizes grow while the level α_s → 0, not at a fixed α.
  - Assumptions 4.1 and 4.2: the p-values are monotone and continuous in β₀.
  - Assumption 4.3: Γ and Δ are "comparable", meaning the two single-study
    intervals converge to the same interval. **This is a tie.**
  - If one study's limit is asymptotically smaller, the combined interval
    converges to that study's interval (p. 22). That is intersection-like
    behaviour away from ties.
- **Finite samples:** shorter intervals at the usual 95% level are shown by
  simulation (Section 5).

### Simulations and application

- **Simulation DGP** (Section 5, p. 23): five observed normal covariates and one
  unobserved one, with true bias levels Γ⋆ and Δ⋆.
- **Comparators** (Supplement S3.2, p. 34):
  - elastic integrative analysis (Yang et al., 2023);
  - the integrative R-learner (Wu and Yang, 2022);
  - Hartman et al. (2015);
  - Stuart et al. (2011).

  These references are as cited by Yu et al.; check them before citing.
- **Application:** lactation and maternal health, using a small primate RCT and
  long-term observational records from the California National Primate Research
  Center.
- **Related multi-parameter sensitivity work they cite** (p. 20): Karmakar and
  Small (2020); Zhao et al. (2022).

### How it relates to our work

- **The setting is different from ours.** The triplet matching design needs
  unit-level data from both studies with shared covariates, plus
  Assumptions 1–2 (selection into the RCT depends only on X, and the RCT's
  support lies inside the OS's). Our methods need no covariates shared across
  datasets:
  - the convex split (M3) needs unit-level data *within* each dataset, to make
    the folds;
  - the CLR methods (M8, M9) need only each dataset's bound estimates and
    standard errors.

  Positioning line: Lanners et al. and Yu et al. both need unit-level data with
  shared covariates, while ours works with summary-level bounds.
- **They combine tests, not bands.** Yu et al. invert a combined test for the
  parameter β⋆; they never intersect identified sets. Their target is a CI for
  the parameter, not for the band, which is the Imbens–Manski distinction.
  - **Their interval is not an intersection of CIs.** At α = 0.05,
    κ ≈ exp(−9.49/2) ≈ 0.0087.
    - It can keep a point that one study's CI excludes: p_o = 0.02 and
      p_r = 0.9 give a product of 0.018, so β₀ stays in.
    - It can drop a point that both CIs keep: p_o = p_r = 0.08 give 0.0064, so
      β₀ is excluded. This is the pooling gain near ties.
  - **No winner's curse.** Fisher's statistic is fixed in advance for each β₀,
    and there is no step that picks a study. It is valid because both p-values
    are valid and the studies are independent.
    - The winner's curse belongs to the *unadjusted* intersection (M1): it picks
      the lowest upper end and still uses 1.96. In the toy at r = 1, g = 0 its
      joint coverage is 0.898.
    - Adjusted intersections (M6 Bonferroni, M8 Šidák) are valid because their
      critical value holds for both datasets at once.
    - The point estimate min(Û_k) is still biased inwards, which is why CLR's
      half-median-unbiased estimate uses its critical value at p = 1/2.
- **Intersection-based in target, not in construction.** Their null requires β₀
  to lie in both studies' identified sets, so their interval converges to the
  intersection of the identified sets (consistent with their remark on p. 22).
  That is the same target as ours and Lanners'.
- **All the valid methods invert a combined test** (my derivation; check before
  using it in the paper). Write ζₖ(u) = (u − Ûₖ)/seₖ for dataset k's
  standardised statistic at a candidate upper end u. Each method keeps u when a
  combined statistic is small enough:

  | Method | Statistic in ζ₁, ζ₂ | Rule |
  |---|---|---|
  | M8, CLR intersection | max(ζ₁, ζ₂) | Tippett (min-p) |
  | M9, convex CLR | max over unit c ≥ 0 of c₁ζ₁ + c₂ζ₂ | sup over weighted-z combinations |
  | M3, convex split | c₁ζ₁ + c₂ζ₂ with c fixed on fold A | one weighted-z combination |
  | Yu et al. | −2 Σ log Φ(−ζₖ) | Fisher |

  - **M8:** its threshold k = Φ⁻¹(√p) pays for the max.
  - **M9:** the weight c(w) = (w·se₁, (1−w)·se₂)/se(w) sweeps the positive
    quadrant, hence the critical value 2.353.
  - **M3:** fixing c on an independent fold avoids any sup, at the cost of the
    split.
  - **Yu et al.:** Fisher has no selection step, so it needs no adjustment.
  - **M1, naive intersection:** Tippett with an unadjusted threshold, which is
    the winner's curse.

  Lit-review framing: Yu et al. share our target but differ in setting (matched
  design, shared covariates, unit-level data) and in combination rule (Fisher).
  This table could also present our menu of methods in the paper.
- **Only the final combination step is comparable** (my reading; not in their
  paper, and not tested). The methods are not similar overall, but intersecting
  the per-dataset CIs is itself a way of combining tests.
  For the upper end, with p_k(u) = Φ((Û_k − u)/se_k), CLR without selection (M8)
  gives
  min_k(Û_k + k·se_k) = sup{u : min_k p_k(u) ≥ 1 − √p}, with k = Φ⁻¹(√p).
  - That is **Tippett's min-p combination** with a Šidák threshold.
  - Yu et al. use **Fisher's product** combination instead.
  - No single combination rule is best against every alternative (Birnbaum,
    1954, from memory; verify).
  - Fisher should gain where both datasets carry similar evidence, near ties at
    similar precision, and lose where one dataset dominates. That matches the
    regime of their Theorem 5.
  - This could soften "not concerned with intersection-based methods" in the
    lit-review outline. Their method is not an intersection and its setting
    differs, but its last step combines two studies' evidence with a different
    rule from ours.
- **A possible toy comparator** (not implemented). Each end at level 1 − α/2:
  - Upper end: sup{u : Φ((Û₁ − u)/se₁)·Φ((Û₂ − u)/se₂) ≥ exp(−χ²_{4;1−α/2}/2)},
    one root find per replicate.
  - Why it is valid: at u = U = min(U₁, U₂), the binding dataset's p-value is
    uniform and the other's is stochastically larger, and the two are
    independent.
### Why we don't compare against it (decided 2026-11-07)

There are already too many comparators, so this is scoped out. Base the reason on
scope and applicability, not on the method being worse, which has not been shown.

- **Defensible reasons:**
  - the matched design needs unit-level data from both studies with shared
    covariates;
  - their inference is tied to their sensitivity models. Only the Fisher step
    would carry over, and transplanting it would be a new method, so name it as
    future work.
- **Reasons a referee could rebut, so avoid them:**
  - "not intersection-based": they target the intersection of identified sets;
  - "winner's curse": Fisher's rule has no selection step;
  - "their CI is for the parameter, not the band": building each side at α/2
    also covers the band's ends;
  - "less efficient": not tested.

Draft paragraph for the lit review (fix the `yu_using_nodate` bib entry first):

> \citet{yu_using_nodate} combine an RCT and an OS through a two-parameter
> sensitivity analysis: a Rosenbaum-type model for unmeasured confounding in
> the OS and an additive bound on generalisability bias for the RCT. After a
> triplet matching design aligns treated and control OS units with RCT units,
> each study yields a sensitivity analysis $p$-value for a hypothesised effect,
> and the two are combined with Fisher's method. Inverting the combined test
> gives a confidence interval that is asymptotically valid and, when the two
> studies' limits coincide, asymptotically shorter than either single-study
> interval. Their target coincides with ours, in that a hypothesised effect is
> retained only if it is compatible with both studies, but their approach
> differs in two respects that place it outside the scope of this work. First,
> the matching design requires unit-level data from both studies on a shared
> set of covariates, whereas our methods require only each study's estimated
> bounds and standard errors, or, for the sample-splitting approach, unit-level
> data within each study. Second, their inference is tied to their sensitivity
> models, using randomisation inference under a constant additive effect in the
> RCT, whereas our methods apply to any pair of asymptotically normal bound
> estimators. Their combination step could nonetheless be applied to bound
> estimators such as ours, and we leave this comparison to future work.

- **Differences in assumptions to note:**
  - their matched design;
  - randomisation inference with a constant additive effect in the RCT;
  - an additive, scale-dependent Δ;
  - the target is the ATOT, not the ATE.
