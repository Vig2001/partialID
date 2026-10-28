# Plan: two-OS simulation based on Dorn, Guo & Kallus (DVDS)

Date: 2026-10-28. Status: **draft**. Open decisions below need answering before any code is written.

## Why
The old ZSB/NIW simulation is archived in `scripts_misc/outdated_help/`. We move to the Dorn, Guo & Kallus (DGK) simulation because their DVDS bound estimators are asymptotically normal, which matches the per-endpoint Gaussian setup ($\hat L_k, \hat U_k$ with variances) that the convex combination method assumes. The plan is to run their DGP twice, giving two independent observational studies (OS): one smaller than the other, possibly observing different covariates so that the amount of unobserved confounding differs.

## Source
Repo: https://github.com/CausalML/DoublyValidDoublySharpSensitivity (paper: https://arxiv.org/abs/2112.11449, Section 7.1). R files not saved locally; re-read them on GitHub when porting.
- `simplified_dvds.R`: estimators only (`dvds.binary`, `dvds.continuous`, `phi.continuous`, ZSB comparison). No simulation.
- `dvds_sim.R`: the simulation (DGP + nuisance fits + bounds).
- `dvds_plotting.r`: true sharp bounds and coverage.

## Their DGP (`dvds_sim.R`)
1000 replications, $n = 1000$, $d = 5$, $\Lambda$ on 20 equally spaced values in $[1, 3]$.
- $X \sim \mathrm{Unif}[-1,1]^5$.
- $P(Z=1 \mid X) = \mathrm{expit}\{-(X_1 + 0.5\,\mathbb{1}\{X_2>0\} + 0.5X_2X_3)\}$, so $e(X) \in [0.12, 0.82]$ (my calculation: good overlap).
- Binary: $Y \sim \mathrm{Bern}(\mathrm{expit}\{-(0.5X_1 + X_2 + 0.25X_2X_3)\})$.
- Continuous: $Y = 2\,\mathrm{sign}(X_1) + X_2 + X_2X_3 + (1+X_4^2)\varepsilon$, $\varepsilon \sim N(0,1)$.
- $X_5$ is pure noise; $X_4$ only affects the outcome variance.

Two features to keep in mind:
1. $Y$ does not depend on $Z$, so the true ATE is 0 and $Y(1) = Y(0)$.
2. There is no unobserved confounding (true $\Lambda = 1$). Their coverage target is the population sharp band $[\mathrm{ATE}^-(\Lambda), \mathrm{ATE}^+(\Lambda)]$ at each hypothesised $\Lambda$, not the ATE.

## Their estimator (`simplified_dvds.R`)
Marginal sensitivity model with parameter $\Lambda$, $\tau = \Lambda/(\Lambda+1)$. Upper bound on $E[Y(1)]$, per unit:
$$\varphi = ZY + (1-Z)\kappa(X) + \tfrac{1-e(X)}{e(X)}\,Z\,(\text{pseudo} - \kappa(X)),$$
$$\kappa(x) = \tfrac{1}{\Lambda}\mu_1(x) + (1-\tfrac{1}{\Lambda})\,\mathrm{CVaR}_\tau(x), \qquad \text{pseudo} = \tfrac{1}{\Lambda}Y + (1-\tfrac{1}{\Lambda})\{Q_\tau(X) + (Y-Q_\tau(X))_+/(1-\tau)\}.$$
- **Treated units** contribute their own $Y$. **Controls' missing $Y(1)$** are replaced by the worst case $\kappa$, a mix of the conditional mean and the upper-tail mean above the $\tau$-quantile. The **augmentation term** makes it doubly robust: valid if either the quantile model or $(\mu, e)$ is correct, sharp if both are.
- **Lower bound:** apply the same thing to $-Y$. **$E[Y(0)]$:** swap $Z \leftrightarrow 1-Z$, $e \leftrightarrow 1-e$. Then $\mathrm{ATE}^+ = E[Y(1)]^+ - E[Y(0)]^-$ and $\mathrm{ATE}^- = E[Y(1)]^- - E[Y(0)]^+$.
- **Binary outcomes** use closed-form $\kappa$ (e.g. $\kappa^+ = \min(1 - 1/\Lambda + \mu/\Lambda,\ \Lambda\mu)$) and a thresholded quantile.
- **Nuisances** are all fitted with `grf`: probability forest for $e$, regression/probability forests for $\mu_0, \mu_1$, a quantile forest for $Q_\tau$, and CVaR from the forest weights of a regression forest fitted at the median $\tau$ (reused for every $\tau$). They use out-of-bag predictions rather than explicit K-fold cross-fitting.
- **CI:** mean $\pm\, z_{1-\alpha/2}\,\mathrm{sd}/\sqrt{n}$ of the influence-function values, separately for each endpoint. Each side misses with probability $\alpha/2$, so the band is covered with probability $\geq 1-\alpha$.

## Their truth and coverage (`dvds_plotting.r`)
- **Continuous**, closed form: $\pm\tfrac{4}{3}\,\tfrac{\Lambda^2-1}{\Lambda}\,\phi(\Phi^{-1}(\tfrac{\Lambda}{\Lambda+1}))$, where $\tfrac43 = E[1+X_4^2]$. This works only because $Y \mid X$ is Gaussian and does not depend on $Z$.
- **Binary:** Monte Carlo with 500,000 draws using the closed-form $\kappa$.
- **Coverage** counts as covered when `lower.CI < true_lower & true_upper < upper.CI`. Reported in their repo (`Simulations/dvds_output/*_coverageRates.tex`):
  - continuous: 90.5% (at $\Lambda = 1$) to 98.2%, average 95.8%;
  - binary: 94.6% to 97.5%, average 96.8%.
- **Caveat:** 90.5% at $\Lambda = 1$ means the CI undercovers at $n = 1000$, so asymptotic normality alone does not give nominal coverage. This matters for the smaller study.

## Implications for the two-OS design
- **Hiding a covariate creates genuine confounding.** If study $k$ drops a covariate that enters both $Z$ and $Y$ (e.g. $X_1$), that covariate becomes the unobserved confounder. The two studies then have different true $\Lambda^\star_k$. Drawing them separately from the same population makes them independent, with the same target ATE.
- **The continuous closed form no longer applies.** With $X_1$ hidden, $Y \mid X_{\text{obs}}, Z$ is a two-component Gaussian mixture (from $2\,\mathrm{sign}(X_1)$) whose weights depend on $Z$. The population bounds $L_k, U_k$ must be computed numerically (large Monte Carlo plus mixture quantiles/CVaR).
- **$\Lambda^\star_k$ must be computed, not chosen.** It is $\sup$ of the odds ratio between $e(X_{\text{full}})$ and $e(X_{\text{obs}})$. Hiding $X_1$ bounds it by $e^2 \approx 7.4$, but the actual value will be smaller and has not been computed yet.
- **Add a treatment effect.** With ATE = 0 and symmetric bounds, bugs can hide (e.g. sign errors that cancel).
- **Variances are estimated here.** This sits outside the known-variance Gaussian section of the paper, so it belongs in a separate section or simulation.

## Open decisions
1. **Language.** Port to Python (quantile-forest package plus forest weights for CVaR; numbers will not match `grf` exactly) or call R `grf` from Python. Needs settling first.
2. **Which covariates each study observes.** E.g. study 1 (large) observes all of $X$ ($\Lambda^\star_1 = 1$) and study 2 (small) hides $X_1$; or each hides a different covariate.
3. **Sample sizes** $n_1, n_2$ and the ratio between them.
4. **Treatment effect:** constant $\beta$ or heterogeneous; continuous outcome only or binary too.
5. **$\Lambda_k$ grid per study**, and whether coverage targets the band at fixed $\Lambda_k$ or the ATE at $\Lambda_k \geq \Lambda^\star_k$.

## Proposed steps (after decisions)
1. **Replicate DGK first.** Port their single-study DGP and estimator, and check that coverage at their settings is close to their reported 90.5–98.2% (continuous) and 94.6–97.5% (binary). This validates the port before changing anything.
2. **DGP module** in `src/helpers/` (e.g. `dvds_dgp.py`): draw the full population with an `rng` argument and a treatment-effect parameter, plus a function that returns the observed version of a study given its covariate list.
3. **DVDS estimator module** in `src/helpers/`: returns $\hat L_k, \hat U_k$ and their standard errors for each $\Lambda$.
4. **Truth:** population bounds $L_k, U_k$ for each study and $\Lambda$ by large Monte Carlo, plus $\Lambda^\star_k$ numerically.
5. **Entry script** `src/dvds_two_os.py` (CLI flags for seed, $n_1$, $n_2$, reps, $\Lambda$ grids). Feed the endpoints into the existing `ci_methods` (intersection, calibrated intersection, convex oracle/split/plug-in) and summarise with `mc_metrics`.
6. **Tests** for the DGP and estimator: band collapses to a point at $\Lambda = 1$ with no hidden covariates; bounds widen monotonically in $\Lambda$; continuous truth matches the DGK closed form in the no-hidden-covariate case.
