# Example of a CLeaR related-work section (style model)

Pasted by the user on 2026-11-09 as an example of a CLeaR literature review.
It is Section 2 of a paper proposing Causal-ICM: multi-task Gaussian
processes for CATE estimation by fusing an RCT with an observational study.
The authors and venue were not given; don't cite it as a source without
finding out.

Use it only as a model for structure, length and register. Don't reuse its
wording in our paper. How it shapes our section is set out in
`quality_reports/plans/2026-11-08_lit-review.md` (section "Model: the CLeaR
example").

Its citations are as the authors wrote them and were **not checked**. In
particular, its "Fawkes et al. (2025)" (multi-task GPs for pseudo-outcome
regression, error bounds on the observational bias) may not be the Fawkes et
al. paper in our `references.bib` (`fawkes_hardness_2025-1`), and its "Lin et
al. (2025)" (Bayesian dynamic borrowing through a power likelihood) is not
in our plan.

---

2. Related Literature

Heterogeneous treatment effect (HTE) estimation has received renewed interest in recent
years. The most common strategy is to focus on the conditional average treatment effect
(CATE) function, which quantifies the expected treatment effect given particular covariate
values. Under suitable identifiability assumptions (Dahabreh and Hernán, 2019), the CATE
can be estimated parametrically - e.g. using linear regression - or non-parametrically, e.g.
using nearest-neighbour matching, kernel methods (Brantner et al., 2023), or tree-based
methods including causal random forests (Wager and Athey, 2018) and their Bayesian
counterpart (Hahn et al., 2020). A related strand of work studies semi-parametric and
debiased machine learning approaches that combine flexible outcome models with orthogonal
score constructions to obtain robust estimates under high-dimensional confounding.

Data fusion across an RCT and an observational study has tended to focus on ATE
estimation (see Lin et al. (2024) and Colnet et al. (2024) for comprehensive reviews). These
methods aim to improve efficiency by integrating information from both sources. Demirel
et al. (2024) focus instead on data fusion for generalisability, proposing a framework that
builds an observational predictor and then treat the difference between the trial outcome
function and that observational predictor as a bias function that captures both confounding
and transportation bias due to the differences between the trial and target populations.
Similar approaches, although outside the immediate scope of our work, have focused on
data fusion for the estimation of long term effects, where randomised data are unconfounded
but contain only short term effects, while observational data, although confounded, contain
information about long term outcomes (Ghassami et al. (2025), Imbens et al. (2025)).

Data fusion methods for CATE estimation, the focus of this paper, follow three broad lines
of work. One strategy is to use experimental data to debias observational estimates by learning
a correction or bias function. Kallus et al. (2018) introduce a two-step approach that corrects
hidden confounding in the absence of covariate overlap, relying on the strong assumption
of linear confounding and the ability to identify this correction term parametrically. Yang
et al. (2025b) relax this linearity assumption and provide identifiability and efficiency results
under more flexible structural models. Hatt et al. (2022) develop a representation-learning
strategy that learns shared features and confounded outcome models from observational data,
using the RCT to estimate a bias function that debiases these models, with finite-sample
bounds that highlight how performance depends on sample size, distribution shift, and
bias complexity. Wu and Yang (2022) propose an R-learner that achieves consistency and
asymptotic efficiency under covariate overlap between the RCT and the observational study.
A second direction consists of "test-then-pool" strategies, in which the trial and observational
datasets are combined only if the observational study appears sufficiently unbiased. For
instance, Yang et al. (2023) implement such a strategy, although their approach is limited
to linear treatment effect models and relies on strong parametric assumptions. A third line
of work uses Bayesian dynamic borrowing to regulate how much information is extracted
from the observational study. For instance, Lin et al. (2025) propose a Bayesian dynamic
borrowing approach through a power likelihood to reduce the effect of confounding and
control the degree of information we borrow from the observational study.

Our method, Causal-ICM, takes a Bayesian dynamic borrowing approach based on multi-
task GPs, whereby each potential outcome is treated as a distinct task. Multi-task GPs
have been successfully employed for causal inference using solely observational data (Alaa
and van der Schaar, 2017) and extended to handle hierarchical hidden confounders (Witty
et al., 2020). Others have built on this framework to develop a general class of counterfactual
multi-task deep kernel models that efficiently estimate causal effects and learn policies by
stacking coregionalized GPs and deep kernels (Caron et al., 2022). GPs have also been
employed as a matching tool for causal inference (Huang et al., 2023). Recent work has also
explored Gaussian-process-based partially linear models for heterogeneous treatment effects,
combining parametric effect components with flexible GP nuisance estimation (Horii and
Chikahara, 2024). However, these approaches focus on single-source observational settings.
In contrast, our work develops a multi-task GP framework specifically designed for causal
data fusion between experimental and observational sources. Our approach is most similar
to concurrent work employing multi-task GPs for pseudo-outcome regression to obtain error
bounds on the observational bias term (Fawkes et al., 2025).
