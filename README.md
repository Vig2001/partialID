# Partial Identificiation of CATEs and Sensitivity Analysis

At the moment, this repository contains Python simulations demonstrating the effects of unmeasured confounding on the Observational Conditional Average Treatment Effect (CATE).

### Aim
* Understand the relation between the correction function $\Delta(X)$ (Fawkes et al. 2025) and $U$ (unobserved confounders)
* What does enforcing a smoothness constraint on $\Delta(X)$ mean on $U$ or the level of confounding?
* Connect the correction function to the marginal sensitivity model
* Does this mean we can retrieve bounds on $\Gamma$ using bounds on $\Delta(X)$?

### Setup
Install dependencies: `pip install -r requirements.txt`
