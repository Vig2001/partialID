"""Convergence of the split convex method's estimated band as n grows.

Toy layout of src/toy_sim.py with equal gaps g_L = g_U = g: dataset 1
attains the true lower end, dataset 2 the true upper end, and the other
dataset's bound at each end is g looser. Three cases: separated (g fixed),
near tie (g = h SE_1(n), shrinking with n) and exact tie (g = 0). The noise
SDs are fixed at the toy's values, so the standard errors shrink like
n^{-1/2}.

For each n_1 (n_2 in the toy's ratio, n_1 / 4), records the split method's
evaluation-fold (fold B) estimates L_hat = w_L L_hat_1 + (1 - w_L) L_hat_2
and U_hat, with the weights from fold A. Fold means and variances are drawn
as sufficient statistics (normal mean, scaled chi-square variance), which
for Gaussian data have exactly the distribution of the toy's estimator.

Usage, from the repo root:
    python -m src.toy_convergence [--seed 7] [--reps 5000] [--gap 0.2]
        [--near-tie-h 2] [--ratio 2] [--frac-a 0.3] [--n1 100 200 ...]

Plan: quality_reports/plans/2026-11-04_toy-bound-convergence.md
Outputs: src/outputs/toy_convergence/ (results.csv, bound_convergence.pdf)
"""

import argparse
import dataclasses
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src import toy_sim
from src.helpers import ci_methods
from src.helpers import mc_metrics
from src.plotting import toy_sim_plots

DEFAULT_SEED = 7
DEFAULT_OUT_DIR = Path(__file__).resolve().parent / "outputs" / "toy_convergence"
DEFAULT_N1 = (100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)
SEPARATED = "separated"
NEAR_TIE = "near tie"
EXACT_TIE = "exact tie"
CASES = (SEPARATED, NEAR_TIE, EXACT_TIE)
CHECK_FAMILY_LEVEL = 0.01  # family-wise level of the Monte Carlo checks


def config_at(base: toy_sim.ToyConfig, n1: int) -> toy_sim.ToyConfig:
    """The base config at sample size n1, with the base noise SDs.

    n_2 keeps the base ratio n_2 / n_1, and SE_1 is sigma_1 / sqrt(n_1)
    with sigma_1 = SE_1 sqrt(n_1) at the base sizes, so that dgp_params
    returns the base sigma at every n.
    """
    n2 = round(n1 * base.n[1] / base.n[0])
    sigma1 = base.se1 * np.sqrt(base.n[0])
    return dataclasses.replace(base, n=(n1, n2), se1=sigma1 / np.sqrt(n1))


def case_gap(case: str, cfg: toy_sim.ToyConfig, gap: float,
             near_tie_h: float) -> float:
    """Gap g_L = g_U of a case at the sample size of cfg."""
    if case == SEPARATED:
        return gap
    if case == NEAR_TIE:
        return near_tie_h * cfg.se1
    if case == EXACT_TIE:
        return 0.0
    raise ValueError(f"Unknown case {case=}; expected one of {CASES}")


def fold_estimates(
    params: toy_sim.DGPParams,
    sizes: Sequence[int],
    var_sizes: Sequence[int],
    z: np.ndarray,
    rng: np.random.Generator,
) -> toy_sim.BoundEstimates:
    """Bounds from one fold of each dataset, drawn as sufficient statistics.

    For Gaussian data the fold mean is N(mu_k, sigma_k^2 / m_k) and the
    sample variance is sigma_k^2 chi^2_{m_k - 1} / (m_k - 1), independently,
    as in toy_sim.estimate_bounds on m_k draws.

    Args:
        params: Population quantities (toy_sim.dgp_params).
        sizes: Fold size m_k of each dataset.
        var_sizes: Sample size dividing each sample variance.
        z: Standard normals for the fold means, shape (R, 2).
        rng: Generator for the chi-square draws.
    """
    sizes = np.asarray(sizes)
    Ybar = params.mu + params.sigma / np.sqrt(sizes) * z
    dof = sizes - 1
    sample_var = params.sigma**2 * rng.chisquare(dof, size=z.shape) / dof
    return toy_sim.BoundEstimates(
        L=Ybar - params.Delta, U=Ybar + params.Delta,
        s2=sample_var / np.asarray(var_sizes),
    )


def split_estimates(
    params: toy_sim.DGPParams,
    cfg: toy_sim.ToyConfig,
    z_a: np.ndarray,
    z_b: np.ndarray,
    rng: np.random.Generator,
) -> ci_methods.CIResult:
    """The split method: weights from fold A, bounds from fold B.

    Folds as in toy_sim.ci_by_method: fold A has frac_a of each dataset, and
    its variances are rescaled to the fold-B size.
    """
    n_a = [int(cfg.frac_a * n_k) for n_k in cfg.n]
    n_b = [n_k - a for n_k, a in zip(cfg.n, n_a)]
    fold_a = fold_estimates(params, n_a, n_b, z_a, rng)
    fold_b = fold_estimates(params, n_b, n_b, z_b, rng)
    w_L, w_U = toy_sim.plug_in_weights(fold_a, cfg.z)
    return ci_methods.convex_ci(fold_b.L, fold_b.U, fold_b.s2, w_L, w_U, cfg.z)


def summarise(est: ci_methods.CIResult, params: toy_sim.DGPParams) -> dict:
    """Mean, MC s.e. and 5% / 95% quantiles of L_hat and U_hat, and weights.

    resid_<end> is the estimate minus the convex combination of the
    population bounds at the estimated weight; its mean is 0 because fold
    B is independent of the weights.
    """
    row = {}
    ends = (("L", est.L_hat, est.w_L, params.L_star),
            ("U", est.U_hat, est.w_U, params.U_star))
    for end, values, w, star in ends:
        row[f"{end}_mean"], row[f"{end}_mcse"] = mc_metrics.mean_se(values)
        row[f"{end}_q05"], row[f"{end}_q95"] = np.quantile(values, [0.05, 0.95])
        row[f"w{end}_mean"], row[f"w{end}_mcse"] = mc_metrics.mean_se(w)
        resid = values - (w * star[0] + (1 - w) * star[1])
        row[f"resid_{end}_mean"], row[f"resid_{end}_mcse"] = (
            mc_metrics.mean_se(resid)
        )
    return row


def simulate(
    base: toy_sim.ToyConfig,
    n1_values: Sequence[int],
    gap: float,
    near_tie_h: float,
    r: float,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """One row per (case, n_1); fold-mean normals shared across cells."""
    z_a, z_b = rng.standard_normal((2, base.reps, 2))
    rows = []
    for n1 in n1_values:
        cfg = config_at(base, n1)
        for case in CASES:
            g = case_gap(case, cfg, gap, near_tie_h)
            params = toy_sim.dgp_params(g, g, r, cfg)
            est = split_estimates(params, cfg, z_a, z_b, rng)
            row = dict(case=case, n1=cfg.n[0], n2=cfg.n[1], gap=g, r=r,
                       true_lo=cfg.true_lo, true_hi=cfg.true_hi,
                       L_2=params.L_star[1], U_1=params.U_star[0])
            row.update(summarise(est, params))
            rows.append(row)
    return pd.DataFrame(rows)


def run_checks(df: pd.DataFrame) -> bool:
    """Prints PASS/FAIL for the theoretical properties; True if all pass.

    Each check is a family of z-tests at Bonferroni level
    CHECK_FAMILY_LEVEL over its comparisons.
    """
    def z_limit(k):
        return stats.norm.ppf(1 - CHECK_FAMILY_LEVEL / (2 * k))

    resid = pd.concat([df.resid_L_mean / df.resid_L_mcse,
                       df.resid_U_mean / df.resid_U_mcse]).abs()
    tie = df[df.case == EXACT_TIE]
    tie_z = pd.concat([(tie.L_mean - tie.true_lo) / tie.L_mcse,
                       (tie.U_mean - tie.true_hi) / tie.U_mcse]).abs()
    largest = df[(df.case == SEPARATED) & (df.n1 == df.n1.max())]
    sep_z = pd.concat([(largest.L_mean - largest.true_lo) / largest.L_mcse,
                       (largest.U_mean - largest.true_hi) / largest.U_mcse]).abs()
    checks = [
        ("fold-B estimates unbiased for the weighted population bounds",
         resid),
        ("exact tie: estimated bounds unbiased at every n", tie_z),
        ("separated: bias zero at the largest n", sep_z),
    ]
    all_ok = True
    for name, z in checks:
        limit = z_limit(len(z))
        ok = bool((z <= limit).all())
        all_ok &= ok
        print(f"{'PASS' if ok else 'FAIL'} {name}: max |z| = {z.max():.2f} "
              f"(limit {limit:.2f}, {len(z)} tests)")
    return all_ok


def case_titles(gap: float, near_tie_h: float) -> list[str]:
    return [f"Separated, $g = {gap:g}$",
            rf"Near tie, $g = {near_tie_h:g}\,\mathrm{{SE}}_1$",
            "Exact tie, $g = 0$"]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    defaults = toy_sim.ToyConfig()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="Seed of the random generator.")
    parser.add_argument("--reps", type=int, default=defaults.reps,
                        help="Monte Carlo replicates per (case, n_1).")
    parser.add_argument("--gap", type=float, default=0.2,
                        help="Gap g_L = g_U of the separated case.")
    parser.add_argument("--near-tie-h", type=float, default=2.0,
                        help="Near-tie gap in units of SE_1(n): g = h SE_1.")
    parser.add_argument("--ratio", type=float, default=2.0,
                        help="r = SE_2 / SE_1.")
    parser.add_argument("--frac-a", type=float, default=defaults.frac_a,
                        help="Fold-A fraction of each dataset.")
    parser.add_argument("--n1", type=int, nargs="+", default=list(DEFAULT_N1),
                        help="Sizes of dataset 1; n_2 keeps the toy's ratio.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="Output directory.")
    args = parser.parse_args(argv)
    if min(args.n1) < 20:
        parser.error(f"--n1 values must be >= 20, got {min(args.n1)}")
    return args


def main() -> None:
    args = parse_args()
    base = toy_sim.ToyConfig(reps=args.reps, frac_a=args.frac_a)
    rng = np.random.default_rng(args.seed)
    df = simulate(base, sorted(args.n1), args.gap, args.near_tie_h,
                  args.ratio, rng)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    results_path = args.out_dir / "results.csv"
    df.to_csv(results_path, index=False)
    run_checks(df)
    fig_path = args.out_dir / "bound_convergence.pdf"
    toy_sim_plots.plot_bound_convergence(
        df, CASES, case_titles(args.gap, args.near_tie_h), fig_path
    )
    print(f"Written to {results_path}\nWritten to {fig_path}")


if __name__ == "__main__":
    main()
