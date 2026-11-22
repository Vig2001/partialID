"""Convergence of the estimated fused band in the toy as n grows.

Toy layout of src/toy_sim.py with equal gaps g_L = g_U = g between the two
datasets' population bounds: dataset 1 attains the true lower end, dataset 2
the true upper end, and the other dataset's bound at each end is g looser.
The population gaps are fixed; the noise SDs are fixed at the toy's values,
so the standard errors shrink like n^{-1/2}.

With equal gaps the two datasets' population means differ by g
(mu_1 - mu_2 = g), and the population bounds by g at each end.

For each n_1 (n_2 in the toy's ratio, n_1 / 4) and gap, records two
estimated bands:
    - convex split: the evaluation-fold (fold B) estimates
      L_hat = w_L L_hat_1 + (1 - w_L) L_hat_2 and U_hat, weights from fold A;
    - intersection: max(L_hat_1, L_hat_2) and min(U_hat_1, U_hat_2) on the
      full sample.
Fold means and variances are drawn as sufficient statistics (normal mean,
scaled chi-square variance), which for Gaussian data have exactly the
distribution of the toy's estimator; the full sample pools the two folds.

Figures: bound_convergence.pdf, the convex split's estimated bounds for a
difference in population means (--gaps) and for equal means (g = 0); and
uniform_convergence.pdf, the standardised outward bias of each end for the
convex split and the intersection at small fixed gaps (--uniform-gaps): it
vanishes for each fixed g, but not uniformly over g.

Usage, from the repo root:
    python -m src.toy_convergence [--seed 7] [--reps 5000]
        [--gaps 0.05 0.2 0.5] [--uniform-gaps 0.01 0.02 0.05 0.1 0.2]
        [--ratio 2] [--frac-a 0.3] [--n1 100 200 ...]

Plan: quality_reports/plans/done/2026-11-04_toy-bound-convergence.md
Outputs: src/outputs/toy_convergence/ (results.csv, bound_convergence.pdf,
    uniform_convergence.pdf)
"""

import argparse
import dataclasses
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

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
DEFAULT_GAPS = (0.05, 0.2, 0.5)
DEFAULT_UNIFORM_GAPS = (0.01, 0.02, 0.05, 0.1, 0.2)
SPLIT = "convex split"
INTERSECTION = "intersection"
CHECK_FAMILY_LEVEL = 0.01  # family-wise level of the Monte Carlo checks
SETTLED_SES = 10  # a gap counts as settled at n if it is this many SE_1


class FoldStats(NamedTuple):
    """Sample means and variances of each dataset in one fold, each (R, 2)."""

    Ybar: np.ndarray
    var: np.ndarray
    sizes: np.ndarray


def config_at(base: toy_sim.ToyConfig, n1: int) -> toy_sim.ToyConfig:
    """The base config at sample size n1, with the base noise SDs.

    n_2 keeps the base ratio n_2 / n_1, and SE_1 is sigma_1 / sqrt(n_1)
    with sigma_1 = SE_1 sqrt(n_1) at the base sizes, so that dgp_params
    returns the base sigma at every n.
    """
    n2 = round(n1 * base.n[1] / base.n[0])
    sigma1 = base.se1 * np.sqrt(base.n[0])
    return dataclasses.replace(base, n=(n1, n2), se1=sigma1 / np.sqrt(n1))


def draw_fold(
    params: toy_sim.DGPParams,
    sizes: Sequence[int],
    z: np.ndarray,
    rng: np.random.Generator,
) -> FoldStats:
    """Fold means and sample variances, drawn as sufficient statistics.

    For Gaussian data the mean of m_k draws is N(mu_k, sigma_k^2 / m_k) and
    the sample variance is sigma_k^2 chi^2_{m_k - 1} / (m_k - 1),
    independently.

    Args:
        params: Population quantities (toy_sim.dgp_params).
        sizes: Fold size m_k of each dataset.
        z: Standard normals for the fold means, shape (R, 2).
        rng: Generator for the chi-square draws.
    """
    sizes = np.asarray(sizes)
    Ybar = params.mu + params.sigma / np.sqrt(sizes) * z
    dof = sizes - 1
    var = params.sigma**2 * rng.chisquare(dof, size=z.shape) / dof
    return FoldStats(Ybar=Ybar, var=var, sizes=sizes)


def pool_folds(a: FoldStats, b: FoldStats) -> FoldStats:
    """Mean and sample variance of the two folds together."""
    n = a.sizes + b.sizes
    Ybar = (a.sizes * a.Ybar + b.sizes * b.Ybar) / n
    within = (a.sizes - 1) * a.var + (b.sizes - 1) * b.var
    between = a.sizes * b.sizes / n * (a.Ybar - b.Ybar) ** 2
    return FoldStats(Ybar=Ybar, var=(within + between) / (n - 1), sizes=n)


def bounds(
    fold: FoldStats, Delta: np.ndarray, var_sizes: Sequence[int]
) -> toy_sim.BoundEstimates:
    """Bounds Ybar_k -/+ Delta_k, as in toy_sim.estimate_bounds."""
    return toy_sim.BoundEstimates(
        L=fold.Ybar - Delta, U=fold.Ybar + Delta,
        s2=fold.var / np.asarray(var_sizes),
    )


def estimated_bands(
    params: toy_sim.DGPParams,
    cfg: toy_sim.ToyConfig,
    z_a: np.ndarray,
    z_b: np.ndarray,
    rng: np.random.Generator,
) -> dict[str, ci_methods.CIResult]:
    """Convex split (weights on fold A, bounds on fold B) and intersection.

    Folds as in toy_sim.ci_by_method: fold A has frac_a of each dataset, and
    its variances are rescaled to the fold-B size.
    """
    n_a = [int(cfg.frac_a * n_k) for n_k in cfg.n]
    n_b = [n_k - a for n_k, a in zip(cfg.n, n_a)]
    stats_a = draw_fold(params, n_a, z_a, rng)
    stats_b = draw_fold(params, n_b, z_b, rng)
    fold_a = bounds(stats_a, params.Delta, n_b)
    fold_b = bounds(stats_b, params.Delta, n_b)
    full = bounds(pool_folds(stats_a, stats_b), params.Delta, cfg.n)
    w_L, w_U = toy_sim.plug_in_weights(fold_a, cfg.z)
    return {
        SPLIT: ci_methods.convex_ci(fold_b.L, fold_b.U, fold_b.s2, w_L, w_U,
                                    cfg.z),
        INTERSECTION: ci_methods.intersection_ci(full.L, full.U, full.s2,
                                                 cfg.z),
    }


def summarise(est: ci_methods.CIResult, params: toy_sim.DGPParams,
              method: str) -> dict:
    """Mean, MC s.e., SD and 5% / 95% quantiles of L_hat and U_hat, weights.

    For the split, resid_<end> is the estimate minus the convex combination
    of the population bounds at the estimated weight; its mean is 0 because
    fold B is independent of the weights. It is not defined for the
    intersection, whose selection uses the same data as its estimate.
    """
    row = {}
    ends = (("L", est.L_hat, est.w_L, params.L_star),
            ("U", est.U_hat, est.w_U, params.U_star))
    for end, values, w, star in ends:
        row[f"{end}_mean"], row[f"{end}_mcse"] = mc_metrics.mean_se(values)
        row[f"{end}_sd"] = values.std(ddof=1)
        row[f"{end}_q05"], row[f"{end}_q95"] = np.quantile(values, [0.05, 0.95])
        row[f"w{end}_mean"], row[f"w{end}_mcse"] = mc_metrics.mean_se(w)
        row[f"resid_{end}_mean"] = row[f"resid_{end}_mcse"] = np.nan
        if method == SPLIT:
            resid = values - (w * star[0] + (1 - w) * star[1])
            row[f"resid_{end}_mean"], row[f"resid_{end}_mcse"] = (
                mc_metrics.mean_se(resid)
            )
    return row


def simulate(
    base: toy_sim.ToyConfig,
    n1_values: Sequence[int],
    gaps: Sequence[float],
    r: float,
    rng: np.random.Generator,
) -> pd.DataFrame:
    """One row per (n_1, gap, method); fold-mean normals shared across cells."""
    z_a, z_b = rng.standard_normal((2, base.reps, 2))
    rows = []
    for n1 in n1_values:
        cfg = config_at(base, n1)
        for g in gaps:
            params = toy_sim.dgp_params(g, g, r, cfg)
            bands = estimated_bands(params, cfg, z_a, z_b, rng)
            for method, est in bands.items():
                row = dict(n1=cfg.n[0], n2=cfg.n[1], se1=cfg.se1, gap=g, r=r,
                           method=method, true_lo=cfg.true_lo,
                           true_hi=cfg.true_hi, L_2=params.L_star[1],
                           U_1=params.U_star[0])
                row.update(summarise(est, params, method))
                rows.append(row)
    return pd.DataFrame(rows)


def _bias_z(df: pd.DataFrame) -> pd.Series:
    """Signed bias of both ends over its MC s.e.; positive means inward."""
    return pd.concat([(df.L_mean - df.true_lo) / df.L_mcse,
                      (df.true_hi - df.U_mean) / df.U_mcse])


def run_checks(df: pd.DataFrame) -> bool:
    """Prints PASS/FAIL for the theoretical properties; True if all pass.

    Each check is a family of z-tests at Bonferroni level
    CHECK_FAMILY_LEVEL over its comparisons.
    """
    def z_limit(k):
        return stats.norm.ppf(1 - CHECK_FAMILY_LEVEL / (2 * k))

    split = df[df.method == SPLIT]
    tie = df[df.gap == 0]
    largest = split[(split.n1 == split.n1.max())
                    & (split.gap >= SETTLED_SES * split.se1)
                    & (split.gap > 0)]
    resid = pd.concat([split.resid_L_mean / split.resid_L_mcse,
                       split.resid_U_mean / split.resid_U_mcse]).abs()
    split_tie = _bias_z(tie[tie.method == SPLIT]).abs()
    settled = _bias_z(largest).abs()
    inward = _bias_z(tie[tie.method == INTERSECTION])
    checks = [
        ("split: fold-B estimates unbiased for the weighted population "
         "bounds", resid, lambda z, lim: z <= lim),
        ("split: unbiased at an exact tie at every n", split_tie,
         lambda z, lim: z <= lim),
        (f"split: bias zero at the largest n for gaps >= {SETTLED_SES} SE_1",
         settled, lambda z, lim: z <= lim),
        ("intersection: biased inwards at an exact tie at every n", inward,
         lambda z, lim: z >= lim),
    ]
    all_ok = True
    for name, z, passes in checks:
        limit = z_limit(len(z))
        ok = bool(passes(z, limit).all())
        all_ok &= ok
        print(f"{'PASS' if ok else 'FAIL'} {name}: |z| in "
              f"[{z.abs().min():.2f}, {z.abs().max():.2f}] "
              f"(limit {limit:.2f}, {len(z)} tests)")
    return all_ok


def standardised_bias(df: pd.DataFrame) -> pd.DataFrame:
    """Adds out_L and out_U: outward bias of each end over its SD.

    out_L = (true_lo - E[L_hat]) / SD(L_hat) and out_U = (E[U_hat] -
    true_hi) / SD(U_hat); positive means the estimated band is too wide.
    """
    return df.assign(out_L=(df.true_lo - df.L_mean) / df.L_sd,
                     out_U=(df.U_mean - df.true_hi) / df.U_sd)


def figure_panels(
    gaps: Sequence[float],
) -> list[toy_sim_plots.ConvergencePanel]:
    """A difference in population means (several g) and equal means."""
    return [
        toy_sim_plots.ConvergencePanel(
            "Difference in population means, $g > 0$", tuple(gaps), (SPLIT,),
            ribbon=False),
        toy_sim_plots.ConvergencePanel(
            "Equal population means, $g = 0$", (0.0,), (SPLIT,), ribbon=True),
    ]


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    defaults = toy_sim.ToyConfig()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="Seed of the random generator.")
    parser.add_argument("--reps", type=int, default=defaults.reps,
                        help="Monte Carlo replicates per (n_1, gap).")
    parser.add_argument("--gaps", type=float, nargs="+",
                        default=list(DEFAULT_GAPS),
                        help="Population gaps g > 0 of the bound figure.")
    parser.add_argument("--uniform-gaps", type=float, nargs="+",
                        default=list(DEFAULT_UNIFORM_GAPS),
                        help="Small population gaps g > 0 of the uniform-"
                             "convergence figure.")
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
    if min(args.gaps) <= 0 or min(args.uniform_gaps) <= 0:
        parser.error("--gaps and --uniform-gaps must be > 0 (g = 0 is "
                     "always included)")
    return args


def main() -> None:
    args = parse_args()
    base = toy_sim.ToyConfig(reps=args.reps, frac_a=args.frac_a)
    rng = np.random.default_rng(args.seed)
    gaps = sorted({0.0, *args.gaps, *args.uniform_gaps})
    df = standardised_bias(simulate(base, sorted(args.n1), gaps, args.ratio,
                                    rng))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    results_path = args.out_dir / "results.csv"
    df.to_csv(results_path, index=False)
    run_checks(df)
    fig_path = args.out_dir / "bound_convergence.pdf"
    toy_sim_plots.plot_bound_convergence(df, figure_panels(args.gaps),
                                         fig_path)
    uniform_path = args.out_dir / "uniform_convergence.pdf"
    toy_sim_plots.plot_standardised_bias(
        df, (SPLIT, INTERSECTION), (0.0, *sorted(args.uniform_gaps)),
        uniform_path,
    )
    print(f"Written to {results_path}\nWritten to {fig_path}\n"
          f"Written to {uniform_path}")


if __name__ == "__main__":
    main()
