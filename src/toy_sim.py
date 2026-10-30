"""Toy example with a known true band [theta - delta, theta + delta].

Two datasets each identify a band that is wider than the truth: dataset 1
attains the true lower end, dataset 2 the true upper end, so their intersection
is exactly the true band. We compare the naive intersection CI with the convex
combination CI (oracle weights, sample-split weights, plug-in weights without
splitting), and with a calibrated intersection CI that uses
z = Phi^{-1}(1 - alpha / 4): Bonferroni over 2 datasets x 2 ends, so it is
valid whichever dataset attains each end.

Usage, from the repo root:
    python -m src.toy_sim [--seed 7] [--reps 5000] [--frac-a 0.3] [--no-plots]

Plans: quality_reports/plans/2026-10-09_band-toy-simulation.md,
       quality_reports/plans/2026-10-11_band-toy-refactor.md,
       quality_reports/plans/2026-10-27_calibrated-intersection.md
Outputs: src/outputs/toy_sim/
"""

import argparse
import dataclasses
import itertools
from collections.abc import Sequence
from pathlib import Path
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats

from src.helpers import ci_methods
from src.helpers import mc_metrics
from src.helpers import weight_optimisers
from src.plotting import toy_sim_plots

DEFAULT_SEED = 7
DEFAULT_OUT_DIR = Path(__file__).resolve().parent / "outputs" / "toy_sim"
PROGRESS_EVERY = 18  # cells between progress messages

M1 = "M1 naive intersection"
M2 = "M2 convex (oracle)"
M3 = "M3 convex (split)"
M4 = "M4 convex (plug-in)"
M5_D1 = "M5 dataset 1 only"
M5_D2 = "M5 dataset 2 only"
M6 = "M6 intersection (calibrated)"
PLOT_METHODS = (M1, M2, M3, M4)


# ----------------------------- configuration ------------------------------
@dataclasses.dataclass(frozen=True)
class ToyConfig:
    """Settings for the toy simulation.

    Attributes:
        theta: True estimand.
        delta: Half-width of the true band.
        alpha: Nominal level. Each end uses z = Phi^{-1}(1 - alpha / 2),
            Bonferroni for the joint band; the calibrated intersection uses
            z_calibrated = Phi^{-1}(1 - alpha / 4).
        n: Dataset sizes (n_1, n_2).
        se1: Standard deviation of Ybar_1.
        frac_a: Fold-A fraction for the sample-split method.
        reps: Monte Carlo replicates per cell.
        gaps: Values for g_L = a_2 and g_U = c_1; must include 0 and 1.
        ratios: Values for r = SE_2 / SE_1; must include 1.
    """

    theta: float = 2.0
    delta: float = 1.0
    alpha: float = 0.05
    n: tuple[int, int] = (2000, 500)
    se1: float = 0.1
    frac_a: float = 0.3
    reps: int = 5000
    gaps: tuple[float, ...] = (0.0, 0.05, 0.1, 0.2, 0.5, 1.0)
    ratios: tuple[float, ...] = (1, 2, 5, 10)

    def __post_init__(self):
        if not 0 < self.alpha < 1:
            raise ValueError(f"alpha must be in (0, 1), got {self.alpha=}")
        if not 0 < self.frac_a < 1:
            raise ValueError(f"frac_a must be in (0, 1), got {self.frac_a=}")
        if self.reps < 2:
            raise ValueError(f"reps must be >= 2, got {self.reps=}")
        if self.delta < 0 or self.se1 <= 0:
            raise ValueError(
                f"need delta >= 0 and se1 > 0, got {self.delta=}, {self.se1=}"
            )

    @property
    def z(self) -> float:
        return stats.norm.ppf(1 - self.alpha / 2)

    @property
    def z_calibrated(self) -> float:
        """Bonferroni over the 4 one-sided (dataset, end) events."""
        return stats.norm.ppf(1 - self.alpha / 4)

    @property
    def true_lo(self) -> float:
        return self.theta - self.delta

    @property
    def true_hi(self) -> float:
        return self.theta + self.delta


# ----------------------------- DGP ------------------------------
class DGPParams(NamedTuple):
    """Population quantities for one cell; arrays are indexed by dataset."""

    mu: np.ndarray
    Delta: np.ndarray
    sigma: np.ndarray
    se: np.ndarray
    L_star: np.ndarray
    U_star: np.ndarray


class BoundEstimates(NamedTuple):
    """Estimated bounds and the variances the CI uses, each shape (R, 2)."""

    L: np.ndarray
    U: np.ndarray
    s2: np.ndarray


def dgp_params(g_L: float, g_U: float, r: float, cfg: ToyConfig) -> DGPParams:
    """Slack a = (0, g_L), c = (g_U, 0) around the true band."""
    a = np.array([0.0, g_L])
    c = np.array([g_U, 0.0])
    mu = cfg.theta + (c - a) / 2
    Delta = cfg.delta + (a + c) / 2
    se = np.array([cfg.se1, r * cfg.se1])
    sigma = se * np.sqrt(cfg.n)
    return DGPParams(mu=mu, Delta=Delta, sigma=sigma, se=se,
                     L_star=mu - Delta, U_star=mu + Delta)


def draw_standard_normals(
    n: Sequence[int], reps: int, rng: np.random.Generator
) -> list[np.ndarray]:
    """Draws the (reps, n_k) standard normals that every cell reuses.

    Common random numbers: differences between cells are then not MC noise.
    """
    return [rng.standard_normal((reps, n_k)) for n_k in n]


def estimate_bounds(
    Y: Sequence[np.ndarray], Delta: np.ndarray, var_sizes: Sequence[int]
) -> BoundEstimates:
    """Bounds Ybar_k -/+ Delta_k, with variance var(Y_k) / var_sizes[k].

    Args:
        Y: Outcomes for each dataset, shapes (R, m_k).
        Delta: Half-width of each dataset's identified band, shape (2,).
        var_sizes: Sample size dividing each sample variance.
    """
    Ybar = np.column_stack([y.mean(axis=1) for y in Y])
    s2 = np.column_stack(
        [y.var(axis=1, ddof=1) / m for y, m in zip(Y, var_sizes)]
    )
    return BoundEstimates(L=Ybar - Delta, U=Ybar + Delta, s2=s2)


# ----------------------------- estimation ------------------------------
def oracle_weights(params: DGPParams, z: float) -> tuple[float, float]:
    """Optimal (w_L, w_U) from the population bounds and standard errors."""
    se2 = params.se**2
    w_L = weight_optimisers.optimal_weight_grid(
        params.L_star[0], params.L_star[1], se2[0], se2[1], lower=True, z=z
    )[0]
    w_U = weight_optimisers.optimal_weight_grid(
        params.U_star[0], params.U_star[1], se2[0], se2[1], lower=False, z=z
    )[0]
    return w_L, w_U


def plug_in_weights(
    est: BoundEstimates, z: float
) -> tuple[np.ndarray, np.ndarray]:
    """Optimal (w_L, w_U) per replicate from estimated bounds, each (R,)."""
    w_L = weight_optimisers.optimal_weight_grid(
        est.L[:, 0], est.L[:, 1], est.s2[:, 0], est.s2[:, 1], lower=True, z=z
    )
    w_U = weight_optimisers.optimal_weight_grid(
        est.U[:, 0], est.U[:, 1], est.s2[:, 0], est.s2[:, 1], lower=False, z=z
    )
    return w_L, w_U


def ci_by_method(
    Y: Sequence[np.ndarray], params: DGPParams, cfg: ToyConfig
) -> tuple[dict[str, ci_methods.CIResult], tuple[float, float]]:
    """CIs for every method, and the oracle weights (w_L, w_U)."""
    z = cfg.z
    n_a = [int(cfg.frac_a * n_k) for n_k in cfg.n]
    n_b = [n_k - a for n_k, a in zip(cfg.n, n_a)]
    full = estimate_bounds(Y, params.Delta, cfg.n)
    # Rows are i.i.d., so the first n_a columns are a random split. Fold-A
    # variance is rescaled to the fold-B size: the variance the CI will use.
    fold_a = estimate_bounds(
        [y[:, :a] for y, a in zip(Y, n_a)], params.Delta, n_b
    )
    fold_b = estimate_bounds(
        [y[:, a:] for y, a in zip(Y, n_a)], params.Delta, n_b
    )
    w_oracle = oracle_weights(params, z)
    cis = {
        M1: ci_methods.intersection_ci(full.L, full.U, full.s2, z),
        M2: ci_methods.convex_ci(full.L, full.U, full.s2, *w_oracle, z),
        M3: ci_methods.convex_ci(
            fold_b.L, fold_b.U, fold_b.s2, *plug_in_weights(fold_a, z), z
        ),
        M4: ci_methods.convex_ci(
            full.L, full.U, full.
            s2, *plug_in_weights(full, z), z
        ),
        M5_D1: ci_methods.convex_ci(full.L, full.U, full.s2, 1.0, 1.0, z),
        M5_D2: ci_methods.convex_ci(full.L, full.U, full.s2, 0.0, 0.0, z),
        M6: ci_methods.intersection_ci(
            full.L, full.U, full.s2, cfg.z_calibrated
        ),
    }
    return cis, w_oracle


def oracle_own_target_coverage(
    ci: ci_methods.CIResult, params: DGPParams, w_L: float, w_U: float
) -> dict[str, float]:
    """Per-end coverage of M2's own target, the population convex bounds."""
    own_lo = w_L * params.L_star[0] + (1 - w_L) * params.L_star[1]
    own_hi = w_U * params.U_star[0] + (1 - w_U) * params.U_star[1]
    row = {}
    row["cov_lo_own"], row["cov_lo_own_mcse"] = mc_metrics.mean_se(
        ci.lo <= own_lo
    )
    row["cov_hi_own"], row["cov_hi_own_mcse"] = mc_metrics.mean_se(
        ci.hi >= own_hi
    )
    return row


def simulate_cell(
    g_L: float,
    g_U: float,
    r: float,
    eps: Sequence[np.ndarray],
    cfg: ToyConfig,
) -> list[dict]:
    """Summary rows, one per method, for one (g_L, g_U, r) cell.

    Args:
        g_L: Lower-end gap of dataset 2.
        g_U: Upper-end gap of dataset 1.
        r: SE_2 / SE_1.
        eps: Standard normal draws per dataset, shapes (reps, n_k).
        cfg: Simulation settings.
    """
    params = dgp_params(g_L, g_U, r, cfg)
    Y = [m + s * e for m, s, e in zip(params.mu, params.sigma, eps)]
    cis, (w_L_oracle, w_U_oracle) = ci_by_method(Y, params, cfg)
    rows = []
    for method, ci in cis.items():
        row = dict(g_L=g_L, g_U=g_U, r=r, method=method,
                   wL_oracle=w_L_oracle, wU_oracle=w_U_oracle)
        row.update(mc_metrics.summarise_ci(ci, cfg.true_lo, cfg.true_hi,
                                           cfg.theta))
        if method == M2:
            row.update(
                oracle_own_target_coverage(ci, params, w_L_oracle, w_U_oracle)
            )
        rows.append(row)
    return rows


def simulate_grid(eps: Sequence[np.ndarray], cfg: ToyConfig) -> pd.DataFrame:
    """Runs every (r, g_L, g_U) cell on the same draws."""
    rows = []
    cells = list(itertools.product(cfg.ratios, cfg.gaps, cfg.gaps))
    for i, (r, g_L, g_U) in enumerate(cells, 1):
        rows.extend(simulate_cell(g_L, g_U, r, eps, cfg))
        if i % PROGRESS_EVERY == 0:
            print(f"cell {i}/{len(cells)}")
    return pd.DataFrame(rows)


# ----------------------------- sanity checks ------------------------------
class Check(NamedTuple):
    name: str
    ok: bool
    detail: str


def _cell(df: pd.DataFrame, method: str, g_L: float, g_U: float,
          r: float) -> pd.Series:
    match = (df.method == method) & (df.g_L == g_L) & (df.g_U == g_U)
    return df[match & (df.r == r)].iloc[0]


def _check_closed_forms(df: pd.DataFrame, cfg: ToyConfig) -> list[Check]:
    """M1 and M2 at g_L = 0, r = 1 against closed-form coverage and bias."""
    z = cfg.z
    m1 = _cell(df, M1, 0.0, 0.0, 1)
    m2 = _cell(df, M2, 0.0, 0.0, 1)
    cov_m1, bias_m1 = stats.norm.cdf(z) ** 2, cfg.se1 / np.sqrt(np.pi)
    cov_m2 = stats.norm.cdf(z)
    return [
        Check("M1 lower-end coverage at g_L=0, r=1 vs Phi(z)^2",
              abs(m1.cov_lo - cov_m1) <= 3 * m1.cov_lo_mcse,
              f"{m1.cov_lo:.4f} vs {cov_m1:.4f} (MC SE {m1.cov_lo_mcse:.4f})"),
        Check("M1 lower bias at g_L=0, r=1 vs SE/sqrt(pi)",
              abs(m1.bias_L - bias_m1) <= 3 * m1.bias_L_mcse,
              f"{m1.bias_L:.4f} vs {bias_m1:.4f} "
              f"(MC SE {m1.bias_L_mcse:.4f})"),
        Check("M2 lower-end coverage at g_L=0, r=1 vs Phi(z)",
              abs(m2.cov_lo - cov_m2) <= 3 * m2.cov_lo_mcse,
              f"{m2.cov_lo:.4f} vs {cov_m2:.4f} (MC SE {m2.cov_lo_mcse:.4f})"),
    ]


def _check_oracle_weights_at_tie(
    df: pd.DataFrame, cfg: ToyConfig
) -> list[Check]:
    """At g_L = g_U = 0 the oracle weights are inverse-variance weights."""
    checks = []
    for r in cfg.ratios:
        target = (r * cfg.se1) ** 2 / (cfg.se1**2 + (r * cfg.se1) ** 2)
        m2 = _cell(df, M2, 0.0, 0.0, r)
        ok = (abs(m2.wL_oracle - target) <= 1e-3
              and abs(m2.wU_oracle - target) <= 1e-3)
        checks.append(Check(
            f"M2 inverse-variance weights at a tie, r={r}", ok,
            f"w_L={m2.wL_oracle:.3f}, w_U={m2.wU_oracle:.3f} vs {target:.3f}"))
    return checks


def _check_oracle_own_target(df: pd.DataFrame, cfg: ToyConfig) -> list[Check]:
    """M2's per-end coverage of its own target is Phi(z) in every cell."""
    # Coverage of the true band is >= Phi(z) when w is interior (the fused
    # population bound is looser than the truth), so check M2's own target.
    m2 = df[df.method == M2]
    nominal = stats.norm.cdf(cfg.z)
    zs = np.concatenate([(m2.cov_lo_own - nominal) / m2.cov_lo_own_mcse,
                         (m2.cov_hi_own - nominal) / m2.cov_hi_own_mcse])
    return [Check(
        "M2 per-end coverage of its own target ~ Phi(z) in every cell",
        np.max(np.abs(zs)) <= 4,
        f"max |z| = {np.max(np.abs(zs)):.2f} over {len(zs)} ends; "
        f"{np.sum(np.abs(zs) > 3)} with |z| > 3 "
        "(cells share draws, so not independent)")]


def _check_m1_matches_m2_at_large_gap(
    df: pd.DataFrame, cfg: ToyConfig
) -> list[Check]:
    """At g_L = g_U = 1, M2 puts all weight on one dataset per end, like M1.

    Only for r with max gap >= z * SE_2, the threshold for w_U* = 0 (paper,
    toy Remark); for larger r the oracle upper weight is still interior.
    """
    checks = []
    for r in cfg.ratios:
        if max(cfg.gaps) < cfg.z * r * cfg.se1:
            continue
        m1, m2 = _cell(df, M1, 1.0, 1.0, r), _cell(df, M2, 1.0, 1.0, r)
        ok = (m2.wL_oracle >= 0.95 and m2.wU_oracle <= 0.05
              and abs(m1.cov_joint - m2.cov_joint) <= 0.01
              and abs(m1.width - m2.width) / m2.width <= 0.02)
        checks.append(Check(
            f"M1 ~ M2 at g_L=g_U=1, r={r}", ok,
            f"w_L={m2.wL_oracle:.3f}, w_U={m2.wU_oracle:.3f}; coverage "
            f"{m1.cov_joint:.4f} vs {m2.cov_joint:.4f}; width "
            f"{m1.width:.4f} vs {m2.width:.4f}"))
    return checks


def _check_coverage_separable(df: pd.DataFrame, cfg: ToyConfig) -> list[Check]:
    """Lower-end coverage cannot depend on g_U, nor upper-end on g_L."""
    worst = 0.0
    for method in PLOT_METHODS:
        sub = df[df.method == method]
        worst = max(worst,
                    sub.groupby(["r", "g_L"]).cov_lo.agg(np.ptp).max(),
                    sub.groupby(["r", "g_U"]).cov_hi.agg(np.ptp).max())
    return [Check("Lower-end coverage constant in g_U (and upper in g_L)",
                  worst <= 2 / cfg.reps,
                  f"max range {worst:.5f} (allowed {2 / cfg.reps:.5f})")]


def _check_calibrated_valid(df: pd.DataFrame, cfg: ToyConfig) -> list[Check]:
    """M6 joint coverage is >= 1 - alpha in every cell (up to MC error)."""
    m6 = df[df.method == M6]
    nominal = 1 - cfg.alpha
    # Cells with coverage 1 have MC SE 0 and give z = +inf, which passes.
    zs = (m6.cov_joint - nominal) / m6.cov_joint_mcse
    return [Check(
        f"M6 joint coverage >= {nominal:g} in every cell",
        zs.min() >= -3,
        f"min coverage {m6.cov_joint.min():.4f}, min z = {zs.min():.2f} "
        f"over {len(zs)} cells")]


def run_sanity_checks(df: pd.DataFrame, cfg: ToyConfig) -> bool:
    """Prints PASS/FAIL for each check and returns True if all pass."""
    checks = []
    for check_fn in (_check_closed_forms, _check_oracle_weights_at_tie,
                     _check_oracle_own_target,
                     _check_m1_matches_m2_at_large_gap,
                     _check_coverage_separable, _check_calibrated_valid):
        checks.extend(check_fn(df, cfg))
    print("\nSanity checks")
    for check in checks:
        print(f"[{'PASS' if check.ok else 'FAIL'}] {check.name}: "
              f"{check.detail}")
    n_passed = sum(check.ok for check in checks)
    print(f"{n_passed}/{len(checks)} checks passed")
    return n_passed == len(checks)


# ----------------------------- output ------------------------------
def print_diagonal(df: pd.DataFrame) -> None:
    """Prints the g_L = g_U cells for the plotted methods."""
    show = df[(df.g_L == df.g_U) & df.method.isin(PLOT_METHODS)]
    columns = ["r", "g_L", "method", "cov_joint", "cov_lo", "cov_hi", "width",
               "wL", "wU"]
    with pd.option_context("display.width", 140,
                           "display.float_format", "{:.3f}".format):
        print("\nDiagonal g_L = g_U (full results in results.csv)")
        print(show[columns].to_string(index=False))


def make_plots(df: pd.DataFrame, cfg: ToyConfig, out_dir: Path) -> None:
    """Writes figures 1-4 to out_dir."""
    styles = toy_sim_plots.method_styles(PLOT_METHODS)
    grid = dict(gaps=cfg.gaps, ratios=cfg.ratios)
    toy_sim_plots.plot_per_end(
        df, "cov_{end_name}", stats.norm.cdf(cfg.z), "Coverage of {end} end",
        PLOT_METHODS, styles, reps=cfg.reps,
        out_path=out_dir / "fig1_per_end_coverage.png", **grid)

    cov_methods = (M1, M3)
    cov = [[toy_sim_plots.pivot_cells(df, m, r, "cov_joint", cfg.gaps)
            for r in cfg.ratios] for m in cov_methods]
    toy_sim_plots.heatmap_grid(
        cov, cov_methods,
        f"Joint coverage of the true band (centre {1 - cfg.alpha:g})",
        1 - cfg.alpha, "{:.3f}",
        out_path=out_dir / "fig2_joint_coverage.png", **grid)

    # Ratios against M6, not M1: M1 undercovers at small gaps, so widths
    # relative to it are not like-for-like.
    def width_ratio(method, r):
        return (toy_sim_plots.pivot_cells(df, method, r, "width", cfg.gaps)
                / toy_sim_plots.pivot_cells(df, M6, r, "width", cfg.gaps))

    ratios = [[width_ratio(m, r) for r in cfg.ratios] for m in (M3, M2)]
    toy_sim_plots.heatmap_grid(
        ratios, ["M3 / M6", "M2 / M6"], "Mean width ratio (centre 1)", 1.0,
        "{:.2f}", out_path=out_dir / "fig3_width_ratio.png", **grid)

    toy_sim_plots.plot_per_end(
        df, "w{end}", None, "Mean weight on dataset 1 ({end} end)",
        PLOT_METHODS[1:], styles, reps=cfg.reps,
        out_path=out_dir / "fig4_weights.png", **grid)


# ----------------------------- main ------------------------------
def run(cfg: ToyConfig, seed: int, out_dir: Path, plots: bool) -> pd.DataFrame:
    """Simulates every cell, saves results.csv and (optionally) figures."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    eps = draw_standard_normals(cfg.n, cfg.reps, rng)
    df = simulate_grid(eps, cfg)
    df.to_csv(out_dir / "results.csv", index=False)
    run_sanity_checks(df, cfg)
    print_diagonal(df)
    if plots:
        make_plots(df, cfg, out_dir)
    print(f"\nWritten to {out_dir}")
    return df


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    defaults = ToyConfig()
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED,
                        help="Seed for the shared standard normal draws.")
    parser.add_argument("--reps", type=int, default=defaults.reps,
                        help="Monte Carlo replicates per cell.")
    parser.add_argument("--frac-a", type=float, default=defaults.frac_a,
                        help="Fold-A fraction for the sample-split method.")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR,
                        help="Directory for results.csv and figures.")
    parser.add_argument("--no-plots", action="store_true",
                        help="Skip the figures.")
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    cfg = ToyConfig(reps=args.reps, frac_a=args.frac_a)
    run(cfg, args.seed, args.out_dir, plots=not args.no_plots)


if __name__ == "__main__":
    main()
