"""Check: analytical against simulated critical values for CLR with selection.

For one cell (Lambda_1, Lambda_2) of a causal run, at both ends of every
replicate, the CLR intersection and convex CLR CIs with adaptive inequality
selection are computed twice from the saved full-sample estimates:

    analytical  as in src/helpers/ci_methods.py: k_gamma and the kept-set
                critical value k_p from closed forms (two bounds: quantiles of
                the max of independent N(0, 1); all convex weights: the arc
                formula of _arc_exceedance);
    simulated   every critical value is instead a quantile of simulated
                sup_{kept w} a(w)'Z over the weight grid, Z ~ N(0, I_2), as in
                CLR's simulation-based algorithm, with fresh draws for every
                replicate and end.

Fresh draws make the simulation errors independent across replicate-ends, so
the average of (simulated - analytical) / MC s.e. over them tests the formulas
(about 0 if they are right, with a spread of about 1); with shared draws every
replicate would carry the same simulation error. Where the simulated k_gamma
puts a bound or weight on the other side of the selection threshold, the two
versions keep different sets and their critical values are not comparable;
those replicate-ends are counted and left out of the k_p and CI end comparisons.

Both use the same weight grid (--grid-size, smaller than the 1001 weights of
the simulation runs, for speed); the summary also compares the analytical CIs
on this grid with those on the 1001-weight grid. The Monte Carlo standard error
of a simulated quantile is sqrt(level (1 - level) / draws) / f(k), with f the
analytical density of the sup at k. Default run: about 2 minutes.

Usage, from the repo root (needs the run's estimates.csv and true_bands.csv):
    python checks/clr_selection_check.py [--run causal_sim_quick_continuous]
        [--outcome continuous] [--lambda-1 5] [--lambda-2 1.25]
        [--draws 200000] [--grid-size 201] [--seed 1]
Outputs, in checks/outputs/:
    clr_selection_<run>_<outcome>_l1-<L1>_l2-<L2>.csv  per replicate and end:
        analytical and simulated k_gamma, kept count, k_p and CI end
    clr_selection_<run>_<outcome>_l1-<L1>_l2-<L2>.txt  summary
"""

import argparse
import re
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from src.helpers import ci_methods

REPO_ROOT = Path(__file__).resolve().parents[1]
OUTPUTS = REPO_ROOT / "src" / "outputs"
CHECK_OUTPUTS = Path(__file__).resolve().parent / "outputs"
LAMBDA_TOL = 1e-9
FULL_GRID_SIZE = 1001  # weight grid of the simulation runs
DENSITY_STEP = 1e-4  # finite-difference step for the density of the sup
DRAW_CHUNK = 25_000  # draws per matrix product, to bound memory


def load_cell(
    run_dir: Path, outcome: str, lambda_1: float, lambda_2: float
) -> tuple[dict[str, tuple[np.ndarray, np.ndarray]], tuple[float, float]]:
    """Each end's (theta, s2), R x 2, and the cell's true band.

    The upper end is min_k U_k; the lower end is handled as the upper end of
    -L, as in ci_methods.
    """
    estimates = pd.read_csv(run_dir / "estimates.csv")
    full = estimates[(estimates.fit == "full") & (estimates.outcome == outcome)]
    truth = pd.read_csv(run_dir / "true_bands.csv")
    studies, bands = [], []
    for study, lam in ((1, lambda_1), (2, lambda_2)):
        rows = full[(full.study == study)
                    & np.isclose(full["lambda"], lam, atol=LAMBDA_TOL)]
        band = truth[(truth.study == study) & (truth.outcome == outcome)
                     & np.isclose(truth["lambda"], lam, atol=LAMBDA_TOL)]
        if rows.empty or band.empty:
            raise ValueError(f"No study {study} rows at Lambda = {lam}")
        studies.append(rows.sort_values("rep"))
        bands.append(band.iloc[0])
    pair = lambda column: np.c_[studies[0][column].to_numpy(),
                                studies[1][column].to_numpy()]
    ends = {"upper": (pair("upper"), pair("s2_upper")),
            "lower": (-pair("lower"), pair("s2_lower"))}
    true_band = (max(b.lower for b in bands), min(b.upper for b in bands))
    return ends, true_band


def run_total_n(run_info: Path) -> int:
    """n_1 + n_2 from a run's run_info.txt ('... | n: 3000, 500 | ...').

    Raises:
        ValueError: If the line is missing.
    """
    match = re.search(r"\| n: (\d+), (\d+)", run_info.read_text())
    if match is None:
        raise ValueError(f"No 'n: n_1, n_2' in {run_info}")
    return int(match[1]) + int(match[2])


def directions(se: np.ndarray, w_grid: np.ndarray) -> np.ndarray:
    """Standardised-error directions a(w) of one replicate, shape (W, 2)."""
    raw = np.c_[w_grid * se[0], (1 - w_grid) * se[1]]
    return raw / np.linalg.norm(raw, axis=1, keepdims=True)


def simulated_sups(
    draws: np.ndarray, dirs: np.ndarray, kept: np.ndarray | None = None
) -> np.ndarray:
    """sup over the (kept) directions of a'Z for every draw, shape (N,)."""
    use = dirs if kept is None else dirs[kept]
    return np.concatenate([
        (draws[i:i + DRAW_CHUNK] @ use.T).max(axis=1)
        for i in range(0, draws.shape[0], DRAW_CHUNK)])


def arc_density(level_k: np.ndarray, arc_length: np.ndarray,
                dropped_arcs: np.ndarray) -> np.ndarray:
    """Density of the sup over the kept arcs at k, by central difference."""
    up = ci_methods._arc_exceedance(level_k + DENSITY_STEP, arc_length,
                                    dropped_arcs)
    down = ci_methods._arc_exceedance(level_k - DENSITY_STEP, arc_length,
                                      dropped_arcs)
    return (down - up) / (2 * DENSITY_STEP)


def convex_end(
    theta: np.ndarray, s2: np.ndarray, w_grid: np.ndarray, gamma: float,
    p: float, rng: np.random.Generator, n_draws: int,
) -> pd.DataFrame:
    """Convex CLR end with selection, analytical and simulated, per replicate.

    Returns:
        One row per replicate: k_gamma, kept count, k_p and CI end for both
        versions, the analytical MC s.e. of the simulated k_p, and the
        analytical CI end on the full 1001-weight grid.
    """
    se = np.sqrt(s2)
    family, s = ci_methods._convex_family(theta, s2, w_grid)
    k_gamma = ci_methods._quadrant_sup_quantile(gamma)
    kept = ci_methods._kept_by_selection(family, s, k_gamma)
    arc_length, dropped_arcs = ci_methods._kept_arcs(kept, se, w_grid)
    k_p = ci_methods._arc_quantile(p, arc_length, dropped_arcs)
    end, _ = ci_methods._corrected_minimum(family, s, k_p)
    mc_se = (np.sqrt(p * (1 - p) / n_draws)
             / arc_density(k_p, arc_length, dropped_arcs))
    full_arc = np.array([np.pi / 2])
    k_gamma_mc_se = (np.sqrt(gamma * (1 - gamma) / n_draws)
                     / arc_density(np.array([k_gamma]), full_arc,
                                   np.zeros((1, 1)))[0])
    sim = {"k_gamma": [], "kept": [], "k_p": [], "end": []}
    for r in range(theta.shape[0]):
        draws = rng.standard_normal((n_draws, 2))
        dirs = directions(se[r], w_grid)
        k_gamma_sim = np.quantile(simulated_sups(draws, dirs), gamma)
        floor = np.min(family[r] + k_gamma_sim * s[r])
        kept_sim = family[r] <= floor + 2 * k_gamma_sim * s[r]
        k_p_sim = np.quantile(simulated_sups(draws, dirs, kept_sim), p)
        sim["k_gamma"].append(k_gamma_sim)
        sim["kept"].append(kept_sim.sum())
        sim["k_p"].append(k_p_sim)
        sim["end"].append(np.min(family[r] + k_p_sim * s[r]))
    full_grid = np.linspace(0, 1, FULL_GRID_SIZE)
    family_full, s_full = ci_methods._convex_family(theta, s2, full_grid)
    end_full, _, _, _ = ci_methods._clr_upper_selected_convex(
        family_full, s_full, se, full_grid, k_gamma, p)
    return pd.DataFrame({
        "k_gamma_analytic": k_gamma, "k_gamma_sim": sim["k_gamma"],
        "k_gamma_mc_se": k_gamma_mc_se,
        "kept_analytic": kept.sum(axis=1), "kept_sim": sim["kept"],
        "k_p_analytic": k_p, "k_p_sim": sim["k_p"], "k_p_mc_se": mc_se,
        "end_analytic": end, "end_sim": sim["end"],
        "end_analytic_full_grid": end_full, "se_min": se.min(axis=1)})


def max_density(k: np.ndarray, m: np.ndarray) -> np.ndarray:
    """Density at k of the max of m independent N(0, 1), m in {1, 2}."""
    return np.where(m == 2, 2 * stats.norm.cdf(k) * stats.norm.pdf(k),
                    stats.norm.pdf(k))


def intersection_end(
    theta: np.ndarray, s2: np.ndarray, gamma: float, p: float,
    rng: np.random.Generator, n_draws: int,
) -> pd.DataFrame:
    """CLR intersection end with selection, analytical and simulated.

    The standardised errors of the two bounds are the coordinates of Z, so
    the simulated sup over the kept bounds is the max of the kept
    coordinates.
    """
    se = np.sqrt(s2)
    k_gamma = ci_methods.max_of_two_quantile(gamma)
    kept = ci_methods._kept_by_selection(theta, se, k_gamma)
    m = kept.sum(axis=1)
    k_p = stats.norm.ppf(p ** (1 / m))
    end, _ = ci_methods._corrected_minimum(theta, se, k_p)
    k_gamma_sim, kept_sim, k_p_sim = [], [], []
    for r in range(theta.shape[0]):
        draws = rng.standard_normal((n_draws, 2))
        k_gamma_r = np.quantile(draws.max(axis=1), gamma)
        floor = np.min(theta[r] + k_gamma_r * se[r])
        kept_r = theta[r] <= floor + 2 * k_gamma_r * se[r]
        k_gamma_sim.append(k_gamma_r)
        kept_sim.append(kept_r.sum())
        k_p_sim.append(np.quantile(draws[:, kept_r].max(axis=1), p))
    k_p_sim = np.array(k_p_sim)
    end_sim, _ = ci_methods._corrected_minimum(theta, se, k_p_sim)
    mc_se = np.sqrt(p * (1 - p) / n_draws) / max_density(k_p, m)
    k_gamma_mc_se = (np.sqrt(gamma * (1 - gamma) / n_draws)
                     / max_density(np.array([k_gamma]), np.array([2]))[0])
    return pd.DataFrame({
        "k_gamma_analytic": k_gamma, "k_gamma_sim": k_gamma_sim,
        "k_gamma_mc_se": k_gamma_mc_se,
        "kept_analytic": m, "kept_sim": kept_sim,
        "k_p_analytic": k_p, "k_p_sim": k_p_sim, "k_p_mc_se": mc_se,
        "end_analytic": end, "end_sim": end_sim,
        "end_analytic_full_grid": end, "se_min": se.min(axis=1)})


def summarise(per_end: pd.DataFrame, true_band: tuple[float, float]) -> str:
    """Summary lines per method: critical values, kept sets, CI ends, and
    joint coverage and mean width with each version."""
    lines = []
    for method, rows in per_end.groupby("method", sort=False):
        z_gamma = ((rows.k_gamma_sim - rows.k_gamma_analytic)
                   / rows.k_gamma_mc_se)
        same = rows.kept_sim == rows.kept_analytic
        agree = rows[same]
        z = (agree.k_p_sim - agree.k_p_analytic) / agree.k_p_mc_se
        end_diff = (agree.end_sim - agree.end_analytic) / agree.se_min
        grid_diff = ((rows.end_analytic - rows.end_analytic_full_grid)
                     / rows.se_min)
        dropped = rows.kept_analytic < rows.kept_analytic.max()
        lines += [
            f"{method}",
            f"  k_gamma: analytic {rows.k_gamma_analytic.iloc[0]:.4f}; "
            f"(simulated - analytic) / MC s.e. over {len(rows)} "
            f"replicate-ends: mean {z_gamma.mean():+.3f} (s.e. "
            f"{1 / np.sqrt(len(rows)):.3f}), sd {z_gamma.std():.3f}",
            f"  some bound or weight dropped (analytic) in {dropped.sum()} "
            f"of {len(rows)} replicate-ends; kept sets differ in "
            f"{(~same).sum()} (selection threshold crossed by the simulated "
            f"k_gamma)",
            f"  k_p where kept sets agree ({same.sum()}): mean analytic "
            f"{agree.k_p_analytic.mean():.4f}, simulated "
            f"{agree.k_p_sim.mean():.4f}; (simulated - analytic) / MC s.e.: "
            f"mean {z.mean():+.3f} (s.e. {1 / np.sqrt(len(z)):.3f}), sd "
            f"{z.std():.3f}, max |.| {np.abs(z).max():.2f}; median MC s.e. "
            f"{agree.k_p_mc_se.median():.4f}",
            f"  the same, replicate-ends where something is dropped "
            f"({(same & dropped).sum()}): mean z "
            f"{z[dropped[same]].mean():+.3f}, sd {z[dropped[same]].std():.3f}",
            f"  CI end where kept sets agree, (simulated - analytic) / "
            f"smaller s.e.: mean {end_diff.mean():+.4f}, max |.| "
            f"{np.abs(end_diff).max():.4f}",
            f"  CI end, analytic on this grid - on 1001 weights, / smaller "
            f"s.e.: max |.| {np.abs(grid_diff).max():.2e}",
        ]
        for version in ("analytic", "sim"):
            upper = rows[rows.end == "upper"][f"end_{version}"].to_numpy()
            lower = -rows[rows.end == "lower"][f"end_{version}"].to_numpy()
            covered = (lower <= true_band[0]) & (upper >= true_band[1])
            lines.append(f"  {version:8s}: joint coverage {covered.mean():.4f},"
                         f" mean width {np.mean(upper - lower):.5f}")
    return "\n".join(lines) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--run", default="causal_sim_quick_continuous",
                        help="Folder in src/outputs/ with estimates.csv.")
    parser.add_argument("--outcome", choices=("continuous", "binary"),
                        default="continuous")
    parser.add_argument("--lambda-1", type=float, default=5.0)
    parser.add_argument("--lambda-2", type=float, default=1.25)
    parser.add_argument("--draws", type=int, default=200_000,
                        help="Simulated N(0, I_2) draws per replicate-end.")
    parser.add_argument("--grid-size", type=int, default=201,
                        help="Number of equally spaced weights in [0, 1].")
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=1)
    return parser.parse_args(argv)


def main() -> None:
    args = parse_args()
    run_dir = OUTPUTS / args.run
    ends, true_band = load_cell(run_dir, args.outcome, args.lambda_1,
                                args.lambda_2)
    n = run_total_n(run_dir / "run_info.txt")
    gamma = ci_methods.clr_selection_level(n)
    p = 1 - args.alpha / 2
    rng = np.random.default_rng(args.seed)
    w_grid = np.linspace(0, 1, args.grid_size)
    frames = []
    for end, (theta, s2) in ends.items():
        for method, frame in (
                ("CLR intersection with selection",
                 intersection_end(theta, s2, gamma, p, rng, args.draws)),
                ("convex CLR with selection",
                 convex_end(theta, s2, w_grid, gamma, p, rng, args.draws))):
            frames.append(frame.assign(method=method, end=end,
                                       rep=np.arange(1, theta.shape[0] + 1)))
    per_end = pd.concat(frames, ignore_index=True)
    CHECK_OUTPUTS.mkdir(exist_ok=True)
    stem = (f"clr_selection_{args.run}_{args.outcome}_l1-{args.lambda_1:g}"
            f"_l2-{args.lambda_2:g}")
    per_end.to_csv(CHECK_OUTPUTS / f"{stem}.csv", index=False)
    header = (f"Run {args.run}, {args.outcome}, Lambda_1 = {args.lambda_1:g}, "
              f"Lambda_2 = {args.lambda_2:g}; true band "
              f"[{true_band[0]:.5f}, {true_band[1]:.5f}]; gamma_n = "
              f"{gamma:.5f} (n = {n}); p = {p}; {args.draws} fresh draws per "
              f"replicate-end (seed {args.seed}); weight grid of "
              f"{args.grid_size}\n")
    summary = header + summarise(per_end, true_band)
    (CHECK_OUTPUTS / f"{stem}.txt").write_text(summary)
    print(summary)
    print(f"Written to {CHECK_OUTPUTS / stem}.csv and .txt")


if __name__ == "__main__":
    main()
