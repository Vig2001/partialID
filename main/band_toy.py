"""Toy example with a known true band [theta - delta, theta + delta].

Two datasets each identify a band that is wider than the truth: dataset 1 attains
the true lower end, dataset 2 the true upper end, so their intersection is exactly
the true band. We compare the naive intersection CI with the convex combination CI
(oracle weights, sample-split weights, plug-in weights without splitting).

Plan: quality_reports/plans/2026-10-09_band-toy-simulation.md
Outputs: main/outputs/band_toy/
"""

import itertools
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from scipy.stats import norm

from helpers.weight_optimisers import optimal_weight_grid

# ----------------------------- configuration ------------------------------
THETA = 2.0                 # true estimand
DELTA = 1.0                 # half-width of the true band
TRUE_LO, TRUE_HI = THETA - DELTA, THETA + DELTA
ALPHA = 0.05
Z = norm.ppf(1 - ALPHA / 2) # per side, Bonferroni for the joint band
N = (2000, 500)             # dataset sizes
FRAC_A = 0.3                # fold-A fraction for the sample-split method
R = 5000                    # MC replicates per cell
SEED = 7
SE1 = 0.1                   # sd of Ybar_1
GAPS = [0.0, 0.05, 0.1, 0.2, 0.5, 1.0]  # g_L = a_2 and g_U = c_1
RATIOS = [1, 2, 5]          # r = SE_2 / SE_1
W_GRID = np.linspace(0.0, 1.0, 1001)

OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs", "band_toy")

METHODS = ["M1 intersection", "M2 convex (oracle)", "M3 convex (split)",
           "M4 convex (plug-in)", "M5 dataset 1 only", "M5 dataset 2 only"]
PLOT_METHODS = METHODS[:4]
COLOURS = dict(zip(PLOT_METHODS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]))
MARKERS = dict(zip(PLOT_METHODS, ["o", "s", "^", "D"]))


# ----------------------------- DGP ------------------------------
def dgp_params(g_L, g_U, r):
    """Slack a = (0, g_L), c = (g_U, 0) around the true band."""
    a = np.array([0.0, g_L])
    c = np.array([g_U, 0.0])
    mu = THETA + (c - a) / 2
    Delta = DELTA + (a + c) / 2
    se = np.array([SE1, r * SE1])
    sigma = se * np.sqrt(N)
    return dict(mu=mu, Delta=Delta, sigma=sigma, se=se,
                L_star=mu - Delta, U_star=mu + Delta)


# ----------------------------- estimation ------------------------------
def optimal_weight(m1, m2, v1, v2, lower):
    """Weight on dataset 1, optimised on W_GRID, vectorised over replicates.

    lower=True maximises  w m1 + (1-w) m2 - z sqrt(w^2 v1 + (1-w)^2 v2);
    lower=False minimises w m1 + (1-w) m2 + z sqrt(w^2 v1 + (1-w)^2 v2).
    """
    m1, m2, v1, v2 = (np.atleast_1d(x)[:, None] for x in (m1, m2, v1, v2))
    W = W_GRID[None, :]
    mean = W * m1 + (1 - W) * m2
    pen = Z * np.sqrt(W**2 * v1 + (1 - W)**2 * v2)
    idx = np.argmax(mean - pen, axis=1) if lower else np.argmin(mean + pen, axis=1)
    return W_GRID[idx]


def fused_ci(L, U, s2, wL, wU):
    """Convex-combination bounds and CI. L, U, s2 have shape (R, 2)."""
    Lf = wL * L[:, 0] + (1 - wL) * L[:, 1]
    Uf = wU * U[:, 0] + (1 - wU) * U[:, 1]
    seL = np.sqrt(wL**2 * s2[:, 0] + (1 - wL)**2 * s2[:, 1])
    seU = np.sqrt(wU**2 * s2[:, 0] + (1 - wU)**2 * s2[:, 1])
    return dict(lo=Lf - Z * seL, hi=Uf + Z * seU, Lhat=Lf, Uhat=Uf,
                wL=np.broadcast_to(wL, Lf.shape), wU=np.broadcast_to(wU, Uf.shape))


def intersection_ci(L, U, s2):
    rows = np.arange(L.shape[0])
    kL, kU = np.argmax(L, axis=1), np.argmin(U, axis=1)
    Lhat, Uhat = L[rows, kL], U[rows, kU]
    return dict(lo=Lhat - Z * np.sqrt(s2[rows, kL]), hi=Uhat + Z * np.sqrt(s2[rows, kU]),
                Lhat=Lhat, Uhat=Uhat,
                wL=(kL == 0).astype(float), wU=(kU == 0).astype(float))


def simulate_cell(g_L, g_U, r):
    p = dgp_params(g_L, g_U, r)
    # Common random numbers: every cell reuses the same standard normal draws,
    # so differences between cells are not MC noise.
    rng = np.random.default_rng(SEED)
    Y = [p["mu"][k] + p["sigma"][k] * rng.standard_normal((R, N[k])) for k in (0, 1)]

    # Full-sample statistics
    Ybar = np.column_stack([y.mean(axis=1) for y in Y])
    s2 = np.column_stack([y.var(axis=1, ddof=1) / N[k] for k, y in enumerate(Y)])
    L, U = Ybar - p["Delta"], Ybar + p["Delta"]

    # Fold statistics. Rows are i.i.d., so the first n_A columns are a random split.
    nA = [int(FRAC_A * n) for n in N]
    nB = [n - a for n, a in zip(N, nA)]
    YbarA = np.column_stack([Y[k][:, :nA[k]].mean(axis=1) for k in (0, 1)])
    YbarB = np.column_stack([Y[k][:, nA[k]:].mean(axis=1) for k in (0, 1)])
    # Fold-A variance rescaled to the fold-B size: the variance the CI will use
    vA = np.column_stack([Y[k][:, :nA[k]].var(axis=1, ddof=1) / nB[k] for k in (0, 1)])
    s2B = np.column_stack([Y[k][:, nA[k]:].var(axis=1, ddof=1) / nB[k] for k in (0, 1)])
    LA, UA = YbarA - p["Delta"], YbarA + p["Delta"]
    LB, UB = YbarB - p["Delta"], YbarB + p["Delta"]

    res = {}
    res["M1 intersection"] = intersection_ci(L, U, s2)

    se2 = p["se"]**2
    wL_or = optimal_weight_grid(p["L_star"][0], p["L_star"][1], se2[0], se2[1], lower=True)[0]
    wU_or = optimal_weight_grid(p["U_star"][0], p["U_star"][1], se2[0], se2[1], lower=False)[0]
    res["M2 convex (oracle)"] = fused_ci(L, U, s2, wL_or, wU_or)

    wL_sp = optimal_weight_grid(LA[:, 0], LA[:, 1], vA[:, 0], vA[:, 1], lower=True)[0]
    wU_sp = optimal_weight_grid(UA[:, 0], UA[:, 1], vA[:, 0], vA[:, 1], lower=False)[0]
    res["M3 convex (split)"] = fused_ci(LB, UB, s2B, wL_sp, wU_sp)

    wL_pi = optimal_weight_grid(L[:, 0], L[:, 1], s2[:, 0], s2[:, 1], lower=True)[0]
    wU_pi = optimal_weight_grid(U[:, 0], U[:, 1], s2[:, 0], s2[:, 1], lower=False)[0]
    res["M4 convex (plug-in)"] = fused_ci(L, U, s2, wL_pi, wU_pi)

    res["M5 dataset 1 only"] = fused_ci(L, U, s2, 1.0, 1.0)
    res["M5 dataset 2 only"] = fused_ci(L, U, s2, 0.0, 0.0)

    # M2's own target: the population convex combination of the bounds
    m2 = res["M2 convex (oracle)"]
    own_lo = wL_or * p["L_star"][0] + (1 - wL_or) * p["L_star"][1]
    own_hi = wU_or * p["U_star"][0] + (1 - wU_or) * p["U_star"][1]
    own = dict(cov_lo_own=mean_se((m2["lo"] <= own_lo).astype(float)),
               cov_hi_own=mean_se((m2["hi"] >= own_hi).astype(float)))

    rows = []
    for method, out in res.items():
        rows.append(dict(g_L=g_L, g_U=g_U, r=r, method=method,
                         wL_oracle=wL_or, wU_oracle=wU_or, **summarise(out)))
        if method == "M2 convex (oracle)":
            for name, (m, se) in own.items():
                rows[-1][name], rows[-1][f"{name}_mcse"] = m, se
    return rows


# ----------------------------- metrics ------------------------------
def mean_se(x):
    return x.mean(), x.std(ddof=1) / np.sqrt(len(x))


def summarise(out):
    lo, hi = out["lo"], out["hi"]
    width = hi - lo
    stats = dict(
        cov_joint=(lo <= TRUE_LO) & (hi >= TRUE_HI),
        cov_lo=lo <= TRUE_LO,
        cov_hi=hi >= TRUE_HI,
        cov_theta=(lo <= THETA) & (hi >= THETA),
        width=width,
        bias_L=out["Lhat"] - TRUE_LO,
        bias_U=out["Uhat"] - TRUE_HI,
        wL=out["wL"],
        wU=out["wU"],
    )
    row = {}
    for name, x in stats.items():
        row[name], row[f"{name}_mcse"] = mean_se(np.asarray(x, dtype=float))
    row["width_median"] = np.median(width)
    row["wL_sd"], row["wU_sd"] = np.std(out["wL"], ddof=1), np.std(out["wU"], ddof=1)
    row["wL_interior"] = np.mean((out["wL"] > 0) & (out["wL"] < 1))
    row["wU_interior"] = np.mean((out["wU"] > 0) & (out["wU"] < 1))
    return row


# ----------------------------- sanity checks ------------------------------
def sanity_checks(df):
    def get(method, g_L, g_U, r):
        return df[(df.method == method) & (df.g_L == g_L) & (df.g_U == g_U) & (df.r == r)].iloc[0]

    results = []

    def check(name, ok, detail):
        results.append(ok)
        print(f"[{'PASS' if ok else 'FAIL'}] {name}: {detail}")

    print("\nSanity checks")
    m1 = get("M1 intersection", 0.0, 0.0, 1)
    target = norm.cdf(Z)**2
    check("M1 lower-end coverage at g_L=0, r=1 vs Phi(z)^2",
          abs(m1.cov_lo - target) <= 3 * m1.cov_lo_mcse,
          f"{m1.cov_lo:.4f} vs {target:.4f} (MC SE {m1.cov_lo_mcse:.4f})")
    target = SE1 / np.sqrt(np.pi)
    check("M1 lower bias at g_L=0, r=1 vs SE/sqrt(pi)",
          abs(m1.bias_L - target) <= 3 * m1.bias_L_mcse,
          f"{m1.bias_L:.4f} vs {target:.4f} (MC SE {m1.bias_L_mcse:.4f})")

    m2 = get("M2 convex (oracle)", 0.0, 0.0, 1)
    check("M2 lower-end coverage at g_L=0, r=1 vs Phi(z)",
          abs(m2.cov_lo - norm.cdf(Z)) <= 3 * m2.cov_lo_mcse,
          f"{m2.cov_lo:.4f} vs {norm.cdf(Z):.4f} (MC SE {m2.cov_lo_mcse:.4f})")

    for r in RATIOS:
        target = (r * SE1)**2 / (SE1**2 + (r * SE1)**2)
        m2 = get("M2 convex (oracle)", 0.0, 0.0, r)
        check(f"M2 inverse-variance weights at a tie, r={r}",
              abs(m2.wL_oracle - target) <= 1e-3 and abs(m2.wU_oracle - target) <= 1e-3,
              f"w_L={m2.wL_oracle:.3f}, w_U={m2.wU_oracle:.3f} vs {target:.3f}")

    m2_all = df[df.method == "M2 convex (oracle)"]
    # Coverage of the true band is >= Phi(z) when w is interior (the fused population
    # bound is looser than the truth), so check coverage of M2's own target instead.
    zs = np.concatenate([(m2_all.cov_lo_own - norm.cdf(Z)) / m2_all.cov_lo_own_mcse,
                         (m2_all.cov_hi_own - norm.cdf(Z)) / m2_all.cov_hi_own_mcse])
    check("M2 per-end coverage of its own target ~ Phi(z) in every cell",
          np.max(np.abs(zs)) <= 4,
          f"max |z| = {np.max(np.abs(zs)):.2f} over {len(zs)} ends; "
          f"{np.sum(np.abs(zs) > 3)} with |z| > 3 (cells share draws, so not independent)")

    for r in RATIOS:
        m1, m2 = get("M1 intersection", 1.0, 1.0, r), get("M2 convex (oracle)", 1.0, 1.0, r)
        ok = (m2.wL_oracle >= 0.95 and m2.wU_oracle <= 0.05
              and abs(m1.cov_joint - m2.cov_joint) <= 0.01
              and abs(m1.width - m2.width) / m2.width <= 0.02)
        check(f"M1 ~ M2 at g_L=g_U=1, r={r}", ok,
              f"w_L={m2.wL_oracle:.3f}, w_U={m2.wU_oracle:.3f}; coverage {m1.cov_joint:.4f} vs "
              f"{m2.cov_joint:.4f}; width {m1.width:.4f} vs {m2.width:.4f}")

    worst = 0.0
    for method in PLOT_METHODS:
        sub = df[df.method == method]
        worst = max(worst,
                    sub.groupby(["r", "g_L"]).cov_lo.agg(np.ptp).max(),
                    sub.groupby(["r", "g_U"]).cov_hi.agg(np.ptp).max())
    check("Lower-end coverage constant in g_U (and upper in g_L)", worst <= 2 / R,
          f"max range {worst:.5f} (allowed {2 / R:.5f})")

    print(f"{sum(results)}/{len(results)} checks passed")
    return all(results)


# ----------------------------- plots ------------------------------
DIVERGING = LinearSegmentedColormap.from_list("div", ["#d03b3b", "#f0efec", "#2a78d6"])


def style(ax):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.grid(axis="y", color="#e5e4e0", lw=0.8)
    ax.set_axisbelow(True)


def plot_per_end(df, column, ref, ylabel, fname, methods):
    """Rows: lower end vs g_L, upper end vs g_U; columns: r. Gaps at equal spacing."""
    fig, axes = plt.subplots(2, len(RATIOS), figsize=(13, 7), sharey="row")
    pos = np.arange(len(GAPS))
    for j, r in enumerate(RATIOS):
        for i, (end, gap, other) in enumerate([("L", "g_L", "g_U"), ("U", "g_U", "g_L")]):
            ax = axes[i, j]
            # Per-end quantities do not depend on the other gap; take it at 0
            sub = df[(df.r == r) & (df[other] == 0.0)].sort_values(gap)
            col = column.format(end=end, end_name="lo" if end == "L" else "hi")
            if ref is not None:
                band = 2 * np.sqrt(ref * (1 - ref) / R)
                ax.axhspan(ref - band, ref + band, color="#e5e4e0", lw=0)
                ax.axhline(ref, color="#73726c", lw=1, ls="--")
            for m in methods:
                s = sub[sub.method == m]
                ax.plot(pos, s[col], color=COLOURS[m], marker=MARKERS[m], ms=6, lw=2, label=m)
            ax.set_xticks(pos, [f"{g:g}" for g in GAPS])
            ax.set_xlabel(f"${gap[0]}_{gap[2]}$ (not to scale)")
            ax.set_title(f"{'Lower' if end == 'L' else 'Upper'} end, r = {r}", fontsize=10)
            style(ax)
        axes[0, j].set_ylabel(ylabel.format(end="lower")) if j == 0 else None
        axes[1, j].set_ylabel(ylabel.format(end="upper")) if j == 0 else None
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=150)
    plt.close(fig)


def heatmap_grid(mats, row_labels, title, fname, centre, fmt):
    vals = np.concatenate([m.ravel() for row in mats for m in row])
    span = max(abs(vals.max() - centre), abs(vals.min() - centre), 1e-3)
    norm_ = TwoSlopeNorm(vcenter=centre, vmin=centre - span, vmax=centre + span)
    fig, axes = plt.subplots(len(mats), len(RATIOS), figsize=(13, 4.2 * len(mats)), squeeze=False)
    for i, row in enumerate(mats):
        for j, (r, M) in enumerate(zip(RATIOS, row)):
            ax = axes[i, j]
            im = ax.imshow(M, origin="lower", cmap=DIVERGING, norm=norm_)
            for (a, b), v in np.ndenumerate(M):
                ax.text(b, a, fmt.format(v), ha="center", va="center", fontsize=7, color="#1a1a19")
            ax.set_xticks(range(len(GAPS)), [f"{g:g}" for g in GAPS])
            ax.set_yticks(range(len(GAPS)), [f"{g:g}" for g in GAPS])
            ax.set_xlabel("$g_U$")
            if j == 0:
                ax.set_ylabel("$g_L$")
            ax.set_title(f"{row_labels[i]}, r = {r}", fontsize=10)
    fig.colorbar(im, ax=axes, shrink=0.8, label=title)
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=150, bbox_inches="tight")
    plt.close(fig)


def pivot(df, method, r, column):
    sub = df[(df.method == method) & (df.r == r)]
    return sub.pivot(index="g_L", columns="g_U", values=column).loc[GAPS, GAPS].to_numpy()


def make_plots(df):
    plot_per_end(df, "cov_{end_name}", norm.cdf(Z), "Coverage of {end} end",
                 "fig1_per_end_coverage.png", PLOT_METHODS)

    cov = [[pivot(df, m, r, "cov_joint") for r in RATIOS]
           for m in ("M1 intersection", "M3 convex (split)")]
    heatmap_grid(cov, ["M1 intersection", "M3 convex (split)"],
                 "Joint coverage of the true band (centre 0.95)",
                 "fig2_joint_coverage.png", 1 - ALPHA, "{:.3f}")

    ratios = [[pivot(df, m, r, "width") / pivot(df, "M1 intersection", r, "width") for r in RATIOS]
              for m in ("M3 convex (split)", "M2 convex (oracle)")]
    heatmap_grid(ratios, ["M3 / M1", "M2 / M1"], "Mean width ratio (centre 1)",
                 "fig3_width_ratio.png", 1.0, "{:.2f}")

    plot_per_end(df, "w{end}", None, "Mean weight on dataset 1 ({end} end)",
                 "fig4_weights.png", PLOT_METHODS[1:])


# ----------------------------- main ------------------------------
def run():
    os.makedirs(OUT_DIR, exist_ok=True)
    rows = []
    cells = list(itertools.product(RATIOS, GAPS, GAPS))
    for i, (r, g_L, g_U) in enumerate(cells, 1):
        rows.extend(simulate_cell(g_L, g_U, r))
        if i % 18 == 0:
            print(f"cell {i}/{len(cells)}")
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT_DIR, "results.csv"), index=False)

    sanity_checks(df)

    show = df[(df.g_L == df.g_U) & df.method.isin(PLOT_METHODS)]
    with pd.option_context("display.width", 140, "display.float_format", "{:.3f}".format):
        print("\nDiagonal g_L = g_U (full results in results.csv)")
        print(show[["r", "g_L", "method", "cov_joint", "cov_lo", "cov_hi", "width",
                    "wL", "wU"]].to_string(index=False))

    make_plots(df)
    print(f"\nWritten to {OUT_DIR}")
    return df


if __name__ == "__main__":
    run()
