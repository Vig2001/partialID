# Plan: band toy simulation (intersection vs convex CIs for a known true band)

Date approved: 2026-10-09
Script: `main/band_toy.py`
Outputs: `main/outputs/band_toy/`

## Goal

Compare the naive intersection CI with the convex-combination CI when the true
(sharp) identified band is known by construction. Two datasets each identify a
band that is wider than the truth; together they identify exactly the true band.

## 1. Fixed quantities

| Symbol | Value | Meaning |
|---|---|---|
| θ | 2 | True estimand |
| δ | 1 | Half-width; true band is [θ − δ, θ + δ] = [1, 3] |
| α, z | 0.05, z = Φ⁻¹(1 − α/2) = 1.96 per side | Bonferroni, joint coverage ≥ 0.95 |
| n_1, n_2 | 2000, 500 | Dataset sizes |
| f | 0.3 | Fold-A fraction: n^A = (600, 150), n^B = (1400, 350) |
| R | 5000 | MC replicates per cell |
| Seed | 7 | Same data used by all methods within a replicate |

## 2. Scenario grid

Slack around the true band: dataset k identifies
L*_k = θ − δ − a_k, U*_k = θ + δ + c_k, with a_k, c_k ≥ 0 and
min(a_1, a_2) = min(c_1, c_2) = 0, so the intersection is exactly [θ − δ, θ + δ].
DGP parameters: μ_k = θ + (c_k − a_k)/2, Δ_k = δ + (a_k + c_k)/2.

Set a_1 = 0, c_2 = 0 (dataset 1 attains the true lower end, dataset 2 the true
upper end) and a_2 = g_L, c_1 = g_U independently:

- g_L, g_U ∈ {0, 0.05, 0.1, 0.2, 0.5, 1.0}
- μ_1 = θ + g_U/2, Δ_1 = δ + g_U/2 → population band [1, 3 + g_U]
- μ_2 = θ − g_L/2, Δ_2 = δ + g_L/2 → population band [1 − g_L, 3]

Precision: SE_1 = σ_1/√n_1 = 0.1 (σ_1 ≈ 4.47); SE_2 = r · SE_1 with
r ∈ {1, 2, 5} (σ_2 = r · 0.1 · √500 ≈ 2.24 r). At r = 2 the datasets have equal
σ and differ only through n.

6 × 6 × 3 = 108 cells.

Per-end coverage depends only on that end's gap (weights are chosen per end);
the 2-D grid matters for joint coverage and total width.

## 3. Data generation (per replicate)

Y_1i ~ N(μ_1, σ_1²), i ≤ 2000; Y_2i ~ N(μ_2, σ_2²), i ≤ 500, independent.
The analyst sees only Y_1, Y_2, Δ_1, Δ_2.

## 4. Methods

- **M1 Intersection (full sample).** Ȳ_k, ŝ_k² = S_k²/n_k, L̂_k = Ȳ_k − Δ_k,
  Û_k = Ȳ_k + Δ_k. k_L = argmax L̂_k, k_U = argmin Û_k.
  CI = [L̂_{k_L} − z ŝ_{k_L}, Û_{k_U} + z ŝ_{k_U}].
- **M2 Convex, oracle weights (full sample).**
  w_L* = argmax_{w∈[0,1]} { w L*_1 + (1−w) L*_2 − z √(w² SE_1² + (1−w)² SE_2²) },
  w_U* analogously (minimise the upper endpoint), using true values. Fuse L̂_k, Û_k
  with these weights; SE = √(w² ŝ_1² + (1−w)² ŝ_2²).
- **M3 Convex, sample split.** Split each dataset into fold A (30%) and B (70%).
  Fold A: L̂_k^A, Û_k^A, v̂_k = S_k^{A2}/n_k^B (rescaled to fold-B size) → ŵ_L, ŵ_U.
  Fold B: fuse with ŵ frozen, ŝ_k^{B2} = S_k^{B2}/n_k^B, build the CI.
- **M4 Convex, plug-in, no split (exploratory).** M3 with both folds = full
  sample, v̂_k = S_k²/n_k. No validity guarantee.
- **M5 Single-dataset CIs (reference).** [L̂_k − z ŝ_k, Û_k + z ŝ_k], k = 1, 2.

Optimisation: concave objective maximised on the grid w ∈ {0, 0.001, …, 1},
vectorised across replicates.

## 5. Metrics (per method, per cell, each with MC SE)

- Joint band coverage 1{CI_lo ≤ 1 and CI_hi ≥ 3}; lower-end and upper-end coverage
- Coverage of θ = 2
- Mean and median width
- Bias of point bounds: mean(L̂) − 1, mean(Û) − 3
- Mean and sd of ŵ_L, ŵ_U; share of weights strictly inside (0, 1)

## 6. Sanity checks (must pass before trusting results)

- g_L = 0, r = 1: M1 lower-end coverage ≈ Φ(z)² ≈ 0.951; M1 lower bias ≈ 0.1/√π ≈ 0.056;
  M2 w_L* = 0.5 and lower-end coverage ≈ 0.975.
- g_L = 0, any r: w_L* = SE_2²/(SE_1² + SE_2²); same for upper end at g_U = 0.
- Every cell: M2 per-end coverage of its own target w L*_1 + (1−w) L*_2 (resp. upper) ≈ 0.975.
  (Corrected 2026-10-09: coverage of the true band is ≥ 0.975 when w is interior.)
- g_L = g_U = 1.0: M1 ≈ M2, w_L* ≈ 1, w_U* ≈ 0.
- Lower-end coverage does not change with g_U (bookkeeping check).

## 7. Outputs

- CSV: one row per (g_L, g_U, r, method).
- Fig 1: per-end coverage vs gap (rows: lower vs g_L, upper vs g_U; columns: r;
  lines: methods; reference line 0.975 with ±2 MC SE band).
- Fig 2: heatmaps of joint coverage over (g_L, g_U) for M1 and M3, per r.
- Fig 3: heatmaps of width ratio M3/M1 and M2/M1, per r.
- Fig 4: weights vs gap (M2, M3, M4).
- Nothing goes into `partialID_Overleaf/` for now.

## 8. Expectations (to check, not assume)

- Small gaps: M1 undercovers at that end, most where the tight dataset is the
  noisy one (upper end, r > 1).
- Large gaps: all methods agree.
- M3 valid but wider than M2. M4 unknown.
