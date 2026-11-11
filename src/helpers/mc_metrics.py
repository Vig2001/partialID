"""Monte Carlo summaries of CI performance against a known true band."""

from collections.abc import Sequence

import numpy as np
import pandas as pd

from src.helpers import ci_methods


def mean_se(x: np.ndarray) -> tuple[float, float]:
    """Returns the mean of x and its Monte Carlo standard error.

    Args:
        x: One value per replicate, shape (R,), R >= 2.

    Raises:
        ValueError: If x has fewer than two entries.
    """
    x = np.asarray(x, dtype=float)
    if x.size < 2:
        raise ValueError(f"x needs at least two entries, got {x.size=}")
    return x.mean(), x.std(ddof=1) / np.sqrt(x.size)


def summarise_ci(
    ci: ci_methods.CIResult, true_lo: float, true_hi: float, theta: float
) -> dict[str, float]:
    """Coverage, CI width, estimated band width, bias and weight summaries.

    One method in one cell. The estimated band width is U_hat - L_hat, the
    width of the fused bound estimates without the sampling margin.

    Each per-replicate quantity q gives columns q (mean) and q_mcse (Monte
    Carlo SE).

    Args:
        ci: CI endpoints, bound estimates and weights over R replicates.
        true_lo: True lower end of the band.
        true_hi: True upper end of the band.
        theta: True estimand, inside [true_lo, true_hi].

    Returns:
        Flat dict of summaries, in a fixed column order.
    """
    width = ci.hi - ci.lo
    per_replicate = dict(
        cov_joint=(ci.lo <= true_lo) & (ci.hi >= true_hi),
        cov_lo=ci.lo <= true_lo,
        cov_hi=ci.hi >= true_hi,
        cov_theta=(ci.lo <= theta) & (ci.hi >= theta),
        width=width,
        band_width=ci.U_hat - ci.L_hat,
        bias_L=ci.L_hat - true_lo,
        bias_U=ci.U_hat - true_hi,
        wL=ci.w_L,
        wU=ci.w_U,
    )
    row = {}
    for name, values in per_replicate.items():
        row[name], row[f"{name}_mcse"] = mean_se(values)
    row["width_median"] = np.median(width)
    row["wL_sd"] = np.std(ci.w_L, ddof=1)
    row["wU_sd"] = np.std(ci.w_U, ddof=1)
    row["wL_interior"] = np.mean((ci.w_L > 0) & (ci.w_L < 1))
    row["wU_interior"] = np.mean((ci.w_U > 0) & (ci.w_U < 1))
    return row


def width_relative_to_best(
    df: pd.DataFrame,
    methods: Sequence[str],
    cell_columns: Sequence[str],
    nominal: float,
    undercover_ses: float,
) -> pd.DataFrame:
    """Each method's mean CI width over the narrowest valid one in its cell.

    A method is valid in a cell when its joint coverage is not clearly below
    nominal, cov_joint >= nominal - undercover_ses * cov_joint_mcse. The
    reference in a cell is the smallest mean width among the methods valid
    there, so every valid method has ratio >= 1 and the narrowest has 1. A
    method that undercovers in a cell gets no ratio there (NaN), since width
    comparisons only make sense between valid CIs; so does every method in a
    cell where none is valid. Each method is a rule fixed in advance, so its
    ratios over cells show how close it stays to the per-cell best.

    Args:
        df: One row per (cell, method) with the cell columns, method, width,
            cov_joint and cov_joint_mcse (e.g. summary_full.csv).
        methods: Methods to compare; rows of other methods are ignored.
        cell_columns: Columns identifying a cell, e.g. lambda_1, lambda_2.
        nominal: Nominal joint coverage, 1 - alpha.
        undercover_ses: Monte Carlo standard errors below nominal at which
            coverage counts as too low.

    Returns:
        One row per (cell, method), with the cell columns, method, width,
        valid and ratio, sorted by cell, then by method in the order given.

    Raises:
        ValueError: If a cell lacks one of the methods or has it twice.
    """
    cell_columns = list(cell_columns)
    out = _valid_widths(df, methods, cell_columns, nominal, undercover_ses)
    best = (out.width.where(out.valid)
            .groupby([out[c] for c in cell_columns]).transform("min"))
    out["ratio"] = (out.width / best).where(out.valid)
    return _sort_by_cell(out, methods, cell_columns)


def width_relative_to_truth(
    df: pd.DataFrame,
    methods: Sequence[str],
    cell_columns: Sequence[str],
    nominal: float,
    undercover_ses: float,
) -> pd.DataFrame:
    """Each method's mean CI width over the true band width of its cell.

    Validity is as in width_relative_to_best: a method that undercovers in a
    cell gets no ratio there (NaN). The reference is the cell's true band
    width (true_band_width), the same for every method, so a ratio of 1
    would be a CI as narrow as the target band itself.

    Args:
        df: One row per (cell, method) with the cell columns, method, width,
            true_band_width, cov_joint and cov_joint_mcse (e.g.
            summary_full.csv).
        methods: Methods to compare; rows of other methods are ignored.
        cell_columns: Columns identifying a cell, e.g. lambda_1, lambda_2.
        nominal: Nominal joint coverage, 1 - alpha.
        undercover_ses: Monte Carlo standard errors below nominal at which
            coverage counts as too low.

    Returns:
        One row per (cell, method), with the cell columns, method, width,
        valid and ratio, sorted by cell, then by method in the order given.

    Raises:
        ValueError: If a cell lacks one of the methods or has it twice, or
            a cell's true band width is not positive (the ratio is undefined
            where the band is a point).
    """
    cell_columns = list(cell_columns)
    rows = df[df.method.isin(methods)]
    if (rows.true_band_width <= 0).any():
        bad = rows.loc[rows.true_band_width <= 0, cell_columns].iloc[0]
        raise ValueError(
            "true band width must be positive in every cell, got "
            f"{rows.true_band_width.min()} at {tuple(bad)}"
        )
    out = _valid_widths(df, methods, cell_columns, nominal, undercover_ses)
    out["ratio"] = (out.width / rows.true_band_width).where(out.valid)
    return _sort_by_cell(out, methods, cell_columns)


def _valid_widths(
    df: pd.DataFrame,
    methods: Sequence[str],
    cell_columns: list[str],
    nominal: float,
    undercover_ses: float,
) -> pd.DataFrame:
    """Cell columns, method, width and valid for the methods' rows.

    Raises:
        ValueError: If a cell lacks one of the methods or has it twice.
    """
    rows = df[df.method.isin(methods)]
    counts = rows.groupby(cell_columns).method.agg(["nunique", "size"])
    incomplete = counts[(counts["nunique"] != len(methods))
                        | (counts["size"] != len(methods))]
    if not incomplete.empty:
        raise ValueError(
            f"{len(incomplete)} cells lack a method or repeat one, e.g. "
            f"{incomplete.index[0]}; methods: {list(methods)}"
        )
    out = rows[[*cell_columns, "method", "width"]].copy()
    out["valid"] = (rows.cov_joint
                    >= nominal - undercover_ses * rows.cov_joint_mcse)
    return out


def _sort_by_cell(
    out: pd.DataFrame, methods: Sequence[str], cell_columns: list[str]
) -> pd.DataFrame:
    """Rows sorted by cell, then by method in the order given."""
    order = out.method.map({m: k for k, m in enumerate(methods)})
    return (out.assign(method_order=order)
            .sort_values([*cell_columns, "method_order"])
            .drop(columns="method_order")
            .reset_index(drop=True))
