"""Paper tables for the causal simulation, read from summary_full.csv.

Does not rerun the simulation: run `Rscript R/causal_sim.R` first.

Rows: for each Lambda_1 in --lambdas-1, the --top values of Lambda_2 at which
the convex split CI is narrowest relative to the calibrated intersection CI,
among the pairs where both have valid joint coverage. Columns: the methods in
--methods (default: each dataset alone, the calibrated intersection and the
convex split).

Usage, from the repo root:
    python -m src.latex.causal_sim_tables [--summary <summary_full.csv>]
        [--outcome continuous] [--lambdas-1 3 3.5 4 4.5 5] [--top 3]
        [--metric ci|band] [--methods dataset1 dataset2 intersection split
        crossfit clr convexclr] (default: the first four)
        [--macro-decimals 3] [--macro-prefix causal] [--alpha 0.05]
Outputs, beside the summary:
    width_coverage_best_<outcome>.tex  rows: Lambda_1, Lambda_2, then
        width (coverage) per method
    band_width_best_<outcome>.tex      rows: Lambda_1, Lambda_2, true band
        width, then mean estimated band width per method
    <table>_macros.tex                 median, min and max MC s.e. in the
        table, as LaTeX macros for the caption, named by outcome:
        \\<prefix><Outcome>WidthMcse{Median,Min,Max} and
        \\<prefix><Outcome>CovMcse{...} (ci) or
        \\<prefix><Outcome>BandMcse{...} (band), e.g.
        \\causalBinaryWidthMcseMedian; a different --macro-prefix (e.g.
        causalForest) lets two runs' tables share a document
"""

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd

from src.latex import toy_sim_slide_figs

REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_SUMMARY = (REPO_ROOT / "src" / "outputs" / "causal_sim_noise"
                   / "summary_full.csv")
DATASET_1 = "dataset 1"
DATASET_2 = "dataset 2"
INTERSECTION = "intersection (calibrated)"
CONVEX_SPLIT = "convex (split)"
# CLI name -> method name in summary_full.csv. The last three exist only in
# runs of R/causal_sim_quick.R.
METHOD_NAMES = {"dataset1": DATASET_1, "dataset2": DATASET_2,
                "intersection": INTERSECTION, "split": CONVEX_SPLIT,
                "crossfit": "convex (cross-fit)",
                "clr": "intersection (CLR)",
                "convexclr": "convex (CLR)"}
# Default columns: the methods of R/causal_sim.R.
DEFAULT_METHODS = ["dataset1", "dataset2", "intersection", "split"]
# Start of every caption-macro name, unless --macro-prefix is given.
DEFAULT_MACRO_PREFIX = "causal"
# Caption macros: name after \<prefix><Outcome> -> column of MC s.e.
CAPTION_MACROS = {
    "ci": {"WidthMcse": "width_mcse", "CovMcse": "cov_joint_mcse"},
    "band": {"BandMcse": "band_width_mcse"},
}
UNDERCOVER_SES = toy_sim_slide_figs.UNDERCOVER_SES
LAMBDA_TOL = 1e-9  # Lambda values are read back from CSV


def valid_coverage(df: pd.DataFrame, nominal: float) -> pd.Series:
    """True where joint coverage is not clearly below nominal."""
    return df.cov_joint >= nominal - UNDERCOVER_SES * df.cov_joint_mcse


def best_pairs(
    df: pd.DataFrame, lambdas_1: Sequence[float], top: int, nominal: float
) -> pd.DataFrame:
    """The (lambda_1, lambda_2) pairs with the smallest split/intersection ratio.

    Args:
        df: summary_full.csv rows for one outcome.
        lambdas_1: Values of Lambda_1 to tabulate; each must be on the grid.
        top: Number of Lambda_2 values kept per Lambda_1.
        nominal: Nominal joint coverage, 1 - alpha.

    Returns:
        Columns lambda_1, lambda_2, ratio (mean width of the convex split CI
        over that of the calibrated intersection CI), grouped by lambda_1 in
        the order given and sorted by lambda_2 within each group.

    Raises:
        ValueError: If a Lambda_1 is not on the grid or has no pair where
            both CIs have valid coverage.
    """
    fused = df[df.method.isin([INTERSECTION, CONVEX_SPLIT])]
    valid = fused.assign(valid=valid_coverage(fused, nominal))
    both_valid = valid.groupby(["lambda_1", "lambda_2"]).valid.all()
    width = fused.pivot_table(index=["lambda_1", "lambda_2"],
                              columns="method", values="width")
    ratio = (width[CONVEX_SPLIT] / width[INTERSECTION])[both_valid]
    ratio = ratio.rename("ratio").reset_index()
    groups = []
    for lambda_1 in lambdas_1:
        rows = ratio[np.isclose(ratio.lambda_1, lambda_1, atol=LAMBDA_TOL)]
        if rows.empty:
            raise ValueError(
                f"No valid pairs at {lambda_1=}; grid has "
                f"{sorted(set(df.lambda_1))}"
            )
        groups.append(rows.nsmallest(top, "ratio").sort_values("lambda_2"))
    return pd.concat(groups, ignore_index=True)


def cells_at(
    df: pd.DataFrame, lambda_1: float, lambda_2: float,
    methods: Sequence[str]
) -> pd.DataFrame:
    """The summary rows of the given methods at one (Lambda_1, Lambda_2)."""
    keep = (df.method.isin(methods)
            & np.isclose(df.lambda_1, lambda_1, atol=LAMBDA_TOL)
            & np.isclose(df.lambda_2, lambda_2, atol=LAMBDA_TOL))
    return df[keep]


def pair_rows(
    df: pd.DataFrame, pairs: pd.DataFrame, methods: Sequence[str]
) -> pd.DataFrame:
    """The summary rows of every tabulated pair and method."""
    return pd.concat([cells_at(df, p.lambda_1, p.lambda_2, methods)
                      for p in pairs.itertuples()])


def table_rows(
    df: pd.DataFrame, pairs: pd.DataFrame, nominal: float, metric: str,
    methods: Sequence[str]
) -> str:
    """LaTeX tabular rows: Lambda_1, Lambda_2, then one cell per method.

    Rows are grouped by Lambda_1, printed on the first row of each group.
    With metric "ci" a cell is width (coverage), with coverage in
    \\undercov{} when clearly below nominal; the document must define
    \\undercov. With metric "band" a column with the true band width comes
    before the methods, and a cell is the mean estimated band width.
    """
    lines = []
    for lambda_1, group in pairs.groupby("lambda_1", sort=False):
        for k, pair in enumerate(group.itertuples()):
            cells = cells_at(df, lambda_1, pair.lambda_2, methods)
            row = [f"{lambda_1:g}" if k == 0 else "", f"{pair.lambda_2:g}"]
            if metric == "band":
                row.append(f"{cells.true_band_width.iloc[0]:.3f}")
            for method in methods:
                result = cells[cells.method == method].iloc[0]
                row.append(toy_sim_slide_figs.band_cell(result)
                           if metric == "band"
                           else toy_sim_slide_figs.table_cell(result, nominal))
            lines.append(" & ".join(row) + r" \\")
        lines.append(r"\addlinespace")
    return "\n".join(lines[:-1]) + "\n"


def caption_macros(
    cells: pd.DataFrame,
    metric: str,
    outcome: str,
    decimals: int,
    macro_prefix: str = DEFAULT_MACRO_PREFIX,
) -> str:
    """LaTeX \\def lines: median, min and max MC s.e. of each tabulated quantity.

    Macros are named \\<macro_prefix><Outcome><quantity><stat>, so the tables
    of both outcomes (and, with different prefixes, of two runs) can share a
    document.
    """
    lines = []
    for quantity, column in CAPTION_MACROS[metric].items():
        prefix = f"{macro_prefix}{outcome.capitalize()}{quantity}"
        stats = {"Median": cells[column].median(), "Min": cells[column].min(),
                 "Max": cells[column].max()}
        lines.extend(rf"\def\{prefix}{name}{{{value:.{decimals}f}}}"
                     for name, value in stats.items())
    return "\n".join(lines) + "\n"


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY,
                        help="summary_full.csv written by R/causal_sim.R.")
    parser.add_argument("--outcome", choices=("continuous", "binary"),
                        default="continuous", help="Outcome to tabulate.")
    parser.add_argument("--lambdas-1", type=float, nargs="+",
                        default=[3.0, 3.5, 4.0, 4.5, 5.0],
                        help="Values of Lambda_1, one row group each.")
    parser.add_argument("--top", type=int, default=3,
                        help="Lambda_2 values per Lambda_1: those with the "
                             "smallest split/intersection width ratio.")
    parser.add_argument("--metric", choices=("ci", "band"), default="ci",
                        help="Table cells: 'ci' for width (coverage), 'band' "
                             "for the mean estimated band width.")
    parser.add_argument("--methods", nargs="+", choices=tuple(METHOD_NAMES),
                        default=DEFAULT_METHODS,
                        help="Method columns, in this order; crossfit, clr "
                             "and convexclr need a R/causal_sim_quick.R run.")
    parser.add_argument("--macro-decimals", type=int, default=3,
                        help="Decimal places of the MC s.e. caption macros.")
    parser.add_argument("--macro-prefix", default=DEFAULT_MACRO_PREFIX,
                        help="Start of the caption-macro names (letters "
                             "only), e.g. causalForest, so that two runs' "
                             "macros can share a document [default "
                             "%(default)s].")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Nominal level of the run (see run_info.txt).")
    parser.add_argument("--out-dir", type=Path, default=None,
                        help="Output directory; defaults to the summary's.")
    args = parser.parse_args(argv)
    if args.top < 1:
        parser.error(f"--top must be >= 1, got {args.top}")
    if len(set(args.methods)) < len(args.methods):
        parser.error(f"--methods has a repeat: {args.methods}")
    if args.macro_decimals < 0:
        parser.error(f"--macro-decimals must be >= 0, got "
                     f"{args.macro_decimals}")
    if not (args.macro_prefix.isascii() and args.macro_prefix.isalpha()):
        parser.error(f"--macro-prefix must be ASCII letters only, got "
                     f"{args.macro_prefix!r}")
    return args


def main() -> None:
    args = parse_args()
    nominal = 1 - args.alpha
    methods = [METHOD_NAMES[name] for name in args.methods]
    summary = pd.read_csv(args.summary)
    df = summary[summary.outcome == args.outcome]
    pairs = best_pairs(df, args.lambdas_1, args.top, nominal)
    out_dir = args.summary.parent if args.out_dir is None else args.out_dir
    prefix = "band_width" if args.metric == "band" else "width_coverage"
    table_path = out_dir / f"{prefix}_best_{args.outcome}.tex"
    header = (f"% Written by python -m src.latex.causal_sim_tables from "
              f"{args.summary.name}\n")
    table_path.write_text(
        header + table_rows(df, pairs, nominal, args.metric, methods)
    )
    macros_path = table_path.with_name(f"{table_path.stem}_macros.tex")
    macros_path.write_text(
        header + caption_macros(pair_rows(df, pairs, methods), args.metric,
                                args.outcome, args.macro_decimals,
                                args.macro_prefix)
    )
    print(pairs.round(3).to_string(index=False))
    print(f"Written to {table_path}\nWritten to {macros_path}")


if __name__ == "__main__":
    main()
