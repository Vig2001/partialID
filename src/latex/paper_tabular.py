"""Self-contained LaTeX tabulars for the paper's comparator-set tables.

Since 2026-11-15 each table file holds the whole tabular (column spec,
comparator-set header, column headers, rules and body), so the paper inputs
it with one \\input and a change of columns needs no edit outside Tables/.
"""

from collections.abc import Sequence

SAMPLE_SPLITTING = "Sample splitting"
SIMULTANEOUS = r"\shortstack{Simultaneous\\ inference}"
SELECTION = r"\shortstack{Adaptive inequality\\ selection}"
# Comparator set and column header of each method, by its name in the causal
# summaries (src/latex/causal_sim_tables.py); toy methods map to these names.
COLUMNS = {
    # Names since 2026-11-16 (supervisor): Intersection/Convex (procedure);
    # Intersection (CLR) is the Sidak-corrected intersection.
    "intersection (split)": (SAMPLE_SPLITTING,
                             r"\shortstack{Intersection\\ (split)}"),
    "convex (split)": (SAMPLE_SPLITTING, r"\shortstack{Convex\\ (split)}"),
    "convex (cross-fit)": (SAMPLE_SPLITTING,
                           r"\shortstack{Convex\\ (cross-fit)}"),
    "intersection (CLR)": (SIMULTANEOUS, r"\shortstack{Intersection\\ (CLR)}"),
    "convex (CLR)": (SIMULTANEOUS, r"\shortstack{Convex\\ (CLR)}"),
    "intersection (CLR, selection)": (SELECTION,
                                      r"\shortstack{Intersection\\ (CLR)}"),
    "convex (CLR, selection)": (SELECTION, r"\shortstack{Convex\\ (CLR)}"),
}


def tabular(lead: Sequence[str], methods: Sequence[str], body: str) -> str:
    """The whole tabular: comparator-set header, column headers and body.

    Args:
        lead: Headers of the leading columns, e.g. ["$r$", "$g$"].
        methods: Method columns, keys of COLUMNS; each comparator set must
            be contiguous.
        body: Table rows, each ending in \\\\.

    Raises:
        ValueError: If a method has no header or a set is not contiguous.
    """
    unknown = [m for m in methods if m not in COLUMNS]
    if unknown:
        raise ValueError(f"No column header for {unknown}")
    sets = [COLUMNS[m][0] for m in methods]
    runs: list[tuple[str, int]] = []  # (set, number of columns), in order
    for name in sets:
        if runs and runs[-1][0] == name:
            runs[-1] = (name, runs[-1][1] + 1)
        elif any(name == done for done, _ in runs):
            raise ValueError(f"Comparator set {name!r} is not contiguous")
        else:
            runs.append((name, 1))
    groups, rules, start = [], [], len(lead) + 1
    for name, width in runs:
        groups.append(rf"\multicolumn{{{width}}}{{c}}{{{name}}}")
        rules.append(rf"\cmidrule(lr){{{start}-{start + width - 1}}}")
        start += width
    headers = [*lead, *(COLUMNS[m][1] for m in methods)]
    lines = [
        rf"\begin{{tabular}}{{{' '.join('c' * len(headers))}}}",
        r"\toprule",
        " & ".join([""] * len(lead) + groups) + r" \\",
        "".join(rules),
        " & ".join(headers) + r" \\",
        r"\midrule",
        body.rstrip("\n"),
        r"\bottomrule",
        r"\end{tabular}",
    ]
    return "\n".join(lines) + "\n"
