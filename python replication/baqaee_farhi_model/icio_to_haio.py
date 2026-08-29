"""
OECD ICIO -> standardized "BF HAIO" dict loader (see main_load_data.py's
module docstring for the exact six-key contract this produces).

This is the sibling to io_reorder.py's WIOD-specific loading logic -- same
arithmetic (transpose to a buyer x seller matrix, row-normalize for
intermediate-input shares, column-normalize for final-consumption shares,
value-added share = VA / (VA + total intermediate cost)), just against
OECD ICIO's plain-CSV row/column layout instead of WIOD's packed binary
block structure, which is most of why this loader is much shorter than
io_reorder.py.

EXPECTED INPUT FORMAT: a single matrix (pandas DataFrame, or a path to a CSV
in this shape) whose row and column labels are '{country}_{sector}' for
every country-sector pair (the same codes on both axes, matching OECD's own
ICIO table layout: https://www.oecd.org/en/data/datasets/inter-country-input-output-tables.html),
PLUS:
  - one extra row per country, its total value added, labeled '{country}_VA'
  - one extra column per country, its total final demand, labeled
    '{country}_FD' -- this is a SIMPLIFICATION of the real OECD file's
    several separate final-demand columns per country (household
    consumption, NPISH, government consumption, gross fixed capital
    formation, inventory change, direct purchases abroad); sum those
    together into one '{country}_FD' column yourself before calling this
    function if you're working from an actual OECD release.

NOT YET RUN AGAINST A REAL DOWNLOADED OECD ICIO RELEASE in this repository
-- see the toy worked example in this file's __main__ block and
baqaee_farhi_model/README.md for exactly what's been verified (the parsing
and normalization logic runs correctly end-to-end, including through
run_scenario(), on a small synthetic table with the same shape) versus what
still needs checking against a real file (ICIO's actual current column
names, whether the "sum the FD columns yourself" simplification above holds
for whatever vintage you download, and the Russia-imputation caveat
discussed separately -- see the top-level README).

Two known ICIO-specific gaps, both documented in baqaee_farhi_model/README.md:
    - No labor/skill factor breakdown (unlike WIOD SEA's 4 categories) --
      `alpha_VA` defaults to a single all-ones column here (no factor
      heterogeneity), matching main_load_data.py's ability to infer any
      factor-category count from alpha_VA's own shape. Pass your own
      alpha_VA if you have a supplementary factor-split source.
    - No trade elasticities of its own -- must be supplied externally
      (e.g. Fontagne-Guimbard-Orefice, already pre-aggregated to ICIO's
      sector classification).
"""
import numpy as np
import pandas as pd


def icio_to_haio(df, countries, sectors, trade_elast, alpha_VA=None,
                  va_suffix='_VA', fd_suffix='_FD'):
    """
    Parameters
    ----------
    df : pandas.DataFrame, the raw ICIO matrix (see module docstring).
    countries : list of country codes to track individually. Any country
        NOT in this list is simply excluded, not automatically folded into
        a "rest of world" composite -- build that composite yourself
        (sum the excluded countries' rows/columns into one synthetic extra
        country's row/column in `df`) before calling this function if you
        want one, matching io_reorder.py's own ROW handling.
    sectors : list of sector codes, in the order they should appear (fixes
        N = len(sectors) and each producer's position within its country's
        block).
    trade_elast : (N,) array-like, sector order must match `sectors`.
    alpha_VA : optional (C*N, F_data) array -- omit for the no-factor-split
        fast path (a single all-ones column, F_data=1).
    va_suffix, fd_suffix : label suffixes identifying each country's VA row
        / FD column.

    Returns
    -------
    dict, the standardized HAIO contract -- pass straight to
    run()/run_scenario() via `haio=`.
    """
    C, N = len(countries), len(sectors)
    cs_labels = [f'{c}_{s}' for c in countries for s in sectors]  # (CN,) fixed row/col order
    va_labels = [f'{c}{va_suffix}' for c in countries]
    fd_labels = [f'{c}{fd_suffix}' for c in countries]

    # Raw ICIO convention (like WIOD/WIOT): row = seller/origin, column =
    # buyer/destination. Transpose so row = buyer, matching Omega's
    # documented (i,j) = "expenditure share of producer i on good j".
    IO = df.loc[cs_labels, cs_labels].to_numpy(dtype=float).T  # (CN, CN), row=buyer
    FD = df.loc[cs_labels, fd_labels].to_numpy(dtype=float)    # (CN, C), row=good, col=household
    VA = df.loc[va_labels, cs_labels].to_numpy(dtype=float).sum(axis=0)  # (CN,)

    row_sums = IO.sum(axis=1, keepdims=True)
    Omega = np.divide(IO, row_sums, out=np.zeros_like(IO), where=row_sums != 0)

    col_sums = FD.sum(axis=0, keepdims=True)
    beta = np.divide(FD, col_sums, out=np.zeros_like(FD), where=col_sums != 0)

    total_intermediate_cost = IO.sum(axis=1)  # buyer i's total intermediate purchases
    alpha = np.divide(VA, VA + total_intermediate_cost,
                       out=np.zeros_like(VA), where=(VA + total_intermediate_cost) != 0)

    GDP_weights = FD.sum(axis=0)
    GDP_weights = GDP_weights / GDP_weights.sum()

    if alpha_VA is None:
        alpha_VA = np.ones((C * N, 1))  # no factor-split source (see module docstring)

    return dict(
        C=C, Omega=Omega, beta=beta, alpha=alpha, alpha_VA=alpha_VA,
        trade_elast=np.asarray(trade_elast, dtype=float), GDP_weights=GDP_weights,
    )


if __name__ == '__main__':
    # --- Toy worked example: 2 countries x 2 sectors, made-up numbers -----
    # NOT real ICIO data -- purely to show (a) the exact input shape
    # icio_to_haio() expects and (b) that the output is a valid HAIO dict
    # that runs cleanly through run_scenario(). See README.md for this
    # printed out alongside the resulting HAIO arrays.
    countries = ['A', 'B']
    sectors = ['S1', 'S2']
    cs = [f'{c}_{s}' for c in countries for s in sectors]  # A_S1, A_S2, B_S1, B_S2

    # Rows/cols: A_S1, A_S2, B_S1, B_S2, A_FD, B_FD (VA rows appended below)
    data = pd.DataFrame(
        [
            [5, 3, 2, 1, 10, 4],   # A_S1 sells to: A_S1 A_S2 B_S1 B_S2 A_FD B_FD
            [2, 6, 1, 2, 8, 3],    # A_S2
            [1, 2, 4, 3, 3, 12],   # B_S1
            [1, 1, 2, 5, 2, 9],    # B_S2
            [9, 7, 0, 0, 0, 0],    # A_VA (value added, only has entries under cs columns)
            [0, 0, 8, 6, 0, 0],    # B_VA
        ],
        index=cs + ['A_VA', 'B_VA'],
        columns=cs + ['A_FD', 'B_FD'],
    )

    print('=== Toy raw ICIO-shaped input ===')
    print(data)
    print()

    trade_elast = [5.0, 8.0]  # made-up, sector order matches `sectors`
    haio = icio_to_haio(data, countries, sectors, trade_elast)

    print('=== Resulting standardized HAIO dict ===')
    for k, v in haio.items():
        print(f'{k}:')
        print(np.round(v, 3) if isinstance(v, np.ndarray) else v)
        print()

    print('=== Feeding it straight into run_scenario() ===')
    from run_model import run_scenario
    result = run_scenario(None, countries, [{'sellers': [1], 'buyers': [2], 'sectors': None, 'intensity': 50}],
                           ngrid=3, haio=haio)
    print(result)
