"""
OECD ICIO -> standardized "BF HAIO" dict loader (see main_load_data.py's
module docstring for the exact six-key contract this produces).

This is the sibling to io_reorder.py's WIOD-specific loading logic -- same
arithmetic (transpose to a buyer x seller matrix, row-normalize for
intermediate-input shares, column-normalize for final-consumption shares,
value-added share = VA / (VA + total intermediate cost)), just against
OECD ICIO's plain-CSV row/column layout instead of WIOD's packed binary
block structure.

EXPECTED INPUT FORMAT: verified against OECD's own ICIO documentation
(the "ReadMe_icio_csv" file and 2025-edition country/sector notes, both
retrieved from stats.oecd.org -- see baqaee_farhi_model/README.md for what
this confirmed/corrected). A single matrix (pandas DataFrame, or a path to
a CSV in this shape):
  - Row and column labels for the intermediate-transactions block are
    '{country}_{sector}' for every country-sector pair (e.g. 'AUS_B06'),
    the same codes on both axes.
  - ONE extra row, value added at basic prices, labeled 'VA' by default --
    NOT one row per country. Its columns are the same country-sector
    labels as the intermediate block.
  - Final demand is 6 named categories per country in the real OECD file:
    HFCE, NPISH, GGFC, GFCF, INVNT, P33 (household consumption, non-profit
    institutions, government consumption, gross fixed capital formation,
    inventory change, direct purchases abroad). This loader sums those 6
    into one total per country itself (`fd_categories` below) -- you do
    not need to pre-sum them.

NOT YET RUN AGAINST A REAL DOWNLOADED OECD ICIO RELEASE in this repository.
The actual current data file (2017-2022_EXT.zip, ~136MB, referenced from
OECD's ICIO page) is hosted on www.oecd.org behind a Cloudflare bot
challenge that a plain HTTP client cannot pass -- unlike the documentation
PDFs/README above, which OECD happens to also serve unprotected from
stats.oecd.org. Getting a real file into this loader therefore needs either
a manual browser download (see baqaee_farhi_model/README.md for exactly
what to do with it once downloaded) or a different access path; the
parsing/normalization logic here is now matched to OECD's own documented
structure (verified, not guessed) but has only been exercised end-to-end on
the synthetic toy table in this file's __main__ block.

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

FD_CATEGORIES = ('HFCE', 'NPISH', 'GGFC', 'GFCF', 'INVNT', 'P33')


def icio_to_haio(df, countries, sectors, trade_elast, alpha_VA=None,
                  va_label='VA', fd_categories=FD_CATEGORIES):
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
        block) -- e.g. OECD's own codes like 'B06' (crude petroleum and
        natural gas extraction) or 'D' (electricity, gas, steam and air
        conditioning supply, i.e. utility-side, NOT the same sector as B06).
    trade_elast : (N,) array-like, sector order must match `sectors`.
    alpha_VA : optional (C*N, F_data) array -- omit for the no-factor-split
        fast path (a single all-ones column, F_data=1).
    va_label : row label for the single value-added row (default 'VA').
    fd_categories : column-label suffixes to sum per country into that
        country's total final demand (default: OECD's own 6 categories).

    Returns
    -------
    dict, the standardized HAIO contract -- pass straight to
    run()/run_scenario() via `haio=`.
    """
    C, N = len(countries), len(sectors)
    cs_labels = [f'{c}_{s}' for c in countries for s in sectors]  # (CN,) fixed row/col order
    fd_labels = {c: [f'{c}_{cat}' for cat in fd_categories] for c in countries}

    # Raw ICIO convention (like WIOD/WIOT): row = seller/origin, column =
    # buyer/destination. Transpose so row = buyer, matching Omega's
    # documented (i,j) = "expenditure share of producer i on good j".
    IO = df.loc[cs_labels, cs_labels].to_numpy(dtype=float).T  # (CN, CN), row=buyer

    # Sum each country's 6 final-demand categories into one column.
    FD = np.column_stack([df.loc[cs_labels, fd_labels[c]].to_numpy(dtype=float).sum(axis=1)
                           for c in countries])  # (CN, C)

    VA = df.loc[va_label, cs_labels].to_numpy(dtype=float)  # (CN,) -- ONE row, not one per country

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
    # icio_to_haio() expects, matched to OECD's own verified structure
    # (one shared VA row, 6 final-demand categories per country), and (b)
    # that the output is a valid HAIO dict that runs cleanly through
    # run_scenario(). See README.md for this printed out alongside the
    # resulting HAIO arrays. Sector codes use OECD's real naming
    # ('B06' = crude petroleum & natural gas extraction) purely for realism.
    countries = ['AUS', 'DEU']
    sectors = ['B06', 'D']  # crude oil/gas extraction; electricity/gas/steam supply
    cs = [f'{c}_{s}' for c in countries for s in sectors]  # AUS_B06, AUS_D, DEU_B06, DEU_D

    fd_cols = [f'{c}_{cat}' for c in countries for cat in FD_CATEGORIES]
    columns = cs + fd_cols

    # Intermediate flows (4x4) -- rows sum arbitrarily, not from real data.
    intermediate = [
        [5, 3, 2, 1],
        [2, 6, 1, 2],
        [1, 2, 4, 3],
        [1, 1, 2, 5],
    ]
    # Final demand: split a made-up total per country evenly-ish across its
    # 6 categories, per row (good).
    fd_totals = {  # (good) -> per-country total FD, before splitting into 6 categories
        'AUS': [10, 8], 'DEU': [3, 12],  # AUS_B06/AUS_D rows -> DEU_B06/DEU_D rows below
    }
    fd_block = [
        [10, 0, 0, 0, 0, 0,  4, 0, 0, 0, 0, 0],   # AUS_B06 row: all 10 in HFCE, all 4 in HFCE
        [0, 8, 0, 0, 0, 0,   0, 3, 0, 0, 0, 0],    # AUS_D row: spread differently, illustrative only
        [0, 0, 3, 0, 0, 0,   0, 0, 12, 0, 0, 0],   # B_S1-equivalent (DEU_B06)
        [0, 0, 0, 2, 0, 0,   0, 0, 0, 9, 0, 0],    # DEU_D row
    ]
    rows_cs = [r + fd for r, fd in zip(intermediate, fd_block)]
    va_row = [9, 7, 8, 6]  # ONE row, value added per country-sector column

    data = pd.DataFrame(
        rows_cs + [va_row],
        index=cs + ['VA'],
        columns=columns,
    )

    print('=== Toy raw ICIO-shaped input (OECD-verified structure) ===')
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
    result = run_scenario(None, countries, [{'sellers': [1], 'buyers': [2], 'sectors': [0], 'intensity': 50}],
                           ngrid=3, haio=haio)
    print(result)
