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
    HFCE, NPISH, GGFC, GFCF, INVNT, DPABR (household consumption, non-profit
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
what to do with it once downloaded) or a different access path; once
obtained, the parsing/normalization logic here has been run end-to-end
against a real 2022 file (81 countries x 50 sectors, the "SML" per-year
release) -- see baqaee_farhi_model/README.md for the real scenario result
and the data quirks that run surfaced.

ALWAYS EXCLUDE SECTOR 'T' when using real ICIO data with this model --
"activities of households as employers... for own use" is, by
national-accounts convention, 100% value-added with zero measured
intermediate cost for every country, which divides by zero in
nested_ces.py's AES formulas. It's a negligible share of GNE everywhere and
irrelevant to anything this project models, so excluding it costs nothing.

Two known ICIO-specific gaps, both documented in baqaee_farhi_model/README.md:
    - No labor/skill factor breakdown (unlike WIOD SEA's 4 categories) --
      `alpha_VA` defaults to a single all-ones column here (no factor
      heterogeneity), matching main_load_data.py's ability to infer any
      factor-category count from alpha_VA's own shape. Pass your own
      alpha_VA if you have a supplementary factor-split source.
    - No trade elasticities of its own -- must be supplied externally.
      `fontagne_icio_trade_elast.csv` in this directory is the real
      Fontagne-Guimbard-Orefice "New ICIO classification" (March 2026
      update, https://sites.google.com/view/product-level-trade-elasticity)
      download -- confirmed to use the exact same sector codes as the real
      ICIO file (`B06`, `C24A`, `C24B`, etc.), so the "New" vs. "Old"
      classification call was correct. See `load_fontagne_trade_elast()`
      below and its SIGN CONVENTION CAVEAT before treating its output as
      ground truth.

fontagne_icio_trade_elast.csv covers only 32 of the ~48 non-'T' ICIO
sectors -- trade elasticities are only estimable for tradable goods, so
most services (F, G, I, K, L, N, O, P, Q, S, ...) have no row at all, which
is expected, not a download error. Of those 32, four have `epsilon_icio ==
NA` even though the sector itself has a row: A02 (forestry), B05 (coal
mining), D (electricity/gas supply), R (recreation) -- Fontagne et al.'s own
estimation didn't converge for these, for whatever reason (thin trade
volume, insufficient tariff variation for identification). B06 (crude
petroleum and natural gas extraction -- what a Hormuz-type scenario would
actually shock) DOES have a valid estimate: -5.44.
"""
import numpy as np
import pandas as pd

FD_CATEGORIES = ('HFCE', 'NPISH', 'GGFC', 'GFCF', 'INVNT', 'DPABR')


def _fold_excluded_into_row(df, countries, row_label, va_label):
    """Aggregate every country present in `df` but NOT in `countries` (and
    not `row_label` itself) into `row_label`'s own rows/columns, by summing
    -- so excluded countries' flows are preserved (as a composite), not
    silently dropped. Mirrors io_reorder.py's WIOD "rest of world" handling.
    Without this, a producer whose real intermediate inputs come mostly
    from an excluded country would show near-zero intermediate cost in the
    reduced view, artificially inflating its value-added share (in the
    worst case to exactly 1.0, which divides by zero in nested_ces.py's
    AES formulas -- this is not a hypothetical, it's what happens on the
    real ICIO file with a small `countries` subset and no row_label)."""
    all_countries = sorted(set(lbl.split('_', 1)[0] for lbl in df.index
                                if lbl not in (va_label, 'TLS', 'OUT')))
    excluded = set(all_countries) - set(countries) - {row_label}
    if not excluded:
        return df

    def remap(label):
        if label in (va_label, 'TLS', 'OUT'):
            return label
        c, _, rest = label.partition('_')
        return f'{row_label}_{rest}' if c in excluded else label

    df = df.rename(index={lbl: remap(lbl) for lbl in df.index},
                    columns={lbl: remap(lbl) for lbl in df.columns})
    df = df.groupby(df.index).sum()
    df = df.T.groupby(df.T.index).sum().T
    return df


def load_fontagne_trade_elast(sectors, csv_path='fontagne_icio_trade_elast.csv'):
    """Load Fontagne-Guimbard-Orefice's ICIO-classification trade
    elasticities and return a (len(sectors),) array in nested_ces.py's
    `trade_elast` convention, ordered to match `sectors`.

    SIGN CONVENTION CAVEAT -- NOT independently verified: the CSV's
    `epsilon_icio` column is NEGATIVE (e.g. -5.44 for B06), consistent with
    the standard trade-literature convention of reporting the elasticity of
    trade flows with respect to (1 + trade cost) -- theoretically
    epsilon = -(sigma - 1) for a CES/Armington elasticity of substitution
    sigma. nested_ces.py's own `trade_elast` variable satisfies
    trade_elast + 1 == sigma (see nested_ces.py's docstring and
    io_reorder.py's AES_func.m-matched formulas), which combined with the
    above gives trade_elast = sigma - 1 = -epsilon_icio -- the transform
    used below. This is a plausible, evidence-consistent derivation (it
    lands in the same 1-15 range as WIOD's own trade_elast_2008.mat, which
    is positive), but it has NOT been empirically cross-checked against
    trade_elast_2008.mat's actual values for overlapping sectors the way
    every other number in this codebase has been -- do that cross-check
    (via Fontagne et al.'s own WIOD-classification download, also on their
    site) before trusting this for anything beyond illustrative use.

    Raises KeyError for any requested sector missing from the CSV entirely
    (~half of ICIO's non-'T' sectors -- see module docstring) OR present but
    with epsilon_icio == NA (A02, B05, D, R at last check) -- both cases
    need a caller decision (exclude the sector, or substitute a literature
    value from ../python-second-order/sigma_literature.py), not a silent
    NaN that would otherwise propagate unnoticed into nested_ces.py's
    arithmetic.
    """
    df = pd.read_csv(csv_path, index_col='icio2025')
    missing = [s for s in sectors if s not in df.index]
    if missing:
        raise KeyError(f'sectors not in {csv_path} (no trade elasticity estimated at all): {missing}')
    na = [s for s in sectors if pd.isna(df.loc[s, 'epsilon_icio'])]
    if na:
        raise KeyError(f'sectors in {csv_path} but epsilon_icio is NA (estimation did not converge): {na}')
    return -df.loc[sectors, 'epsilon_icio'].to_numpy(dtype=float)


def load_combined_trade_elast(sectors, csv_path='combined_trade_elast.csv', include_approximate=True):
    """Load trade_elast from combined_trade_elast.csv (built by
    build_combined_trade_elast.py) -- goods sectors from Fontagne-Guimbard-Orefice
    plus services sectors from Ahmad & Schreiber (2024)'s GTAP-level
    estimates, each already converted to nested_ces.py's convention, with
    provenance tracked in that file's `source`/`mapping_quality`/`note`
    columns (see build_combined_trade_elast.py's module docstring for the
    exact GTAP->ICIO mapping judgment calls this made). Covers 44 of
    ICIO's 49 non-'T' sectors -- still missing: A02, B05, B09, D, O (no
    estimate from either source).

    include_approximate : if False, raise for any sector whose
        mapping_quality is 'approximate' (the many-to-one GTAP mappings --
        K, J61, J62_63, N, R, S) instead of silently including it, for
        callers who want only the more directly-matched estimates.

    Raises KeyError for any requested sector missing from the CSV
    entirely -- same reasoning as load_fontagne_trade_elast(): a caller
    decision (exclude the sector, or substitute something else), not a
    silent default.
    """
    df = pd.read_csv(csv_path, index_col='icio_sector')
    missing = [s for s in sectors if s not in df.index]
    if missing:
        raise KeyError(f'sectors not in {csv_path} (no trade elasticity from either source): {missing}')
    if not include_approximate:
        approx = [s for s in sectors if df.loc[s, 'mapping_quality'] == 'approximate']
        if approx:
            raise KeyError(f'sectors only have an approximate (many-to-one GTAP) mapping: {approx}')
    return df.loc[sectors, 'trade_elast'].to_numpy(dtype=float)


def icio_to_haio(df, countries, sectors, trade_elast, alpha_VA=None,
                  va_label='VA', fd_categories=FD_CATEGORIES, row_label=None):
    """
    Parameters
    ----------
    df : pandas.DataFrame, the raw ICIO matrix (see module docstring).
    countries : list of country codes to track individually. Any country
        NOT in this list, and not equal to `row_label`, is folded into
        `row_label`'s own composite (summed in, not dropped) if `row_label`
        is given; otherwise it is simply excluded -- see
        `_fold_excluded_into_row()`'s docstring for why that silent
        exclusion is dangerous with a small `countries` subset. Include
        `row_label` itself in `countries` if you want it tracked as one of
        the model's countries (matching io_reorder.py's WIOD ROW, which is
        always tracked).
    row_label : optional country code (e.g. 'ROW') to fold every excluded
        country into. Omit for the old behavior (excluded countries
        silently dropped, not recommended except when `countries` already
        covers every country in `df`).
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
    if row_label is not None:
        df = _fold_excluded_into_row(df, countries, row_label, va_label)

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
