# baqaee_farhi_model

A self-contained, importable Python package implementing the linearized
multi-country, multi-sector nested-CES network model of Baqaee and Farhi
("Networks, Barriers, and Trade"), as used in Section 2 / Appendix A.5 of
Bachmann, Baqaee, Bayer, Kuhn, Löschel, Moll, Peichl, Pittel and Schularick
(2024), "What if? The macroeconomic and distributional effects for Germany
of a stop of energy imports from Russia," *Economica* 91(364).

This is a translation of the original authors' MATLAB code
(`../../replication/baqaee_farhi_model/`), validated element-wise against
its actual output under GNU Octave 8.4.0 -- see "Validation" below -- and
checked formula-by-formula against the theory paper's own Computational
Appendix (Appendix D) and general Allen-Uzawa elasticity results (Appendix
E). It is not a re-implementation from the paper's description; every
closed-form step traces back to either the original code's own output or
the theory paper's own derivation.

## Installing / importing this as a standalone package

There's no `pip install` step -- just put the parent directory on
`sys.path` and import the folder as a package:

```python
import sys
sys.path.insert(0, "/path/to/python replication")   # the directory *containing* baqaee_farhi_model/
import baqaee_farhi_model as bf

result = bf.run_scenario(
    keep_c=[2, 10, 34, 41],             # AUT, DEU, RUS, USA (ROW appended automatically)
    countries=['AUT', 'DEU', 'RUS', 'USA', 'ROW'],
    shocks=[{'sellers': [3], 'buyers': [1, 2], 'sectors': None, 'intensity': 150}],
    ngrid=5,
)
print(result['dlogW'])   # {'AUT': -0.0024, 'DEU': -0.0041, 'RUS': -0.0075, ...} -- LOG GNE change, not %
```

`data_dir` defaults to this package's own directory (where the WIOD `.mat`
files live), so this works from any working directory without extra setup.

## Using a non-WIOD dataset: the HAIO contract

`run()`/`run_scenario()` don't actually require WIOD -- `keep_c` and
`data_dir` are just the WIOD-specific path (via `io_reorder()`) to build a
standardized six-key "BF HAIO" (Homothetic Aggregate Input-Output) dict,
named after the theory paper's own term for this object (Appendix D/E of
Baqaee & Farhi's "Networks, Barriers, and Trade"). Everything downstream of
that dict -- `main_load_data()`, `nested_ces.py`, the discretization loop in
`run()` -- is already source-agnostic: it infers the number of sectors and
factor categories from the arrays' own shapes rather than hardcoding WIOD's
30 sectors / 4 factor categories.

You can build and pass that dict yourself, bypassing `io_reorder()`/WIOD
entirely, via the `haio=` keyword (pass anything, e.g. `None`, for `keep_c`
when doing this -- it's ignored):

```python
haio = dict(
    C=...,             # int, number of countries (incl. any aggregate "rest of world")
    Omega=...,         # (C*N, C*N): producer i's expenditure share on good j
    beta=...,          # (C*N, C): household c's expenditure share on good i
    alpha=...,         # (C*N,): value-added share of producer i
    alpha_VA=...,      # (C*N, F_data): producer i's share spent on factor f, any F_data
    trade_elast=...,   # (N,): cross-country trade elasticity by sector
    GDP_weights=...,   # (C,): each country's share of world GDP/GNE, sums to 1
)
result = bf.run_scenario(None, countries, shocks, ngrid=5, haio=haio)
```

This is the plug-in point for a GTAP or OECD-ICIO-based build: write a
loader that turns that source's raw data into this six-key dict (its own
sector/country classification, its own factor breakdown or none at all --
`alpha_VA` can be a single all-ones column if the source has no factor
split), and `run_scenario()` runs it completely unchanged. See
`main_load_data.py`'s module docstring for the exact per-key shapes and
semantics -- that docstring is the authoritative contract definition.
`io_reorder.py` is simply the WIOD implementation of it.

### `icio_to_haio.py`: an OECD ICIO loader, verified against OECD's own docs

`icio_to_haio.py` implements this contract for OECD ICIO's plain-CSV layout
-- same normalization arithmetic as `io_reorder.py` (transpose to a buyer x
seller matrix, row-normalize for intermediate-input shares, column-normalize
for final-consumption shares, value-added share = VA / (VA + intermediate
cost)), just against ICIO's row/column bookkeeping instead of WIOD's packed
binary blocks.

**The expected input format is now checked against OECD's own published
documentation** (their "ReadMe_icio_csv" file and 2025-edition country/sector
notes, both retrieved from `stats.oecd.org`), not just a guess -- this
corrected two assumptions the first version of this loader made:
value added is **one shared row** for the whole table (not one row per
country), and final demand is **6 named categories per country**
(`HFCE, NPISH, GGFC, GFCF, INVNT, P33`), which the loader now sums itself
rather than asking the caller to pre-sum. It also surfaced a correction to
something said earlier in this project: crude oil/gas *extraction*
(`B06`) is its own ICIO sector, separate from `D` (electricity/gas/steam
*supply*) -- `B06` is the right target for a Hormuz/crude-specific shock.

**Still not run against a real downloaded OECD ICIO release.** The actual
data file (`2017-2022_EXT.zip`, ~136MB, linked from OECD's ICIO page) sits
behind a Cloudflare bot challenge on `www.oecd.org` that no plain HTTP
client can pass -- unlike the documentation PDFs, which OECD happens to also
serve unprotected from `stats.oecd.org`. Getting a real file into this
loader needs a manual browser download; once you have it:

```python
import pandas as pd
from icio_to_haio import icio_to_haio

df = pd.read_csv('2020.csv', index_col=0)   # one file per year in the OECD zip
haio = icio_to_haio(df, countries=['USA', 'SAU', 'IRQ', ...], sectors=['B06', ...],
                     trade_elast=[...])
result = run_scenario(None, countries, shocks, haio=haio)
```

Running `python icio_to_haio.py` prints a small **synthetic** (not real
data) 2-country x 2-sector example showing exactly what "the standard data
should look like after cleaning" means in practice, using OECD's real sector
codes for realism (`B06`, `D`) even though the numbers are made up. The toy
input -- country-sector labels on both axes, ONE value-added row, 6
final-demand columns per country:

```
         AUS_B06  AUS_D  DEU_B06  DEU_D  AUS_HFCE ... AUS_P33  DEU_HFCE ... DEU_P33
AUS_B06        5      3        2      1        10 ...       0         4 ...       0
AUS_D          2      6        1      2         0 ...       0         3 ...       0
DEU_B06        1      2        4      3         0 ...       0        12 ...       0
DEU_D          1      1        2      5         0 ...       0         0 ...       9
VA             9      7        8      6       NaN ...     NaN       NaN ...     NaN
```

(`VA`'s row is only ever read at the 4 country-sector columns, so the `NaN`s
under the final-demand columns are irrelevant and never touched.) This
produces the standardized HAIO dict:

```
Omega (4x4, row-normalized intermediate-input shares):
[[0.556 0.222 0.111 0.111]
 [0.25  0.5   0.167 0.083]
 [0.222 0.111 0.444 0.222]
 [0.091 0.182 0.273 0.455]]

beta (4x2, column-normalized final-consumption shares):
[[0.435 0.143]
 [0.348 0.107]
 [0.13  0.429]
 [0.087 0.321]]

alpha (value-added share of gross output, per producer):
[0.5   0.368 0.471 0.353]

alpha_VA (no factor-split source -> single all-ones column):
[[1.] [1.] [1.] [1.]]

GDP_weights (each country's share of total final demand):
[0.451 0.549]
```

which then runs cleanly through `run_scenario(None, ['AUS', 'DEU'], shocks,
haio=haio)` and returns real (if economically meaningless, since the input
numbers are made up) `dlogW` output -- confirming the contract and the
loader logic both work against OECD's actual documented structure,
independent of the still-open question of getting the real file past
Cloudflare.

If you're working from inside this directory instead (the original usage
pattern -- `cd baqaee_farhi_model && python run_paper_scenario.py`), nothing
changes: every module still works as a plain script exactly as before. The
package wrapper is purely additive.

## Quick reference

| Function | What it does |
|---|---|
| `run_scenario(keep_c, countries, shocks, ngrid=20, sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05, data_dir=None, haio=None)` | The one you want. Runs the model and returns a labeled, self-describing result. |
| `run(keep_c, shocks, ngrid=20, ..., haio=None)` | Lower-level: same computation, returns a raw `(C,)` numpy array instead of a labeled dict. |
| `main_load_data(haio, initial_tariff_index, factor_index)` | Source-agnostic: turns a standardized HAIO dict into the model's standard-form inputs. Normally called for you by `run()`. |
| `io_reorder(keep_c, data_dir)` | The WIOD-specific loader: WIOD 2008 `.mat` files -> a HAIO dict. |
| `icio_to_haio(df, countries, sectors, trade_elast, ...)` | A first OECD-ICIO loader: a plain-CSV-shaped DataFrame -> a HAIO dict. Not yet run against a real ICIO download -- see below. |
| `value_added_shares`, `response`, `solve_dlambda_F_all` | The model internals (Allen elasticities, one discretization step's equilibrium response, the linear-system solve). You won't normally call these directly. |

### `run_scenario()`'s return value

```python
{
    'dlogW': {country_code: log_GNE_change, ...},   # NOT percent -- multiply by 100 yourself
    'World': log_GNE_change,                         # GNE-weighted world aggregate, same units
    'elasticities': {'sigma': ..., 'theta': ..., 'gamma': ..., 'epsilon': ...},
    'ngrid': ...,
    'shocks': [...],   # echoed back, so the result is self-describing
}
```

### The `shocks` parameter

A list of iceberg-cost shock legs (multiple legs are summed, so independent
shocks can be combined in one run):

```python
{
    'sellers':   [3],           # 1-indexed positions *within keep_c* (not WIOD country codes)
    'buyers':    [1, 2],        # same indexing; None = all countries
    'sectors':   None,          # 0-indexed WIOD sector indices (0..29); None = all 30 sectors
    'intensity': 150,           # see "Interpreting intensity" below
}
```

The paper's own scenario (EU vs. Russia, all sectors) is exactly:
`shocks=[{'sellers': [RUS], 'buyers': EU, 'sectors': None, 'intensity': 150}]`
(see `run_paper_scenario.py`).

### Interpreting `intensity`

`intensity` sets the cumulative log iceberg-cost wedge applied between the
shocked sellers and buyers: `Δlogτ_total = ln(1 + intensity/100)`, split
evenly across `ngrid` discretization steps. This is a **technology/
productivity shock** (real resources lost in transit, the classic iceberg
formulation), not a tariff or markup -- the paper explicitly distinguishes
"productivity shocks, which nest iceberg shocks" from "wedge shocks, which
nest tariff changes," and only the iceberg/technology channel is wired up
here (the tariff channel, `dlogt`, exists in the model's data structures but
is always zero in this implementation, matching the paper's own main driver
script).

`intensity = 150` means the effective delivered price is 2.5x higher
(`ln(2.5)`) due to the friction. Whether that's enough to send trade to
*near zero* depends on the shocked sector's own trade elasticity
(`trade_elast_2008.mat`, median ~5 across the 30 WIOD sectors): at a 2.5x
price wedge and an elasticity of ~6 (median + 1, per `AES_func.m`'s
convention), quantity traded falls to roughly `2.5^-6 ≈ 0.26%` of its
original level. 150 was tuned by the original authors for an *aggregate,
all-sectors* shock -- a single, more price-inelastic sector may need a
larger intensity to fully collapse trade, and a highly substitutable sector
may need much less.

## Why this is slow, and why that's inherent

The full 41-country run (`run_paper_scenario.py`) solves a ~1271-dimensional
linear system once per discretization step, and needs roughly
`(C+CF+1) ≈ 1272` evaluations of the equilibrium-response function per step
-- expect on the order of **tens of minutes** on a single core at
`ngrid=20`. This isn't a Python inefficiency: the original MATLAB code faces
the exact same system, just solved via an explicit Jacobian construction
that needs ~30-60GB of RAM (confirmed to OOM under GNU Octave 8.4.0 with
15GB available in this environment). See `nested_ces.py`'s module docstring
for the full explanation of the closed-form-AES + probing strategy that
makes this tractable in ~50MB instead, at the cost of wall-clock time. For
iterating on scenario design, use a small `keep_c` subset and a low `ngrid`
(3-5) -- see the example above -- and only scale up to the full run for
final numbers.

## Data included

Three WIOD 2013-release, benchmark-year-2008 files (see the top-level
`python replication/README.md` for the full data provenance and the other
three `.mat` files present but unused by this code path):

| File | Contents |
|---|---|
| `wiott2008.mat` | World Input-Output Table: the 41-country x 31-industry bilateral flow matrix |
| `wiodsea2008.mat` | Socio-Economic Accounts: labor/capital compensation by country-industry |
| `trade_elast_2008.mat` | Industry-level cross-country (Armington-style) trade elasticities |

## Validation

Every layer checked element-wise against the actual MATLAB/Octave output at
a scale small enough for Octave to run without exhausting memory (`< 1e-14`
to `< 3e-17` depending on layer -- see the top-level `python
replication/README.md`'s "Validation" section for the full breakdown by
file), plus a full formula-by-formula check against the underlying theory
paper's Computational Appendix (Appendix D) and general Allen-Uzawa
elasticity result (Appendix E) -- every closed-form AES expression in
`nested_ces.py` was independently re-derived from the theory paper's own
notation and matched term-for-term, not just checked numerically. The full
41-country run itself (this being the one case Octave can't execute here to
cross-check bit-for-bit) reproduces the paper's own stated Table 2 column
(2) figure (Germany GNE loss 0.2-0.3%) almost exactly: see
`paper_scenario_result.json`.

## See also

- `../README.md` -- the Python replication package this lives inside.
- `../../README.md` -- the top-level overview of all four packages
  (this one, the rest of the Python replication, and the two closed-form
  second-order-approximation packages in Python and R).
