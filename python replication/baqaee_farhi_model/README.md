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

## The pipeline

Two loaders converge on the same standardized "HAIO" contract; everything
to the right of `main_load_data.py` doesn't know or care which loader fed
it (see "Using a non-WIOD dataset" below):

```
                       io_reorder.py   ---\
DATA (.mat / .csv)                          --->  main_load_data.py  --->  nested_ces.py  --->  run_model.py  --->  run_paper_scenario.py / run_paper_scenario_icio.py
                       icio_to_haio.py ---/
                       (or your own loader)              ^
                       build_combined_trade_elast.py -----+ (trade_elast for icio_to_haio.py)
```

| File | Role |
|---|---|
| `io_reorder.py` | WIOD `.mat` files -> HAIO dict |
| `icio_to_haio.py` | OECD ICIO CSV -> HAIO dict |
| `build_combined_trade_elast.py` | builds `combined_trade_elast.csv` (the trade elasticities `icio_to_haio.py` needs, ICIO has none of its own) |
| `main_load_data.py` | HAIO dict -> the model's standard-form inputs (source-agnostic) |
| `nested_ces.py` | the model itself: Allen elasticities, one discretization step's equilibrium response, the linear-system solve |
| `run_model.py` | the outer discretization loop (`run()`) and the labeled convenience wrapper (`run_scenario()`) |
| `run_paper_scenario.py` | a ready-to-run driver with the paper's exact settings on WIOD 2008 (all 41 WIOD countries, EU-vs-Russia) |
| `run_paper_scenario_icio.py` | the same EU-vs-Russia scenario, rebuilt on real OECD ICIO 2022 data + `combined_trade_elast.csv` instead of WIOD -- see "Running the paper scenario on ICIO instead of WIOD" below. Much slower: expect hours, not minutes (see that section). |
| `__init__.py` | makes this directory importable as a package from anywhere (see below) |

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

**The expected input format was checked against OECD's own published
documentation** (their "ReadMe_icio_csv" file and 2025-edition country/sector
notes, both retrieved from `stats.oecd.org`) and then confirmed against a
**real downloaded 2022 ICIO file** (the "SML" per-year release, 81 countries
x 50 sectors) -- both steps corrected assumptions the first version of this
loader made: value added is **one shared row** for the whole table (not one
row per country), final demand is **6 named categories per country**
(`HFCE, NPISH, GGFC, GFCF, INVNT, DPABR` -- the real file uses `DPABR`, not
`P33` as OECD's older documentation suggested), which the loader sums
itself, and crude oil/gas *extraction* (`B06`) is its own ICIO sector,
separate from `D` (electricity/gas/steam *supply*) -- `B06` is the right
target for a Hormuz/crude-specific shock.

**Real-data run, end to end, with real trade elasticities**: Saudi Arabia +
UAE (sellers) facing a 150% iceberg cost on `B06` exports to the US,
tracking `SAU`, `ARE`, `USA` and a `ROW` composite for every other country,
restricted to the 28 sectors Fontagne-Guimbard-Orefice actually estimated a
usable elasticity for:

```python
import pandas as pd
from icio_to_haio import icio_to_haio, load_fontagne_trade_elast
from run_model import run_scenario

df = pd.read_csv('2022_SML.csv', index_col=0)   # not committed -- see "Getting the real data" below
countries = ['SAU', 'ARE', 'USA', 'ROW']
fontagne = pd.read_csv('fontagne_icio_trade_elast.csv', index_col='icio2025')
sectors = [s for s in fontagne[fontagne['epsilon_icio'].notna()].index if s != 'T']
trade_elast = load_fontagne_trade_elast(sectors)   # see its SIGN CONVENTION CAVEAT below

haio = icio_to_haio(df, countries, sectors, trade_elast, row_label='ROW')
result = run_scenario(None, countries, [{'sellers': [1, 2], 'buyers': [3],
                                          'sectors': [sectors.index('B06')], 'intensity': 150}],
                      ngrid=5, haio=haio)
# {'SAU': -0.43%, 'ARE': +0.14%, 'USA': -0.015%, 'ROW': +0.001%, 'World': -0.005%} (log GNE)
```

Real, sensible general-equilibrium output: Saudi Arabia (the more
US-dependent of the two sellers) loses noticeably more than with a
placeholder elasticity; the US barely notices (two sellers, one sector, out
of its whole economy); UAE's small *gain* is a genuine substitution effect
worth digging into further, not a bug.

**`fontagne_icio_trade_elast.csv`** (committed, 2.4KB -- unlike the 90MB
ICIO file itself) is the real Fontagne-Guimbard-Orefice "New ICIO
classification" download; `load_fontagne_trade_elast(sectors)` reads it and
converts to `nested_ces.py`'s convention. **Two things worth knowing before
trusting its numbers**:
1. It only covers **32 of ICIO's 49 non-`T` sectors** -- trade elasticities
   are only estimable for tradable goods, so most services (`F`, `G`, `I`,
   `K`, `L`, `N`, `O`, `P`, `Q`, `S`, ...) have no row at all, and it's
   missing `B09` too. Of the 32 it does have, 4 are `NA` anyway (`A02`,
   `B05`, `D`, `R`) -- Fontagne et al.'s own estimation didn't converge for
   those. `B06` (crude petroleum and natural gas -- what a Hormuz scenario
   actually needs) does have a valid value: -5.44. (`combined_trade_elast.csv`,
   introduced further below, fills all of these gaps except `O` -- use that
   instead of this raw file unless you specifically want goods-only,
   tariff-based estimates.)
2. **The sign/offset conversion is a plausible derivation, not an
   independently verified fact.** The CSV reports a negative number
   (`epsilon_icio`); `nested_ces.py` needs a positive `trade_elast` where
   `trade_elast + 1` is a CES elasticity of substitution. `load_fontagne_trade_elast()`
   uses `trade_elast = -epsilon_icio`, based on the standard trade-literature
   convention that a reported elasticity like this equals `-(sigma-1)` --
   it lands in the same 1-15 range as WIOD's own (positive) `trade_elast_2008.mat`,
   which is a good sign, but this has **not** been cross-checked the way
   every other number in this codebase has been. The rigorous check: download
   Fontagne et al.'s WIOD-classification file from the same page, apply the
   same transform, and compare against the real `trade_elast_2008.mat`
   values already in this repo for overlapping sectors.

**Country coverage caveat found along the way**: ICIO does **not** track
Iran, Iraq, or Kuwait individually -- only Saudi Arabia and the UAE get
their own rows among Gulf states; everything else Gulf-related is inside
`ROW`. A serious Hormuz scenario naming Iran specifically (as in the
back-of-envelope work earlier in this project) can't isolate Iran's own
export collapse with ICIO alone.

**Real data quirks in ICIO** (all inherent to the data, not bugs):
1. **Still exclude sector `T`** ("activities of households as employers...
   for own use") when using ICIO with this model -- not because it crashes
   anymore (the `alpha==1.0` divide it triggers in `nested_ces.py`'s AES
   formulas is now guarded, same fix as item 2 below), but because it's a
   national-accounts imputation with no real trade behavior (always 100%
   value-added, zero measured intermediate cost, for every country) and a
   negligible share of GNE everywhere -- excluding it costs nothing and
   keeps the sector list to economically meaningful categories:
   `sectors = [s for s in ALL_ICIO_SECTORS if s != 'T']`.
2. **A country can have literally zero recorded activity in a niche
   sector** (Saudi Arabia has no coal mining, `B05`; Luxembourg has no
   mining or oil extraction at all) -- ICIO's 50-sector classification is
   finer than WIOD's 30, so this shows up more than it would have with
   WIOD, and gets more common, not less, as more countries are tracked
   individually (a 41-country run hits ~20 such cases). This used to
   produce a `LinAlgError: Singular matrix` (a very small `countries`
   subset) or silent NaN in every country's result (a larger one, via a
   `dlambda_F / lambda_F` divide-by-zero in `nested_ces.py` that NaN-poisons
   the shared Leontief inverse) -- both are now fixed by guarding every
   divide-by-a-zero-share in `nested_ces.py`/`run_model.py` to return 0
   instead. A zero-output producer genuinely has nothing to say about its
   own price in this linearization, so treating it as contributing zero is
   the correct fix, not a workaround -- no sector needs to be excluded for
   this reason anymore, at any `countries` size.
3. **A handful of real cells are negative** -- a small country's value
   added at basic prices in a niche sector (e.g. Cyprus, air transport),
   or inventory drawdown in final demand. `icio_to_haio.py` clamps these to
   0 before computing shares, exactly like `io_reorder.py` already does for
   WIOD -- otherwise a share (`alpha`, `beta`) could come out negative.

**`row_label`, added after finding a real numerical failure mode**: the
first version of this loader simply *dropped* any country not in
`countries`. Testing against the real file with only 4 countries kept
immediately produced silently-wrong `alpha` values (falsely inflated toward
1.0, since a producer's real intermediate purchases from an excluded
country vanished instead of being counted) and, in the worst case, an exact
`alpha=1.0` division-by-zero. Passing `row_label='ROW'` now folds every
excluded country into that composite (summing, not dropping) before
building the HAIO dict -- see `_fold_excluded_into_row()`'s docstring.
Always pass this unless your `countries` list already covers every country
in the file.

**Getting `trade_elast`**: ICIO has none of its own (see `icio_to_haio.py`'s
module docstring, "Two known ICIO-specific gaps"). Three layers exist here:

1. `fontagne_icio_trade_elast.csv` / `load_fontagne_trade_elast()` -- the
   raw Fontagne, Guimbard and Orefice download (product-level trade
   elasticities, https://sites.google.com/view/product-level-trade-elasticity,
   "New ICIO classification" -- confirmed to match `2022_SML.csv`'s sector
   codes exactly). Covers 28 of 49 non-`T` sectors -- **goods only**;
   tariff-based identification doesn't work for services, which have no
   tariffs to vary.
2. Ahmad and Schreiber (2024, USITC Working Paper), "Estimating
   Elasticities for Tradable Services in Policy Simulations," estimates
   services elasticities directly (a markup/monopolistic-competition
   identification, since there's no tariff variation to exploit) at the
   GTAP sector level -- adds 16 more sectors.
3. A small fallback layer sourced from the **original Bachmann-Baqaee
   et al. paper's own `trade_elast_2008.mat`** (already in this repo),
   found by matching the original MATLAB driver script's hardcoded
   30-sector label list against that file's values by index -- adds 3
   more real sectors (`A02`, `B05`, `B09`) plus one flagged placeholder
   (`D`, see below).

**`combined_trade_elast.csv` / `load_combined_trade_elast()`** is all
three combined -- **48 of ICIO's 49 non-`T` sectors**, every row carrying
a `source` column (`fontagne_guimbard_orefice_2022`, `ahmad_schreiber_2024`,
or `bachmann_baqaee_2024_wiod`) and a three-tier `mapping_quality` column,
weakest last:
   - `direct` -- a clean 1:1 sector match.
   - `approximate` -- a many-to-one mapping (GTAP or WIOD bundles several
     ICIO sectors into one broader category) -- e.g. `A02` (forestry) and
     `B05`/`B09` (mining) come from WIOD's own coarser "Agriculture,
     Hunting, Forestry and Fishing" and "Mining and Quarrying" categories.
   - `placeholder` -- not a real estimate at all: `D` (electricity/gas
     supply)'s value is **exactly 5.00, the paper's own admitted flat
     default** -- checking `trade_elast_2008.mat`'s actual values revealed
     that *every* WIOD sector from electricity/gas onward (all services)
     is exactly 5.00 with zero variation, meaning the original authors used
     one placeholder number for everything past manufacturing rather than
     sourcing real values individually. Using it here just adopts that
     same admitted guess, not a genuine resolution.

**Only `O` (public administration) is left out on purpose**, even though
the paper's own value is there too (also the same flat 5.00 placeholder):
government services are essentially non-traded by nature, so having *no*
trade elasticity for `O` is arguably correct rather than a gap -- filling
it with a number already known to be an unsourced placeholder, for a
sector where the concept barely applies, would be strictly worse than
leaving it out.

**The three sources report different objects and need different
sign/offset transforms** -- getting this wrong is the easiest way to
silently corrupt every number downstream:
   - Fontagne's `epsilon_icio` is a *negative* trade-cost elasticity
     (`epsilon = -(sigma-1)` is the standard convention), so
     `trade_elast = sigma - 1 = -epsilon_icio`.
   - Ahmad & Schreiber's EOS is *already* a directly-estimated elasticity of
     substitution (`sigma` itself, positive), so `trade_elast = EOS - 1`.
   - `trade_elast_2008.mat`'s values are already in this project's own
     `trade_elast` convention (it's the real WIOD data file `io_reorder.py`
     reads directly) -- used as-is, no transform.

**GTAP/WIOD -> ICIO mapping judgment calls** (all made explicitly in
`build_combined_trade_elast.py`, not silently): most GTAP services sectors
match an ICIO code 1:1 (`cns`->`F`, `trd`->`G`, `edu`->`P`, `hht`->`Q`,
`rsa`->`L`, etc. -- marked `direct`). A few don't split as cleanly, marked
`approximate`: GTAP's one broad `cmn` (ICT) sector is duplicated across
ICIO's `J61` and `J62_63` (but *not* `J58T60`, which Fontagne already
covers directly); ICIO's single `K` (financial+insurance) averages GTAP's
separate `ins` and `ofi`; `N` uses GTAP's `obs` (its NAICS 561
admin/support component is the closer match, while `M` is left to
Fontagne's own direct estimate since NAICS 541 professional services maps
better there); GTAP's `ros` is duplicated across ICIO's `R` and `S` (a
weaker match for `S` specifically); WIOD's "Mining and Quarrying" is used
for both `B05` and `B09`. See `build_combined_trade_elast.py`'s module
docstring for the full reasoning behind each call -- re-run that script if
any upstream source is revised.

**Result with broader coverage**: rerunning the Saudi Arabia + UAE
scenario with all 48 covered sectors (`B05` no longer needs excluding for
Saudi Arabia's zero coal-mining activity -- see "Real data quirks in ICIO"
above) changes the numbers again -- `SAU -0.201%`, `ARE +0.179%`,
`USA -0.007%` log GNE (vs. `-0.19%`/`+0.12%`/`-0.006%` with the earlier
44-sector coverage, and `-0.43%`/`+0.14%`/`-0.015%` with Fontagne's
coverage) -- a concrete demonstration that sector coverage isn't a detail,
it changes the answer.

**Getting the real data into your own copy**: the actual ICIO zip
(`2017-2022_EXT.zip` or similar, ~136MB, linked from OECD's ICIO page) sits
behind a Cloudflare bot challenge on `www.oecd.org` that no plain HTTP
client can pass -- unlike the documentation PDFs, which OECD happens to also
serve unprotected from `stats.oecd.org`. Download it manually in a browser,
unzip to get one CSV per year, and drop the year(s) you want into this
directory -- `.gitignore`d (via the year-prefixed pattern
`baqaee_farhi_model/[0-9][0-9][0-9][0-9]*.csv`, e.g. `2022_SML.csv`, which
deliberately does *not* match the small committed elasticity CSVs like
`combined_trade_elast.csv`) since a single year's file is far too large for
a normal GitHub push.

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
| `icio_to_haio(df, countries, sectors, trade_elast, row_label=None, ...)` | An OECD-ICIO loader, verified against a real downloaded file: a plain-CSV-shaped DataFrame -> a HAIO dict. Pass `row_label` (e.g. `'ROW'`) to fold excluded countries in rather than dropping them -- see below. |
| `load_fontagne_trade_elast(sectors, csv_path=...)` | Loads real Fontagne-Guimbard-Orefice trade elasticities for ICIO (goods, 28 sectors), converted to `nested_ces.py`'s convention -- see its sign-convention caveat below. |
| `load_combined_trade_elast(sectors, csv_path=..., include_approximate=True, include_placeholder=True)` | Loads Fontagne + Ahmad-Schreiber (2024) + a WIOD fallback combined (48 of 49 sectors), with `source`/`mapping_quality` tracked per sector -- see below. |
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
    'sectors':   None,          # 0-indexed positions into whatever sector list the data uses
                                 # (WIOD: 0..29; a directly-supplied haio: however many it has);
                                 # None = all sectors
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

## Running the paper scenario on ICIO instead of WIOD

`run_paper_scenario_icio.py` recreates the same EU-vs-Russia iceberg-cost
scenario as `run_paper_scenario.py`, but built entirely on real OECD ICIO
2022 data (`icio_to_haio.py`) and `combined_trade_elast.csv` instead of
WIOD 2008 and `trade_elast_2008.mat` -- same 40 individually-tracked
countries (mapped to their ICIO codes) plus a `ROW` composite, same EU
member list, same `ngrid=20`/`intensity=150` defaults, but all 48
non-`T`/non-`O` ICIO sectors instead of WIOD's 30.

This is **much slower** than the WIOD run: ICIO's 48-sector classification
makes the linear system roughly 1.6x larger per country (`C*N` grows from
1230 to 1968), and since the per-step solve cost scales with the cube of
the system size, expect **hours, not minutes**, at the default `ngrid=20`
-- run it in the background. It is also a materially different exercise,
not a like-for-like regression check against the WIOD numbers: different
base year (2022 vs. 2008), different sector classification, and
elasticities from three combined literature sources instead of the paper's
own (mostly-placeholder) `trade_elast_2008.mat` -- expect the magnitude to
differ from `run_paper_scenario.py`'s output for these reasons, not because
either run is wrong.

Since this script's sector set includes most of the goods sectors
`load_combined_trade_elast()` sources from Fontagné-Guimbard-Orefice, the
sign-convention caveat on that source (see "Getting `trade_elast`" above --
`trade_elast = -epsilon_icio` is a plausible but not independently verified
derivation) applies to roughly half of the 48 sectors this script uses, not
just to an isolated illustrative case.

## Data included

Three WIOD 2013-release, benchmark-year-2008 files -- everything `io_reorder.py`
actually reads (see the top-level `python replication/README.md` for the
full data provenance). The original MATLAB replication package
(`../../replication/baqaee_farhi_model/`) ships three additional `.mat`
files (`WIOD_SEA_14.mat`, `ahs_all.mat`, `wiot2008_row_apr12.mat`) used only
by `IO_reorder_init_tariff.m`, the "with initial tariffs" variant that was
never ported (the paper's own driver script doesn't use it either) -- those
are intentionally not duplicated here. Plus two small committed CSVs
`icio_to_haio.py` reads (real OECD ICIO files themselves, like
`2022_SML.csv`, are `.gitignore`d -- see "Getting the real data into your
own copy" above -- and must be downloaded separately).

| File | Contents |
|---|---|
| `wiott2008.mat` | World Input-Output Table: the 41-country x 31-industry bilateral flow matrix |
| `wiodsea2008.mat` | Socio-Economic Accounts: labor/capital compensation by country-industry |
| `trade_elast_2008.mat` | Industry-level cross-country (Armington-style) trade elasticities |
| `fontagne_icio_trade_elast.csv` | Raw Fontagné-Guimbard-Orefice ICIO-classification trade elasticities (28 usable sectors, goods only) |
| `combined_trade_elast.csv` | Fontagné + Ahmad-Schreiber (2024) + WIOD-fallback trade elasticities combined (48 of 49 non-`T` sectors), built by `build_combined_trade_elast.py` |

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
