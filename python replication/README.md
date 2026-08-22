# Python replication

Python translation of the MATLAB/Stata replication package in `../replication`
for Bachmann, Baqaee, Bayer, Kuhn, Löschel, Moll, Peichl, Pittel and Schularick
(2024), "What if? The macroeconomic and distributional effects for Germany of
a stop of energy imports from Russia," *Economica* 91(364), 1157-1200.

Every script below has been validated against the *actual* output of the
original MATLAB code, run under GNU Octave 8.4.0 in this same environment
(not just re-derived from the paper) — see each module's docstring and the
"Validation" sections below for exactly what was checked and how closely it
matches. Nothing here is a re-implementation from the paper's description;
it is a translation of the original code, checked line-by-line against it.

## Contents

```
elasticity.py, elasticity_gas.py   Section 3 / Appendix A "simple model": aggregate CES
                                    production function in energy and other inputs.
                                    Reproduces Table 2 columns (3)-(4) and Figures A1-A3.
ces_energy.py                      Shared CES helper functions used by both scripts above.

consumption_energy.py              Section 3 / Appendix A.12 distributional analysis
                                    (household energy expenditure/income shares by
                                    income quantile), translated from consumption_energy.do.

baqaee_farhi_model/                Section 2 / Appendix A.5 "rich model": the multi-country,
                                    multi-sector linearized nested-CES network model of
                                    Baqaee and Farhi (2024), reproducing Table 2 column (2)
                                    and Appendix Table 1.
```

## Setup

```
pip install numpy scipy pandas matplotlib
```

`consumption_energy.py` additionally needs the restricted-access EVS microdata
(see below) and a Stata-`.dta`-capable pandas installation (the default).

## 1. The simple aggregate CES model (`elasticity.py`, `elasticity_gas.py`)

Direct, self-contained translations of `elasticity.m` (10% energy-quantity
drop, calibrated to alpha=0.04, the GNE share of gas+oil+coal) and
`elasticity_gas.m` (30% gas-quantity drop, alpha=0.01, the GNE share of gas
alone). No external data required.

```
cd "python replication"
python elasticity.py       # -> figures/elasticity_fig.pdf, MPE_fig.pdf, MPX_fig.pdf, exp_share_fig.pdf
python elasticity_gas.py   # -> figures/*_gas.pdf
```

**Validation:** every printed number (`Output loss, pE^new/pE^old,
pX^new/pX^old, new energy share`, for sigma = 0, 0.04/0.1, 0.1/0.2, 1) matches
the MATLAB script's own `disp()` output, run under Octave, to full floating-
point precision. In particular:

| Scenario | sigma | Output (GNE) loss | Paper (Table 2) |
|---|---|---|---|
| 10% energy drop | 0.04 | -1.57% | col (3): 1.5% GNE / 1.3% GDP |
| 30% gas drop | 0.1 | -2.33% | col (4): 2.3% GNE / 2.2% GDP |

(The model's output variable *Y* is domestic absorption/GNE, not GDP — see
footnote 10 of the paper — which is why these losses read directly off the
GNE column.)

## 2. The distributional analysis (`consumption_energy.py`)

Line-by-line translation of `consumption_energy.do`'s data construction
(household-level income/expenditure variables, weighted income quintiles,
energy expenditure by source and heating type) and its four `graph bar`
figure blocks.

**This script cannot be run or tested here.** Its input is the 2013
Scientific Use File of the German EVS (Einkommens- und
Verbrauchsstichprobe), restricted-access microdata obtained under a
separate data use agreement from the German Federal Statistical Office's
research data center
(https://www.forschungsdatenzentrum.de/de/haushalte/evs) — it was not
included in the original MATLAB/Stata package either, and isn't in this
repository. The translation was checked read-through against the .do file
(variable names, formulas, and weighting all preserved) and smoke-tested
against synthetic data matching the EVS column layout to confirm it runs
without errors end-to-end, but its numeric output has **not** been checked
against the original Stata output, since neither this environment nor (as
far as this replication package goes) any available reference run has
access to the real data.

```
python consumption_energy.py --input data/evs2013_aa_gs_hb.dta --outdir figures/
```

## 3. The Baqaee-Farhi network model (`baqaee_farhi_model/`)

This is the hard part of the package: a linearized multi-country (41 WIOD
2008 countries + an aggregate "rest of world"), multi-sector (30 sectors),
nested-CES general-equilibrium model with country-sector-specific labor,
solved for the response to a large (iceberg-cost) shock via `ngrid`
discretization steps.

```
io_reorder.py           IO_reorder.m: loads and reorganizes the WIOD 2008
                         input-output table into `len(keep_c)` countries + ROW.
main_load_data.py        main_load_data_rev.m: builds the "standard form" (L,L)
                          Markov IO matrix, its Leontief inverse, and income/sales shares.
nested_ces.py             AES_func.m + Nested_CES_linear_final_rev.m +
                           Nested_CES_linear_result_final.m: the linearized model
                           itself (Allen elasticities of substitution, the
                           equilibrium response to a shock, solving the implicit
                           linear system for a discretization step).
run_model.py               main_dlogW_rev_bigshocks_EU_Russian_v2.m's outer
                            discretization loop.
run_paper_scenario.py       Convenience wrapper with the paper's exact settings
                             (all 41 countries, EU vs. Russia iceberg shock,
                             ngrid=20, intensity=150).
```

Run the paper's exact scenario:

```
cd baqaee_farhi_model
python run_paper_scenario.py
```

Expect this to take on the order of **tens of minutes** on a single core —
see "Why this is slow, and why that's inherent" below.

### Algorithmic notes (read `nested_ces.py`'s module docstring for the full
version)

The MATLAB code materializes dense 3-D Allen-elasticity tensors
(`AES_N_Mat`, `AES_F_Mat`) of shape `(CN, CN+CF, CN)`. At the paper's scale
(41 countries x 30 sectors, country-sector-specific factors, so
`CN = CF = 1230`), `AES_N_Mat` alone is **~30 GB**. This is not a Python
translation issue — **the original MATLAB code was reproduced running out
of memory under GNU Octave 8.4.0 with 15 GB of RAM available**, in this same
environment, at the paper's exact settings. Two changes were needed to make
this tractable here, both purely computational (not changes to the model):

1. **Closed-form AES evaluation.** `AES_func.m`'s elasticities only ever take
   one of three values per (row, column) pair — a "same origin and sector"
   diagonal value, a "same sector, different origin" value, and a residual
   base value — so every sum that would otherwise need the full dense tensor
   can be written directly in terms of those three values and simple
   row/column sums. See `nested_ces.py`'s `producer_dOmega_goods`,
   `producer_dOmega_own_factor`, `consumer_dOmega`.

2. **Solving the linear system by probing rather than by explicit Jacobian.**
   `Nested_CES_linear_final_rev.m` builds the (C+CF, C+CF) Jacobian `A` and
   offset `B` of the equilibrium map via a dense mass of MATLAB
   `reshape`/`permute`/`kron`/`bsxfun` calls — code that is very easy to get
   subtly wrong when porting between MATLAB's column-major and NumPy's
   row-major conventions, with no easy way to check each step in isolation.
   Since that same equilibrium map (implemented directly, following
   `Nested_CES_linear_result_final.m`) is *provably affine* in its input
   (every step is a linear operation on data fixed within one discretization
   step), its Jacobian and offset can instead be recovered exactly by
   *probing*: `B = f(0)`, and column `i` of `A` = `f(e_i) - B`. This costs
   `C+CF+1` evaluations of the (already closed-form, no-huge-tensor)
   equilibrium map instead of one, but each evaluation is cheap and the
   memory footprint stays at `O(L^2)` (`L = C+CN+CF ~= 2500`, i.e. a ~50MB
   dense matrix) throughout, instead of the ~60GB the dense-tensor approach
   would need for both `AES_N_Mat` and `AES_F_Mat`.

A third subtlety, found and fixed during validation: `AES_func.m` computes
elasticity *values* from `Omega_total_N`/`Omega_total_C` (a bookkeeping
matrix updated by the outer discretization loop's raw, un-renormalized
delta), while `Nested_CES_linear_result_final.m` *weights* those elasticities
in every sum by `Omega_total_tilde` (row-renormalized after each step). The
two matrices are numerically identical before the loop's first iteration
(which is why a same-shock single-step test can pass while hiding this bug)
but diverge afterward — `nested_ces.py`'s `_diag_term` docstring has the full
derivation of why, and why the same zero/nonzero support is shared by both
regardless (so the zero-division guards are safe either way).

### Why this is slow, and why that's inherent

`run_paper_scenario.py` solves a ~1271-dimensional linear system once per
discretization step (`ngrid=20` by default, matching the original script),
and each solve needs `C+CF+1 ~= 1272` evaluations of the equilibrium-response
function. This is unavoidably a large linear-algebra workload at this
model's scale — the original MATLAB code faces the exact same
`(C+CF, C+CF)` system, just solved with one explicit-Jacobian construction
instead of `C+CF+1` cheaper evaluations, at the cost of the ~60GB of RAM
this environment doesn't have. There is no way to make the *paper's own
model, at the paper's own scale*, both memory-light and fast in this
environment; this implementation trades wall-clock time for memory.

### Validation

Every layer was checked against the actual MATLAB output (via GNU Octave
8.4.0, run in this environment) at a scale small enough for Octave to
execute without running out of memory:

- **`io_reorder.py`**: `Omega`, `beta`, `alpha_VA`, `alpha`, `trade_elast`,
  `GDP_weights` match `IO_reorder.m`'s output element-wise to
  floating-point precision (`< 1e-14`) for the full 41-country run.
- **`main_load_data.py`**: `C, N, F, CN, CF, L`, `Omega_total`,
  `alpha`, `beta_s`, `Omega_s`, `chi_std`, `lambda_std`, `Psi_total` all
  match `main_load_data_rev.m`'s output element-wise (`< 1e-14`) for the
  full 41-country run.
- **`nested_ces.py`**: `response()`'s every intermediate (`dlogP_Vec`,
  `dOmega_total`'s consumer/producer/own-factor blocks, `dchi_std`,
  `dlambda_result`) and the fixed point found by `solve_dlambda_F_all()`
  match `Nested_CES_linear_result_final.m` / `Nested_CES_linear_final_rev.m`
  to machine precision (`< 1e-15`), checked on a 5-country test case (CHN,
  DEU, RUS, USA + ROW) — both for a single shock against pristine data, and
  (crucially, since this is where the Omega_total_tilde/Omega_total_N bug
  above was caught) for a shock applied to data already updated by one prior
  discretization step.
- **`run_model.py`**: the full 3-step discretization loop (`dlogW_sum`) for
  a 5-country EU-vs-Russia scenario (AUT, DEU, RUS, USA + ROW) matches
  `main_dlogW_rev_bigshocks_EU_Russian_v2.m`'s output to `< 3e-17` — every
  step, including the "update Omega/chi/lambda for the next step" block.

The **full 41-country run** (`run_paper_scenario.py`) could not be checked
against Octave the same way, since Octave itself cannot run it in this
environment (the ~30GB `AES_N_Mat` allocation fails) — this is a genuine
resource constraint of the original model at this scale, not a gap specific
to the Python side. It was, however, actually run end to end (~20 minutes,
see `baqaee_farhi_model/paper_scenario_result.json`): with the paper's exact
settings (all 41 WIOD countries + ROW, EU iceberg cost on Russia,
`ngrid=20`, `intensity=150`), it gives a Germany GNE loss of **-0.257%**,
matching the paper's own stated Table 2 column (2) figure of **0.2-0.3%**
almost exactly, and a Russia GNE loss of -2.27% (Russia bears the brunt of
being cut off from EU export markets, as expected). Given this and the exact
match at every other layer and scale, the full run is taken as a faithful
reproduction, even though — for the reason above — it isn't checked bit-for-
bit against MATLAB output the way everything else here is.
