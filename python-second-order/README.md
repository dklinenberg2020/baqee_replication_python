# Second-order approximation and calibration (Python)

A lightweight, self-contained toolkit for the closed-form parts of Appendix
A of Bachmann, Baqaee, Bayer, Kuhn, Löschel, Moll, Peichl, Pittel and
Schularick (2024), "What if? The macroeconomic and distributional effects
for Germany of a stop of energy imports from Russia," *Economica* 91(364):

- **Full calibration** (`ces.py`, Subsection A.9): calibrate the simple CES
  production function's share parameter alpha to match ANY target
  expenditure share, at ANY relative prices and ANY elasticity of
  substitution sigma — generalizing the paper's own two hard-coded
  calibrations (`../python replication/elasticity.py`'s alpha=0.04/0.01).
- **Second-order approximation** (`second_order.py`, Lemmas 1 and 2): the
  paper's closed-form Taylor approximations for the welfare/output loss from
  an energy shock — both the simple two-input model (equations A4/A5) and,
  strikingly, a version that applies to the FULL Baqaee-Farhi multi-country
  network model (equations A6/A7) using only the change in one number (the
  energy import share of GNE) — no simulation required.
- **Empirical sigma calibration** (`sigma_literature.py`, Subsection A.4): a
  small reference table of the literature elasticity estimates the paper
  reviews (Labandeira et al. 2017; Auffhammer and Rubin 2018), alongside the
  paper's own (deliberately more conservative) choices.

Depends only on the Python standard library — no numpy/scipy needed for
this package, since every formula here is closed-form scalar algebra. Run
`python examples.py` to see every number below reproduced and checked.

## What this package is *for*

This is not a replacement for the exact computation in
`../python replication/elasticity.py` or the full network simulation in
`../python replication/baqaee_farhi_model/` — it's the fast, closed-form
complement to both: a few lines of algebra instead of solving a CES
optimization or a ~1300-dimensional linear system. Use it to explore
calibrations quickly, or to sanity-check the exact/simulated results against
a simple formula — **not** as a substitute for either, especially at low
sigma (see "Where this breaks down" below, and the paper's own footnote 12).

## Quick start

```python
from ces import calibrate_alpha, ces_expenditure_share
from second_order import second_order_simple, second_order_from_share
from sigma_literature import sigma_from_literature

# Calibrate alpha to match a 4% GNE share of gas+oil+coal (Subsection A.9.1)
alpha = calibrate_alpha(target_share=0.04)          # -> 0.04, for any sigma
sigma = sigma_from_literature('paper.aggregate_energy')  # -> 0.04

# Second-order approximation to a 10% energy-supply drop (equation A4)
dlogY = second_order_simple(alpha_tilde=alpha, dlogE=-0.10, sigma=sigma)
print(f'{100*dlogY:.2f}%')   # -0.86% (see "Where this breaks down": the exact
                              # answer, from elasticity.py, is -1.57%)

# Sufficient-statistic (network-model) approximation, equation (A7):
# gas-only scenario, import share assumed to triple from 1.2% to 3.6%
dlogW = second_order_from_share(share=0.012, dlogE=-0.30, dshare=0.024)
print(f'{100*dlogW:.2f}%')   # -0.72%, the paper's own equation (A8)
```

## Files

| File | Paper section | Contents |
|---|---|---|
| `ces.py` | A.2, A.9 | CES production function, cost-minimizing demand (A9)/(A10), and `calibrate_alpha()` (the general inverse of A9/A10 for alpha) |
| `second_order.py` | Lemma 1 (A4/A5), Lemma 2 (A6/A7) | The Taylor approximations, for both the simple model and the sufficient-statistic network-model version |
| `sigma_literature.py` | A.4, A.9.2 | Reference elasticity estimates and the paper's own calibration choices |
| `examples.py` | — | Validation script: reproduces every number checked below |

## Validated against the paper

Running `python examples.py` checks, and prints:

1. **Equations (A4) and (A5) are the same approximation.** The paper states
   this but its own printed formula for Δ(share) (right under equation A5)
   turns out to be missing a factor of α̃ — see the "apparent erratum" note
   in `second_order.py`'s docstring for the derivation and a numerical check
   against a finite-difference derivative of the true CES expenditure share.
   This package uses the corrected formula throughout, which is the only
   version making (A4) and (A5) actually equal (as the paper claims) *and*
   which matches the independent re-derivation.
2. **`calibrate_alpha()` reproduces the target share exactly**, for every
   sigma, at pE=pX=1 (Subsection A.9.1's stated calibration property) and at
   non-unit relative prices (the general form of equations A9/A10).
3. **The three back-of-envelope Δlog W numbers of Subsection A.5.3**
   (equation A7 applied to three scenarios): -1.5%, -0.63%, and the paper's
   own preferred estimate, equation (A8)'s -0.72% — all reproduced to within
   the paper's own rounding.
4. **Where the approximation breaks down** (footnote 12): comparing equation
   (A4) against the *exact* computation
   (`../python replication/elasticity.py`, itself validated against the
   original MATLAB code) at the paper's own two calibrations shows a
   divergence of 0.7-1.6 percentage points — this is exactly why the paper
   computes its headline numbers exactly rather than via (A4), and why this
   package is a complement to, not a replacement for, the exact/simulated
   results.

## See also

- `../r-second-order/`: the same package in R.
- `../python replication/`: the exact computation (`elasticity.py`,
  `elasticity_gas.py`) and the full Baqaee-Farhi network simulation.
- `../README.md`: the top-level overview of all four packages.
