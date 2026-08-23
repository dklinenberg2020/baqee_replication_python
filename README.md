# Replication packages for Bachmann et al. (2024)

Replication materials for Bachmann, Baqaee, Bayer, Kuhn, Löschel, Moll,
Peichl, Pittel and Schularick (2024), "What if? The macroeconomic and
distributional effects for Germany of a stop of energy imports from
Russia," *Economica* 91(364), 1157-1200
([paper](Economica%20-%202024%20-%20Bachmann%20-%20What%20if%20%20The%20macroeconomic%20and%20distributional%20effects%20for%20Germany%20of%20a%20stop%20of%20energy%20imports.pdf)).

## What's here

| Folder | Language | What it is |
|---|---|---|
| `replication/` | MATLAB / Stata | The **original** replication package released with the paper. Unmodified. |
| `python replication/` | Python | A full translation of `replication/`: the exact simple CES model, the full Baqaee-Farhi multi-country network simulation, and the (untestable, restricted-data) distributional analysis. See its own README for what's validated and how. |
| `python-second-order/` | Python | A lightweight, closed-form companion package: the paper's *second-order approximation* (Appendix Lemmas 1-2) and *full calibration* procedure (Appendix A.9) — algebra, not simulation. |
| `r-second-order/` | R | The same closed-form package as `python-second-order/`, ported line-for-line for R users. |

The two "second-order" packages are **not** replacements for
`python replication/` — they're the fast, closed-form complement to it: a
few lines of algebra instead of solving a CES optimization or a
~1300-dimensional linear system. See "Which package do I want?" below.

## Which package do I want?

- **Want the paper's exact headline numbers** (Table 2, Figures A1-A3, the
  full multi-country simulation, Table 2 column 2 / Appendix Table 1)? →
  `python replication/`. This is the one that was actually run end-to-end
  and checked against the original MATLAB output (see its README for
  exactly how closely, and what wasn't checkable).
- **Want to quickly try a different calibration** (a different target
  expenditure share, a different sigma from the empirical literature, a
  different assumed change in that share) without running any simulation? →
  `python-second-order/` or `r-second-order/`, whichever language you work
  in. These reproduce the paper's own back-of-envelope numbers
  (Subsection A.5.3, equation A8) and let you swap in your own inputs.
- **Want the original MATLAB/Stata code as the authors released it?** →
  `replication/`.

## A note on scope: why isn't there an R port of the full model?

`python replication/baqaee_farhi_model/` is a genuine linearized
general-equilibrium simulation — a ~1300-dimensional linear system solved
once per discretization step, needing real validation effort (every layer
was checked against the original MATLAB output element-wise; see that
README's "Algorithmic notes" for a subtle bug that validation actually
caught). Porting that to R as well would mean redoing all of that
validation a second time in a second language, for a model that already has
a validated Python implementation sitting right next to where an R version
would go. The two *second-order* packages, by contrast, are pure closed-form
algebra — porting them to both languages costs little and both are
validated the same way (`examples.py` / `examples.R` produce identical
output, checked side by side while building them). If you need the full
network model from R, the practical path is to shell out to
`python replication/baqaee_farhi_model/run_paper_scenario.py` (e.g. via
`reticulate` or `system2()`) rather than a from-scratch R port.

## An erratum found while building the second-order packages

The published paper's formula for Δ(energy expenditure share) directly
under equation (A5) is missing a factor of α̃ (the paper states equations A4
and A5 are the same approximation written two ways, but the printed Δ-share
formula doesn't actually make that true). `python-second-order/second_order.py`
and `r-second-order/second_order.R` both use the corrected formula, with the
full derivation — including a numerical check against a finite-difference
derivative of the exact CES expenditure share — in `second_order.py`'s /
`second_order.R`'s header comment.

## Setup

```
# python replication/ and python-second-order/
pip install -r "python replication/requirements.txt"   # numpy, scipy, pandas, matplotlib
cd python-second-order && python examples.py            # no extra dependencies

# r-second-order/
cd r-second-order && Rscript examples.R                 # base R only, no packages required
```
