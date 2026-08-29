"""
Standalone-package entry point for the Baqaee-Farhi linearized nested-CES
network model (Section 2 / Appendix A.5 of Bachmann et al. 2024; see
README.md in this directory).

The individual modules (io_reorder.py, main_load_data.py, nested_ces.py,
run_model.py) use plain flat imports (`from io_reorder import io_reorder`,
etc.) rather than package-relative ones, so that they keep working exactly
as validated -- run directly as scripts from inside this directory, per the
original usage pattern -- without touching any of the numerically-validated
files. This __init__.py instead makes the *directory* importable as a
package from anywhere by putting it on sys.path once, then re-exporting the
public functions.

It also fixes the one thing that would otherwise break when this package is
imported from a different working directory: `run`/`run_scenario` default
`data_dir` to this directory (where the WIOD .mat files actually live)
instead of the caller's current working directory.
"""
import os
import sys

_PKG_DIR = os.path.dirname(os.path.abspath(__file__))
if _PKG_DIR not in sys.path:
    sys.path.insert(0, _PKG_DIR)

import run_model as _run_model
import main_load_data as _main_load_data
import io_reorder as _io_reorder
import icio_to_haio as _icio_to_haio
import nested_ces as _nested_ces


def run(keep_c, shocks, ngrid=20, sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05, data_dir=None, haio=None):
    """See run_model.run()'s docstring for the `shocks` format and for what
    passing `haio` directly (a non-WIOD dataset) does. `data_dir` defaults
    to this package's own directory (where the .mat files live) rather than
    the caller's cwd; irrelevant when `haio` is given."""
    return _run_model.run(keep_c, shocks, ngrid=ngrid, sigma=sigma, theta=theta,
                           gamma=gamma, epsilon=epsilon, data_dir=data_dir or _PKG_DIR, haio=haio)


def run_scenario(keep_c, countries, shocks, ngrid=20, sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05,
                  data_dir=None, haio=None):
    """See run_model.run_scenario()'s docstring. `data_dir` defaults to this
    package's own directory rather than the caller's cwd; irrelevant when
    `haio` is given."""
    return _run_model.run_scenario(keep_c, countries, shocks, ngrid=ngrid, sigma=sigma, theta=theta,
                                    gamma=gamma, epsilon=epsilon, data_dir=data_dir or _PKG_DIR, haio=haio)


main_load_data = _main_load_data.main_load_data
io_reorder = _io_reorder.io_reorder
icio_to_haio = _icio_to_haio.icio_to_haio
value_added_shares = _nested_ces.value_added_shares
response = _nested_ces.response
solve_dlambda_F_all = _nested_ces.solve_dlambda_F_all

__all__ = [
    'run', 'run_scenario', 'main_load_data', 'io_reorder', 'icio_to_haio',
    'value_added_shares', 'response', 'solve_dlambda_F_all',
]
