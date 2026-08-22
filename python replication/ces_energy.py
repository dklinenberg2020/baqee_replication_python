"""
Shared CES production-function helpers for the two-input model
Y = (alpha^(1/sigma) E^((sigma-1)/sigma) + (1-alpha)^(1/sigma) X^((sigma-1)/sigma))^(sigma/(sigma-1))

used by elasticity.py and elasticity_gas.py, translated from the MATLAB
originals replication/elasticity.m and replication/elasticity_gas.m.

E is the energy input, X the composite of all other inputs (held fixed at its
initial endowment level, X = (1-alpha)/alpha, so that E=1 is the pre-shock
benchmark with price pE = pX = 1). alpha is the initial expenditure share on
energy. See Appendix Subsections A.2-A.4 of Bachmann et al. (2024) for the
economics.
"""
import numpy as np


def ces_quantities(E, alpha, sigma, pX=1.0):
    """Output, marginal products and energy expenditure share for a CES(sigma)
    aggregator, given a grid of energy quantities E and the endowment of the
    other input X implied by alpha (see module docstring)."""
    E = np.asarray(E, dtype=float)
    X = (1 - alpha) / alpha

    Y = (alpha ** (1 / sigma) * E ** ((sigma - 1) / sigma)
         + (1 - alpha) ** (1 / sigma) * X ** ((sigma - 1) / sigma)) ** (sigma / (sigma - 1))
    Ybench = (alpha ** (1 / sigma) + (1 - alpha) ** (1 / sigma) * X ** ((sigma - 1) / sigma)) ** (sigma / (sigma - 1))

    MPE = (alpha ** (1 / sigma) * E ** ((sigma - 1) / sigma) + (1 - alpha) ** (1 / sigma) * X ** ((sigma - 1) / sigma)) \
        ** (sigma / (sigma - 1) - 1) * alpha ** (1 / sigma) * E ** (-1 / sigma)
    MPX = (alpha ** (1 / sigma) * E ** ((sigma - 1) / sigma) + (1 - alpha) ** (1 / sigma) * X ** ((sigma - 1) / sigma)) \
        ** (sigma / (sigma - 1) - 1) * (1 - alpha) ** (1 / sigma) * X ** (-1 / sigma)

    P = (alpha * MPE ** (1 - sigma) + (1 - alpha) * pX ** (1 - sigma)) ** (1 / (1 - sigma))
    exp_share = (MPE * E) / (P * Y)

    return Y, Ybench, MPE, MPX, exp_share


def ces_cobb_douglas(E, alpha, X=None):
    """sigma -> 1 (Cobb-Douglas) limit, evaluated directly to avoid 0/0."""
    E = np.asarray(E, dtype=float)
    if X is None:
        X = (1 - alpha) / alpha
    Y = E ** alpha * X ** (1 - alpha)
    Ybench = 1.0 ** alpha * X ** (1 - alpha)
    MPE = alpha * E ** (alpha - 1) * X ** (1 - alpha)
    MPX = (1 - alpha) * E ** alpha * X ** (-alpha)
    exp_share = np.full_like(E, alpha)
    return Y, Ybench, MPE, MPX, exp_share


def ces_leontief(E, alpha, X=None):
    """sigma -> 0 (Leontief) limit: Y = min(E/alpha, X/(1-alpha))."""
    E = np.asarray(E, dtype=float)
    if X is None:
        X = (1 - alpha) / alpha
    Y = np.minimum(E / alpha, X / (1 - alpha))
    Ybench = min(1 / alpha, X / (1 - alpha))

    MPE = np.where(E < 1, 1 / alpha, np.where(E == 1, 1.0, 0.0))
    MPX = np.where(E < 1, 0.0, np.where(E == 1, 1.0, 1 / (1 - alpha)))
    exp_share = np.where(E < 1, 1.0, np.where(E == 1, alpha, 0.0))

    return Y, Ybench, MPE, MPX, exp_share


def nearest_index(E, target):
    """Index of the grid point in E closest to target (mirrors MATLAB's
    `[value, index] = min(abs(E - target))`)."""
    return int(np.argmin(np.abs(np.asarray(E) - target)))


def print_row(label, output_loss, pE_rel, pX_rel, share):
    print(label)
    print('Output loss, pE^new/pE^old, pX^new/pX^old, new energy share')
    print(f'  {output_loss: .6f}  {pE_rel}  {pX_rel}  {share: .6f}')
