#' The simple two-input CES production function of Appendix Subsection A.2 of
#' Bachmann, Baqaee, Bayer, Kuhn, Loeschel, Moll, Peichl, Pittel and
#' Schularick (2024), "What if? The macroeconomic and distributional effects
#' for Germany of a stop of energy imports from Russia", Economica 91(364):
#'
#'   Y = (alpha^(1/sigma) E^((sigma-1)/sigma) + (1-alpha)^(1/sigma) X^((sigma-1)/sigma))^(sigma/(sigma-1))   (A1)
#'
#' and its "full calibration" (Appendix Subsection A.9): given a target
#' expenditure share for energy (and, for the general case, energy and
#' non-energy prices pE, pX), solve for the share parameter alpha so that the
#' CES cost-minimizing demand system (A9)/(A10) reproduces that share
#' exactly, for ANY value of the elasticity of substitution sigma -- this is
#' the paper's key calibration property (Subsection A.9.1): "the model fits
#' the share of energy imports in German GNE for any value of the elasticity
#' of substitution sigma."
#'
#' This module supplements (not replaces)
#' `../python replication/ces_energy.py` and
#' `elasticity.py`/`elasticity_gas.py`, which hard-code alpha and sigma to
#' the paper's own two headline calibrations (0.04/0.04 and 0.01/0.1). Here
#' alpha is calibrated from an arbitrary target share and sigma from an
#' arbitrary elasticity, generalizing (A1)/(A9)/(A10) and their inversion so
#' they can be reused for the empirical elasticities in `sigma_literature.R`
#' or any other calibration target.

#' Energy's expenditure share pE*E/(P*Y) implied by cost minimization,
#' equation (A9)/(A10)'s share form. With pE == pX (the paper's own
#' normalization), this reduces to exactly alpha, independent of sigma.
ces_expenditure_share <- function(alpha, sigma, pE = 1.0, pX = 1.0) {
  if (sigma == 1.0) {
    return(alpha)
  }
  (alpha * pE^(1 - sigma)) / (alpha * pE^(1 - sigma) + (1 - alpha) * pX^(1 - sigma))
}

#' Price index P, equation (A10).
ces_price_index <- function(alpha, sigma, pE = 1.0, pX = 1.0) {
  if (sigma == 1.0) {
    return(pE^alpha * pX^(1 - alpha))
  }
  (alpha * pE^(1 - sigma) + (1 - alpha) * pX^(1 - sigma))^(1 / (1 - sigma))
}

#' Invert ces_expenditure_share() for alpha: the value of alpha such that a
#' CES(sigma) cost function facing prices (pE, pX) generates exactly
#' `target_share` as energy's expenditure share (Subsection A.9.1).
#'
#' With pE == pX, this is just alpha = target_share for every sigma -- the
#' paper's own calibration (e.g. alpha=0.04 to match a 4% GNE share of
#' gas+oil+coal, alpha=0.012 to match a 1.2% GNE share of gas alone). The
#' general form below also supports calibrating against non-unit relative
#' prices, which the paper's own scripts never need (they always normalize
#' pE = pX = 1) but which equation (A9) supports in general.
calibrate_alpha <- function(target_share, sigma = 1.0, pE = 1.0, pX = 1.0) {
  if (target_share <= 0 || target_share >= 1) {
    stop('target_share must be strictly between 0 and 1')
  }
  if (sigma == 1.0) {
    return(target_share)
  }
  ratio <- (target_share / (1 - target_share)) * (pX^(1 - sigma) / pE^(1 - sigma))
  ratio / (1 + ratio)
}

#' Optimal factor demands E, X for total expenditure PY, equation (A9).
#' Returns a list with elements E, X, P.
cost_minimizing_demand <- function(alpha, sigma, PY, pE = 1.0, pX = 1.0) {
  P <- ces_price_index(alpha, sigma, pE, pX)
  if (sigma == 1.0) {
    E <- alpha * PY / pE
    X <- (1 - alpha) * PY / pX
  } else {
    denom <- alpha * pE^(1 - sigma) + (1 - alpha) * pX^(1 - sigma)
    E <- alpha * pE^(-sigma) / denom * PY
    X <- (1 - alpha) * pX^(-sigma) / denom * PY
  }
  list(E = E, X = X, P = P)
}

#' Y from the production function itself, equation (A1). Vectorized in E.
ces_output <- function(E, X, alpha, sigma) {
  if (sigma == 1.0) {
    return(E^alpha * X^(1 - alpha))
  }
  if (sigma == 0.0) {
    return(pmin(E / alpha, X / (1 - alpha)))
  }
  (alpha^(1 / sigma) * E^((sigma - 1) / sigma) + (1 - alpha)^(1 / sigma) * X^((sigma - 1) / sigma))^(sigma / (sigma - 1))
}
