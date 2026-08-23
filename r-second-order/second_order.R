#' The second-order (Taylor) approximations of Appendix Lemma 1 and Lemma 2
#' of Bachmann et al. (2024), giving the output/welfare loss from an
#' energy-supply shock WITHOUT solving the full CES optimization (`ces.R`)
#' or, for the network-model version, WITHOUT running the full Baqaee-Farhi
#' simulation in `../python replication/baqaee_farhi_model/` at all (that
#' simulation exists only in the Python package -- see this folder's
#' README.md).
#'
#' Two versions, both derived (see the paper's Subsections A.10/A.11) from a
#' second-order Taylor expansion of a homothetic aggregator around the
#' initial equilibrium:
#'
#' - The SIMPLE model (Lemma 1, equations A4/A5): a single aggregate CES
#'   production function in energy E and everything else X (`ces.R`).
#' - The NETWORK/sufficient-statistic model (Lemma 2, equations A6/A7):
#'   applies to the FULL Baqaee-Farhi multi-country, multi-sector model (or,
#'   per the paper, "a much wider class of general equilibrium models" than
#'   just that one) -- Lemma 2's striking implication is that the entire
#'   input-output structure collapses into a single sufficient statistic,
#'   the change in energy's expenditure share, so this needs no simulation
#'   at all.
#'
#' IMPORTANT CAVEAT (paper's own footnote 12, and see `examples.R`): the
#' second-order approximation can be inaccurate for sigma very close to
#' zero -- this is exactly why the paper's own headline numbers (Table 2,
#' Figures A1-A3, `../python replication/elasticity.py`) are computed
#' EXACTLY rather than via this approximation. Use this module for quick
#' back-of-envelope estimates and for cross-checking the exact/simulated
#' results, not as a replacement for either.
#'
#' APPARENT ERRATUM IN THE PUBLISHED PAPER: the printed formula directly
#' under equation (A5), "Delta(pE*E/PY) ~= (1-1/sigma)*(1-alpha_tilde)*Delta
#' log E", is missing a factor of alpha_tilde. The paper states that (A4)
#' and (A5) are the same approximation written two ways, but substituting
#' the printed Delta-share formula into (A5) does NOT reproduce (A4) unless
#' alpha_tilde=1 (the two second-order coefficients are
#' (1-1/sigma)*(1-alpha_tilde) vs. the required
#' (1-1/sigma)*alpha_tilde*(1-alpha_tilde)). `implied_dshare()` below uses
#' the alpha_tilde*(1-alpha_tilde) form instead, for two independent
#' reasons, both checked in examples.R:
#'   1. It is the only version that makes (A4) and (A5) algebraically
#'      equivalent, as the paper's own text claims they are.
#'   2. It matches a direct re-derivation from the CES optimality
#'      conditions: since
#'      share/(1-share) = (alpha/(1-alpha))^(1/sigma) * (E/X)^((sigma-1)/sigma)
#'      at fixed X, d(log(share/(1-share)))/d(log E) = 1-1/sigma, and the
#'      chain rule (d(share) = share*(1-share)*d(log(share/(1-share))))
#'      gives d(share)/d(log E) = alpha_tilde*(1-alpha_tilde)*(1-1/sigma)
#'      exactly -- confirmed numerically against a finite-difference
#'      derivative of ces_expenditure_share() to 3+ significant figures (the
#'      printed formula is off by orders of magnitude for small
#'      alpha_tilde, since it omits that factor).

#' Equation (A4): second-order approximation to Delta log Y for the simple
#' CES(sigma) model, in terms of the CURRENT energy expenditure share
#' alpha_tilde (== alpha under the paper's pE=pX=1 normalization, see
#' ces.R) and the log change in energy dlogE.
#'
#' Deliberately does NOT special-case sigma == 0 (Leontief): the
#' (1 - 1/sigma) term blows up there by construction, which is exactly the
#' paper's own point (footnote 12) that this second-order approximation is
#' unreliable for sigma close to zero -- the paper computes its own headline
#' numbers exactly instead, never via this formula. See examples.R for a
#' demonstration of the breakdown.
second_order_simple <- function(alpha_tilde, dlogE, sigma) {
  alpha_tilde * dlogE + 0.5 * (1 - 1 / sigma) * alpha_tilde * (1 - alpha_tilde) * dlogE^2
}

#' The first-order-accurate change in the energy expenditure share,
#' Delta(pE*E/PY) ~= (1 - 1/sigma)*alpha_tilde*(1-alpha_tilde)*dlogE, used to
#' build equation (A5) from equation (A4) (see the paragraph following
#' Lemma 1). NOTE: this includes an alpha_tilde factor the paper's own
#' printed formula omits -- see this file's header comment's "APPARENT
#' ERRATUM" note for why. Same caveat about sigma == 0 as
#' second_order_simple() above.
implied_dshare <- function(alpha_tilde, dlogE, sigma) {
  (1 - 1 / sigma) * alpha_tilde * (1 - alpha_tilde) * dlogE
}

#' The first-order (Lemma 2's first equation) term alone: share * dlogE.
first_order_from_share <- function(share, dlogE) {
  share * dlogE
}

#' Generic sufficient-statistic form used by BOTH equation (A5) (simple
#' model, share = alpha_tilde, dshare = implied_dshare(...)) and equations
#' (A6)/(A7) (network model / Lemma 2, share = the energy import share of
#' GNE, dshare = an assumed or model-implied change in that share):
#'
#'   Delta log W ~= share * dlogE + 0.5 * dshare * dlogE
#'
#' This is also exactly equation (A7)'s form when `share`/`dlogE` refer to
#' energy IMPORTS specifically (pE*mE/GNE, Delta log mE) rather than total
#' energy use -- the paper's whole point in Lemma 2 is that this same
#' formula holds "in a much richer open-economy model with a complex
#' production network" (see this file's header comment).
second_order_from_share <- function(share, dlogE, dshare) {
  first_order_from_share(share, dlogE) + 0.5 * dshare * dlogE
}

#' Equation (A5): algebraically identical to second_order_simple(...)
#' (equation A4), just re-expressed via implied_dshare(...) -- see the
#' paper's remark right after Lemma 1 that this is "an advantage ... it is
#' likely easier to decide on what is a reasonable change in the
#' expenditure share than what is a reasonable elasticity of substitution."
#' Included mainly to demonstrate/verify that equivalence -- see examples.R.
second_order_share_implied <- function(alpha_tilde, dlogE, sigma) {
  second_order_from_share(alpha_tilde, dlogE, implied_dshare(alpha_tilde, dlogE, sigma))
}
