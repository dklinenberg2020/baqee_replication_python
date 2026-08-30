"""
Python translation of replication/baqaee_farhi_model/main_dlogW_rev_bigshocks_EU_Russian_v2.m
(the outer discretization loop), restricted to shock_index == 3 (EU iceberg
trade cost on Russia), factor_index == 2 and initial_tariff_index == 1 --
what the paper's own driver script uses. See nested_ces.py's module
docstring for the memory/algorithmic notes that make this tractable.

The nonlinear equilibrium is approximated, as in the original code, by
discretizing the shock into `ngrid` small linear steps and updating the
economy's shares (Omega, chi, lambda, ...) after each step -- this is what
lets the *linearized* CES model (exact only for infinitesimal shocks)
approximate the effect of a *large* shock (here, an iceberg trade cost large
enough to send EU-Russia trade to roughly zero).
"""
import numpy as np
from io_reorder import io_reorder
from main_load_data import main_load_data
from nested_ces import value_added_shares, response, solve_dlambda_F_all


def run(keep_c, shocks, ngrid=20, sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05,
        data_dir='.', haio=None):
    """
    Parameters
    ----------
    keep_c : 1-indexed WIOD country numbers to track individually (ROW=35
        excluded, appended automatically as the last country). Ignored if
        `haio` is given directly.
    haio : optional pre-built standardized "BF HAIO" dict (see
        main_load_data.py's module docstring for the exact six-key
        contract: C, Omega, beta, alpha, alpha_VA, trade_elast,
        GDP_weights). Pass this instead of relying on `keep_c` to run the
        model on a non-WIOD dataset (GTAP, OECD ICIO, ...) without touching
        this function at all -- when given, `keep_c` and `data_dir` are
        ignored and `io_reorder()` (the WIOD-specific loader) is never
        called. A source-specific loader only needs to produce this dict;
        everything downstream of it (this function, nested_ces.py) is
        already source-agnostic.
    shocks : list of dicts, each describing one iceberg-cost shock leg. All
        country positions below are 1-indexed positions within whichever
        country ordering the data uses (keep_c's order for the WIOD path,
        or `haio['GDP_weights']`'s/`haio['Omega']`'s implied order for a
        directly-supplied haio) -- e.g. sellers=[3] means the 3rd country in
        that ordering, not a WIOD country code. Each leg:
            'sellers'   : positions of the seller country/countries whose
                          exports become more costly to buy (required).
            'buyers'    : positions of the buyer country/countries facing
                          the higher cost, or None for all countries
                          (default: None).
            'sectors'   : 0-indexed positions into whatever sector list the
                          data uses (0..N-1 -- WIOD's N=30, or however many
                          sectors a directly-supplied `haio` has), or None
                          for all sectors (default: None) -- e.g. a single
                          energy sector instead of an economy-wide trade
                          cost.
            'intensity' : shock size in percent -- see run_scenario()'s
                          module-level docs / README for how to interpret
                          this number. Applied evenly across `ngrid` steps.
        Multiple legs are summed, so independent shocks (e.g. a Gulf-wide
        cutoff and a separate bilateral effect) can be combined in one run.
        The paper's own EU-vs-Russia scenario is
        `shocks=[{'sellers': [RUS], 'buyers': EU, 'sectors': None, 'intensity': 150}]`.
    ngrid : number of discretization steps for the shock.

    Returns
    -------
    dlogW_sum : (C,) array, cumulative log change in real income by country
        (last entry is the appended ROW country).
    """
    if haio is None:
        Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights = io_reorder(keep_c, data_dir=data_dir)
        haio = dict(C=len(keep_c) + 1, Omega=Omega, beta=beta, alpha=alpha,
                    alpha_VA=alpha_VA, trade_elast=trade_elast, GDP_weights=GDP_weights)
    data, shock = main_load_data(haio, 1, 2)
    C, N, CN, CF = data['C'], data['N'], data['CN'], data['CF']
    data['sigma'], data['theta'], data['gamma'], data['epsilon'] = sigma, theta, gamma, epsilon
    # Pre-shock GNE weights, for a GNE-weighted world aggregate. Must be
    # captured here, before the loop below starts mutating data['chi_std']
    # step by step -- matching main_dlogW_rev_bigshocks_EU_Russian_v2.m's own
    # GNE_weights (captured once, before its loop), not silently drifting
    # into a post-shock weighting the way reading data['chi_std'] after the
    # loop would.
    data['GNE_weights_initial'] = data['chi_std'][:C].copy()

    F = data['F']
    Phi_F = np.zeros((C, CF))
    for c in range(C):
        Phi_F[c, c * F:(c + 1) * F] = 1
    data['Phi_F'] = Phi_F

    # Precompute each leg's 0-indexed positions and per-step intensity once.
    legs = []
    for leg in shocks:
        sellers = np.asarray(leg['sellers'], dtype=int) - 1
        buyers = leg.get('buyers')
        buyers = np.arange(C) if buyers is None else np.asarray(buyers, dtype=int) - 1
        sectors = leg.get('sectors')
        sectors = np.arange(N) if sectors is None else np.asarray(sectors, dtype=int)
        intensity_grid = np.log(1 + leg['intensity'] / 100) / ngrid
        legs.append((sellers, buyers, sectors, intensity_grid))

    dlogW = np.zeros((C, ngrid))

    for i in range(ngrid):
        dlogt = np.zeros((C + CN, CN + CF))
        dlogtau = np.zeros((C + CN, CN + CF))
        for sellers, buyers, sectors, intensity_grid in legs:
            for s in sellers:
                cols = s * N + sectors  # this seller's shocked-sector columns
                dlogtau[np.ix_(buyers, cols)] += intensity_grid  # buyers' households
                for b in buyers:
                    rows = C + b * N + np.arange(N)  # buyer b's producers, all its own sectors
                    dlogtau[np.ix_(rows, cols)] += intensity_grid
        dX = (data['Omega_total_tilde'][:C + CN, C:] * (dlogt + dlogtau)).sum(axis=1)
        shock = dict(dlogt=dlogt, dlogtau=dlogtau, dX=dX)

        blocks = value_added_shares(data)
        x = solve_dlambda_F_all(data, shock, blocks)
        g, dlogP_Vec, dOmega_N_goods, dOmega_N_own, dOmega_C, dchi_std, dlambda_result = response(x, data, shock, blocks)

        # Only the household (0..C-1) block of chi_std is ever nonzero, so
        # only it is needed here; avoid a spurious divide-by-zero warning
        # from the (always-zero) producer/factor block.
        dlogW[:, i] = dchi_std[:C] / data['chi_std'][:C] - dlogP_Vec[:C]

        # --- Update Omega_total_tilde and all quantities derived from it,
        # mirroring the MATLAB "Update variables" block -----------------
        L = data['L']
        dOmega_total_tilde = np.zeros((L, L))
        dOmega_total_tilde[:C, C:C + CN] = dOmega_C
        dOmega_total_tilde[C:C + CN, C:C + CN] = dOmega_N_goods
        dOmega_total_tilde[C + np.arange(CN), C + CN + np.arange(CN)] = dOmega_N_own

        Omega_total_tilde_new = np.clip(data['Omega_total_tilde'] + dOmega_total_tilde, 0, None)
        dOmega_total_tilde = Omega_total_tilde_new - data['Omega_total_tilde']
        row_sums = Omega_total_tilde_new.sum(axis=1, keepdims=True)
        Omega_total_tilde_new = np.divide(Omega_total_tilde_new, row_sums,
                                           out=np.zeros_like(Omega_total_tilde_new), where=row_sums != 0)
        data['Omega_total_tilde'] = Omega_total_tilde_new
        data['Omega_total'] = Omega_total_tilde_new

        Psi_total_tilde = np.linalg.solve(np.eye(L) - Omega_total_tilde_new, np.eye(L))
        data['Psi_total_tilde'] = Psi_total_tilde
        data['Psi_total'] = Psi_total_tilde

        data['chi_std'] = data['chi_std'] + dchi_std
        lambda_total = data['chi_std'] @ Psi_total_tilde
        data['lambda_CN'] = lambda_total[:C + CN]
        data['lambda_F'] = lambda_total[C + CN:]
        data['lambda_std'] = np.concatenate([data['lambda_CN'], data['lambda_F']])

        dbeta_temp = dOmega_total_tilde[:C, C:C + CN]
        dOmega_temp = dOmega_total_tilde[C:C + CN, C:]
        # Same alpha==1.0 edge case guarded in nested_ces.py's
        # value_added_shares(): a producer with zero measured intermediate
        # cost has an all-zero row here anyway, so the filler value is inert.
        one_minus_alpha_safe = np.where(data['alpha'] == 1, 1, 1 - data['alpha'])
        dOmega_temp2 = dOmega_total_tilde[C:C + CN, C:C + CN] * (1 / one_minus_alpha_safe)[:, None]
        dbeta_s = np.zeros((C, N))
        dOmega_s = np.zeros((CN, N))
        for c in range(C):
            dbeta_s += dbeta_temp[:, c * N:(c + 1) * N]
            dOmega_s += dOmega_temp2[:, c * N:(c + 1) * N]
        data['beta_s'] = data['beta_s'] + dbeta_s
        data['Omega_s'] = data['Omega_s'] + dOmega_s
        data['Omega_total_C'] = data['Omega_total_C'] + dbeta_temp
        data['Omega_total_N'] = data['Omega_total_N'] + dOmega_temp

    dlogW_sum = dlogW.sum(axis=1)
    return dlogW_sum, data


def run_scenario(keep_c, countries, shocks, ngrid=20,
                  sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05, data_dir='.', haio=None):
    """Convenience wrapper around `run()` for one-off scenario calls: runs the
    model and returns log GNE changes (dlogW, NOT percent -- multiply by 100
    yourself for a percentage) labeled by country code, the GNE-weighted
    world aggregate (same units), and the elasticity parameters that produced
    them, so a result is self-describing without needing to track down which
    run() call it came from. See `run()`'s docstring for the `shocks` format
    and for what passing `haio` directly (a non-WIOD dataset) does.

    Returns
    -------
    dict with:
        'dlogW' : {country code -> log GNE change}, one entry per `countries`
        'World' : GNE-weighted world aggregate log GNE change
        'elasticities' : {'sigma', 'theta', 'gamma', 'epsilon'} used for this run
        'ngrid', 'shocks' : shock discretization settings used
    """
    dlogW_sum, data = run(keep_c, shocks, ngrid=ngrid,
                           sigma=sigma, theta=theta, gamma=gamma, epsilon=epsilon,
                           data_dir=data_dir, haio=haio)
    world = data['GNE_weights_initial'] @ dlogW_sum
    return dict(
        dlogW=dict(zip(countries, dlogW_sum)),
        World=world,
        elasticities=dict(sigma=sigma, theta=theta, gamma=gamma, epsilon=epsilon),
        ngrid=ngrid,
        shocks=shocks,
    )
