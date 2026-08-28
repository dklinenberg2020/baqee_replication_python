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
from main_load_data import main_load_data
from nested_ces import value_added_shares, response, solve_dlambda_F_all


def run(keep_c, shocks, ngrid=20, sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05,
        data_dir='.'):
    """
    Parameters
    ----------
    keep_c : 1-indexed WIOD country numbers to track individually (ROW=35
        excluded, appended automatically as the last country).
    shocks : list of dicts, each describing one iceberg-cost shock leg. All
        country positions below are 1-indexed *positions within keep_c* (not
        WIOD country numbers), matching keep_c itself -- e.g. sellers=[3]
        means the 3rd entry of keep_c, not WIOD country 3. Each leg:
            'sellers'   : positions of the seller country/countries whose
                          exports become more costly to buy (required).
            'buyers'    : positions of the buyer country/countries facing
                          the higher cost, or None for all countries
                          (default: None).
            'sectors'   : 0-indexed WIOD sector indices (0..N-1, N=30)
                          affected, or None for all sectors (default: None)
                          -- e.g. a single energy sector instead of an
                          economy-wide trade cost.
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
    data, shock = main_load_data(keep_c, 1, 2, data_dir=data_dir)
    C, N, CN, CF = data['C'], data['N'], data['CN'], data['CF']
    data['sigma'], data['theta'], data['gamma'], data['epsilon'] = sigma, theta, gamma, epsilon

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
                dlogtau[np.ix_(buyers, cols)] = intensity_grid  # buyers' households
                for b in buyers:
                    rows = C + b * N + np.arange(N)  # buyer b's producers, all its own sectors
                    dlogtau[np.ix_(rows, cols)] = intensity_grid
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
        dOmega_temp2 = dOmega_total_tilde[C:C + CN, C:C + CN] * (1 / (1 - data['alpha']))[:, None]
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
                  sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05, data_dir='.'):
    """Convenience wrapper around `run()` for one-off scenario calls: runs the
    model and returns log GNE changes (dlogW, NOT percent -- multiply by 100
    yourself for a percentage) labeled by country code, the GNE-weighted
    world aggregate (same units), and the elasticity parameters that produced
    them, so a result is self-describing without needing to track down which
    run() call it came from. See `run()`'s docstring for the `shocks` format.

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
                           data_dir=data_dir)
    GNE_weights = data['chi_std'][:data['C']]
    world = GNE_weights @ dlogW_sum
    return dict(
        dlogW=dict(zip(countries, dlogW_sum)),
        World=world,
        elasticities=dict(sigma=sigma, theta=theta, gamma=gamma, epsilon=epsilon),
        ngrid=ngrid,
        shocks=shocks,
    )
