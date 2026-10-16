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
            'attribute_to' : 'seller' (default) or 'buyer' -- which side's
                          domestic technology absorbs the iceberg-cost
                          shock for GDP-accounting purposes (see `dlogY`
                          below). The theory paper models an iceberg cost
                          as a Hicks-neutral productivity shock to a
                          fictional specialized importer/exporter, and is
                          explicit that *which* country that fictional
                          producer belongs to is a free modeling choice
                          that does not affect welfare or real allocations
                          at all (confirmed below: `dlogW` is identical
                          regardless of this key) -- it only changes how
                          the SAME real event gets attributed across
                          countries' GDP. 'seller' matches standard FOB
                          export accounting (the exporter's measured
                          output absorbs the friction); 'buyer' instead
                          attributes it to the importer's absorption.
        Multiple legs are summed, so independent shocks (e.g. a Gulf-wide
        cutoff and a separate bilateral effect) can be combined in one run.
        The paper's own EU-vs-Russia scenario is
        `shocks=[{'sellers': [RUS], 'buyers': EU, 'sectors': None, 'intensity': 150}]`.
    ngrid : number of discretization steps for the shock.
    sigma, theta, gamma, epsilon : the model's four nested-CES elasticities
        of substitution, applied uniformly to every country/sector (not
        estimated per-sector like `trade_elast`). Defaults are the paper's
        own calibration, straight from its driver script's own comments
        (main_dlogW_rev_bigshocks_EU_Russian_v2.m lines 83-86) -- each
        governs a different tier of the nesting:
            sigma   (0.9)  : across household Consumption -- substitution
                             ACROSS SECTORS in what a household buys (food
                             vs. clothing, say). The WITHIN-a-sector,
                             across-origin-country dimension (German vs.
                             French wine) is instead governed by
                             `trade_elast`+1, per sector -- sigma is only
                             the cross-sector base rate households fall
                             back to (see nested_ces.py's sameSector_C).
            theta   (0.05) : across a producer's Composite Value-added and
                             Intermediates -- how easily a producer
                             substitutes its whole VA bundle for its whole
                             intermediate-input bundle, the top tier of the
                             producer's own nest.
            gamma   (0.5)  : across Primary Factors -- substitution within
                             the value-added bundle itself, across factor
                             types; irrelevant in the single-factor-per-
                             producer path this codebase uses
                             (factor_index==2, see main_load_data.py) since
                             there is nothing to substitute across when a
                             producer buys from only one factor, but kept
                             as a parameter for a future factor_index==1
                             caller.
            epsilon (0.05) : across a producer's Intermediate Inputs --
                             the producer-side analogue of sigma: cross-
                             SECTOR substitution in what a producer buys as
                             inputs (steel vs. plastic), with the within-
                             sector-across-origin-country dimension again
                             handled separately by `trade_elast`+1.
        All four are unrelated to `trade_elast` (a per-sector Armington
        elasticity from real data, see icio_to_haio.py/io_reorder.py) --
        these four are structural CES parameters the paper calibrates once
        and holds fixed across every sector and country, while `trade_elast`
        varies by sector and is the only one of the five actually estimated
        from data.

    Returns
    -------
    dlogW_sum : (C,) array, cumulative log change in real income (welfare)
        by country (last entry is the appended ROW country).
    dlogY_sum : (C,) array, cumulative log change in real GDP by country,
        via Theorem 1 of the theory paper under no initial wedges
        (Corollary 1): a Domar-weighted sum of the technology shocks
        attributed to each country's own producers, using each step's
        CURRENT sales shares (`lambda_CN`) to both weight and to derive
        each country's own GDP share (summing its producers' world-sales
        shares) -- GDP is genuinely a different object from GNE/welfare at
        the country level (only equal in aggregate, at the world level),
        and under this no-wedge scope a country's GDP only moves to the
        extent ITS OWN producers were directly handed a technology shock
        (via `attribute_to` above) -- a pure *welfare* loss from paying
        more for someone else's good, with nothing about your own
        producers' technology changing, does not move your own GDP to
        first order (see Corollary 1's "without wedges... foreign shocks
        do not affect real GDP").
    dlogy_sum : (CN,) array, cumulative log change in each individual
        producer's own real output/quantity (Corollary 3:
        `dlogy_i = dlambda_i/lambda_i - dlogp_i`, using the same per-step
        sales-share and price changes `response()` already computes for
        the converged step, just never previously returned). Producer i's
        position is `c*N + n` for country c (0-indexed within whichever
        `countries`/`keep_c` ordering the data uses) and sector n
        (0-indexed, matching `shocks[...]['sectors']`'s own indexing) --
        unlike `dlogY_sum`, this needs no `attribute_to` choice: it is
        unambiguous and unaffected by it.
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
    # Also precompute dlogA (CN,): the fixed per-step technology shock
    # attributed to whichever side's producers `attribute_to` names, for
    # the GDP calculation below (see Returns -- this is entirely separate
    # bookkeeping from dlogtau/dX above, which is unaffected by
    # `attribute_to` and drives the already-validated dlogW/equilibrium
    # machinery exactly as before).
    #
    # The naive version of this -- hitting the attributed country's WHOLE
    # sector with the full shock intensity -- massively overstates the GDP
    # effect: e.g. a Russia-to-Germany-only gas shock would otherwise imply
    # ALL of Russia's gas output (including what it sells elsewhere, or
    # uses domestically) lost productivity, not just the fraction actually
    # exported to Germany. So each attributed producer's dlogA is scaled by
    # the SHARE of its relevant trade flow that this leg's buyer/seller set
    # actually represents, using baseline (pre-shock) data -- a producer
    # whose sales are 10% to the shocked buyer set gets 10% of the per-step
    # intensity, not the full amount.
    chi_std_initial = data['chi_std'][:C].copy()
    lambda_CN_initial = data['lambda_CN'].copy()
    Omega_total_C_initial = data['Omega_total_C'].copy()
    Omega_total_tilde_initial = data['Omega_total_tilde'].copy()
    Omega_s_initial = data['Omega_s'].copy()

    legs = []
    dlogA = np.zeros(CN)
    for leg in shocks:
        sellers = np.asarray(leg['sellers'], dtype=int) - 1
        buyers = leg.get('buyers')
        buyers = np.arange(C) if buyers is None else np.asarray(buyers, dtype=int) - 1
        sectors = leg.get('sectors')
        sectors = np.arange(N) if sectors is None else np.asarray(sectors, dtype=int)
        intensity_grid = np.log(1 + leg['intensity'] / 100) / ngrid
        legs.append((sellers, buyers, sectors, intensity_grid))

        buyer_producer_rows = C + (buyers[:, None] * N + np.arange(N)[None, :]).ravel()

        attribute_to = leg.get('attribute_to', 'seller')
        if attribute_to == 'seller':
            # For each shocked (seller, sector): what share of that
            # producer's total sales (lambda_CN_initial) goes to the
            # shocked buyer set -- summing household purchases
            # (Omega_total_C) and every shocked buyer's own producers'
            # purchases (Omega_total_tilde), each weighted by the buying
            # party's own size (chi_std/lambda_CN), since Omega is a
            # share, not a level.
            for s in sellers:
                cols = s * N + sectors
                total_sales = lambda_CN_initial[C + cols]
                hh = (chi_std_initial[buyers, None] * Omega_total_C_initial[buyers][:, cols]).sum(axis=0)
                prod = (lambda_CN_initial[buyer_producer_rows, None]
                        * Omega_total_tilde_initial[buyer_producer_rows][:, C + cols]).sum(axis=0)
                share = np.divide(hh + prod, total_sales, out=np.zeros_like(total_sales), where=total_sales != 0)
                dlogA[cols] -= intensity_grid * share
        elif attribute_to == 'buyer':
            # Symmetric logic from the buyer's side: what share of each
            # shocked buyer-producer's total spending on this SECTOR
            # (Omega_s, summed across every origin country) is specifically
            # these shocked sellers -- i.e. how reliant is this buyer's
            # sector-n input bundle on these particular sellers, versus
            # other origins (including domestic) for the same good. Unlike
            # the seller case, no household term: GDP (Theorem 1) sums over
            # producers, and the "specialized importer" is itself a
            # producer -- household purchases don't have one.
            for b in buyers:
                prod_rows = C + b * N + np.arange(N)  # L-indexed: for Omega_total_tilde's rows
                cn_rows = b * N + sectors  # CN-indexed (no +C): for Omega_s and dlogA
                sector_total_prod = Omega_s_initial[cn_rows, sectors]
                from_sellers_prod = Omega_total_tilde_initial[prod_rows[sectors][:, None],
                                                                C + (sellers[None, :] * N + sectors[:, None])].sum(axis=1)
                share_prod = np.divide(from_sellers_prod, sector_total_prod,
                                        out=np.zeros_like(sector_total_prod), where=sector_total_prod != 0)
                dlogA[cn_rows] -= intensity_grid * share_prod
        else:
            raise ValueError(f"attribute_to must be 'seller' or 'buyer', got {attribute_to!r}")

    dlogW = np.zeros((C, ngrid))
    dlogY = np.zeros((C, ngrid))
    dlogy = np.zeros((CN, ngrid))

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
        # from the (always-zero) producer/factor block. The divide below is
        # itself still guarded (not just sliced to :C) in case some tracked
        # "country" genuinely has zero GNE weight (e.g. an empty ROW
        # composite for a `countries` list that already covers everyone).
        chi_std_C = data['chi_std'][:C]
        dlogW[:, i] = np.divide(dchi_std[:C], chi_std_C, out=np.zeros(C), where=chi_std_C != 0) - dlogP_Vec[:C]

        # --- GDP (Theorem 1 / Corollary 1, no initial wedges) and producer
        # output (Corollary 3), using THIS step's current sales shares.
        # A country's GDP share is derived on the fly as the sum of its own
        # producers' world sales shares (lambda_CN's producer block), the
        # same Domar-weight object Hulten's theorem always uses -- not
        # drawn from the (GNE-based, generally different) GDP_weights/
        # GNE_weights_initial baseline.
        lambda_producers = data['lambda_CN'][C:]  # (CN,), current step
        GDP_share_c = lambda_producers.reshape(C, N).sum(axis=1)  # (C,)
        GDP_share_safe = np.where(GDP_share_c == 0, 1, GDP_share_c)
        lambdaY = lambda_producers / np.repeat(GDP_share_safe, N)
        dlogY[:, i] = (lambdaY * dlogA).reshape(C, N).sum(axis=1)

        dlambda_producers = dlambda_result[C:C + CN]
        dlogP_producers = dlogP_Vec[C:C + CN]
        dlogy[:, i] = np.divide(dlambda_producers, lambda_producers,
                                 out=np.zeros(CN), where=lambda_producers != 0) - dlogP_producers

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
    dlogY_sum = dlogY.sum(axis=1)
    dlogy_sum = dlogy.sum(axis=1)
    return dlogW_sum, dlogY_sum, dlogy_sum, data


def run_scenario(keep_c, countries, shocks, ngrid=20,
                  sigma=0.9, theta=0.05, gamma=0.5, epsilon=0.05, data_dir='.', haio=None):
    """Convenience wrapper around `run()` for one-off scenario calls: runs the
    model and returns log GNE changes (dlogW, NOT percent -- multiply by 100
    yourself for a percentage) labeled by country code, the GNE-weighted
    world aggregate (same units), and the elasticity parameters that produced
    them, so a result is self-describing without needing to track down which
    run() call it came from. See `run()`'s docstring for the `shocks` format
    (including `attribute_to`, which only affects `dlogY`/`World_Y` below,
    never `dlogW`/`World`) and for what passing `haio` directly (a non-WIOD
    dataset) does.

    Returns
    -------
    dict with:
        'dlogW' : {country code -> log GNE/welfare change}, one entry per `countries`
        'World' : GNE-weighted world aggregate log GNE change
        'dlogY' : {country code -> log real GDP change} -- a genuinely
            different object from `dlogW` at the country level (only equal
            in aggregate); see `run()`'s docstring for exactly what it
            captures and why it depends on `attribute_to`.
        'World_Y' : GDP-weighted world aggregate log GDP change (uses the
            same baseline weights as `World`, since this codebase's HAIO
            contract provides a single baseline GDP/GNE share array -- see
            main_load_data.py's module docstring)
        'dlogy' : {country code -> list of per-sector log real output
            changes}, 0-indexed in whatever sector order the data uses
            (same convention as `shocks[...]['sectors']`) -- Corollary 3's
            producer-level quantity change, unrelated to `attribute_to`.
        'elasticities' : {'sigma', 'theta', 'gamma', 'epsilon'} used for this run
        'ngrid', 'shocks' : shock discretization settings used
    """
    dlogW_sum, dlogY_sum, dlogy_sum, data = run(
        keep_c, shocks, ngrid=ngrid, sigma=sigma, theta=theta, gamma=gamma, epsilon=epsilon,
        data_dir=data_dir, haio=haio)
    world_weights = data['GNE_weights_initial']
    world = world_weights @ dlogW_sum
    world_y = world_weights @ dlogY_sum
    N = data['N']
    return dict(
        dlogW=dict(zip(countries, dlogW_sum)),
        World=world,
        dlogY=dict(zip(countries, dlogY_sum)),
        World_Y=world_y,
        dlogy={c: dlogy_sum[i * N:(i + 1) * N].tolist() for i, c in enumerate(countries)},
        elasticities=dict(sigma=sigma, theta=theta, gamma=gamma, epsilon=epsilon),
        ngrid=ngrid,
        shocks=shocks,
    )
