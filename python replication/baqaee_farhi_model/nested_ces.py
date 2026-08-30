"""
Python translation of the linearized nested-CES network model of Baqaee and
Farhi (2024), i.e. replication/baqaee_farhi_model/AES_func.m,
Nested_CES_linear_final_rev.m and Nested_CES_linear_result_final.m.

IMPORTANT ALGORITHMIC NOTE (read before modifying):
The MATLAB code materializes full 3-D Allen-elasticity-of-substitution (AES)
tensors: AES_N_Mat has shape (CN, CN+CF, CN) and AES_F_Mat has shape
(CF, CN+CF, CN). At the paper's actual scale (C=41 countries,
country-sector-specific factors so CN=CF=1230), AES_N_Mat alone is
1230*2460*1230*8 bytes =~ 30 GB (confirmed: the original MATLAB code runs out
of memory under GNU Octave 8.4.0 with 15GB RAM available). This module
instead computes the *effect* of those elasticities directly in closed form,
which is mathematically identical -- AES_func.m's elasticities only ever
take one of three values per (row, column) pair (a same-origin-and-sector
"diagonal" value, a same-sector-different-origin value, and a residual base
value), so the row sums needed downstream can be written directly in terms
of those three values without ever forming the (CN, CN+CF, CN) array. This
is NOT a different model, just a much cheaper way of evaluating the same
formulas -- validated element-wise against AES_func.m's own output and
against Nested_CES_linear_result_final.m's full output on a 5-country test
case, see baqaee_farhi_model/README.md.

Additionally, `factor_index == 2` (country-sector-specific factors, what the
paper's main driver uses) means every producer buys value-added from exactly
one factor -- its own (see main_load_data.py) -- so the factor-substitution
block of the AES machinery collapses to a single "own factor" term per
producer instead of the full (CF, CN+CF) block AES_F_Mat carries per
producer in the original code.

For the same memory reason, solving the implicit linear system for the
equilibrium response (what Nested_CES_linear_final_rev.m does via an
explicit, heavily reshaped/permuted Jacobian) is instead done here by
*probing*: the map g(x) = response(x)'s [factor-income-share ;
redistribution] block is provably affine in x (every step from x to that
output is a linear operation on fixed data), so its Jacobian A and offset B
can be recovered exactly as B = g(0) and column i of A = g(e_i) - g(0). This
avoids hand-translating Nested_CES_linear_final_rev.m's reshape/permute/kron
machinery (very high risk of an unnoticed row/column-major bug) at the cost
of C+CF evaluations of `response()` instead of one -- affordable because
`response()` no longer needs the huge dense tensors either. See
`solve_dlambda_F_all` below.

SCOPE: `factor_index == 2`, `initial_tariff_index == 1` (no initial
tariffs) only -- what the paper's main driver script uses. It further
assumes `dlogt` (tariff shocks) is identically zero, true for the
iceberg-cost shocks (shock_index 1 and 3) the paper actually uses: under no
initial tariffs and no dlogt, the row-1 regularization in
Nested_CES_linear_final_rev.m reduces to the simple global adding-up
constraint sum(dlambda_F) + sum(dlambda_F_star) = 0, which is imposed
directly instead of replicated via its G/H bookkeeping. Both
simplifications are verified against the MATLAB output for the exact shock
the paper uses (EU iceberg cost on Russia); see baqaee_farhi_model/README.md.
"""
import numpy as np


def value_added_shares(data):
    """Closed-form AES building blocks that depend only on `data` (not on
    the shock), reused across every `response()` call within one iteration
    of the outer discretization loop."""
    sigma, theta, epsilon = data['sigma'], data['theta'], data['epsilon']
    alpha = data['alpha']
    beta_s = data['beta_s']
    Omega_s = data['Omega_s']
    trade_elast = data['trade_elast']

    one_minus_alpha = 1 - alpha
    # alpha==1.0 (all value-added, zero measured intermediate cost) is not
    # only sector T -- any real producer with no measured intermediate
    # purchases hits it too (ICIO's finer sector grid makes this more likely
    # than WIOD's). Guard the same way Omega_s_safe/beta_s_safe below
    # already guard their own zero-denominator case: a producer with
    # one_minus_alpha==0 has zero weight everywhere this feeds into
    # (Omega_total_N's off-diagonal terms are separately scaled by
    # (1-alpha)==0), so the filler value used here for the division itself
    # is inert -- it only exists to avoid inf/NaN poisoning the matmuls
    # downstream via Psi_total, exactly like the lambda_F fix in response().
    one_minus_alpha_safe = np.where(one_minus_alpha == 0, 1, one_minus_alpha)
    base_kc = epsilon / one_minus_alpha_safe + theta * (1 - 1 / one_minus_alpha_safe)  # (CN,)

    Omega_s_safe = np.where(Omega_s == 0, 1, Omega_s)  # (CN, N)
    sameSector_minus_base = ((trade_elast[None, :] + 1 - epsilon)
                              / (one_minus_alpha_safe[:, None] * Omega_s_safe))  # (CN, N)

    beta_s_safe = np.where(beta_s == 0, 1, beta_s)  # (C, N)
    sameSector_C = (trade_elast[None, :] + 1) / beta_s_safe + sigma * (1 - 1 / beta_s_safe)  # (C, N)

    return dict(base_kc=base_kc, one_minus_alpha=one_minus_alpha,
                sameSector_minus_base=sameSector_minus_base, sameSector_C=sameSector_C)


def _diag_term(Om_denom, Om_weight, dlogP, trade_elast_of_col):
    """(AES_diag(row,col) - AES_sameSector(row,col)) * Om_weight(row,col) * dlogP(row,col).
    AES_diag and AES_sameSector are identical except for AES_diag's extra
    -(trade_elast+1)/Om_denom(row,col) term (every other term is shared and
    cancels exactly), so the difference is just -(trade_elast+1)/Om_denom,
    guarded against 0/0 where Om_denom(row,col) == 0 -- there the whole term
    is exactly zero, matching AES_func.m's explicit `if Om==0, AES=1` branch
    (which only ever gets multiplied by that same Om==0 term).

    Om_denom (Omega_total_N / Omega_total_C, the AES formula's own
    denominator per AES_func.m) and Om_weight (the weight actually
    multiplying (AES-1) in the sum, which Nested_CES_linear_result_final.m
    always takes from Omega_total_tilde) are the SAME array initially, but
    diverge after the outer discretization loop's first iteration: the two
    are updated by different, non-interchangeable formulas -- see
    run_model.py's "Update variables" block and AES_func.m vs.
    Nested_CES_linear_result_final.m's choice of which Omega_total_* they
    read. Both share the same zero/nonzero support throughout (a cell that
    starts at exactly zero has an exactly-zero update in both, by induction
    -- see nested_ces.py's module docstring), so the 0/0 guard is safe to
    base on either; only the nonzero *values* must come from the right one."""
    nz = Om_denom != 0
    Om_safe = np.where(nz, Om_denom, 1.0)
    corr = np.where(nz, -(trade_elast_of_col + 1) / Om_safe, 0.0)
    return corr * Om_weight * dlogP


def producer_dOmega_goods(data, blocks, Og, Pg, own_factor_Om_weight, own_factor_P):
    """dOmega_CN_tilde[kc, i] for i=1..CN, the intermediate-goods block of
    every producer kc's cost function (Nested_CES_linear_result_final.m's
    `temp1`). Og, Pg are (CN, CN): the WEIGHT (Omega_total_tilde's producer
    goods-goods block -- see _diag_term's docstring for why this is not
    Omega_total_N once the outer loop has run) and the bilateral price
    change producer kc faces buying good i.

    `temp1(i)`'s underlying sum in the MATLAB code actually runs over ALL
    CN+CF columns (goods AND factors) of producer kc's cost function, not
    just the CN goods -- see AES_N_Mat(:,CN+1:end,kc) = theta. Since
    Omega(kc, factor j) is zero except at j = kc's own factor (factor_index
    == 2, see main_load_data.py), that adds one extra CONSTANT term (same
    for every candidate good i) to `default_term`:
    own_factor_Om_weight(kc)*(theta-1)*dlogP_own_factor(kc)."""
    C, N, CN = data['C'], data['N'], data['CN']
    theta = data['theta']
    trade_elast = data['trade_elast']
    Omega_total_N = data['Omega_total_N']
    base_kc = blocks['base_kc']
    sameSector_minus_base = blocks['sameSector_minus_base']

    w = Og * Pg
    S = w.sum(axis=1)  # (CN,): goods-only weighted price sum
    S_sector = w.reshape(CN, C, N).sum(axis=1)  # (CN, N)

    sector_of_col = np.tile(np.arange(N), C)  # sector(i), i=0..CN-1
    own_factor_term = own_factor_Om_weight * (theta - 1) * own_factor_P
    default_term = (base_kc - 1) * S + own_factor_term
    sector_term = (sameSector_minus_base * S_sector)[:, sector_of_col]

    Om_denom = Omega_total_N[:, :CN]
    te_col = trade_elast[sector_of_col][None, :]
    diag = _diag_term(Om_denom, Og, Pg, te_col)  # one power of Og (=Om_weight)

    return Og * Pg + Og * (default_term[:, None] + sector_term) + Og * diag, S  # outer Og gives the second power


def producer_dOmega_own_factor(data, S, own_factor_Om_denom, own_factor_Om_weight, own_factor_P):
    """dOmega_CN_tilde[kc, CN+kc] -- change in producer kc's expenditure
    share on its own country-sector-specific factor (the only factor it
    ever buys when factor_index == 2). `S` is the goods-weighted price sum
    from producer_dOmega_goods (reused, not recomputed).

    AES_F_Mat's own-factor diagonal value is
    -gamma/Omega_total_N(kc,CN+kc) + gamma/alpha(kc) + theta*(1-1/alpha(kc))
    (own_factor_Om_denom = Omega_total_N(kc,CN+kc), never Omega_total_tilde
    -- see AES_func.m). Initially Omega_total_N(kc,CN+kc) == alpha(kc)
    exactly (see main_load_data.py), so the gamma terms cancel to
    theta*(1-1/alpha(kc)) -- but `alpha` is a fixed calibration parameter
    that is never updated across the outer discretization loop's iterations
    (run_model.py), while Omega_total_N's own-factor share *does* evolve
    with the shock, so after the first iteration they are no longer equal
    and the full formula must be used. Separately, the WEIGHT multiplying
    (AES-1) in the sum (own_factor_Om_weight, in the Om*P terms below) is
    always Omega_total_tilde's own-factor diagonal, per
    Nested_CES_linear_result_final.m -- see _diag_term's docstring."""
    theta, gamma = data['theta'], data['gamma']
    alpha = data['alpha']

    Om_nz = own_factor_Om_denom != 0
    Om_safe = np.where(Om_nz, own_factor_Om_denom, 1.0)
    alpha_nz = alpha != 0
    alpha_safe = np.where(alpha_nz, alpha, 1.0)
    own_factor_AES = np.where(Om_nz, -gamma / Om_safe + gamma / alpha_safe + theta * (1 - 1 / alpha_safe), 1.0)
    own_factor_AES_minus_1 = own_factor_AES - 1

    temp2_own = (theta - 1) * S + own_factor_AES_minus_1 * own_factor_Om_weight * own_factor_P
    return own_factor_Om_weight * own_factor_P + own_factor_Om_weight * temp2_own


def consumer_dOmega(data, blocks, Og, Pg):
    """dOmega_CN_tilde[c, i] for i=1..CN -- change in household c's
    expenditure share on sector-country good i, mirroring
    producer_dOmega_goods but using AES_C_Mat (sigma-based)."""
    C, N = data['C'], data['N']
    sigma = data['sigma']
    trade_elast = data['trade_elast']
    Omega_total_C = data['Omega_total_C']
    sameSector_C = blocks['sameSector_C']

    w = Og * Pg
    S = w.sum(axis=1)  # (C,)
    S_sector = w.reshape(C, C, N).sum(axis=1)  # (C, N)

    sector_of_col = np.tile(np.arange(N), C)
    default_term = (sigma - 1) * S
    sector_term = ((sameSector_C - sigma) * S_sector)[:, sector_of_col]

    te_col = trade_elast[sector_of_col][None, :]
    diag = _diag_term(Omega_total_C, Og, Pg, te_col)

    return Og * Pg + Og * (default_term[:, None] + sector_term) + Og * diag


def response(x, data, shock, blocks):
    """Translation of Nested_CES_linear_result_final.m, restricted to the
    factor_index==2 / no-initial-tariff / dlogt==0 scope documented in the
    module docstring. Returns the [dlambda_F ; dlambda_F_star]-shaped vector
    g(x) whose fixed point is the equilibrium response to the shock, plus
    dlogP_Vec and dOmega_total for downstream use (welfare, reallocation)."""
    C, CN, CF = data['C'], data['CN'], data['CF']
    lambda_F = data['lambda_F']
    lambda_CN = data['lambda_CN']
    Psi_total = data['Psi_total']
    Phi_F = data['Phi_F']

    dlogtau = shock['dlogtau']
    dX = shock['dX']

    dlambda_F = x[:CF]
    dlambda_F_star = x[CF:]

    # A producer/factor with exactly zero baseline income (real for some
    # ICIO country-sector pairs at fine sector resolution, e.g. Luxembourg's
    # zero mining/oil output -- never arose in WIOD's coarser 30 sectors) has
    # dlambda_F == lambda_F == 0 here (its share of a shrinking/growing
    # world total starts and stays at 0 in this linearization). Guard the
    # divide as zero rather than NaN: the alternative propagates NaN through
    # every country's dlogP_CN below via the dense Psi_total matmul.
    dlogP_F = np.divide(dlambda_F, lambda_F, out=np.zeros_like(dlambda_F), where=lambda_F != 0)
    dlogP_CN = Psi_total[:C + CN, :C + CN] @ dX + Psi_total[:C + CN, C + CN:] @ dlogP_F  # (C+CN,)
    dlogP_Vec = np.concatenate([dlogP_CN, dlogP_F])

    dlogP_goods = dlogP_CN[C:]  # seller's own price, seller j -> position C+j in dlogP_CN
    dlogtau_goods = dlogtau[:, :CN]  # (C+CN, CN); dlogt assumed 0 (see docstring)
    dlogP_goods_mat = dlogtau_goods + dlogP_goods[None, :]  # (C+CN, CN)

    dlogtau_owndiag = np.diagonal(dlogtau[C:, CN:CN + CN]).copy()
    dlogP_own_factor = dlogtau_owndiag + dlogP_F  # (CN,)

    # Weights multiplying (AES-1) in every sum below always come from
    # Omega_total_tilde (== Omega_total here, no initial tariffs), NOT from
    # the separately-bookkept Omega_total_C/Omega_total_N -- see
    # _diag_term's docstring for why these differ after the outer
    # discretization loop's first iteration.
    Omega_total_tilde = data['Omega_total_tilde']
    Omega_total_N = data['Omega_total_N']  # only used below for the AES-formula denominator

    Og_C_weight = Omega_total_tilde[:C, C:C + CN]
    dOmega_C = consumer_dOmega(data, blocks, Og_C_weight, dlogP_goods_mat[:C, :])  # (C, CN)

    own_factor_Om_denom = Omega_total_N[np.arange(CN), CN + np.arange(CN)]
    own_factor_Om_weight = Omega_total_tilde[C + np.arange(CN), C + CN + np.arange(CN)]

    Og_N_weight = Omega_total_tilde[C:C + CN, C:C + CN]
    Pg_N = dlogP_goods_mat[C:, :]
    dOmega_N_goods, S = producer_dOmega_goods(data, blocks, Og_N_weight, Pg_N, own_factor_Om_weight, dlogP_own_factor)  # (CN, CN)

    dOmega_N_ownfactor = producer_dOmega_own_factor(data, S, own_factor_Om_denom, own_factor_Om_weight, dlogP_own_factor)  # (CN,)

    # --- Equation (2): dOmega_tilde -> dOmega (no initial tariffs => equal) --
    dOmega_C_total = dOmega_C
    dOmega_N_goods_total = dOmega_N_goods
    dOmega_N_ownfactor_total = dOmega_N_ownfactor

    # --- Equation (5): dchi = Phi_F * dlambda_F + dlambda_F_star (no tariff
    # revenue redistribution terms since init_t == 1 everywhere) -----------
    dchi = Phi_F @ dlambda_F + dlambda_F_star  # (C,)

    # --- dlambda_result = (dchi_std' Psi_total + chi_std' dPsi_total)' ------
    # chi_std' dPsi_total = chi_std' Psi_total dOmega_total Psi_total
    #                     = lambda_std' dOmega_total Psi_total   (avoids
    # forming the (L,L) matrix dPsi_total = Psi_total dOmega_total Psi_total
    # explicitly, see module docstring).
    dOmega_total_row_C = dOmega_C_total  # (C, CN): rows 0..C-1 of dOmega_total's (C+CN, CN) goods block
    dOmega_total_row_N = dOmega_N_goods_total  # (CN, CN)

    # lambda_std' @ dOmega_total : only the (C+CN, CN) goods block and the
    # (CN,) own-factor diagonal of dOmega_total are ever nonzero.
    v = np.zeros(data['L'])
    # households contribute lambda_CN[0:C] against dOmega_C rows
    v[C:C + CN] += lambda_CN[:C] @ dOmega_total_row_C  # -> goods columns (C..C+CN)
    v[C:C + CN] += lambda_CN[C:C + CN] @ dOmega_total_row_N  # producers' goods purchases
    v[C + CN:] += lambda_CN[C:C + CN] * dOmega_N_ownfactor_total  # producers' own-factor purchases (diagonal)

    dchi_std = np.concatenate([dchi, np.zeros(CN + CF)])
    dlambda_result = dchi_std @ Psi_total + v @ Psi_total

    # dlambda_F_star_result = dlambda_CN_result' * sum_j Phi_T .* Omega_total .* (init_t-1)
    # (init_t == 1 everywhere in the no-initial-tariff scope) => identically zero
    dlambda_F_star_result = np.zeros(C)

    g = np.concatenate([dlambda_result[C + CN:], dlambda_F_star_result])
    return g, dlogP_Vec, dOmega_total_row_N, dOmega_N_ownfactor_total, dOmega_C_total, dchi_std, dlambda_result


def solve_dlambda_F_all(data, shock, blocks):
    """Solve the fixed point x = g(x) by probing (see module docstring):
    since g is affine, B = g(0) and A[:, i] = g(e_i) - B."""
    C, CF = data['C'], data['CF']
    n = CF + C
    B, *_ = response(np.zeros(n), data, shock, blocks)
    A = np.zeros((n, n))
    for i in range(n):
        e = np.zeros(n)
        e[i] = 1.0
        gi, *_ = response(e, data, shock, blocks)
        A[:, i] = gi - B

    # Row-1 regularization: under no initial tariffs and dlogt==0, the
    # redundant first equation is replaced by the global adding-up
    # constraint sum(dlambda_F) + sum(dlambda_F_star) = 0 (verified exactly
    # against Octave's A[0,:], B[0] output -- see baqaee_farhi_model/README.md).
    A[0, :] = 0.0
    A[0, 1:] = -1.0
    B[0] = 0.0

    x = np.linalg.solve(np.eye(n) - A, B)
    return x
