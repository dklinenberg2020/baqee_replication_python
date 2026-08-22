"""
Python translation of replication/baqaee_farhi_model/main_load_data_rev.m
(the `initial_tariff_index == 1` path only -- the paper's main driver script
never uses the with-tariff path, so IO_reorder_init_tariff.m is not ported;
see baqaee_farhi_model/README.md).

Turns the country x sector input-output table from io_reorder() into the
"standard form" world economy used by the linearized nested-CES model: a
single (L, L) Markov-like matrix Omega_total_tilde stacking household
consumption, intermediate/factor purchases by producers, and a zero block
for (endowed) factors, where L = C + C*N + C*F, plus its Leontief inverse
Psi_total and the implied sales/income shares chi_std, lambda_std.

Validated element-wise against the original MATLAB code (GNU Octave 8.4.0)
for the full 41-country, country-sector-specific-factor run -- see
baqaee_farhi_model/README.md.
"""
import numpy as np
from io_reorder import io_reorder


def main_load_data(keep_c, initial_tariff_index, factor_index, data_dir='.'):
    if initial_tariff_index != 1:
        raise NotImplementedError(
            'initial_tariff_index == 2 (with initial tariffs) is not ported; '
            'the paper\'s main driver script only uses initial_tariff_index == 1.')

    Omega, beta, alpha_VA, alpha, trade_elast, GDP_weights = io_reorder(keep_c, data_dir=data_dir)
    Omega_tilde = Omega
    GNE_weights = GDP_weights

    # MATLAB: alpha_VA(:,3) <-> alpha_VA(:,4) swap (medium-skill/high-skill columns)
    alpha_VA = alpha_VA.copy()
    alpha_VA[:, [2, 3]] = alpha_VA[:, [3, 2]]

    C = len(keep_c) + 1  # kept countries + aggregate ROW
    N = 30
    F = N if factor_index == 2 else 4
    CN = C * N
    CF = C * F

    init_t = np.ones((C + CN, CN + CF))

    Omega_tilde = np.nan_to_num(Omega_tilde[:CN, :CN], posinf=1.0, neginf=1.0)
    beta = np.nan_to_num(beta[:CN, :C], posinf=1.0, neginf=1.0)
    alpha = np.nan_to_num(alpha[:CN], posinf=1.0, neginf=1.0)
    alpha_VA = np.nan_to_num(alpha_VA[:CN, :], posinf=1.0, neginf=1.0)
    GNE_weights = GNE_weights[:C]

    # --- Transform raw shares into "who buys from whom" disaggregated form --
    beta_s = np.zeros((C, N))
    Omega_s = np.zeros((CN, N))
    for i in range(C):
        beta_s += beta[i * N:(i + 1) * N, :].T
        Omega_s += Omega_tilde[:, i * N:(i + 1) * N]

    # beta_disagg[c, n] = share of household c's spending on sector-n good
    # that goes to each origin country (length-C vector)
    beta_s_safe = np.where(beta_s == 0, 1, beta_s)
    beta_disagg_temp = beta.T / np.tile(beta_s_safe, (1, C))  # (C, CN)
    beta_disagg_temp = beta_disagg_temp.reshape(C, C, N)  # [c, origin, n]

    Omega_s_safe = np.where(Omega_s == 0, 1, Omega_s)
    Omega_disagg_temp = Omega_tilde / np.tile(Omega_s_safe, (1, C))  # (CN, CN)
    Omega_disagg_temp = Omega_disagg_temp.reshape(CN, C, N)  # [buyer, origin, n]

    # Omega_total_C[c, i]: how much household c spends on sector-country i
    Omega_total_C = np.zeros((C, CN))
    for i in range(CN):
        ic, in_ = divmod(i, N)
        Omega_total_C[:, i] = beta_s[:, in_] * beta_disagg_temp[:, ic, in_]

    # Omega_total_N[i, j]: how much producer i spends on sector-country /
    # factor j
    ic_of = np.arange(CN) // N
    Omega_total_N = np.zeros((CN, CN + CF))
    for j in range(CN):
        jc, jn = divmod(j, N)
        Omega_total_N[:, j] = (1 - alpha) * Omega_s[:, jn] * Omega_disagg_temp[:, jc, jn]
    if factor_index == 1:
        for i in range(CN):
            ic = ic_of[i]
            Omega_total_N[i, CN + ic * F: CN + (ic + 1) * F] = alpha[i] * alpha_VA[i, :]
    else:  # factor_index == 2: country-sector-specific labor, one factor per producer
        Omega_total_N[np.arange(CN), CN + np.arange(CN)] = alpha

    # --- Standard-form Omega_total_tilde / Omega_total (L, L) ---------------
    L = C + CN + CF
    Omega_total_tilde = np.zeros((L, L))
    Omega_total_tilde[:C, C:C + CN] = Omega_total_C
    Omega_total_tilde[C:C + CN, C:] = Omega_total_N
    # last CF rows (factors) are endowed: stay zero

    Omega_total = Omega_total_tilde  # initial_tariff_index == 1

    Psi_total_tilde = np.linalg.solve(np.eye(L) - Omega_total_tilde, np.eye(L))
    Psi_total = Psi_total_tilde  # same matrix when there is no initial tariff wedge

    chi_std = np.concatenate([GNE_weights, np.zeros(CN + CF)])
    lambda_std = chi_std @ Psi_total

    lambda_CN = lambda_std[:C + CN]
    lambda_F = lambda_std[C + CN:]

    data = dict(
        C=C, N=N, F=F, CN=CN, CF=CF, L=L,
        trade_elast=trade_elast,
        alpha=alpha,
        beta_s=beta_s,
        Omega_s=Omega_s,
        Omega_total_C=Omega_total_C,
        Omega_total_N=Omega_total_N,
        Omega_total_tilde=Omega_total_tilde,
        Omega_total=Omega_total,
        Psi_total_tilde=Psi_total_tilde,
        Psi_total=Psi_total,
        lambda_std=lambda_std,
        lambda_CN=lambda_CN,
        lambda_F=lambda_F,
        chi_std=chi_std,
    )
    shock = dict(init_t=init_t)
    return data, shock
