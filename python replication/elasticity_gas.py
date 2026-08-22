"""
Python translation of replication/elasticity_gas.m.

Same aggregate CES(sigma) model as elasticity.py, but calibrated to the
"30% gas drop" scenario of Table 2, column (4): initial gas expenditure
share alpha = 0.01 of GNE, subjected to a 30% reduction in gas quantity.

Validated against the original MATLAB output (GNU Octave 8.4.0):
    Elasticity = 0 (Leontief):    output loss -29.96%
    Elasticity = 0.1:              output loss  -2.33%  (Table 2, col 4: 2.3% GNE / 2.2% GDP)
    Elasticity = 0.2:              output loss  -0.77%
    Elasticity = 1 (Cobb-Douglas): output loss  -0.36%
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from ces_energy import ces_quantities, ces_cobb_douglas, ces_leontief, nearest_index, print_row

pE, pX = 1.0, 1.0
alpha = 0.01  # GNE share of gas

E_drop = 0.3
E_after_drop = 1 - E_drop

X = (1 - alpha) / alpha

Emin, Emax, N = E_after_drop - 0.05, 1.05, 501
E = np.linspace(Emin, Emax, N)

sigmas = [0.1, 0.2, 0.999]
Y = np.zeros((N, len(sigmas)))
Ybench = np.zeros(len(sigmas))
MPE = np.zeros((N, len(sigmas)))
MPX = np.zeros((N, len(sigmas)))
exp_share = np.zeros((N, len(sigmas)))

for i, sigma in enumerate(sigmas):
    Y[:, i], Ybench[i], MPE[:, i], MPX[:, i], exp_share[:, i] = ces_quantities(E, alpha, sigma, pX)

Y_CD, Ybench_CD, MPE_CD, MPX_CD, exp_share_CD = ces_cobb_douglas(E, alpha, X)
Y_Leontief, Ybench_Leontief, MPE_Leontief, MPX_Leontief, exp_share_Leontief = ces_leontief(E, alpha, X)

index_bench = nearest_index(E, 1.0)

# --- Figure 1: output ---------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(E, Y_Leontief / Ybench_Leontief, '--', lw=2, label=r'$\sigma=0$ (Leontief)')
ax.plot(E, Y[:, 0] / Ybench[0], '-.', lw=2, label=r'$\sigma=0.1$')
ax.plot(E, Y[:, 1] / Ybench[1], '-.', lw=2, label=r'$\sigma=0.16$')
ax.plot(E, Y_CD / Ybench_CD, lw=2, label=r'$\sigma=1$ (Cobb-Douglas)')
ax.axvline(E_after_drop, color='k', lw=1.5)
ax.set_xlim(Emin, Emax)
ax.set_ylim(Emin, 1.02)
ax.set_xlabel('Energy, $E$')
ax.set_ylabel('Production, $Y$')
ax.legend(loc='lower right')
ax.grid(True)
fig.savefig('figures/elasticity_fig_gas.pdf')
fig.savefig('figures/elasticity_fig_gas.eps')
plt.close(fig)

# --- Figure 2: price of energy relative to baseline ----------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(E, MPE_Leontief / MPE_Leontief[index_bench], '--', lw=2, label=r'$\sigma=0$ (Leontief)')
ax.plot(E, MPE[:, 0] / MPE[index_bench, 0], '-.', lw=2, label=r'$\sigma=0.1$')
ax.plot(E, MPE[:, 1] / MPE[index_bench, 0], '-.', lw=2, label=r'$\sigma=0.16$')
ax.plot(E, MPE_CD / MPE_CD[index_bench], lw=2, label=r'$\sigma=1$ (Cobb-Douglas)')
ax.axvline(E_after_drop, color='k', lw=1.5)
ax.set_xlim(Emin, Emax)
ax.set_ylim(0.5, 10)
ax.set_xlabel('Energy, $E$')
ax.set_ylabel('Price $p_E$ ($=MPE$) rel. to baseline')
ax.legend(loc='upper right')
ax.grid(True)
fig.savefig('figures/MPE_fig_gas.pdf')
fig.savefig('figures/MPE_fig_gas.eps')
plt.close(fig)

# --- Figure 3: price of other inputs relative to baseline -----------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(E, MPX_Leontief / MPX_Leontief[index_bench], '--', lw=2, label=r'$\sigma=0$ (Leontief)')
ax.plot(E, MPX[:, 0] / MPX[index_bench, 0], '-.', lw=2, label=r'$\sigma=0.1$')
ax.plot(E, MPX[:, 1] / MPX[index_bench, 0], '-.', lw=2, label=r'$\sigma=0.16$')
ax.plot(E, MPX_CD / MPX_CD[index_bench], lw=2, label=r'$\sigma=1$ (Cobb-Douglas)')
ax.axvline(E_after_drop, color='k', lw=1.5)
ax.set_xlim(Emin, Emax)
ax.set_ylim(0, 1.05)
ax.set_xlabel('Energy, $E$')
ax.set_ylabel('Price $p_X$ ($=MPX$) rel. to baseline')
ax.legend(loc='lower right')
ax.grid(True)
fig.savefig('figures/MPX_fig_gas.pdf')
fig.savefig('figures/MPX_fig_gas.eps')
plt.close(fig)

# --- Figure 4: energy expenditure share -----------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(E, exp_share_Leontief, '--', lw=2, label=r'$\sigma=0$ (Leontief)')
ax.plot(E, exp_share[:, 0], '-.', lw=2, label=r'$\sigma=0.1$')
ax.plot(E, exp_share[:, 1], '-.', lw=2, label=r'$\sigma=0.16$')
ax.plot(E, exp_share_CD, lw=2, label=r'$\sigma=1$ (Cobb-Douglas)')
ax.axvline(E_after_drop, color='k', lw=1.5)
ax.set_xlim(Emin, Emax)
ax.set_ylim(0, 1)
ax.set_xlabel('Energy, $E$')
ax.set_ylabel('Expenditure share on energy $E$')
ax.legend(loc='upper right')
ax.grid(True)
fig.savefig('figures/exp_share_fig_gas.pdf')
fig.savefig('figures/exp_share_fig_gas.eps')
plt.close(fig)

# --- Numeric summary, matching the MATLAB `disp` blocks -------------------
index = nearest_index(E, E_after_drop)

print_row('Elasticity = 0 (Leontief)', Y_Leontief[index] / Ybench_Leontief - 1, 'Inf', 0, 1)

print_row('Elasticity = 0.1', Y[index, 0] / Ybench[0] - 1,
          MPE[index, 0] / MPE[index_bench, 0], MPX[index, 0] / MPX[index_bench, 0], exp_share[index, 0])

print_row('Elasticity = 0.2', Y[index, 1] / Ybench[1] - 1,
          MPE[index, 1] / MPE[index_bench, 1], MPX[index, 1] / MPX[index_bench, 1], exp_share[index, 1])

print_row('Elasticity = 1 (Cobb-Douglas)', Y_CD[index] / Ybench_CD - 1,
          MPE_CD[index] / MPE_CD[index_bench], MPX_CD[index] / MPX_CD[index_bench], alpha)
