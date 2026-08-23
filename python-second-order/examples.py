"""
Worked examples reproducing every second-order-approximation and
calibration number printed in the paper, run as a lightweight validation
script: `python examples.py`.
"""
from ces import calibrate_alpha, ces_expenditure_share
from second_order import second_order_simple, second_order_share_implied, second_order_from_share
from sigma_literature import sigma_from_literature


def check(label, got, want, tol):
    ok = abs(got - want) < tol
    status = 'OK' if ok else 'MISMATCH'
    print(f'[{status}] {label}: got {got:+.4f}%, paper says {want:+.4f}%')
    assert ok, f'{label}: {got} vs {want}'


print('=== 1. Equations (A4) and (A5) are algebraically identical ===')
for sigma in (0.04, 0.1, 0.5, 0.9, 2.0):
    for alpha, dlogE in ((0.04, -0.10), (0.01, -0.30), (0.25, 0.05)):
        a4 = second_order_simple(alpha, dlogE, sigma)
        a5 = second_order_share_implied(alpha, dlogE, sigma)
        assert abs(a4 - a5) < 1e-12, (sigma, alpha, dlogE, a4, a5)
print('OK: second_order_simple (A4) == second_order_share_implied (A5) for all tested (alpha, dlogE, sigma).\n')


print('=== 2. Calibration (Subsection A.9.1): alpha reproduces the target share ===')
for target in (0.04, 0.012, 0.10):
    for sigma in (0.04, 0.5, 1.0, 2.0):
        alpha = calibrate_alpha(target, sigma=sigma)
        share = ces_expenditure_share(alpha, sigma)
        assert abs(share - target) < 1e-12
        assert abs(alpha - target) < 1e-12  # with pE=pX=1, alpha == target share for every sigma
print('OK: calibrate_alpha() reproduces the target share exactly, for every sigma (pE=pX=1).')
print(f"  paper's own calibration: alpha={calibrate_alpha(0.04):.3f} (aggregate energy, 4% of GNE)")
print(f"  paper's own calibration: alpha={calibrate_alpha(0.012):.3f} (gas alone, 1.2% of GNE)\n")

print('=== 2b. Calibration with non-unit relative prices (general form of A9/A10) ===')
alpha_general = calibrate_alpha(0.04, sigma=0.5, pE=1.5, pX=1.0)
share_check = ces_expenditure_share(alpha_general, sigma=0.5, pE=1.5, pX=1.0)
check('share reproduced at pE=1.5', 100 * share_check, 100 * 0.04, 1e-9)
print()

print('=== 3. Literature-based sigma calibration (Subsection A.4) ===')
print(f"  Labandeira et al. (2017), natural gas, short run: sigma = {sigma_from_literature('labandeira2017.natural_gas_short_run')}")
print(f"  Auffhammer & Rubin (2018), natural gas, short run: sigma in "
      f"[{sigma_from_literature('auffhammer_rubin2018.natural_gas_short_run_low')}, "
      f"{sigma_from_literature('auffhammer_rubin2018.natural_gas_short_run_high')}]")
print(f"  Paper's own (conservative) choice: sigma = {sigma_from_literature('paper.aggregate_energy')} "
      f"(aggregate energy) vs. literature short-run ~0.18-0.24 -- deliberately far more conservative.\n")

print('=== 4. Back-of-envelope Delta log W calculations, Subsection A.5.3 ===')
# (i) All Russian energy imports cut, no substitution, extreme tripling of the import share.
check('Extreme scenario (eq. before A8)',
      100 * second_order_from_share(share=0.025, dlogE=-0.30, dshare=0.05),
      -1.5, 0.02)
# (ii) Substitute oil/coal but not gas; import share doubles.
check('Main scenario (eq. before A8)',
      100 * second_order_from_share(share=0.025, dlogE=-0.17, dshare=0.025),
      -0.63, 0.02)
# (iii) Gas as a separate input; import share triples -- the paper's preferred back-of-envelope number.
check('Preferred back-of-envelope, gas-only (eq. A8)',
      100 * second_order_from_share(share=0.012, dlogE=-0.30, dshare=0.024),
      -0.72, 0.02)
print()

print('=== 5. Where the second-order approximation breaks down (footnote 12) ===')
print('The paper explicitly does NOT use equation (A4) for its own headline numbers')
print('because it is unreliable for sigma close to zero. Comparing against the exact')
print('computation (../python replication/elasticity.py, validated against the original')
print('MATLAB code) at the paper\'s own two calibrations:\n')

scenarios = [
    ('10% energy drop (alpha=0.04, sigma=0.04)', 0.04, 0.04, -0.10, -1.57),
    ('30% gas drop (alpha=0.01, sigma=0.1)', 0.01, 0.10, -0.30, -2.33),
]
for label, alpha, sigma, dlogE, exact_pct in scenarios:
    approx_pct = 100 * second_order_simple(alpha, dlogE, sigma)
    print(f'  {label}:')
    print(f'    second-order approximation (A4): {approx_pct:+.2f}%')
    print(f'    exact computation (elasticity.py): {exact_pct:+.2f}%')
    print(f'    -> {"large" if abs(approx_pct - exact_pct) > 0.3 else "small"} divergence '
          f'({approx_pct - exact_pct:+.2f} pp), as the paper warns.\n')

print('All checks passed.')
