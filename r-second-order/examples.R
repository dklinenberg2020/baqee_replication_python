#' Worked examples reproducing every second-order-approximation and
#' calibration number printed in the paper, run as a lightweight validation
#' script: `Rscript examples.R`.

source('ces.R')
source('second_order.R')
source('sigma_literature.R')

check <- function(label, got, want, tol) {
  ok <- abs(got - want) < tol
  status <- if (ok) 'OK' else 'MISMATCH'
  cat(sprintf('[%s] %s: got %+.4f%%, paper says %+.4f%%\n', status, label, got, want))
  stopifnot(ok)
}

cat('=== 1. Equations (A4) and (A5) are algebraically identical ===\n')
for (sigma in c(0.04, 0.1, 0.5, 0.9, 2.0)) {
  for (case in list(c(0.04, -0.10), c(0.01, -0.30), c(0.25, 0.05))) {
    alpha <- case[1]; dlogE <- case[2]
    a4 <- second_order_simple(alpha, dlogE, sigma)
    a5 <- second_order_share_implied(alpha, dlogE, sigma)
    stopifnot(abs(a4 - a5) < 1e-12)
  }
}
cat('OK: second_order_simple (A4) == second_order_share_implied (A5) for all tested (alpha, dlogE, sigma).\n\n')

cat('=== 2. Calibration (Subsection A.9.1): alpha reproduces the target share ===\n')
for (target in c(0.04, 0.012, 0.10)) {
  for (sigma in c(0.04, 0.5, 1.0, 2.0)) {
    alpha <- calibrate_alpha(target, sigma = sigma)
    share <- ces_expenditure_share(alpha, sigma)
    stopifnot(abs(share - target) < 1e-12)
    stopifnot(abs(alpha - target) < 1e-12)  # with pE=pX=1, alpha == target share for every sigma
  }
}
cat('OK: calibrate_alpha() reproduces the target share exactly, for every sigma (pE=pX=1).\n')
cat(sprintf("  paper's own calibration: alpha=%.3f (aggregate energy, 4%% of GNE)\n", calibrate_alpha(0.04)))
cat(sprintf("  paper's own calibration: alpha=%.3f (gas alone, 1.2%% of GNE)\n\n", calibrate_alpha(0.012)))

cat('=== 2b. Calibration with non-unit relative prices (general form of A9/A10) ===\n')
alpha_general <- calibrate_alpha(0.04, sigma = 0.5, pE = 1.5, pX = 1.0)
share_check <- ces_expenditure_share(alpha_general, sigma = 0.5, pE = 1.5, pX = 1.0)
check('share reproduced at pE=1.5', 100 * share_check, 100 * 0.04, 1e-9)
cat('\n')

cat('=== 3. Literature-based sigma calibration (Subsection A.4) ===\n')
cat(sprintf('  Labandeira et al. (2017), natural gas, short run: sigma = %s\n',
            sigma_from_literature('labandeira2017.natural_gas_short_run')))
cat(sprintf('  Auffhammer & Rubin (2018), natural gas, short run: sigma in [%s, %s]\n',
            sigma_from_literature('auffhammer_rubin2018.natural_gas_short_run_low'),
            sigma_from_literature('auffhammer_rubin2018.natural_gas_short_run_high')))
cat(sprintf("  Paper's own (conservative) choice: sigma = %s (aggregate energy) vs. literature short-run ~0.18-0.24 -- deliberately far more conservative.\n\n",
            sigma_from_literature('paper.aggregate_energy')))

cat('=== 4. Back-of-envelope Delta log W calculations, Subsection A.5.3 ===\n')
check('Extreme scenario (eq. before A8)',
      100 * second_order_from_share(share = 0.025, dlogE = -0.30, dshare = 0.05),
      -1.5, 0.02)
check('Main scenario (eq. before A8)',
      100 * second_order_from_share(share = 0.025, dlogE = -0.17, dshare = 0.025),
      -0.63, 0.02)
check('Preferred back-of-envelope, gas-only (eq. A8)',
      100 * second_order_from_share(share = 0.012, dlogE = -0.30, dshare = 0.024),
      -0.72, 0.02)
cat('\n')

cat('=== 5. Where the second-order approximation breaks down (footnote 12) ===\n')
cat('The paper explicitly does NOT use equation (A4) for its own headline numbers\n')
cat('because it is unreliable for sigma close to zero. Comparing against the exact\n')
cat("computation (../python replication/elasticity.py, validated against the original\n")
cat("MATLAB code) at the paper's own two calibrations:\n\n")

scenarios <- list(
  list('10% energy drop (alpha=0.04, sigma=0.04)', 0.04, 0.04, -0.10, -1.57),
  list('30% gas drop (alpha=0.01, sigma=0.1)', 0.01, 0.10, -0.30, -2.33)
)
for (s in scenarios) {
  label <- s[[1]]; alpha <- s[[2]]; sigma <- s[[3]]; dlogE <- s[[4]]; exact_pct <- s[[5]]
  approx_pct <- 100 * second_order_simple(alpha, dlogE, sigma)
  cat(sprintf('  %s:\n', label))
  cat(sprintf('    second-order approximation (A4): %+.2f%%\n', approx_pct))
  cat(sprintf('    exact computation (elasticity.py): %+.2f%%\n', exact_pct))
  divergence <- approx_pct - exact_pct
  size <- if (abs(divergence) > 0.3) 'large' else 'small'
  cat(sprintf('    -> %s divergence (%+.2f pp), as the paper warns.\n\n', size, divergence))
}

cat('All checks passed.\n')
