#' Empirical estimates of the elasticity of substitution between energy and
#' other inputs, reviewed in Appendix Subsection A.4 of Bachmann et al.
#' (2024), and the paper's own conservative calibration choices (Subsection
#' A.9.2).
#'
#' Under the standard identification the paper relies on, the (absolute
#' value of the) own-price elasticity of energy demand equals the elasticity
#' of substitution sigma between energy and other inputs, so these two are
#' used interchangeably here, matching the paper's own text.

# Labandeira, Labeaga and Lopez-Otero (2017), "A meta-analysis on the price
# elasticity of energy demand", Energy Policy 102, 549-68 -- averages across
# their full sample of estimates, by energy type and horizon (short run =
# within one year, long run = more than one year).
LABANDEIRA_ET_AL_2017 <- list(
  energy_short_run = -0.221,        # 376 studies, energy in general
  energy_long_run = -0.584,
  natural_gas_short_run = -0.18,    # 230 studies, natural gas only
  natural_gas_long_run = -0.684,
  heating_oil_short_run = -0.017,   # 44 studies, heating oil
  heating_oil_long_run = -0.185,
  aggregate_short_run = -0.236,     # full sample average, before dropping outliers
  aggregate_long_run = -0.596
)

# Auffhammer and Rubin (2018), "Natural gas price elasticities and optimal
# cost recovery under consumer heterogeneity: evidence from 300 million
# natural gas bills", NBER WP 24295 -- residential natural gas, short run.
AUFFHAMMER_RUBIN_2018 <- list(
  natural_gas_short_run_low = -0.17,
  natural_gas_short_run_high = -0.20
)

# The paper's own calibration (Subsection A.9.2): deliberately conservative,
# well below the empirical range above, "to build in a dose of caution"
# (Section 3 of the main text).
PAPER_CALIBRATION <- list(
  aggregate_energy = 0.04,               # used for the 10% energy-drop scenario (elasticity.R)
  natural_gas = 0.1,                     # used for the 30% gas-drop scenario (elasticity_gas.R)
  natural_gas_implausible_floor = 0.04   # footnote 4: below this, results become implausible
)

#' Look up |own-price elasticity of energy demand| == elasticity of
#' substitution sigma from the tables above. `source` is a dotted key, e.g.
#' 'labandeira2017.natural_gas_short_run', 'paper.aggregate_energy'.
sigma_from_literature <- function(source) {
  tables <- list(labandeira2017 = LABANDEIRA_ET_AL_2017,
                 auffhammer_rubin2018 = AUFFHAMMER_RUBIN_2018,
                 paper = PAPER_CALIBRATION)
  parts <- strsplit(source, '.', fixed = TRUE)[[1]]
  table_name <- parts[1]
  key <- parts[2]
  if (is.null(tables[[table_name]]) || is.null(tables[[table_name]][[key]])) {
    available <- unlist(lapply(names(tables), function(t) paste0(t, '.', names(tables[[t]]))))
    stop(paste0("Unknown source '", source, "'. Available: ", paste(available, collapse = ', ')))
  }
  value <- tables[[table_name]][[key]]
  if (table_name != 'paper') abs(value) else value
}
