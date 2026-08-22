"""
Convenience wrapper reproducing the exact scenario run by
replication/baqaee_farhi_model/main_dlogW_rev_bigshocks_EU_Russian_v2.m:
all 41 WIOD 2008 countries (minus the standalone "rest of world" country,
which is appended automatically), an EU iceberg trade cost on Russian goods
large enough to send EU-Russia trade to roughly zero, factor_index == 2
(country-sector-specific labor), initial_tariff_index == 1 (no initial
tariffs).

Usage: python run_paper_scenario.py [--ngrid 20] [--intensity 150]

Runtime note: this solves a (C+CF) ~= 1271-dimensional linear system by
probing (see nested_ces.py's module docstring) once per discretization step,
so it is CPU-heavy -- expect on the order of tens of minutes on a single
core (each of the ngrid steps costs roughly (C+CF+1) evaluations of
`response()`, and `response()` itself is O(CN^2 + L^2) per call). This is
unavoidably a large linear-algebra workload at this scale; the original
MATLAB code faces the same C+CF-dimensional system, just solved via an
explicit Jacobian construction that additionally requires ~30-60 GB of RAM
this environment does not have (confirmed: it OOMs under GNU Octave 8.4.0
with 15 GB available) -- see nested_ces.py's docstring for why probing
avoids that memory cost.
"""
import argparse
import time
from run_model import run

# WIOD 2008 country order used throughout the replication package (see
# main_dlogW_rev_bigshocks_EU_Russian_v2.m); ROW is appended as the 42nd.
COUNTRIES = ['AUS', 'AUT', 'BEL', 'BGR', 'BRA', 'CAN', 'CHN', 'CYP', 'CZE', 'DEU', 'DNK', 'ESP',
             'EST', 'FIN', 'FRA', 'GBR', 'GRC', 'HUN', 'IDN', 'IND', 'IRL', 'ITA', 'JPN', 'KOR',
             'LTU', 'LUX', 'LVA', 'MEX', 'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'RUS', 'SVK', 'SVN',
             'SWE', 'TUR', 'TWN', 'USA', 'ROW']

KEEP_C = list(range(1, 35)) + list(range(36, 42))  # WIOD 1..41 minus 35 (ROW)
# EU positions *within KEEP_C* (not WIOD country codes -- see run_model.py's
# docstring); copied verbatim from the original MATLAB script's EU/RUS.
EU = [2, 3, 4, 8, 9, 10, 11, 12, 13, 14, 15, 17, 18, 21, 22, 25, 26, 27, 29, 30, 31, 32, 33, 35, 36, 37]
RUS = 34


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ngrid', type=int, default=20)
    parser.add_argument('--intensity', type=float, default=150)
    args = parser.parse_args()

    t0 = time.time()
    dlogW_sum, data = run(KEEP_C, EU=EU, RUS=RUS, ngrid=args.ngrid, intensity=args.intensity, data_dir='.')
    elapsed = time.time() - t0

    GNE_weights = data['chi_std'][:data['C']]
    world = GNE_weights @ dlogW_sum

    print(f'\nSolved in {elapsed:.1f}s ({args.ngrid} discretization steps)\n')
    print(f'{"Country":8s} {"dlogW (%)":>10s}')
    for name, w in zip(COUNTRIES, dlogW_sum):
        print(f'{name:8s} {100 * w:10.3f}')
    print(f'{"World":8s} {100 * world:10.3f}')


if __name__ == '__main__':
    main()
