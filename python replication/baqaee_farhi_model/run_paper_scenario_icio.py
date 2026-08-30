"""
Recreates the paper's "What if" scenario (EU iceberg trade cost on Russia)
on the real OECD ICIO 2022 "SML" release instead of WIOD 2008, using
combined_trade_elast.csv (Fontagne-Guimbard-Orefice + Ahmad & Schreiber 2024
GTAP services + WIOD fallback) instead of trade_elast_2008.mat.

Country set: the same 40 individually-tracked countries as
run_paper_scenario.py's WIOD scenario (all present in ICIO under the same
ISO3 codes), plus every other ICIO country -- including ICIO's own
pre-existing 'ROW' aggregate -- folded into a composite 'ROW' via
icio_to_haio.py's `row_label` mechanism. Sectors: all 48 non-'T', non-'O'
ICIO sectors (T is always excluded -- see icio_to_haio.py; O/public
administration is excluded because it has no trade elasticity, by the same
non-traded-by-nature reasoning already applied on the WIOD side).

This is a materially different exercise from run_paper_scenario.py, not a
bugfix regression check: different base year (2022 vs 2008), different
sector classification (48 ICIO sectors vs 30 WIOD sectors), different
country coverage details, and elasticities from three combined literature
sources instead of the paper's own (mostly-placeholder) trade_elast_2008.mat.
Expect the magnitude to differ from the WIOD run for these reasons, not
because either run is wrong.
"""
import argparse
import time
import numpy as np
import pandas as pd

from icio_to_haio import icio_to_haio, load_combined_trade_elast
from run_model import run_scenario

COUNTRIES = ['AUS', 'AUT', 'BEL', 'BGR', 'BRA', 'CAN', 'CHN', 'CYP', 'CZE', 'DEU', 'DNK', 'ESP',
             'EST', 'FIN', 'FRA', 'GBR', 'GRC', 'HUN', 'IDN', 'IND', 'IRL', 'ITA', 'JPN', 'KOR',
             'LTU', 'LUX', 'LVA', 'MEX', 'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'RUS', 'SVK', 'SVN',
             'SWE', 'TUR', 'TWN', 'USA', 'ROW']

EU = ['AUT', 'BEL', 'BGR', 'CYP', 'CZE', 'DEU', 'DNK', 'ESP', 'EST', 'FIN', 'FRA', 'GRC', 'HUN',
      'IRL', 'ITA', 'LTU', 'LUX', 'LVA', 'MLT', 'NLD', 'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'SWE']
RUS = 'RUS'

ALL_SECTORS = ['A01', 'A02', 'A03', 'B05', 'B06', 'B07', 'B08', 'B09', 'C10T12', 'C13T15', 'C16',
               'C17_18', 'C19', 'C20', 'C21', 'C22', 'C23', 'C24A', 'C24B', 'C25', 'C26', 'C27',
               'C28', 'C29', 'C301', 'C302T309', 'C31T33', 'D', 'E', 'F', 'G', 'H49', 'H50', 'H51',
               'H52', 'H53', 'I', 'J58T60', 'J61', 'J62_63', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R',
               'S', 'T']
SECTORS = [s for s in ALL_SECTORS if s not in ('T', 'O')]  # T: div-by-zero; O: no trade elasticity


def build_haio(icio_path='2022_SML.csv'):
    df = pd.read_csv(icio_path, index_col=0)
    trade_elast = load_combined_trade_elast(SECTORS)
    return icio_to_haio(df, COUNTRIES, SECTORS, trade_elast, row_label='ROW')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ngrid', type=int, default=20)
    parser.add_argument('--intensity', type=float, default=150)
    parser.add_argument('--icio-path', default='2022_SML.csv')
    args = parser.parse_args()

    print('Building HAIO from real ICIO data + combined_trade_elast.csv...')
    t0 = time.time()
    haio = build_haio(args.icio_path)
    print(f'HAIO built in {time.time() - t0:.1f}s (C={haio["C"]}, N={len(SECTORS)})\n')

    eu_pos = [COUNTRIES.index(c) + 1 for c in EU]
    rus_pos = [COUNTRIES.index(RUS) + 1]
    shocks = [{'sellers': rus_pos, 'buyers': eu_pos, 'sectors': None, 'intensity': args.intensity}]

    t0 = time.time()
    result = run_scenario(None, COUNTRIES, shocks, ngrid=args.ngrid, haio=haio)
    elapsed = time.time() - t0
    print(f'\nSolved in {elapsed:.1f}s ({args.ngrid} discretization steps)\n')
    print(f'Elasticities used: {result["elasticities"]}\n')
    print(f'{"Country":8s} {"dlogW (%)":>10s}')
    for name, w in result['dlogW'].items():
        print(f'{name:8s} {100 * w:10.3f}')
    print(f'{"World":8s} {100 * result["World"]:10.3f}')


if __name__ == '__main__':
    main()
