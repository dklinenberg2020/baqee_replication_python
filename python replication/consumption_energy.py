"""
Python translation of replication/consumption_energy.do.

Constructs household energy expenditure/income shares from the 2013 SUF of
the German EVS (Einkommens- und Verbrauchsstichprobe) and reproduces the
paper's Section 3 / Appendix A.12 figures on the distribution of energy
expenditure and income shares.

DATA ACCESS: the EVS Scientific Use File is restricted-access microdata from
the German Federal Statistical Office (Forschungsdatenzentrum), obtained
under a separate data use agreement -- see
https://www.forschungsdatenzentrum.de/de/haushalte/evs. It is NOT included
in this replication package (nor was it in the original MATLAB/Stata
package) and this script cannot be run or tested without it. This is a
line-by-line logic translation of consumption_energy.do only; it has NOT
been executed against real data.

Usage (once you have the data):
    python consumption_energy.py --input data/evs2013_aa_gs_hb.dta --outdir figures/

The input .dta file needs pandas' `read_stata`, which itself requires the
original EVS `.dta` extract (Stata's own DTA format) to be readable by
pandas -- this is normally the case for a data extract exported from Stata
by the research data center.
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def weighted_quantile_group(values, weights, n_quantiles=5):
    """Equivalent of Stata's `xtile ... [pw=weight], nq(n)`: assigns each
    observation to one of `n_quantiles` groups of (approximately) equal
    total weight, ordered by `values`. NaN values are left as NaN."""
    result = pd.Series(np.nan, index=values.index)
    mask = values.notna() & weights.notna()
    v = values[mask]
    w = weights[mask]
    order = v.sort_values().index
    cum_w = w.loc[order].cumsum()
    total_w = w.loc[order].sum()
    group = np.ceil(cum_w / total_w * n_quantiles).clip(1, n_quantiles)
    result.loc[order] = group.values
    return result


def weighted_mean(values, weights):
    mask = values.notna() & weights.notna()
    return np.average(values[mask], weights=weights[mask])


def build_consumption_dataset(df):
    """Section: CONSTRUCT GENERAL VARIABLES / INCOME AND EXPENDITURE VARIABLES.
    `df` is the raw EVS 2013 SUF extract (one row per household)."""
    d = df.copy()

    # --- weights, quarter, age, household composition ---------------------
    d['wgt'] = d['EF107']
    d['quarter'] = d['EF6']
    d['age_group'] = d['EF40']
    d['age_group_large'] = np.select(
        [d['age_group'] < 5, d['age_group'].between(5, 6), d['age_group'].between(7, 9),
         d['age_group'].between(10, 12), d['age_group'] > 12],
        [1, 2, 3, 4, 5], default=np.nan)

    d['age'] = 2013 - d['EF8U3']
    d['age_group_own'] = pd.cut(
        d['age'], bins=[-np.inf, 30, 40, 45, 50, 55, 60, 65, 70, np.inf],
        labels=[1, 2, 3, 4, 5, 6, 7, 8, 9], right=False).astype(float)

    d['renter'] = (d['EF20'] == 3).astype(int)
    d['home_owner'] = d['EF20'].isin([1, 2]).astype(int)

    d['size'] = d['EF7']
    children_map = {**{k: 0 for k in [1, 2, 9, 10, 21, 22]}, **{k: 1 for k in [3, 11, 12, 23, 24]},
                     **{k: 2 for k in [5, 13, 14, 25, 26]}, **{k: 3 for k in [7, 15, 16, 27, 28]},
                     **{k: 4 for k in [17, 18]}, **{k: 5 for k in [19, 20]}}
    d['children'] = d['EF39'].map(children_map)
    d['size_3'] = d['size'].clip(upper=3)

    d['city_size'] = np.select(
        [d['EF4'].isin([1, 2, 6]), d['EF4'].isin([3, 8]), d['EF4'].isin([4, 9]), d['EF4'] == 5],
        [1, 2, 3, 4], default=np.nan)

    d['main_earner_male'] = (d['EF8U2'] == 1).astype(int)

    for i in range(1, 9):
        d[f'social_pos_{i}'] = d[f'EF{i + 7}U8']

    d['occupation'] = d['EF8U19']
    d.loc[d['EF8U8'] == 9, 'occupation'] = 23
    d.loc[d['EF8U8'] == 10, 'occupation'] = 24
    d.loc[d['EF8U8'] == 11, 'occupation'] = 25

    d['education'] = d['EF8U7']
    d['education_school'] = np.select(
        [d['EF8U7'] < 5, d['EF8U7'] > 4], [1, 2], default=np.nan)
    d.loc[d['education'] > 10, 'education_school'] = 3

    for i in range(1, 9):
        t = i + 7
        emp = pd.Series(np.nan, index=d.index)
        emp[(d[f'social_pos_{i}'] != 10) & (d[f'social_pos_{i}'] != 11)] = 0
        emp[(d[f'EF{t}U13'] > 0) & d[f'EF{t}U13'].notna()] = 1
        d[f'emp_{i}'] = emp
    d['number_employed'] = sum(d[f'emp_{i}'] for i in range(1, 9))
    d['one_earner_HH'] = np.where(d['emp_1'].notna(), (d['number_employed'] == 1).astype(int), np.nan)

    d['unemployed'] = (d['social_pos_1'] == 9).astype(int)

    d['living_space'] = d['EF21']
    d['living_space_pc'] = d['living_space'] / d['size']

    # --- consumption ---------------------------------------------------
    d['c_insurance'] = d['EF98']
    d['c_further'] = d['EF530'] + d['EF472'] + d['EF473'] + d['EF474'] + d['EF476'] + d['EF531']
    d['c_all'] = d['EF89'] - d['EF77'] + d['c_insurance'] + d['c_further']
    d['c_all_imp_rents'] = d['c_all'] + d['EF77']
    d['c_all_no_insurance'] = d['c_all'] - d['c_insurance']
    d['c_all_no_ins_no_exp'] = d['EF89'] - d['EF77']

    d['c_food'] = d['EF73'] + d['EF74']
    d['c_house'] = d['EF76'] + d['EF78'] + d['EF79'] + d['EF530']
    d['c_h_imp_rents'] = d['c_house'] + d['EF77']
    d['c_h_credit'] = d['c_house'] + d['EF102']
    d['c_cloth'] = d['EF75']
    d['c_transp'] = d['EF82']
    d['c_comm'] = d['EF84']
    d['c_leisure'] = d['EF85'] + d['EF87'] + d['EF472'] + d['EF476']

    d['c_other'] = (d['c_all'] - d['c_food'] - d['c_house'] - d['c_cloth'] - d['c_transp']
                     - d['c_comm'] - d['c_leisure'] - d['c_insurance'])
    d['c_other_no_insurance'] = (d['c_all_no_insurance'] - d['c_food'] - d['c_house'] - d['c_cloth']
                                  - d['c_transp'] - d['c_comm'] - d['c_leisure'])
    d['c_all_credit'] = d['c_all'] + d['EF102']

    # --- expenditures ---------------------------------------------------
    d['exp_further_taxes'] = d['EF96']
    d['exp_voluntary_insurance'] = d['EF98']
    d['exp_further_transfers'] = d['EF100']
    d['exp_accumulate_wealth'] = d['EF101']
    d['exp_repay_loans'] = d['EF102']
    d['exp_further'] = d['EF103']
    d['expenditures_all'] = d['c_all'] + d['exp_further_taxes'] + d['exp_voluntary_insurance'] + d['exp_further_transfers']

    # --- income ---------------------------------------------------------
    d['income_total'] = d['EF72']
    d['income_net'] = d['EF62']
    d['income_net_disp_predefined'] = d['EF65']
    d['income_net_disp'] = (d['EF65'] - d['exp_further_taxes']
                             - d[['EF237U1', 'EF237U2', 'EF237U3', 'EF237U4', 'EF237U5', 'EF237U6']].sum(axis=1)
                             - d['EF475'] - d['EF477'] - d['EF529']
                             - d[['EF238U1', 'EF238U2', 'EF238U3', 'EF238U4', 'EF238U5', 'EF238U6']].sum(axis=1))
    d.loc[d['income_net_disp'] < 1000, 'income_net_disp'] = np.nan

    return d


def compute_income_quantiles(d):
    """Section: COMPUTE QUINTILES OF INCOME."""
    d = d.copy()
    d['income_total_no_size'] = weighted_quantile_group(d['income_total'], d['wgt'])
    d['income_net_no_size'] = weighted_quantile_group(d['income_net'], d['wgt'])
    d['income_dec_no_size'] = weighted_quantile_group(d['income_net_disp'], d['wgt'])

    for hhsize in (1, 2, 3):
        m = d['size_3'] == hhsize
        d.loc[m, f'income_total_{hhsize}'] = weighted_quantile_group(d.loc[m, 'income_total'], d.loc[m, 'wgt'])
        d.loc[m, f'income_net_{hhsize}'] = weighted_quantile_group(d.loc[m, 'income_net'], d.loc[m, 'wgt'])
        d.loc[m, f'income_dec_{hhsize}'] = weighted_quantile_group(d.loc[m, 'income_net_disp'], d.loc[m, 'wgt'])

    return d


def compute_energy_sources(d):
    """Section: ENERGY SOURCES. Drops households heating with electricity
    (EF23==1), constructs energy expenditure by source and heating type."""
    d = d[d['EF23'] != 1].copy()

    d['energy_gas'] = d['EF317']
    d['energy_oil'] = d['EF320']
    d['energy_coal_wood'] = d['EF321']
    d['energy_warm_water'] = d['EF323']
    d['energy_dist_heat'] = d['EF324']
    d['energy_fuel'] = d['EF383']
    d['energy_total'] = (d['energy_gas'] + d['energy_oil'] + d['energy_coal_wood']
                          + d['energy_warm_water'] + d['energy_dist_heat'] + d['energy_fuel'])

    d['energy_main_gas'] = (d['EF23'] == 2).astype(int)
    d['energy_main_oil'] = (d['EF23'] == 3).astype(int)
    d['energy_main_coal_wood'] = (d['EF23'] == 4).astype(int)
    d['energy_main_other'] = (d['EF23'] == 5).astype(int)
    d['energy_main_no_info'] = (d['EF23'] == 0).astype(int)

    d['umlagen_dist_heat'] = d['EF327']
    d['umlagen_gas'] = d['EF329']
    d['umlagen_oil'] = d['EF330']
    d['umlagen_other'] = d['EF331'] + d['EF332']
    d['umlagen_total'] = d['umlagen_dist_heat'] + d['umlagen_gas'] + d['umlagen_oil'] + d['umlagen_other']

    d['expenditure_gas'] = d['energy_gas'] + d['umlagen_gas']
    d['expenditure_oil'] = d['energy_oil'] + d['umlagen_oil']
    d['expenditure_dist_heat'] = d['energy_dist_heat'] + d['umlagen_dist_heat']
    d['expenditure_coal_wood'] = d['energy_coal_wood']

    d.loc[d['energy_main_gas'] == 1, 'expenditure_gas'] += d.loc[d['energy_main_gas'] == 1, 'energy_warm_water']
    d.loc[d['energy_main_oil'] == 1, 'expenditure_oil'] += d.loc[d['energy_main_oil'] == 1, 'energy_warm_water']
    d.loc[d['energy_main_coal_wood'] == 1, 'expenditure_coal_wood'] += d.loc[d['energy_main_coal_wood'] == 1, 'energy_warm_water']
    d.loc[d['EF23'].isin([0, 5]), 'expenditure_dist_heat'] += d.loc[d['EF23'].isin([0, 5]), 'energy_warm_water']

    d['expenditure_total'] = (d['expenditure_gas'] + d['expenditure_oil']
                               + d['expenditure_coal_wood'] + d['expenditure_dist_heat'])
    d['expenditure_fuel'] = d['energy_fuel']

    for cat in ('gas', 'oil', 'coal_wood', 'dist_heat', 'total', 'fuel'):
        d[f'energy_share_c_{cat}'] = d[f'expenditure_{cat}'] / d['c_all'] * 100
        d[f'energy_share_y_{cat}'] = d[f'expenditure_{cat}'] / d['income_net'] * 100
        cons_share = weighted_mean(d[f'energy_share_c_{cat}'], d['wgt'])
        income_share = weighted_mean(d[f'energy_share_y_{cat}'], d['wgt'])
        print(f'{cat:>15} (cons/inc): {cons_share:6.1f} % {income_share:6.1f} %')

    d['heatingtype'] = d['EF23'].replace({2: 1, 3: 2, 4: 3, 5: 4, 0: 4})
    return d


HEATING_LABELS = {1: 'Gas', 2: 'Oil', 3: 'Coal & Wood', 4: 'District heating & other'}


def make_figures(d, outdir):
    """Section: bar-chart figures by heating type / income quantile /
    household size, mirroring the four `graph bar` blocks in the .do file."""
    for shi, estr in (('c', 'expenditure'), ('y', 'net income')):
        # All households by type of heating (heating + fuel shares)
        g = d.groupby('heatingtype').apply(
            lambda x: pd.Series({
                'heating': weighted_mean(x[f'energy_share_{shi}_total'], x['wgt']),
                'fuel': weighted_mean(x[f'energy_share_{shi}_fuel'], x['wgt']),
            }))
        g.index = g.index.map(HEATING_LABELS)
        ax = g.plot(kind='bar', stacked=False, figsize=(7, 5))
        ax.set_ylabel(f'{estr} share (in %)')
        plt.tight_layout()
        plt.savefig(f'{outdir}/{shi}_household_type_all.png')
        plt.close()

        # All households by income quantile, stacked by source
        g2 = d.groupby('income_net_no_size').apply(
            lambda x: pd.Series({src: weighted_mean(x[f'energy_share_{shi}_{src}'], x['wgt'])
                                  for src in ('gas', 'oil', 'coal_wood', 'dist_heat')}))
        ax = g2.plot(kind='bar', stacked=True, figsize=(7, 5))
        ax.legend(['Gas', 'Oil', 'Coal & Wood', 'District Heating'])
        plt.tight_layout()
        plt.savefig(f'{outdir}/{shi}_household_income_all.png')
        plt.close()

        # By heating type and income group (gas/oil/district heating only)
        sub = d[d['heatingtype'].isin([1, 2, 4])]
        g3 = sub.groupby(['heatingtype', 'income_net_no_size']).apply(
            lambda x: weighted_mean(x[f'energy_share_{shi}_total'], x['wgt'])).unstack(0)
        ax = g3.plot(kind='bar', figsize=(8, 5))
        ax.set_ylabel(f'{estr} share (in %)')
        plt.tight_layout()
        plt.savefig(f'{outdir}/{shi}_household_income_type_all.png')
        plt.close()

        # Same, split by household size
        for hhsize in (1, 2, 3):
            subs = d[(d['size_3'] == hhsize) & d['heatingtype'].isin([1, 2, 4])]
            g4 = subs.groupby(['heatingtype', f'income_net_{hhsize}']).apply(
                lambda x: weighted_mean(x[f'energy_share_{shi}_total'], x['wgt'])).unstack(0)
            ax = g4.plot(kind='bar', figsize=(8, 5))
            ax.set_ylabel(f'{estr} share (in %)')
            plt.tight_layout()
            plt.savefig(f'{outdir}/{shi}_household_income_type_size{hhsize}.png')
            plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', default='data/evs2013_aa_gs_hb.dta',
                         help='Path to the restricted-access EVS 2013 SUF extract (Stata .dta)')
    parser.add_argument('--outdir', default='figures')
    args = parser.parse_args()

    df = pd.read_stata(args.input)
    d = build_consumption_dataset(df)
    d = compute_income_quantiles(d)
    d = compute_energy_sources(d)
    make_figures(d, args.outdir)


if __name__ == '__main__':
    main()
