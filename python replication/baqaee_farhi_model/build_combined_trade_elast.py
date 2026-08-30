"""
Builds combined_trade_elast.csv: one trade_elast value per ICIO sector,
combining Fontagne-Guimbard-Orefice (goods, tariff-based), Ahmad &
Schreiber 2024 (services, markup-based, GTAP sector level), and a small
third layer sourced from the original Bachmann-Baqaee et al. paper's own
WIOD trade_elast_2008.mat for the handful of sectors neither of the first
two cover, each converted to nested_ces.py's convention, with a `source`
column tracking provenance. Re-run this script
(`python build_combined_trade_elast.py`) if any upstream source is
revised; the output CSV is what `load_combined_trade_elast()` in
icio_to_haio.py actually reads.

CONVENTION: nested_ces.py needs trade_elast such that trade_elast + 1 equals
the CES elasticity of substitution across origin countries (see
icio_to_haio.py's module docstring). The three sources report different
objects and need DIFFERENT transforms:
  - Fontagne's epsilon_icio is a NEGATIVE trade-cost elasticity
    (epsilon = -(sigma-1) is the standard convention), so
    trade_elast = sigma - 1 = -epsilon_icio.
  - Ahmad & Schreiber's EOS is already a directly-estimated CES elasticity
    of substitution (sigma itself, positive, from a markup/monopolistic-
    competition identification), so trade_elast = EOS - 1.
  - The original paper's trade_elast_2008.mat is already in exactly this
    project's own `trade_elast` convention (it's WIOD's real data file,
    read as-is by io_reorder.py) -- no transform needed, use its values
    directly.
Mixing these up (e.g. negating an EOS value, or not negating an
epsilon_icio value) would silently produce a wrong-signed or wrong-scale
number, so this distinction is the main thing to get right when extending
this script.

GTAP -> ICIO sector mapping judgment calls (all explicit here, not buried):
  - Direct 1:1 matches (same real-world activity, different code): cns->F,
    trd->G, otp->H49, wtp->H50, atp->H51, whs->H52, afs->I, edu->P, hht->Q.
  - cmn (GTAP's one broad ICT/communications bucket) is duplicated across
    ICIO's J61 (telecom) and J62_63 (computer/IT services) -- NOT J58T60
    (publishing/broadcasting), which Fontagne already covers directly with
    its own tariff-based estimate, so that one is left alone.
  - K (ICIO's single "Financial and insurance activities") averages GTAP's
    separate ins (insurance) and ofi (banking/other financial
    intermediation) -- ICIO doesn't split these, GTAP does.
  - obs (GTAP's "other business services", NAICS 532/533/541/561) is used
    for ICIO's N (administrative and support services, closest to NAICS
    561) specifically, NOT M (professional/scientific/technical, NAICS
    541) -- M is left to Fontagne's own direct estimate instead, since
    NAICS 541 is the better match for M and Fontagne already has M covered.
  - ros (GTAP's "recreation and other services") is duplicated across
    ICIO's R (arts/entertainment/recreation) and S (other service
    activities) -- an imperfect match for S specifically (broader in ICIO
    than GTAP's NAICS 81 personal-services scope), flagged 'approximate'.

WIOD fallback layer (A02, B05, B09, D) -- found by matching the original
MATLAB driver script's hardcoded 30-sector label list
(main_dlogW_rev_bigshocks_EU_Russian_v2.m) against trade_elast_2008.mat's
actual values by index:
  - A02 (forestry) <- WIOD's "Agriculture, Hunting, Forestry and Fishing"
    (8.11, a REAL non-default estimate, just bundling forestry with
    agriculture and fishing at WIOD's coarser aggregation) -- quality
    'approximate'.
  - B05 (coal mining) and B09 (mining support) <- WIOD's "Mining and
    Quarrying" (15.72, also a real bundled estimate; ICIO's B06/B07/B08
    already have their own finer Fontagne estimates, so this only needs
    to cover the two mining sub-sectors that fell through the cracks) --
    quality 'approximate'.
  - D (electricity/gas/steam supply) <- WIOD's "Electricity, Gas and Water
    Supply" (5.00) -- BUT this is NOT a real estimate: every WIOD sector
    from index 14 onward (electricity through health/social work, all
    services) is exactly 5.00 with zero variation, meaning the paper's own
    authors used a single flat placeholder for everything past
    manufacturing rather than sourcing real values. Using it here is
    adopting the paper's own admitted guess, not resolving the gap --
    quality 'placeholder', tracked as a genuinely different (weaker) tier
    than 'direct'/'approximate'.

O (public administration) is deliberately NOT filled from this WIOD
layer, even though the paper's own value is there too (also 5.00, the
same flat placeholder as D) -- government services are essentially
non-traded by nature, so leaving O with no trade elasticity at all is
arguably correct, not a gap, and filling it with a number already known
to be an un-sourced placeholder for a sector where the concept barely
applies would be strictly worse than leaving it out. Sector T is excluded
everywhere in this project regardless (see icio_to_haio.py).
"""
import pandas as pd
import scipy.io as sio
import numpy as np

FONTAGNE_CSV = 'fontagne_icio_trade_elast.csv'
WIOD_TRADE_ELAST_MAT = 'trade_elast_2008.mat'
OUTPUT_CSV = 'combined_trade_elast.csv'

# icio_sector -> (WIOD 30-sector label-list index (0-indexed, matching
# io_reorder.py's post-merge trade_elast array), mapping_quality, note)
WIOD_TO_ICIO = {
    'A02': (0, 'approximate', "WIOD 'Agriculture, Hunting, Forestry and Fishing' -- real estimate, bundles forestry with agriculture/fishing"),
    'B05': (1, 'approximate', "WIOD 'Mining and Quarrying' -- real estimate, bundles coal mining with all other mining"),
    'B09': (1, 'approximate', "WIOD 'Mining and Quarrying' -- real estimate, bundles mining support with all other mining"),
    'D':   (14, 'placeholder', "WIOD 'Electricity, Gas and Water Supply' -- NOT a real estimate, this is the paper's own flat 5.00 default used for every services/utility sector past manufacturing"),
}

# Ahmad & Schreiber (2024), USITC Working Paper, Table 10, EOS (2013-2022) column.
GTAP_EOS = {
    'afs': 7.371, 'atp': 5.534, 'cmn': 2.914, 'cns': 5.418, 'edu': 8.150,
    'hht': 9.691, 'ins': 4.128, 'obs': 4.827, 'ofi': 5.299, 'otp': 4.696,
    'ros': 5.843, 'rsa': 1.684, 'trd': 5.194, 'whs': 18.28, 'wtp': 7.264,
}

# icio_sector -> (gtap_sector_or_list, mapping_quality, note)
GTAP_TO_ICIO = {
    'F':       (['cns'], 'direct',      'construction'),
    'G':       (['trd'], 'direct',      'wholesale/retail trade'),
    'H49':     (['otp'], 'direct',      'land transport'),
    'H50':     (['wtp'], 'direct',      'water transport'),
    'H51':     (['atp'], 'direct',      'air transport'),
    'H52':     (['whs'], 'direct',      'warehousing/support activities'),
    'I':       (['afs'], 'direct',      'accommodation/food services'),
    'P':       (['edu'], 'direct',      'education'),
    'Q':       (['hht'], 'direct',      'human health/social work'),
    'L':       (['rsa'], 'direct',      'real estate'),
    'J61':     (['cmn'], 'approximate', 'GTAP cmn is one broad ICT bucket duplicated across ICIO J61/J62_63'),
    'J62_63':  (['cmn'], 'approximate', 'GTAP cmn is one broad ICT bucket duplicated across ICIO J61/J62_63'),
    'K':       (['ins', 'ofi'], 'approximate', 'ICIO K combines insurance+banking, GTAP splits them -- averaged'),
    'N':       (['obs'], 'approximate', 'GTAP obs (NAICS 561 admin/support closest) used for N, not M (Fontagne covers M directly)'),
    'R':       (['ros'], 'approximate', 'GTAP ros used for both ICIO R and S'),
    'S':       (['ros'], 'approximate', 'GTAP ros used for both ICIO R and S; weaker match for S specifically'),
}


def build():
    fontagne = pd.read_csv(FONTAGNE_CSV, index_col='icio2025')
    rows = []

    for sector, row in fontagne.iterrows():
        if pd.notna(row['epsilon_icio']):
            rows.append(dict(
                icio_sector=sector,
                trade_elast=-row['epsilon_icio'],
                source='fontagne_guimbard_orefice_2022',
                mapping_quality='direct',
                note=f"epsilon_icio={row['epsilon_icio']:.3f}, trade_elast = -epsilon_icio",
            ))

    for sector, (gtap_sectors, quality, note) in GTAP_TO_ICIO.items():
        eos_values = [GTAP_EOS[g] for g in gtap_sectors]
        eos_mean = sum(eos_values) / len(eos_values)
        rows.append(dict(
            icio_sector=sector,
            trade_elast=eos_mean - 1,
            source='ahmad_schreiber_2024',
            mapping_quality=quality,
            note=f"GTAP {'+'.join(gtap_sectors)} EOS={eos_mean:.3f}, trade_elast = EOS - 1 ({note})",
        ))

    wiod_raw = sio.loadmat(WIOD_TRADE_ELAST_MAT)['trade_elast_2008'].ravel()
    wiod_trade_elast = np.concatenate([wiod_raw[:8], wiod_raw[9:]])  # exact io_reorder.py processing
    for sector, (wiod_idx, quality, note) in WIOD_TO_ICIO.items():
        value = wiod_trade_elast[wiod_idx]
        rows.append(dict(
            icio_sector=sector,
            trade_elast=value,
            source='bachmann_baqaee_2024_wiod',
            mapping_quality=quality,
            note=f"trade_elast_2008.mat[{wiod_idx}]={value:.3f}, used as-is ({note})",
        ))

    out = pd.DataFrame(rows).sort_values('icio_sector').reset_index(drop=True)
    out.to_csv(OUTPUT_CSV, index=False)
    print(f'Wrote {len(out)} sectors to {OUTPUT_CSV}')
    print(f"  from Fontagne: {(out['source']=='fontagne_guimbard_orefice_2022').sum()}")
    print(f"  from Ahmad-Schreiber/GTAP: {(out['source']=='ahmad_schreiber_2024').sum()}")
    print(f"  from paper's own WIOD data: {(out['source']=='bachmann_baqaee_2024_wiod').sum()}")
    return out


if __name__ == '__main__':
    build()
