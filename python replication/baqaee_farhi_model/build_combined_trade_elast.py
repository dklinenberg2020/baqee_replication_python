"""
Builds combined_trade_elast.csv: one trade_elast value per ICIO sector,
combining Fontagne-Guimbard-Orefice (goods, tariff-based) and Ahmad &
Schreiber 2024 (services, markup-based, GTAP sector level), each converted
to nested_ces.py's convention, with a `source` column tracking provenance.
Re-run this script (`python build_combined_trade_elast.py`) if either
upstream source is revised; the output CSV is what `load_combined_trade_elast()`
in icio_to_haio.py actually reads.

CONVENTION: nested_ces.py needs trade_elast such that trade_elast + 1 equals
the CES elasticity of substitution across origin countries (see
icio_to_haio.py's module docstring). The two sources report fundamentally
different objects and need DIFFERENT transforms:
  - Fontagne's epsilon_icio is a NEGATIVE trade-cost elasticity
    (epsilon = -(sigma-1) is the standard convention), so
    trade_elast = sigma - 1 = -epsilon_icio.
  - Ahmad & Schreiber's EOS is already a directly-estimated CES elasticity
    of substitution (sigma itself, positive, from a markup/monopolistic-
    competition identification), so trade_elast = EOS - 1.
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

Sectors with NO estimate from either source, left out entirely (not
zero-filled, not defaulted): A02 (forestry), B05 (coal mining -- also
happens to be Saudi Arabia's real zero-activity sector from the README's
real-data test), B09 (mining support services), D (electricity/gas/steam
supply), O (public administration -- arguably non-traded by nature, so
this omission may simply be correct rather than a gap). Sector T is
excluded everywhere in this project regardless (see icio_to_haio.py).
"""
import pandas as pd

FONTAGNE_CSV = 'fontagne_icio_trade_elast.csv'
OUTPUT_CSV = 'combined_trade_elast.csv'

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

    out = pd.DataFrame(rows).sort_values('icio_sector').reset_index(drop=True)
    out.to_csv(OUTPUT_CSV, index=False)
    print(f'Wrote {len(out)} sectors to {OUTPUT_CSV}')
    print(f"  from Fontagne: {(out['source']=='fontagne_guimbard_orefice_2022').sum()}")
    print(f"  from Ahmad-Schreiber/GTAP: {(out['source']=='ahmad_schreiber_2024').sum()}")
    return out


if __name__ == '__main__':
    build()
