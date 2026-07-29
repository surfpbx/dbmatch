"""
CdSe epitaxial-match candidate scoring
=======================================

Scores every material in a best-matches database by how well its low-index
facets can form a commensurate epitaxial interface with CdSe. The score has
three components:

    total_score = w_hits * hits_score + w_sg * sg_score + w_chem_flag * chem_flag_score

  - hits_score        reps-aware, 2D facet-tier reward. For each distinct CdSe
                      facet a material hits, the single best (highest
                      tier / interface-area) match is kept; these per-facet
                      winners are summed and capped at 1.0.
  - sg_score          binary bonus for a wurtzite-compatible space group.
  - chem_flag_score   binary bonus for containing >=1 CdSe-compatible element.

Final default weights: 0.75 / 0.15 / 0.10.

Usage
-----
    from cdse_scoring_clean import score_materials
    df = score_materials("best_matches_4_4.db")           # defaults == final scoring
    final_set = df[df.total_score > 0.4]

Or from the command line (writes the scored CSV):
    python cdse_scoring_clean.py best_matches_4_4.db scores.csv

Requires: ase, pandas, numpy   (no pymatgen -- no valence chemistry here)
"""

import re
import numpy as np
import pandas as pd
from ase.db import connect


# ---------------------------------------------------------------------
# Facet-hit component: 2D tier (CdSe-facet major/minor) x (film-facet major/minor)
# ---------------------------------------------------------------------
DEFAULT_MAJOR_CDSE_FACETS = {'0,0,1', '0,1,0', '1,1,0'}
DEFAULT_MINOR_CDSE_FACETS = {'1,0,1', '1,1,1'}


def build_tier_values(base_major_major=1 / 3, cdse_major_minor_ratio=4.0,
                      film_major_minor_ratio=2.0):
    """
    Construct the 4-tier per-facet score table.

    base_major_major: score for one CdSe MAJOR facet hit via a film MAJOR facet.
        Default 1/3 so that hitting all 3 major CdSe facets, each via a major
        film facet at 1:1 repetition, gives hits_score == 1.0 exactly.
    cdse_major_minor_ratio: how much more a CdSe MAJOR facet is worth than a
        CdSe MINOR facet (film tier fixed). Default 4.0.
    film_major_minor_ratio: how much more a film MAJOR facet is worth than a
        film MINOR facet (CdSe tier fixed). Default 2.0.

    Ordering (major-major > major-minor > minor-major > minor-minor) holds as
    long as cdse_major_minor_ratio > film_major_minor_ratio.
    """
    major_minor_cdse = base_major_major / cdse_major_minor_ratio
    return {
        ('major', 'major'): base_major_major,
        ('major', 'minor'): base_major_major / film_major_minor_ratio,
        ('minor', 'major'): major_minor_cdse,
        ('minor', 'minor'): major_minor_cdse / film_major_minor_ratio,
    }


def hits_score_for_material(rows, tier_values,
                            major_cdse_facets=DEFAULT_MAJOR_CDSE_FACETS,
                            cap=1.0, reps_alpha=1.0, use_reps=True):
    """
    Score one material and report per-facet diagnostics.

    rows: iterable of dicts with 'sub_hkl', 'with_major', 'n_sub_reps',
          'n_flm_reps' for a single material (one entry per facet-pairing match).

    For each distinct CdSe facet the material hits, the single match row
    maximizing
        tier_value / (n_sub_reps * n_flm_reps) ** reps_alpha
    is kept ("best per facet"). The kept values are summed and capped.

    Returns a dict:
        hits_score              capped sum of the per-facet best values
        n_facets_scored         number of distinct facets kept
        n_major_facets_scored   how many of those are CdSe major facets
        major_facets_scored     sorted list of the major facet hkl strings kept
        mean_n_sub_reps         mean n_sub_reps over the kept (best-per-facet) rows
        mean_n_flm_reps         mean n_flm_reps over the kept rows
    """
    # per facet -> (score_value, n_sub_reps, n_flm_reps, is_major_facet)
    best = {}
    for row in rows:
        is_major = row['sub_hkl'] in major_cdse_facets
        cdse_tier = 'major' if is_major else 'minor'
        film_tier = 'major' if row['with_major'] else 'minor'
        val = tier_values[(cdse_tier, film_tier)]
        nsub = float(row['n_sub_reps'])
        nflm = float(row['n_flm_reps'])
        if use_reps:
            val = val / ((nsub * nflm) ** reps_alpha)
        facet = row['sub_hkl']
        if facet not in best or val > best[facet][0]:
            best[facet] = (val, nsub, nflm, is_major)

    hits = min(cap, sum(v[0] for v in best.values()))

    n_facets = len(best)
    major_facets = sorted(f for f, v in best.items() if v[3])
    subs = [v[1] for v in best.values()]
    flms = [v[2] for v in best.values()]

    return {
        'hits_score': hits,
        'n_facets_scored': n_facets,
        'n_major_facets_scored': len(major_facets),
        'major_facets_scored': major_facets,
        'mean_n_sub_reps': (sum(subs) / len(subs)) if subs else np.nan,
        'mean_n_flm_reps': (sum(flms) / len(flms)) if flms else np.nan,
    }


# ---------------------------------------------------------------------
# Space-group component
# ---------------------------------------------------------------------
DEFAULT_WURTZITE_COMPATIBLE_SG = {186, 194, 216, 225}


def sg_score(sg_number, wurtzite_compatible=DEFAULT_WURTZITE_COMPATIBLE_SG):
    """1.0 if the space group is wurtzite-compatible, else 0.0."""
    return 1.0 if sg_number in wurtzite_compatible else 0.0


# ---------------------------------------------------------------------
# Chemistry-compatibility flag (binary presence test, no valence modelling)
# ---------------------------------------------------------------------
#   group 12 : Zn, Cd, Hg            (the "II" in II-VI)
#   group 13 : Al, Ga, In           (B, Tl excluded)
#   group 15 : N, P, As, Sb, Bi     (the "V" in III-V; N kept -- nitrides are wurtzite)
#   group 16 : S, Se, Te            (chalcogen anions; O EXCLUDED -- oxides match worst)
# Group 14 (C, Si, Ge, Sn) is deliberately EXCLUDED: it credited rare-earth
# silicide/germanide intermetallics with no genuine CdSe chemistry.
DEFAULT_CDSE_COMPATIBLE_ELEMENTS = {
    'Zn', 'Cd', 'Hg',
    'Al', 'Ga', 'In',
    'N', 'P', 'As', 'Sb', 'Bi',
    'S', 'Se', 'Te',
}


def _formula_elements(formula):
    """Set of element symbols in a formula string."""
    return set(re.findall(r'[A-Z][a-z]?', formula))


def chem_flag_score(formula, compatible_elements=DEFAULT_CDSE_COMPATIBLE_ELEMENTS):
    """1.0 if the material contains >=1 CdSe-compatible element, else 0.0."""
    return 1.0 if _formula_elements(formula) & compatible_elements else 0.0


# ---------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------
def load_best_matches(db_path):
    """Load an ASE .db of matches into a DataFrame (one row per facet-pairing)."""
    db = connect(db_path)
    rows = []
    for r in db.select():
        d = dict(r.key_value_pairs)
        d['id'] = r.id
        d['formula'] = r.formula
        rows.append(d)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------
MAJOR_FACET_JOIN = ';'   # separator for the major_facets_scored string
                         # (hkl values themselves contain commas, so not comma)


def score_materials(
    best_matches,
    # final weights
    w_hits=0.75, w_sg=0.15, w_chem_flag=0.10,
    # facet-tier params
    base_major_major=1 / 3, cdse_major_minor_ratio=4.0, film_major_minor_ratio=2.0,
    major_cdse_facets=DEFAULT_MAJOR_CDSE_FACETS,
    # repetition-compactness params
    reps_alpha=1.0, use_reps=True,
    # space-group params
    wurtzite_compatible=DEFAULT_WURTZITE_COMPATIBLE_SG,
    # chemistry-flag params
    compatible_elements=DEFAULT_CDSE_COMPATIBLE_ELEMENTS,
):
    """
    best_matches: path to a matches .db, OR an already-loaded DataFrame from
                  load_best_matches().

    Returns one row per unique material (cod_id) with the score components,
    the combined total_score, and the per-facet scoring diagnostics.
    """
    df_bm = load_best_matches(best_matches) if isinstance(best_matches, str) else best_matches.copy()
    tier_values = build_tier_values(base_major_major, cdse_major_minor_ratio, film_major_minor_ratio)

    # static per-material properties
    mat = df_bm.drop_duplicates('cod_id').set_index('cod_id').copy()

    # facet-hit component + diagnostics (one call per material)
    diag = (
        df_bm.groupby('cod_id')
        .apply(lambda g: hits_score_for_material(
            g.to_dict('records'), tier_values, major_cdse_facets,
            reps_alpha=reps_alpha, use_reps=use_reps))
        .apply(pd.Series)
    )
    mat = mat.join(diag)

    # space-group component
    mat['sg_score'] = mat['sg_number'].apply(lambda sg: sg_score(sg, wurtzite_compatible))

    # chemistry-flag component
    mat['chem_flag_score'] = mat['reduced_formula'].apply(
        lambda f: chem_flag_score(f, compatible_elements))

    # combined score
    mat['total_score'] = (w_hits * mat['hits_score']
                          + w_sg * mat['sg_score']
                          + w_chem_flag * mat['chem_flag_score'])

    return mat.reset_index().sort_values('total_score', ascending=False)


# columns written to CSV, in order
CSV_COLUMNS = [
    'cod_id', 'reduced_formula',
    'hits_score', 'sg_score', 'chem_flag_score', 'total_score',
    'n_major_facets_scored', 'major_facets_scored', 'n_facets_scored',
    'mean_n_sub_reps', 'mean_n_flm_reps',
]


def write_scores_csv(scored, path):
    """Write the requested scoring columns to CSV.

    The major_facets_scored list is serialised as a MAJOR_FACET_JOIN-separated
    string (e.g. '0,0,1;1,1,0') so the comma-bearing hkl labels don't collide
    with the CSV delimiter.
    """
    out = scored.copy()
    out['major_facets_scored'] = out['major_facets_scored'].apply(
        lambda lst: MAJOR_FACET_JOIN.join(lst) if isinstance(lst, (list, tuple)) else '')
    out[CSV_COLUMNS].to_csv(path, index=False)
    return path


if __name__ == '__main__':
    import sys
    db_path = sys.argv[1] if len(sys.argv) > 1 else 'best_matches_4_4.db'
    out_path = sys.argv[2] if len(sys.argv) > 2 else 'material_scores.csv'
    scored = score_materials(db_path)
    write_scores_csv(scored, out_path)
    print(f"Scored {len(scored)} materials -> {out_path}")
    print(scored[['reduced_formula', 'hits_score', 'sg_score', 'chem_flag_score',
                  'total_score', 'n_major_facets_scored', 'major_facets_scored']].head(15).to_string(index=False))
