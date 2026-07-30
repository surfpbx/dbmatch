"""
Assign a per-material score based on the geometric matches found in the input database.

The score has three components:

    total_score = w_geom * geom_score + w_sg * sg_score + w_comp * comp_score

  - geom_score        reps-aware, 2D facet-tier reward. For each distinct
                      substrate facet a material hits, the single best
                      (highest tier / interface-area) match is kept; these
                      per-facet winners are summed and capped at 1.0.
  - sg_score          binary bonus for a substrate-compatible space group.
  - comp_score        binary bonus for containing >=1 substrate-compatible
                      element.
"""

import argparse
import math
import os
import re
import warnings
from ase.db import connect
from db_ogre_match import config
from db_ogre_match.match import CsvWriter
from tqdm import tqdm


def build_tier_values(base_major_major=config.score.base_major_major,
                      sub_major_minor_ratio=config.score.sub_major_minor_ratio,
                      flm_major_minor_ratio=config.score.flm_major_minor_ratio):
    """
    Construct the 4-tier per-facet score table.

    base_major_major: score for one substrate MAJOR facet hit via a film
        MAJOR facet. Default 1/3 so that hitting all 3 major substrate
        facets, each via a major film facet at 1:1 repetition, gives
        geom_score == 1.0 exactly.
    sub_major_minor_ratio: how much more a substrate MAJOR facet is worth
        than a substrate MINOR facet (film tier fixed). Default 4.0.
    flm_major_minor_ratio: how much more a film MAJOR facet is worth than a
        film MINOR facet (substrate tier fixed). Default 2.0.

    Ordering (major-major > major-minor > minor-major > minor-minor) holds as
    long as sub_major_minor_ratio > flm_major_minor_ratio.
    """
    major_minor_sub = base_major_major / sub_major_minor_ratio
    return {
        ('major', 'major'): base_major_major,
        ('major', 'minor'): base_major_major / flm_major_minor_ratio,
        ('minor', 'major'): major_minor_sub,
        ('minor', 'minor'): major_minor_sub / flm_major_minor_ratio,
    }


def geom_score_for_material(rows, tier_values,
                            major_sub_facets=config.score.major_sub_facets):
    """
    Score one material and report per-facet diagnostics.

    rows: iterable of dicts with 'sub_hkl', 'flm_hkl', 'bravais', 'n_sub_reps',
          'n_flm_reps' for a single material (one entry per facet-pairing match).
          The film tier (major/minor) is derived from flm_hkl against
          config.score.major_facets_by_bravais[row['bravais']] -- the film
          material's own bravais class. An 'id' key, if present (as set by
          group_db), is carried through into selected_match_ids.

    For each distinct substrate facet the material hits, the single match row
    maximizing
        tier_value / (n_sub_reps * n_flm_reps)
    is kept ("best per facet"). The kept values are summed and capped.

    Returns a dict:
        geom_score              capped sum of the per-facet best values
        n_facets_scored         number of distinct facets kept
        n_major_facets_scored   how many of those are substrate major facets
        major_facets_scored     sorted list of the major facet hkl strings kept
        mean_n_sub_reps         mean n_sub_reps over the kept (best-per-facet) rows
        mean_n_flm_reps         mean n_flm_reps over the kept rows
        selected_match_ids      sorted list of the matches-db ids of the kept
                                 (best-per-facet) rows
    """
    best = {}
    for row in rows:
        is_major = row['sub_hkl'] in major_sub_facets
        bravais = row['bravais']
        major_flm_facets = config.score.major_facets_by_bravais[bravais]

        sub_tier = 'major' if is_major else 'minor'
        flm_tier = 'major' if row['flm_hkl'] in major_flm_facets else 'minor'
        val = tier_values[(sub_tier, flm_tier)]

        nsub = float(row['n_sub_reps'])
        nflm = float(row['n_flm_reps'])
        val /= (nsub * nflm)

        facet = row['sub_hkl']
        if facet not in best or val > best[facet][0]:
            best[facet] = (val, nsub, nflm, is_major, row.get('id'))

    hits = min(1.0, sum(v[0] for v in best.values()))

    n_facets = len(best)
    major_facets = sorted(f for f, v in best.items() if v[3])
    subs = [v[1] for v in best.values()]
    flms = [v[2] for v in best.values()]
    match_ids = sorted(v[4] for v in best.values() if v[4] is not None)

    return {
        'geom_score': hits,
        'n_facets_scored': n_facets,
        'n_major_facets_scored': len(major_facets),
        'major_facets_scored': major_facets,
        'mean_n_sub_reps': (sum(subs) / len(subs)) if subs else None,
        'mean_n_flm_reps': (sum(flms) / len(flms)) if flms else None,
        'selected_match_ids': match_ids,
    }


def sg_score(sg_number, compatible_sg=config.score.sub_compatible_sg):
    """1.0 if the space group is substrate-compatible, else 0.0."""
    return 1.0 if sg_number in compatible_sg else 0.0


def _formula_elements(formula):
    """Set of element symbols in a formula string."""
    return set(re.findall(r'[A-Z][a-z]?', formula))


def comp_score(formula, compatible_elements=config.score.sub_compatible_elements):
    """1.0 if the material contains >=1 substrate-compatible element, else 0.0."""
    return 1.0 if _formula_elements(formula) & compatible_elements else 0.0


def group_db(db):
    """
    Group an ASE .db of matches by cod_id.

    Returns a dict cod_id -> list of row dicts (key_value_pairs plus id and
    formula), one dict per facet-pairing match, in db order.
    """
    grouped = {}
    for r in db.select():
        d = dict(r.key_value_pairs)
        d['id'] = r.id
        d['formula'] = r.formula
        grouped.setdefault(d['cod_id'], []).append(d)
    return grouped


def _csv_row(mat):
    """mat, with major_facets_scored and selected_match_ids serialised as
    ';'-separated strings (e.g. '0,0,1;1,1,0') so the comma-bearing hkl
    labels don't collide with the CSV delimiter."""
    return {
        **mat,
        'major_facets_scored': ';'.join(mat['major_facets_scored']),
        'selected_match_ids': ';'.join(str(i) for i in mat['selected_match_ids']),
    }


def score_materials(
    matches_db,
    output_csv,
    output_db,
    # final weights
    w_geom=config.score.w_geom, w_sg=config.score.w_sg, w_comp=config.score.w_comp,
    # facet-tier params
    base_major_major=config.score.base_major_major,
    sub_major_minor_ratio=config.score.sub_major_minor_ratio,
    flm_major_minor_ratio=config.score.flm_major_minor_ratio,
    major_sub_facets=config.score.major_sub_facets,
    # space-group params
    compatible_sg=config.score.sub_compatible_sg,
    # chemistry-flag params
    compatible_elements=config.score.sub_compatible_elements,
):
    """
    matches_db:   path to a matches .db.
    output_csv:   CSV filename, written under 'csv/'.
    output_db:    db filename, written under 'db/'; one row per
                  scored material, reusing that material's own atoms.

    Each material is scored and written (to csv/output_csv and
    db/output_db) as soon as it's scored -- rows land in matches_db's own
    cod_id order, not sorted by total_score. Writes only; returns nothing.

    matches_db's absolute path is recorded in output_db's metadata (as
    'matches_db'), so that a material's selected_match_ids can later be
    resolved back to their source rows without having to separately track
    which matches_db a given scores_db came from.
    """
    print(f'\nScoring materials from {matches_db}...')

    if not math.isclose(w_geom + w_sg + w_comp, 1.0, abs_tol=1e-9):
        warnings.warn(
            f'score weights do not sum to 1: w_geom={w_geom}, w_sg={w_sg}, '
            f'w_comp={w_comp} (sum={w_geom + w_sg + w_comp})'
        )

    src_db = connect(matches_db)
    grouped = group_db(src_db)

    tier_values = build_tier_values(
        base_major_major,
        sub_major_minor_ratio,
        flm_major_minor_ratio
    )
    os.makedirs('csv', exist_ok=True)
    os.makedirs('db', exist_ok=True)
    csv_writer = CsvWriter(os.path.join('csv', output_csv), append=False)
    newdb = connect(os.path.join('db', output_db), append=False)
    newdb.metadata = {'matches_db': os.path.abspath(matches_db)}

    for rows in tqdm(grouped.values()):
        # static per-material properties, taken from the first match row
        mat = dict(rows[0])

        mat.update(
            geom_score_for_material(
                rows, tier_values, major_sub_facets
            )
        )
        mat['sg_score'] = sg_score(mat['sg_number'], compatible_sg)
        mat['comp_score'] = comp_score(
            mat['reduced_formula'], compatible_elements
        )
        mat['total_score'] = (w_geom * mat['geom_score']
                              + w_sg * mat['sg_score']
                              + w_comp * mat['comp_score'])
        row = _csv_row(mat)
        csv_writer.writerow(row)

        kvp = {k: v for k, v in row.items() if k not in ('id', 'formula')}
        newdb.write(src_db.get(id=mat['id']), **kvp)

    csv_writer.file.close()

    print(
        f'Scoring finished. Scored {len(grouped)} materials, '
        f"written to db/{output_db} and csv/{output_csv}."
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'matches_db',
        help='path to the matches database'
    )
    parser.add_argument(
        '-o', '--output-db', default=config.score.output_db,
        help="name of the output database file, written under 'db/' (default: %(default)s)"
    )
    parser.add_argument(
        '--output-csv', default=config.score.output_csv,
        help="name of the output CSV file, written under 'csv/' (default: %(default)s)"
    )
    parser.add_argument(
        '--w-geom', type=float, default=config.score.w_geom,
        help='weight of the geometric-match score component (default: %(default)s)'
    )
    parser.add_argument(
        '--w-sg', type=float, default=config.score.w_sg,
        help='weight of the space-group score component (default: %(default)s)'
    )
    parser.add_argument(
        '--w-comp', type=float, default=config.score.w_comp,
        help='weight of the chemistry-compatibility score component (default: %(default)s)'
    )
    args = parser.parse_args()

    score_materials(
        args.matches_db,
        args.output_csv,
        args.output_db,
        w_geom=args.w_geom,
        w_sg=args.w_sg,
        w_comp=args.w_comp,
    )
