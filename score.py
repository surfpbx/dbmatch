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
from db_ogre_match.match import CsvWriter, read_csv_rows, read_match_metadata
from db_ogre_match.ogre_custom import MATCH_RESULT_KEYS
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
          group_csv), is carried through into selected_match_ids.

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
        selected_match_ids      sorted list of the matches.csv row ids of the kept
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


def group_csv(csv_path):
    """
    Group matches.csv by cod_id.

    Returns a dict cod_id -> list of row dicts (each carrying its own
    'id', its row position in csv_path -- see read_csv_rows), one dict
    per facet-pairing match, in file order.
    """
    grouped = {}
    for row in read_csv_rows(csv_path):
        grouped.setdefault(row['cod_id'], []).append(row)
    return grouped


def _csv_row(mat):
    """mat, with every list-valued field (e.g. geom_score_for_material's own
    major_facets_scored/selected_match_ids, or any list a custom scoring_fn
    returns) serialised as a ';'-terminated string -- so comma-bearing values
    (e.g. hkl labels like '0,0,1') don't collide with the CSV delimiter, and
    so ase db's key_value_pairs (str/int/float/bool only) can store it too,
    via the same dict reused for newdb.write(**kvp). The trailing ';' (even
    for a single-element list) matters: without it, a one-element list of
    ids like [20007] would serialise to the bare string '20007', which
    ase.db's key_value_pairs check rejects as "a string but can be
    interpreted as int"."""
    return {
        k: (';'.join(str(x) for x in v) + ';' if v else '') if isinstance(v, list) else v
        for k, v in mat.items()
    }


def default_scoring_function(
    w_geom=config.score.w_geom, w_sg=config.score.w_sg, w_comp=config.score.w_comp,
    base_major_major=config.score.base_major_major,
    sub_major_minor_ratio=config.score.sub_major_minor_ratio,
    flm_major_minor_ratio=config.score.flm_major_minor_ratio,
    major_sub_facets=config.score.major_sub_facets,
    compatible_sg=config.score.sub_compatible_sg,
    compatible_elements=config.score.sub_compatible_elements,
):
    """Build score()'s own built-in scoring function: a closure over the
    tunables above, matching the scoring_fn contract score() expects --
    rows (one cod_id's matches group, as group_csv groups matches_db) in,
    a dict with at least total_score out. sg_number/reduced_formula are
    read off rows[0] since every match row already carries the material's
    own mother-db fields (match.py's store_matches copies them into every
    match row it writes), not just its own match-specific ones -- so this
    closure is self-sufficient given just rows, same as any custom
    scoring_fn would be."""
    tier_values = build_tier_values(
        base_major_major, sub_major_minor_ratio, flm_major_minor_ratio
    )

    def scoring_fn(rows):
        result = geom_score_for_material(rows, tier_values, major_sub_facets)
        result['sg_score'] = sg_score(rows[0]['sg_number'], compatible_sg)
        result['comp_score'] = comp_score(rows[0]['reduced_formula'], compatible_elements)
        result['total_score'] = (w_geom * result['geom_score']
                                 + w_sg * result['sg_score']
                                 + w_comp * result['comp_score'])
        return result

    return scoring_fn


def score(
    matches_db,
    output_basename,
    scoring_fn=None,
    mother_db=None,
    # only used to build the default scoring_fn (default_scoring_function)
    # when scoring_fn isn't given -- ignored otherwise
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
    matches_db:       path to matches.csv.
    output_basename:  basename for the output files, written to the
                       current directory -- results go to
                       <output_basename>.db (one row per scored material,
                       reusing that material's own atoms, re-read from
                       mother_db by cod_id) and <output_basename>.csv.
    mother_db:    path to the mother database, used to re-read each
                  material's atoms/formula by cod_id (one lookup per
                  material). Defaults to None, which reads it from
                  matches_db's own metadata sidecar (recorded there by
                  match -- see read_match_metadata); pass this explicitly
                  to avoid relying on that sidecar (e.g. if it could be
                  lost or moved) or to point at a different mother_db than
                  the one match originally used.
    scoring_fn:   rows (one cod_id's matches group, as group_csv groups
                  matches_db) -> a dict with at least a total_score key;
                  any other keys are carried into the output db/CSV
                  alongside the material's own mother-db fields. Defaults
                  to None, which builds default_scoring_function from this
                  call's own w_geom/w_sg/.../compatible_elements kwargs
                  (the CLI always uses this default -- an arbitrary Python
                  callable has no CLI encoding).

    Every material is scored first, then all of them are written (to the
    output db/CSV) sorted by total_score descending, best match first.
    Writes only; returns nothing.

    matches_db's absolute path is recorded in the output db's metadata (as
    'matches_db'), so that a material's selected_match_ids can later be
    resolved back to their source rows without having to separately track
    which matches_db a given scores_db came from. mother_db (whether given
    explicitly or read from matches_db's sidecar) and matches_db's own
    'substrate' metadata, if present, are forwarded to the output db's
    metadata too, so refine.py can find the substrate and re-read film
    atoms by cod_id directly from a scores db alone.
    """
    print(f'\nScoring materials from {matches_db}...')

    output_db = f'{output_basename}.db'
    output_csv = f'{output_basename}.csv'

    if scoring_fn is None:
        if not math.isclose(w_geom + w_sg + w_comp, 1.0, abs_tol=1e-9):
            warnings.warn(
                f'score weights do not sum to 1: w_geom={w_geom}, w_sg={w_sg}, '
                f'w_comp={w_comp} (sum={w_geom + w_sg + w_comp})'
            )
        scoring_fn = default_scoring_function(
            w_geom, w_sg, w_comp,
            base_major_major, sub_major_minor_ratio, flm_major_minor_ratio,
            major_sub_facets, compatible_sg, compatible_elements,
        )

    match_metadata = read_match_metadata(matches_db)
    if mother_db is None:
        mother_db = match_metadata.get('mother_db')
    mother_db_conn = connect(mother_db) if mother_db is not None else None

    grouped = group_csv(matches_db)

    csv_writer = CsvWriter(output_csv, append=False)
    newdb = connect(output_db, append=False)
    metadata = {'matches_db': os.path.abspath(matches_db)}
    if 'substrate' in match_metadata:
        metadata['substrate'] = match_metadata['substrate']
    if mother_db is not None:
        metadata['mother_db'] = os.path.abspath(mother_db)
    newdb.metadata = metadata

    materials = []
    for rows in tqdm(grouped.values()):
        # static per-material properties, taken from the first match row --
        # excluding that row's own per-match fields (MATCH_RESULT_KEYS),
        # which describe one particular facet pairing, not the material
        # itself, and belong only in matches_db
        mat = {k: v for k, v in rows[0].items() if k not in MATCH_RESULT_KEYS}

        result = scoring_fn(rows)
        if 'total_score' not in result:
            raise ValueError(
                f"scoring_fn must return a dict with a 'total_score' key (got: {sorted(result)})"
            )
        mat.update(result)
        materials.append(mat)

    materials.sort(key=lambda mat: mat['total_score'], reverse=True)

    for mat in materials:
        mother_row = mother_db_conn.get(cod_id=mat['cod_id'])

        row = _csv_row(mat)
        row['formula'] = mother_row.formula
        csv_writer.writerow(row)

        kvp = {k: v for k, v in row.items() if k not in ('id', 'formula')}
        newdb.write(mother_row.toatoms(), **kvp)

    csv_writer.file.close()

    print(
        f'Scoring finished. Scored {len(grouped)} materials, '
        f'written to {output_db} and {output_csv}.'
    )


def add_arguments(parser):
    """Add score's CLI arguments to parser (shared by this file's own
    __main__ block and by cli.py's `dbm score` subcommand)."""
    parser.add_argument(
        'matches_db',
        help='path to the matches CSV (written by match.py)'
    )
    parser.add_argument(
        '-o', '--output', default=config.score.output_basename,
        help=(
            'basename for the output files -- results are written to '
            '<output>.db and <output>.csv (default: %(default)s)'
        )
    )
    parser.add_argument(
        '--mother-db', default=None,
        help=(
            "path to the mother database, used to re-read each material's "
            "atoms/formula by cod_id (default: read from matches_db's own "
            "metadata sidecar, as recorded there by match.py)"
        )
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
    return parser


def main(args):
    """Run score from a parsed add_arguments() namespace -- shared by this
    file's own __main__ block and by cli.py's `dbm score`."""
    score(
        args.matches_db,
        args.output,
        mother_db=args.mother_db,
        w_geom=args.w_geom,
        w_sg=args.w_sg,
        w_comp=args.w_comp,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    main(parser.parse_args())
