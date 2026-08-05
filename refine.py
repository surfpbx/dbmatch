"""Energetically refine every material picked out of the scores database by
a selection string (e.g. 'total_score>0.7'), running OgreInterface's ionic
surface-matching optimization (in-plane scan, then an interfacial-distance
scan) on every termination combination of every match score.py selected for
it, without cluttering the working directory: each material lands under its
own

    <reduced_formula>-<cod_id>/
        interfaces.db                      # one row per termination combo:
                                            #   the optimized structure, plus
                                            #   its energies, identifying
                                            #   keys (cod_id, match_id,
                                            #   sub_hkl/flm_hkl, termination
                                            #   indices), and subfolder (the
                                            #   combo's own plot directory,
                                            #   relative to this root)
        <match_id>_sub<sub_hkl>_flm<flm_hkl>/
            interface_view.png             # one plot per match -- the
                                            # lattice-registry geometry is
                                            # identical across every
                                            # termination combo of the same
                                            # match, only exposed atomic
                                            # layers differ
            sub<i>_flm<j>/
                PES.png                     # in-plane 2D scan
                z_shift.png                 # interfacial-distance scan

This deliberately never writes a POSCAR/plot per combination outside of
PES.png/z_shift.png -- since every combination is already reconstructible on
demand from (cod_id, match_id, sub_termination, film_termination) via
ogre_custom.interface_from_row, the optimized structure landing in
interfaces.db is enough to regenerate anything else later.
"""

import argparse
import os

import numpy as np
from ase.db import connect
from ase.io import read
from pymatgen.io.ase import AseAtomsAdaptor
from OgreInterface.generate import InterfaceGenerator, SurfaceGenerator
from OgreInterface.surface_matching import IonicSurfaceMatcher

from db_ogre_match import config
from db_ogre_match.ogre_custom import hkl_from_str, interface_from_row
from tqdm import tqdm


def _sanitize_hkl(hkl):
    return hkl.replace(',', '')


def _z_shift_range(matcher, lower=config.refine.z_shift_min, n_points=config.refine.z_shift_n_points):
    """Interfacial distances to scan. Wider than optimizePSO's own default
    z_bounds (see BaseSurfaceMatcher._get_max_z/optimizePSO) -- several
    combos were landing right on that narrower upper edge, meaning their
    true optimum was being cut off."""
    max_z = matcher._get_max_z()
    upper = max(5.0, 1.5 * max_z)
    return np.linspace(lower, upper, n_points)


def selected_matches(scores_db, cod_id, matches_db=None):
    """The matches-db rows score.py picked as best-per-facet for cod_id.

    If matches_db isn't given, it's read from scores_db's own metadata
    (recorded there by score_materials), so a scores db is enough on its
    own to find its way back to the matches it came from."""
    scores = connect(scores_db)
    scored = scores.get(cod_id=cod_id)
    match_ids = [int(i) for i in scored.selected_match_ids.split(';')]

    if matches_db is None:
        matches_db = scores.metadata['matches_db']

    matches = connect(matches_db)
    return [matches.get(id=i) for i in match_ids]


def substrate_from_scores_db(scores_db):
    """The fixed substrate structure path recorded in scores_db's own
    metadata -- forwarded there by score_materials from matches_db's
    metadata, itself recorded by match_database."""
    scores = connect(scores_db)
    scores.count()  # metadata is only readable after some query
    return scores.metadata['substrate']


def cod_ids_for_selection(scores_db, selection):
    """The cod_id of every row in scores_db matching selection -- any query
    string db.select() accepts (e.g. 'total_score>0.7', 'bravais=HEX',
    'cod_id=2300704') -- in the order scores_db yields them. scores_db has
    exactly one row per material (score_materials groups matches_db by
    cod_id before scoring), so no cod_id can repeat here."""
    scores = connect(scores_db)
    return [row.cod_id for row in scores.select(selection)]


def refine_material(
    selection,
    scores_db,
    substrate=None,
    matches_db=None,
    layers=config.refine.layers,
    vacuum=config.refine.vacuum,
    interfacial_distance=config.refine.interfacial_distance,
):
    """selection is any ase db query string db.select() accepts (e.g.
    'total_score>0.7', 'cod_id=2300704'); every cod_id it resolves to in
    scores_db is refined in turn. If substrate isn't given, it's read from
    scores_db's own metadata (see substrate_from_scores_db), the same way
    matches_db defaults from scores_db's metadata in selected_matches.

    Returns the list of root folders written, one per refined cod_id."""
    cod_ids = cod_ids_for_selection(scores_db, selection)

    if substrate is None:
        substrate = substrate_from_scores_db(scores_db)
    substrate_atoms = read(substrate)

    roots = []
    for cod_id in cod_ids:
        roots.append(
            _refine_one_material(
                cod_id, scores_db, substrate_atoms, matches_db,
                layers, vacuum, interfacial_distance,
            )
        )
    return roots


def _refine_one_material(
    cod_id, scores_db, substrate_atoms, matches_db,
    layers, vacuum, interfacial_distance,
):
    matches = selected_matches(scores_db, cod_id, matches_db)
    reduced_formula = matches[0].reduced_formula

    root = f'{reduced_formula}-{cod_id}'
    os.makedirs(root, exist_ok=True)
    # append=False: re-running refine on the same material should replace
    # interfaces.db, not accumulate duplicate rows alongside the old ones
    results_db = connect(os.path.join(root, 'interfaces.db'), append=False)

    for match in tqdm(matches):
        film_atoms = match.toatoms()

        # match.py's own Matcher scans with refine_structure=False -- must
        # match here, or sub_hkl/flm_hkl silently mean different facets.
        subs = SurfaceGenerator(
            bulk=substrate_atoms,
            miller_index=hkl_from_str(match.sub_hkl),
            layers=layers,
            vacuum=vacuum,
            refine_structure=False,
        )
        films = SurfaceGenerator(
            bulk=film_atoms,
            miller_index=hkl_from_str(match.flm_hkl),
            layers=layers,
            vacuum=vacuum,
            refine_structure=False,
        )

        match_dir = os.path.join(
            root,
            f'{match.id}_sub{_sanitize_hkl(match.sub_hkl)}_flm{_sanitize_hkl(match.flm_hkl)}',
        )
        os.makedirs(match_dir, exist_ok=True)

        print(
            f'match id={match.id}: {len(subs)} substrate termination(s), '
            f'{len(films)} film termination(s)'
        )

        view_saved = False

        for si, sub in enumerate(subs):
            for fi, film in enumerate(films):
                interface_generator = InterfaceGenerator(
                    substrate=sub,
                    film=film,
                    max_strain=config.match.max_strain,
                    max_area_mismatch=config.match.max_area_mismatch,
                    max_area=config.match.max_area,
                    interfacial_distance=interfacial_distance,
                    vacuum=40,
                    verbose=False,
                )
                interface = interface_from_row(interface_generator, match, sub, film)

                if not view_saved:
                    interface.plot_interface(
                        output=os.path.join(match_dir, 'interface_view.png')
                    )
                    view_saved = True

                combo_dir = os.path.join(match_dir, f'sub{si}_flm{fi}')
                os.makedirs(combo_dir, exist_ok=True)

                matcher = IonicSurfaceMatcher(interface=interface, verbose=False)

                matcher.run_surface_matching(output=os.path.join(combo_dir, 'PES.png'))
                matcher.get_optimized_structure()

                matcher.run_z_shift(
                    interfacial_distances=_z_shift_range(matcher),
                    output=os.path.join(combo_dir, 'z_shift.png'),
                )
                matcher.get_optimized_structure()

                adhesion_energy, interface_energy = matcher.get_current_energy()

                print(
                    f'  sub{si}/flm{fi}: adhesion={adhesion_energy:.4f} eV/A^2, '
                    f'interface={interface_energy:.4f} eV/A^2'
                )

                atoms = AseAtomsAdaptor.get_atoms(matcher.iface)
                results_db.write(
                    atoms,
                    cod_id=cod_id,
                    reduced_formula=reduced_formula,
                    match_id=match.id,
                    sub_hkl=match.sub_hkl,
                    flm_hkl=match.flm_hkl,
                    sub_termination=si,
                    film_termination=fi,
                    adhesion_energy=float(adhesion_energy),
                    interface_energy=float(interface_energy),
                    interfacial_distance=float(interface.interfacial_distance),
                    subfolder=os.path.join(os.path.basename(match_dir), f'sub{si}_flm{fi}'),
                )

    return root


def add_arguments(parser):
    """Add refine_material's CLI arguments to parser (shared by this file's
    own __main__ block and by cli.py's `dbm refine` subcommand)."""
    parser.add_argument(
        'selection',
        help=(
            "an ase db selection picking which cod_ids to refine, e.g. "
            "'cod_id=2300704' or 'total_score>0.7' -- any query "
            "db.select() accepts, resolved against the scores database"
        )
    )
    parser.add_argument(
        '--scores-db', default=os.path.join('db', config.score.output_db),
        help='path to the scores database (default: %(default)s)'
    )
    parser.add_argument(
        '--substrate', default=None,
        help=(
            "path to the substrate structure file (default: read from the "
            "scores database's own metadata, as recorded there by score.py)"
        )
    )
    parser.add_argument(
        '--matches-db', default=None,
        help=(
            "path to the matches database (default: read from the scores "
            "database's own metadata, as recorded there by score.py)"
        )
    )
    parser.add_argument('--layers', type=int, default=config.refine.layers)
    parser.add_argument('--vacuum', type=float, default=config.refine.vacuum)
    parser.add_argument(
        '--interfacial-distance', type=float, default=config.refine.interfacial_distance
    )
    return parser


def main(args):
    """Run refine_material from a parsed add_arguments() namespace -- shared
    by this file's own __main__ block and by cli.py's `dbm refine`."""
    refine_material(
        args.selection,
        args.scores_db,
        args.substrate,
        args.matches_db,
        layers=args.layers,
        vacuum=args.vacuum,
        interfacial_distance=args.interfacial_distance,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    main(parser.parse_args())
