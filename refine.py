"""For one cod_id, run OgreInterface's ionic surface-matching (in-plane PES +
z-shift) energy optimization for every substrate/film termination
combination of every match score.py selected for that material, without
cluttering the working directory: everything lands under a single

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


def selected_matches(matches_db, scores_db, cod_id):
    """The matches-db rows score.py picked as best-per-facet for cod_id."""
    scores = connect(scores_db)
    scored = scores.get(cod_id=cod_id)
    match_ids = [int(i) for i in scored.selected_match_ids.split(';')]

    matches = connect(matches_db)
    return [matches.get(id=i) for i in match_ids]


def refine_material(
    substrate,
    cod_id,
    matches_db,
    scores_db,
    layers=config.refine.layers,
    vacuum=config.refine.vacuum,
    interfacial_distance=config.refine.interfacial_distance,
):
    matches = selected_matches(matches_db, scores_db, cod_id)
    reduced_formula = matches[0].reduced_formula

    root = f'{reduced_formula}-{cod_id}'
    os.makedirs(root, exist_ok=True)
    results_db = connect(os.path.join(root, 'interfaces.db'))

    substrate_atoms = read(substrate)

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('substrate', help='path to the substrate structure file')
    parser.add_argument('cod_id', type=int, help="a material's cod_id, as found in the scores database")
    parser.add_argument(
        '--matches-db', default=os.path.join('db', config.match.output_db),
        help='path to the matches database (default: %(default)s)'
    )
    parser.add_argument(
        '--scores-db', default=os.path.join('db', config.score.output_db),
        help='path to the scores database (default: %(default)s)'
    )
    parser.add_argument('--layers', type=int, default=config.refine.layers)
    parser.add_argument('--vacuum', type=float, default=config.refine.vacuum)
    parser.add_argument(
        '--interfacial-distance', type=float, default=config.refine.interfacial_distance
    )
    args = parser.parse_args()

    refine_material(
        args.substrate,
        args.cod_id,
        args.matches_db,
        args.scores_db,
        layers=args.layers,
        vacuum=args.vacuum,
        interfacial_distance=args.interfacial_distance,
    )
