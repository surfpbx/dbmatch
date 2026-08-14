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
                                            #   indices, sub_species/
                                            #   film_species -- the element
                                            #   actually facing the
                                            #   interface on each side for
                                            #   this termination), and
                                            #   subfolder (the combo's own
                                            #   plot directory, relative to
                                            #   this root)
        results.txt                        # the same per-combo fields as
                                            # interfaces.db, as a plain-text
                                            # table for a quick human read
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
                interface.cif               # optimized structure, same as
                                             # this combo's interfaces.db row

interface.cif is written purely for convenience (a quick look with any CIF
viewer, no ase db query needed) -- the same structure is already stored in
interfaces.db, so it's not needed to regenerate anything.
"""

import argparse
import os
import subprocess
from unittest.mock import patch

import matplotlib.colors
import numpy as np
from ase.db import connect
from ase.io import read
from ase.visualize import view as ase_view
from pymatgen.io.ase import AseAtomsAdaptor
from OgreInterface.generate import InterfaceGenerator, SurfaceGenerator
from OgreInterface.surface_matching import IonicSurfaceMatcher
from OgreInterface.surface_matching import base_surface_matcher

from db_ogre_match import config
from db_ogre_match.ogre_custom import hkl_from_str, interface_from_row
from tqdm import tqdm


def _sanitize_hkl(hkl):
    return hkl.replace(',', '')


def _terminating_atoms(sub_structure, film_structure):
    """The two atoms that will actually face each other across the
    interface for this termination combo: the substrate's topmost atom and
    the film's bottommost atom (each structure in its own z-coordinate
    frame -- substrate slabs and film slabs are both generated
    independently before being merged into an interface, but merging only
    translates each rigid slab along z, so the atom picked out here is the
    same one that ends up facing the interface after merging). Returned as
    pymatgen sites, so callers can pull out whatever they need (ionic
    radius for _contact_distance, element symbol for display, ...) without
    re-deriving which atom is "at the surface" more than once."""
    sub_z = sub_structure.cart_coords[:, -1]
    film_z = film_structure.cart_coords[:, -1]

    return sub_structure[sub_z.argmax()], film_structure[film_z.argmin()]


def _contact_distance(sub_structure, film_structure, factor=config.refine.contact_distance_factor):
    """The natural ionic contact distance between the two atoms that will
    actually face each other across the interface (see _terminating_atoms).
    D is factor times the sum of their pymatgen Element.average_ionic_radius
    -- used both as the interface's starting interfacial_distance and as
    the basis for the z-shift scan range (see _z_shift_range). factor
    defaults to config.refine.contact_distance_factor, > 1 to push the raw
    ionic-radii sum (which tends to sit a bit closer than a real relaxed
    contact distance) outward."""
    sub_atom, film_atom = _terminating_atoms(sub_structure, film_structure)
    return factor * (float(sub_atom.specie.average_ionic_radius) + float(film_atom.specie.average_ionic_radius))


def _z_shift_range(D, n_points=config.refine.z_shift_n_points):
    """Interfacial distances to scan around the natural ionic contact
    distance D (see _contact_distance): (0.5*D, 2*D).

    Replaces OgreInterface's own _get_max_z bound (2x the largest covalent
    radius over every element in the whole interface, substrate and film
    alike) -- that heuristic ignored which atoms were actually at the
    interface and was cutting several combos' true PES optimum off at its
    narrow upper edge."""
    return np.linspace(0.5 * D, 2 * D, n_points)


def _run_surface_matching(matcher, output, bound=config.refine.pes_colormap_bound):
    """matcher.run_surface_matching(output=output), but with the PES plot's
    colorbar pinned to a fixed +/-bound eV/A^2 range instead of
    OgreInterface's own data-min/max auto-scaling -- a few outlier grid
    points (e.g. atoms overlapping at a bad registry) would otherwise wash
    out the color contrast of the physically relevant part of the scan.
    run_surface_matching's own public signature has no vmin/vmax parameter;
    matplotlib.colors.Normalize is the only thing base_surface_matcher calls
    to build that colorbar, and only inside _plot_heatmap (confirmed: it's
    not reused by run_z_shift's own line plot), so patching that one name
    for the duration of this call reproduces a vmin/vmax kwarg without
    touching OgreInterface's source."""
    def _fixed_normalize(*args, **kwargs):
        return matplotlib.colors.Normalize(vmin=-bound, vmax=bound, clip=True)

    with patch.object(base_surface_matcher, 'Normalize', _fixed_normalize):
        return matcher.run_surface_matching(output=output)


_RESULTS_TABLE_COLUMNS = [
    # (header, row key, alignment)
    ('Match ID', 'match_id', '>'),
    ('Substrate (hkl)', 'sub_hkl', '>'),
    ('Film (hkl)', 'flm_hkl', '>'),
    ('Termination', 'termination', '>'),
    ('Eadh (eV/Å²)', 'adhesion_energy', '>'),
    ('Eint (eV/Å²)', 'interface_energy', '>'),
    ('d_int (Å)', 'interfacial_distance', '>'),
    ('Subfolder', 'subfolder', '<'),
]


def _results_table(rows):
    """rows: a list of dicts, one per termination combo, with a string value
    for every key in _RESULTS_TABLE_COLUMNS -- rendered as a column-aligned
    plain-text table (written to <root>/results.txt by
    _refine_one_material)."""
    widths = [
        max([len(header)] + [len(row[key]) for row in rows])
        for header, key, _ in _RESULTS_TABLE_COLUMNS
    ]

    def fmt_row(values):
        return '  '.join(
            f'{value:{align}{width}}'
            for value, (_, _, align), width in zip(values, _RESULTS_TABLE_COLUMNS, widths)
        )

    lines = [
        fmt_row([header for header, _, _ in _RESULTS_TABLE_COLUMNS]),
        '  '.join('-' * width for width in widths),
    ]
    lines += [fmt_row([row[key] for _, key, _ in _RESULTS_TABLE_COLUMNS]) for row in rows]

    return '\n'.join(lines) + '\n'


def selected_matches(scores_db, cod_id, matches_db=None):
    """The matches-db rows score.py picked as best-per-facet for cod_id.

    If matches_db isn't given, it's read from scores_db's own metadata
    (recorded there by score), so a scores db is enough on its own to find
    its way back to the matches it came from."""
    scores = connect(scores_db)
    scored = scores.get(cod_id=cod_id)
    match_ids = [int(i) for i in scored.selected_match_ids.split(';')]

    if matches_db is None:
        matches_db = scores.metadata['matches_db']

    matches = connect(matches_db)
    return [matches.get(id=i) for i in match_ids]


def substrate_from_scores_db(scores_db):
    """The fixed substrate structure path recorded in scores_db's own
    metadata -- forwarded there by score from matches_db's metadata,
    itself recorded by match."""
    scores = connect(scores_db)
    scores.count()  # metadata is only readable after some query
    return scores.metadata['substrate']


def cod_ids_for_selection(scores_db, selection):
    """The cod_id of every row in scores_db matching selection -- any query
    string db.select() accepts (e.g. 'total_score>0.7', 'bravais=HEX',
    'cod_id=2300704') -- in the order scores_db yields them. scores_db has
    exactly one row per material (score groups matches_db by cod_id before
    scoring), so no cod_id can repeat here."""
    scores = connect(scores_db)
    return [row.cod_id for row in scores.select(selection)]


def refine(
    selection,
    scores_db,
    substrate=None,
    matches_db=None,
    layers=config.refine.layers,
    vacuum=config.refine.vacuum,
):
    """selection is any ase db query string db.select() accepts (e.g.
    'total_score>0.7', 'cod_id=2300704'); every cod_id it resolves to in
    scores_db is refined in turn. If substrate isn't given, it's read from
    scores_db's own metadata (see substrate_from_scores_db), the same way
    matches_db defaults from scores_db's metadata in selected_matches.

    Returns the list of root folders written, one per refined cod_id."""
    cod_ids = cod_ids_for_selection(scores_db, selection)
    print(f"\nRefining {len(cod_ids)} material(s) matching '{selection}' from {scores_db}...")

    if substrate is None:
        substrate = substrate_from_scores_db(scores_db)
    substrate_atoms = read(substrate)

    roots = []
    for cod_id in cod_ids:
        roots.append(
            _refine_one_material(
                cod_id, scores_db, substrate_atoms, matches_db, layers, vacuum,
            )
        )
    return roots


def _refine_one_material(
    cod_id, scores_db, substrate_atoms, matches_db, layers, vacuum,
):
    matches = selected_matches(scores_db, cod_id, matches_db)
    reduced_formula = matches[0].reduced_formula
    nmatches = len(matches)

    print(f'\n{reduced_formula}-{cod_id}: found {nmatches} selected match(es)')

    root = f'{reduced_formula}-{cod_id}'
    os.makedirs(root, exist_ok=True)
    # append=False: re-running refine on the same material should replace
    # interfaces.db, not accumulate duplicate rows alongside the old ones
    results_db = connect(os.path.join(root, 'interfaces.db'), append=False)
    results_rows = []

    for i, match in enumerate(matches):
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
            f'\n({i+1}/{nmatches}) ({_sanitize_hkl(match.sub_hkl)})║({_sanitize_hkl(match.flm_hkl)}): '
            f'{len(subs)} substrate termination(s), {len(films)} film termination(s)'
        )

        view_saved = False

        for si, sub in enumerate(subs):
            for fi, film in enumerate(films):
                # starting interfacial_distance from the same ionic-radii
                # criterion as the z-shift scan below -- computed per
                # termination combo, since different terminations expose
                # different atoms at the surface
                sub_surface = sub.get_surface(orthogonal=True)
                film_surface = film.get_surface(orthogonal=True)
                D = _contact_distance(sub_surface, film_surface)
                sub_atom, film_atom = _terminating_atoms(sub_surface, film_surface)
                sub_species = sub_atom.specie.symbol
                film_species = film_atom.specie.symbol

                interface_generator = InterfaceGenerator(
                    substrate=sub,
                    film=film,
                    max_strain=config.match.max_strain,
                    max_area_mismatch=config.match.max_area_mismatch,
                    max_area=config.match.max_area,
                    interfacial_distance=D,
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

                _run_surface_matching(matcher, output=os.path.join(combo_dir, 'PES.png'))
                matcher.get_optimized_structure()

                matcher.run_z_shift(
                    interfacial_distances=_z_shift_range(D),
                    output=os.path.join(combo_dir, 'z_shift.png'),
                )
                matcher.get_optimized_structure()

                adhesion_energy, interface_energy = matcher.get_current_energy()
                subfolder = os.path.join(os.path.basename(match_dir), f'sub{si}_flm{fi}')

                print(
                    f'  sub{si}({sub_species})/flm{fi}({film_species}):  '
                    f'Eadh= {adhesion_energy:7.4f} eV/Å²  '
                    f'Eint= {interface_energy:7.4f} eV/Å²  '
                    f'{os.path.join(root, subfolder)}'
                )

                atoms = AseAtomsAdaptor.get_atoms(matcher.iface)
                atoms.write(os.path.join(combo_dir, 'interface.cif'))

                results_db.write(
                    atoms,
                    cod_id=cod_id,
                    reduced_formula=reduced_formula,
                    match_id=match.id,
                    sub_hkl=match.sub_hkl,
                    flm_hkl=match.flm_hkl,
                    sub_termination=si,
                    film_termination=fi,
                    sub_species=sub_species,
                    film_species=film_species,
                    adhesion_energy=float(adhesion_energy),
                    interface_energy=float(interface_energy),
                    interfacial_distance=float(interface.interfacial_distance),
                    subfolder=subfolder,
                )

                results_rows.append({
                    'match_id': str(match.id),
                    'sub_hkl': f'({_sanitize_hkl(match.sub_hkl)})',
                    'flm_hkl': f'({_sanitize_hkl(match.flm_hkl)})',
                    'termination': f'{sub_species}/{film_species}',
                    'adhesion_energy': f'{adhesion_energy:.4f}',
                    'interface_energy': f'{interface_energy:.4f}',
                    'interfacial_distance': f'{interface.interfacial_distance:.4f}',
                    'subfolder': subfolder,
                })

    with open(os.path.join(root, 'results.txt'), 'w') as f:
        f.write(_results_table(results_rows))

    print(
        f'\n{reduced_formula}-{cod_id}: {len(results_rows)} interface(s) saved to '
        f'{os.path.join(root, "interfaces.db")}, summary written to '
        f'{os.path.join(root, "results.txt")}'
    )

    return root


def view_interface(interfaces_db, row_id):
    """Visualize one interfaces.db row's optimized structure (ase.visualize.view)
    and open its PES.png/z_shift.png plots (xdg-open) -- both found via that
    row's own subfolder, relative to interfaces_db's own directory (see
    _refine_one_material's subfolder=... kwarg)."""
    db = connect(interfaces_db)
    row = db.get(id=row_id)

    ase_view(row.toatoms())

    combo_dir = os.path.join(os.path.dirname(os.path.abspath(interfaces_db)), row.subfolder)
    for plot in ('PES.png', 'z_shift.png'):
        subprocess.Popen(['xdg-open', os.path.join(combo_dir, plot)])


def add_arguments(parser):
    """Add refine's CLI arguments to parser (shared by this file's own
    __main__ block and by cli.py's `dbm refine` subcommand)."""
    parser.add_argument(
        'selection', nargs='?', default=None,
        help=(
            "an ase db selection picking which cod_ids to refine, e.g. "
            "'cod_id=2300704' or 'total_score>0.7' -- any query "
            "db.select() accepts, resolved against the scores database. "
            "Not needed (and ignored) if --view is given."
        )
    )
    parser.add_argument(
        '-d', '--scores-db', default=config.score.output_db,
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
        '--view', type=int, default=None, metavar='ID',
        help=(
            "visualize one interfaces.db row instead of running refine: ID "
            "is the row's own ase db id. Opens its optimized structure via "
            "ase.visualize.view and its PES.png/z_shift.png via xdg-open."
        )
    )
    parser.add_argument(
        '--interfaces-db', default='interfaces.db',
        help=(
            "path to the interfaces database to read from with --view "
            "(default: %(default)s, i.e. run from inside the material's "
            "own <reduced_formula>-<cod_id>/ folder)"
        )
    )
    return parser


def main(args):
    """Run refine from a parsed add_arguments() namespace -- shared by this
    file's own __main__ block and by cli.py's `dbm refine`. If --view is
    given, visualizes that interfaces.db row instead and returns without
    running refine at all."""
    if args.view is not None:
        view_interface(args.interfaces_db, args.view)
        return

    if args.selection is None:
        raise SystemExit('refine: error: selection is required unless --view is given')

    refine(
        args.selection,
        args.scores_db,
        args.substrate,
        args.matches_db,
        layers=args.layers,
        vacuum=args.vacuum,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    main(parser.parse_args())
