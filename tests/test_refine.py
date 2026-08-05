import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from pymatgen.core import Lattice, Structure
from pymatgen.core.periodic_table import Element

from db_ogre_match.refine import (
    _contact_distance,
    _sanitize_hkl,
    _z_shift_range,
    cod_ids_for_selection,
    selected_matches,
    substrate_from_scores_db,
)


def test_sanitize_hkl_strips_commas():
    assert _sanitize_hkl('1,1,0') == '110'
    assert _sanitize_hkl('0,0,1') == '001'


def _structure_at_z(species_zs):
    """A pymatgen Structure with one atom per (symbol, z) pair, all at
    (x=5, y=5) in a large cubic box, for testing z-coordinate-picking logic
    in isolation from real slab geometry."""
    lattice = Lattice.cubic(20.0)
    species = [symbol for symbol, _ in species_zs]
    coords = [[5.0, 5.0, z] for _, z in species_zs]
    return Structure(lattice, species, coords, coords_are_cartesian=True)


def test_contact_distance_picks_substrate_top_atom_and_film_bottom_atom():
    """D should come from the substrate's topmost atom (Se, not Cd -- Cd is
    lower) and the film's bottommost atom (Te, not Mn -- Mn is higher),
    regardless of each atom's position in its own species list."""
    sub_structure = _structure_at_z([('Se', 2.0), ('Cd', 0.0)])
    film_structure = _structure_at_z([('Mn', 7.0), ('Te', 5.0)])

    D = _contact_distance(sub_structure, film_structure)

    expected = float(Element('Se').average_ionic_radius) + float(Element('Te').average_ionic_radius)
    assert D == pytest.approx(expected)


def test_z_shift_range_scans_half_to_twice_the_contact_distance():
    result = _z_shift_range(D=2.0, n_points=5)
    assert np.array_equal(result, np.linspace(1.0, 4.0, 5))


def test_z_shift_range_respects_custom_n_points():
    result = _z_shift_range(D=2.0, n_points=7)
    assert len(result) == 7


def test_selected_matches_returns_rows_in_selected_match_ids_order(tmp_path):
    matches_db = connect(str(tmp_path / 'matches.db'))
    matches_db.write(Atoms('H'), tag='a')  # id=1
    matches_db.write(Atoms('H'), tag='b')  # id=2
    matches_db.write(Atoms('H'), tag='c')  # id=3

    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='3;1')

    result = selected_matches(str(tmp_path / 'scores.db'), 42, str(tmp_path / 'matches.db'))

    assert [row.tag for row in result] == ['c', 'a']


def test_selected_matches_reads_matches_db_from_scores_db_metadata_by_default(tmp_path):
    matches_db = connect(str(tmp_path / 'matches.db'))
    matches_db.write(Atoms('H'), tag='a')  # id=1
    matches_db.write(Atoms('H'), tag='b')  # id=2

    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='1;2')
    scores_db.metadata = {'matches_db': str(tmp_path / 'matches.db')}

    result = selected_matches(str(tmp_path / 'scores.db'), 42)

    assert [row.tag for row in result] == ['a', 'b']


def test_substrate_from_scores_db_reads_metadata(tmp_path):
    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='1;2')
    scores_db.metadata = {'matches_db': 'irrelevant.db', 'substrate': '/some/substrate.cif'}

    assert substrate_from_scores_db(str(tmp_path / 'scores.db')) == '/some/substrate.cif'


def test_cod_ids_for_selection_filters_by_query(tmp_path):
    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=1, total_score=1.0)
    scores_db.write(Atoms('H'), cod_id=2, total_score=0.5)
    scores_db.write(Atoms('H'), cod_id=3, total_score=0.2)

    assert cod_ids_for_selection(str(tmp_path / 'scores.db'), 'total_score>0.4') == [1, 2]


def test_cod_ids_for_selection_with_single_cod_id_query(tmp_path):
    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=1, total_score=1.0)
    scores_db.write(Atoms('H'), cod_id=2, total_score=0.5)

    assert cod_ids_for_selection(str(tmp_path / 'scores.db'), 'cod_id=2') == [2]


def test_selected_matches_propagates_keyerror_for_unknown_cod_id(tmp_path):
    matches_db = connect(str(tmp_path / 'matches.db'))
    matches_db.write(Atoms('H'))

    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='1;2')

    with pytest.raises(KeyError):
        selected_matches(str(tmp_path / 'scores.db'), 999, str(tmp_path / 'matches.db'))
