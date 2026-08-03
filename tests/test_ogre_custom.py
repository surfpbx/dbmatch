from pathlib import Path

import numpy as np
import pytest
from ase.db import connect
from ase.io import read
from OgreInterface.generate import InterfaceGenerator, SurfaceGenerator

from db_ogre_match import config
from db_ogre_match.ogre_custom import (
    _align_to_x_axis,
    _embed_3x3,
    _inv2x2,
    hkl_from_str,
    interface_from_row,
    transform_from_str,
    transform_to_str,
)

DATA_DIR = Path(__file__).resolve().parent / "data"


def test_transform_to_str_and_from_str_round_trip():
    array = np.array([[1, 0], [-2, 3]])
    assert np.array_equal(transform_from_str(transform_to_str(array)), array)


def test_hkl_from_str():
    assert hkl_from_str('1,0,1') == [1, 0, 1]


def test_embed_3x3_places_2x2_block_and_identity_elsewhere():
    block = np.array([[2, 1], [0, 3]])
    full = _embed_3x3(block)
    assert full.shape == (3, 3)
    assert np.array_equal(full[:2, :2], block)
    assert np.array_equal(full[2], [0, 0, 1])
    assert np.array_equal(full[:2, 2], [0, 0])


def test_align_to_x_axis_rotates_first_vector_onto_x_axis():
    vectors = np.array([[3.0, 4.0, 0.0], [1.0, 1.0, 0.0]])
    rotation = _align_to_x_axis(vectors)
    aligned = vectors @ rotation

    assert aligned[0, 1] == pytest.approx(0.0, abs=1e-9)
    assert np.linalg.norm(aligned[0]) == pytest.approx(np.linalg.norm(vectors[0]))
    assert aligned[0, 2] == pytest.approx(vectors[0, 2])


def test_inv2x2_matches_numpy_inverse():
    for m in [np.array([[1.0, 2.0], [3.0, 4.0]]), np.array([[2.0, 0.0], [0.0, 5.0]])]:
        assert np.allclose(_inv2x2(m), np.linalg.inv(m))


def test_interface_from_row_matches_ogreinterfaces_own_search():
    """`refine.py` reconstructs each selected match's Interface via
    interface_from_row instead of re-running InterfaceGenerator's own
    ZurMcGill search -- this is a permanent regression test for the one
    assumption that whole approach depends on: that the closed-form
    reconstruction reproduces exactly what generate_interface(index=0)
    would have found."""
    substrate_atoms = read(str(DATA_DIR / 'CdSe.cif'))
    mother_db = connect(str(DATA_DIR / 'test-mother.db'))
    golden_db = connect(str(DATA_DIR / 'test-matches.db'))

    row = golden_db.get(cod_id=2300704, sub_hkl='0,1,0', flm_hkl='0,1,0')
    film_atoms = mother_db.get(cod_id=2300704).toatoms()

    subs = SurfaceGenerator(
        bulk=substrate_atoms, miller_index=hkl_from_str(row.sub_hkl),
        layers=2, vacuum=10.0, refine_structure=False,
    )
    films = SurfaceGenerator(
        bulk=film_atoms, miller_index=hkl_from_str(row.flm_hkl),
        layers=2, vacuum=10.0, refine_structure=False,
    )
    sub, film = subs[0], films[0]

    def build_generator():
        return InterfaceGenerator(
            substrate=sub, film=film,
            max_strain=0.05, max_area_mismatch=config.match.max_area_mismatch, max_area=100,
            interfacial_distance=2.0, vacuum=40, verbose=False,
        )

    reconstructed = interface_from_row(build_generator(), row, sub, film)
    searched = build_generator().generate_interface(interface_index=0)

    reconstructed_structure = reconstructed.get_interface(orthogonal=True)
    searched_structure = searched.get_interface(orthogonal=True)

    assert np.allclose(reconstructed_structure.lattice.matrix, searched_structure.lattice.matrix)
    assert [str(s) for s in reconstructed_structure.species] == [str(s) for s in searched_structure.species]
    assert np.allclose(
        reconstructed_structure.frac_coords, searched_structure.frac_coords, atol=1e-6
    )
