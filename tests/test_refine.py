import numpy as np
import pytest
from ase import Atoms
from ase.db import connect
from unittest.mock import MagicMock

from db_ogre_match import config
from db_ogre_match.refine import _sanitize_hkl, _z_shift_range, selected_matches


def test_sanitize_hkl_strips_commas():
    assert _sanitize_hkl('1,1,0') == '110'
    assert _sanitize_hkl('0,0,1') == '001'


def test_z_shift_range_floors_upper_bound_at_5_when_max_z_is_small():
    matcher = MagicMock(_get_max_z=MagicMock(return_value=1.0))
    result = _z_shift_range(matcher)
    expected = np.linspace(config.refine.z_shift_min, 5.0, config.refine.z_shift_n_points)
    assert np.array_equal(result, expected)


def test_z_shift_range_scales_with_max_z_when_larger_than_the_floor():
    matcher = MagicMock(_get_max_z=MagicMock(return_value=10.0))
    result = _z_shift_range(matcher)
    expected = np.linspace(config.refine.z_shift_min, 15.0, config.refine.z_shift_n_points)
    assert np.array_equal(result, expected)


def test_z_shift_range_respects_custom_lower_and_n_points():
    matcher = MagicMock(_get_max_z=MagicMock(return_value=1.0))
    result = _z_shift_range(matcher, lower=0.5, n_points=5)
    expected = np.linspace(0.5, 5.0, 5)
    assert np.array_equal(result, expected)


def test_selected_matches_returns_rows_in_selected_match_ids_order(tmp_path):
    matches_db = connect(str(tmp_path / 'matches.db'))
    matches_db.write(Atoms('H'), tag='a')  # id=1
    matches_db.write(Atoms('H'), tag='b')  # id=2
    matches_db.write(Atoms('H'), tag='c')  # id=3

    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='3;1')

    result = selected_matches(str(tmp_path / 'matches.db'), str(tmp_path / 'scores.db'), 42)

    assert [row.tag for row in result] == ['c', 'a']


def test_selected_matches_propagates_keyerror_for_unknown_cod_id(tmp_path):
    matches_db = connect(str(tmp_path / 'matches.db'))
    matches_db.write(Atoms('H'))

    scores_db = connect(str(tmp_path / 'scores.db'))
    scores_db.write(Atoms('H'), cod_id=42, selected_match_ids='1;2')

    with pytest.raises(KeyError):
        selected_matches(str(tmp_path / 'matches.db'), str(tmp_path / 'scores.db'), 999)
