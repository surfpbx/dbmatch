from ase import Atoms
from ase.db import connect

from db_ogre_match.match import read_csv_rows
from db_ogre_match.utils import csv_to_db, db_to_csv


def _rows_as_sets(csv_path):
    return {frozenset(row.items()) for row in read_csv_rows(csv_path)}


def test_db_to_csv_writes_every_rows_key_value_pairs(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('mother.db')
    db.write(Atoms('H'), cod_id=1, reduced_formula='H')
    db.write(Atoms('He'), cod_id=2, reduced_formula='He')

    db_to_csv('mother.db', 'mother.csv')

    rows = read_csv_rows('mother.csv')
    assert sorted(row['cod_id'] for row in rows) == [1, 2]
    assert {row['cod_id']: row['reduced_formula'] for row in rows} == {1: 'H', 2: 'He'}


def test_db_to_csv_writes_no_metadata_sidecar_when_db_has_none(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('mother.db')
    db.write(Atoms('H'), cod_id=1)

    db_to_csv('mother.db', 'mother.csv')

    assert not (tmp_path / '.mother.csv.meta.json').exists()


def test_db_to_csv_carries_over_a_non_empty_metadata(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('scores.db')
    db.write(Atoms('H'), cod_id=1)
    db.metadata = {'matches_db': 'matches.csv', 'substrate': '/some/substrate.cif'}

    db_to_csv('scores.db', 'scores.csv')

    from db_ogre_match.match import read_match_metadata
    assert read_match_metadata('scores.csv') == {
        'matches_db': 'matches.csv', 'substrate': '/some/substrate.cif',
    }


def test_csv_to_db_round_trips_key_value_pairs_with_a_placeholder_atom(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('mother.db')
    db.write(Atoms('H'), cod_id=1, reduced_formula='H')
    db.write(Atoms('He'), cod_id=2, reduced_formula='He')
    db_to_csv('mother.db', 'mother.csv')

    csv_to_db('mother.csv', 'mother2.db')

    fresh = connect('mother2.db')
    rows = {row.cod_id: row for row in fresh.select()}
    assert set(rows) == {1, 2}
    assert rows[1].reduced_formula == 'H'
    assert rows[1].symbols == ['H']  # placeholder Atoms('H'), not real structure


def test_csv_to_db_restores_a_metadata_sidecar(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from db_ogre_match.match import write_match_metadata
    db = connect('mother.db')
    db.write(Atoms('H'), cod_id=1)
    db_to_csv('mother.db', 'mother.csv')
    write_match_metadata('mother.csv', mother_db='mother.db', substrate='sub.cif')

    csv_to_db('mother.csv', 'mother2.db')

    newdb = connect('mother2.db')
    newdb.count()  # metadata is only readable after some query
    assert newdb.metadata == {
        'mother_db': str((tmp_path / 'mother.db').resolve()),
        'substrate': str((tmp_path / 'sub.cif').resolve()),
    }


def test_csv_to_db_drops_reserved_id_and_formula_keys(tmp_path, monkeypatch):
    """read_csv_rows injects its own 'id' (row position); a real matches.csv
    materialized through score.py can also carry a 'formula' column. Both
    are ase db reserved keys -- csv_to_db must strip them before writing,
    or newdb.write(**kvp) raises."""
    monkeypatch.chdir(tmp_path)
    with open('rows.csv', 'w') as f:
        f.write('cod_id,formula\n1,H2O\n')

    csv_to_db('rows.csv', 'rows.db')  # must not raise

    row = connect('rows.db').get(cod_id=1)
    assert row.cod_id == 1


def test_csv_to_db_and_db_to_csv_round_trip_end_to_end(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('mother.db')
    db.write(Atoms('H'), cod_id=1, reduced_formula='H', sg_number=225)
    db.write(Atoms('He'), cod_id=2, reduced_formula='He', sg_number=194)

    db_to_csv('mother.db', 'round1.csv')
    csv_to_db('round1.csv', 'round1.db')
    db_to_csv('round1.db', 'round2.csv')

    assert _rows_as_sets('round1.csv') == _rows_as_sets('round2.csv')
