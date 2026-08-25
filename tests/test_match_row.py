import csv
from unittest.mock import MagicMock

from ase import Atoms
from ase.db import connect

from db_ogre_match.match import CsvWriter, clear_stale_matches, match_row


def test_match_row_calls_matcher_with_atoms_and_returns_triplet():
    atoms = object()
    row = MagicMock()
    row.toatoms.return_value = atoms
    row.key_value_pairs = {'cod_id': 123}

    results = ['result1', 'result2']
    matcher = MagicMock(return_value=results)

    out_atoms, out_kwp, out_results = match_row(row, matcher)

    assert out_atoms is atoms
    assert out_kwp == {'cod_id': 123}
    assert out_results is results
    matcher.assert_called_once_with(atoms)


def test_clear_stale_matches_deletes_only_rows_for_that_cod_id(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    newdb = connect('matches.db')
    newdb.write(Atoms('H'), cod_id=1)
    newdb.write(Atoms('H'), cod_id=1)
    newdb.write(Atoms('H'), cod_id=2)

    clear_stale_matches(newdb, cod_id=1)

    assert [row.cod_id for row in newdb.select()] == [2]


def test_clear_stale_matches_is_a_noop_when_no_matches_exist_for_cod_id(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    newdb = connect('matches.db')
    newdb.write(Atoms('H'), cod_id=2)

    clear_stale_matches(newdb, cod_id=1)

    assert [row.cod_id for row in newdb.select()] == [2]


def _read_csv(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


def test_csv_writer_drops_a_later_rows_unknown_key(tmp_path):
    path = tmp_path / 'out.csv'
    writer = CsvWriter(str(path), append=False)
    writer.writerow({'a': 1, 'b': 2})
    writer.writerow({'a': 3, 'b': 4, 'c': 5})
    writer.file.close()

    rows = _read_csv(path)
    assert rows == [
        {'a': '1', 'b': '2'},
        {'a': '3', 'b': '4'},
    ]


def test_csv_writer_handles_a_missing_key_on_a_later_row(tmp_path):
    path = tmp_path / 'out.csv'
    writer = CsvWriter(str(path), append=False)
    writer.writerow({'a': 1, 'b': 2})
    writer.writerow({'a': 3})
    writer.file.close()

    rows = _read_csv(path)
    assert rows == [
        {'a': '1', 'b': '2'},
        {'a': '3', 'b': ''},
    ]
