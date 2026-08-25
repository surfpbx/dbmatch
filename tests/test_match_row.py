import csv
from unittest.mock import MagicMock

from db_ogre_match.match import (
    CsvWriter,
    clear_stale_matches,
    match_row,
    read_csv_rows,
    read_match_metadata,
    write_match_metadata,
)


def test_match_row_calls_matcher_with_atoms_and_returns_pair():
    atoms = object()
    row = MagicMock()
    row.toatoms.return_value = atoms
    row.key_value_pairs = {'cod_id': 123}

    results = ['result1', 'result2']
    matcher = MagicMock(return_value=results)

    out_kwp, out_results = match_row(row, matcher)

    assert out_kwp == {'cod_id': 123}
    assert out_results is results
    matcher.assert_called_once_with(atoms)


def _write_csv(path, rows):
    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def test_clear_stale_matches_deletes_only_rows_for_that_cod_id(tmp_path):
    path = tmp_path / 'matches.csv'
    _write_csv(path, [
        {'cod_id': '1', 'tag': 'a'},
        {'cod_id': '1', 'tag': 'b'},
        {'cod_id': '2', 'tag': 'c'},
    ])

    clear_stale_matches(str(path), cod_id=1)

    assert [row['cod_id'] for row in read_csv_rows(str(path))] == [2]


def test_clear_stale_matches_is_a_noop_when_no_matches_exist_for_cod_id(tmp_path):
    path = tmp_path / 'matches.csv'
    _write_csv(path, [{'cod_id': '2', 'tag': 'c'}])

    clear_stale_matches(str(path), cod_id=1)

    assert [row['cod_id'] for row in read_csv_rows(str(path))] == [2]


def test_clear_stale_matches_is_a_noop_when_the_file_does_not_exist(tmp_path):
    path = tmp_path / 'matches.csv'

    clear_stale_matches(str(path), cod_id=1)  # must not raise

    assert not path.exists()


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


def test_read_csv_rows_coerces_types_and_assigns_1_based_ids(tmp_path):
    path = tmp_path / 'matches.csv'
    _write_csv(path, [
        {'cod_id': '1000050', 'strain': '0.0366', 'sub_hkl': '0,1,0', 'has_temp': 'False'},
        {'cod_id': '1000051', 'strain': '0.01', 'sub_hkl': '1,0,0', 'has_temp': 'True'},
    ])

    rows = read_csv_rows(str(path))

    assert rows[0] == {
        'cod_id': 1000050, 'strain': 0.0366, 'sub_hkl': '0,1,0', 'has_temp': False, 'id': 1,
    }
    assert rows[1]['id'] == 2
    assert rows[1]['has_temp'] is True


def test_write_and_read_match_metadata_round_trip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_match_metadata('matches.csv', mother_db='mother.db', substrate='sub.cif')

    metadata = read_match_metadata('matches.csv')

    assert metadata == {
        'mother_db': str((tmp_path / 'mother.db').resolve()),
        'substrate': str((tmp_path / 'sub.cif').resolve()),
    }


def test_read_match_metadata_is_empty_when_no_sidecar_exists(tmp_path):
    assert read_match_metadata(str(tmp_path / 'matches.csv')) == {}
