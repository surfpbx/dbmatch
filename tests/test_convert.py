import pytest
from ase import Atoms
from ase.db import connect

from db_ogre_match.convert import convert
from db_ogre_match.match import read_csv_rows


def test_convert_csv_to_db(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    with open('rows.csv', 'w') as f:
        f.write('cod_id\n1\n2\n')

    convert('rows.csv', 'rows.db')

    assert sorted(row.cod_id for row in connect('rows.db').select()) == [1, 2]


def test_convert_db_to_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    db = connect('rows.db')
    db.write(Atoms('H'), cod_id=1)
    db.write(Atoms('H'), cod_id=2)

    convert('rows.db', 'rows.csv')

    assert sorted(row['cod_id'] for row in read_csv_rows('rows.csv')) == [1, 2]


@pytest.mark.parametrize('src, dst', [
    ('a.csv', 'b.csv'),   # same extension
    ('a.db', 'b.db'),     # same extension
    ('a.txt', 'b.db'),    # unrecognized source extension
    ('a.csv', 'b.txt'),   # unrecognized destination extension
    ('a', 'b'),           # no extensions at all
])
def test_convert_raises_for_unsupported_extension_pairs(tmp_path, monkeypatch, src, dst):
    monkeypatch.chdir(tmp_path)
    with pytest.raises(ValueError):
        convert(src, dst)


def test_convert_raises_for_a_mismatched_case_extension(tmp_path, monkeypatch):
    """Extensions are matched exactly, not case-insensitively -- ase.db's
    own connect() only recognizes a lowercase '.db', so silently accepting
    '.DB' here would just fail one step down instead."""
    monkeypatch.chdir(tmp_path)
    with open('rows.csv', 'w') as f:
        f.write('cod_id\n1\n')

    with pytest.raises(ValueError):
        convert('rows.csv', 'rows.DB')
