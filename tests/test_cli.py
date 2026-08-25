import sys
from pathlib import Path

import pytest
from ase.db import connect

from db_ogre_match import cli

DATA_DIR = Path(__file__).resolve().parent / "data"
MATCH_KEYS = [
    'area', 'strain', 'sub_hkl', 'flm_hkl',
    'sub_transform', 'flm_transform', 'n_sub_reps', 'n_flm_reps',
]


def _match_tuple(row):
    return tuple(row.key_value_pairs[k] for k in MATCH_KEYS)


def test_dbm_help_lists_all_three_subcommands(capsys):
    with pytest.raises(SystemExit) as excinfo:
        cli.main(['--help'])
    assert excinfo.value.code == 0
    out = capsys.readouterr().out
    assert 'match' in out and 'score' in out and 'refine' in out


def test_dbm_with_no_args_exits_2():
    with pytest.raises(SystemExit) as excinfo:
        cli.main([])
    assert excinfo.value.code == 2


def test_dbm_match_help_does_not_import_refine(monkeypatch, capsys):
    # other test modules already import db_ogre_match.refine in this same
    # pytest session -- clear it first, or "not in sys.modules" would pass
    # trivially regardless of whether cli.py's lazy dispatch actually works.
    monkeypatch.delitem(sys.modules, 'db_ogre_match.refine', raising=False)

    with pytest.raises(SystemExit) as excinfo:
        cli.main(['match', '--help'])
    assert excinfo.value.code == 0
    assert 'db_ogre_match.refine' not in sys.modules
    assert 'substrate' in capsys.readouterr().out


def test_dbm_score_help_does_not_import_refine(monkeypatch, capsys):
    monkeypatch.delitem(sys.modules, 'db_ogre_match.refine', raising=False)

    with pytest.raises(SystemExit) as excinfo:
        cli.main(['score', '--help'])
    assert excinfo.value.code == 0
    assert 'db_ogre_match.refine' not in sys.modules
    assert 'matches_db' in capsys.readouterr().out


def test_dbm_match_cli_takes_substrate_then_mother_db(tmp_path, monkeypatch):
    """Runs `dbm match <substrate> <mother_db> ...` through cli.main's real
    argv parsing (not match(**kwargs)), so a positional-order mistake in
    add_arguments/main would actually be caught here."""
    monkeypatch.chdir(tmp_path)

    cli.main([
        'match',
        str(DATA_DIR / 'CdSe.cif'),
        str(DATA_DIR / 'test-mother.db'),
        '--max-substrate-index', '1',
        '--max-film-index', '1',
        '--max-strain', '0.05',
        '--max-area', '100',
        '-o', 'matches',
    ])

    fresh_db = connect('matches.db')
    golden_db = connect(str(DATA_DIR / 'test-matches.db'))

    for cod_id in {row.cod_id for row in golden_db.select()}:
        actual = sorted(_match_tuple(r) for r in fresh_db.select(cod_id=cod_id))
        expected = sorted(_match_tuple(r) for r in golden_db.select(cod_id=cod_id))
        assert actual == expected
