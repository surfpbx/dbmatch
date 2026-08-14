from pathlib import Path

from ase.db import connect

from db_ogre_match.match import match_database
from db_ogre_match.score import score_materials

DATA_DIR = Path(__file__).resolve().parent / "data"
SCORE_KEYS = [
    'geom_score', 'n_facets_scored', 'n_major_facets_scored',
    'major_facets_scored', 'mean_n_sub_reps', 'mean_n_flm_reps',
    'selected_match_ids', 'sg_score', 'comp_score', 'total_score',
]


def result_tuple(row):
    return tuple(row.key_value_pairs[k] for k in SCORE_KEYS)


def test_match_then_score_reproduces_golden_scores(tmp_path, monkeypatch):
    """Runs the real match_database -> score_materials workflow (the same
    two calls example/run_example.py makes) over tests/data/test-mother.db
    against CdSe.cif, and diffs the resulting scores db against
    tests/data/test-scores.db, material by material."""
    monkeypatch.chdir(tmp_path)

    match_database(
        mother_db=str(DATA_DIR / 'test-mother.db'),
        substrate=str(DATA_DIR / 'CdSe.cif'),
        max_substrate_index=1,
        max_film_index=1,
        max_strain=0.05,
        max_area=100,
        output_db='matches.db',
        output_csv='matches.csv',
        restart=False,
    )
    score_materials(
        matches_db='matches.db',
        output_csv='scores.csv',
        output_db='scores.db',
    )

    fresh_db = connect('scores.db')
    golden_db = connect(str(DATA_DIR / 'test-scores.db'))

    golden_rows = {row.cod_id: row for row in golden_db.select()}
    assert {row.cod_id for row in fresh_db.select()} == set(golden_rows)

    for row in fresh_db.select():
        assert result_tuple(row) == result_tuple(golden_rows[row.cod_id])

    # match_database's substrate path is forwarded through score_materials
    # into the scores db's own metadata, so refine.py can find it without
    # being told separately (see refine.substrate_from_scores_db).
    assert fresh_db.metadata['substrate'] == str((DATA_DIR / 'CdSe.cif').resolve())
