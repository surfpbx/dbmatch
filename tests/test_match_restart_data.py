from pathlib import Path

from ase.db import connect

from db_ogre_match.match import checkpoint_path_for, match, write_checkpoint

DATA_DIR = Path(__file__).resolve().parent / "data"

MATCH_COMPARISON_KEYS = [
    'cod_id', 'sub_hkl', 'flm_hkl', 'sub_transform', 'flm_transform',
    'area', 'strain', 'n_sub_reps', 'n_flm_reps',
]

MATCH_KWARGS = dict(
    mother_db=str(DATA_DIR / 'test-mother.db'),
    substrate=str(DATA_DIR / 'CdSe.cif'),
    max_substrate_index=1,
    max_film_index=1,
    max_strain=0.05,
    max_area=100,
)


def match_rows(db_path):
    return sorted(
        tuple(row.key_value_pairs[k] for k in MATCH_COMPARISON_KEYS)
        for row in connect(db_path).select()
    )


def test_restart_redoes_checkpointed_row_and_reproduces_a_full_run(tmp_path, monkeypatch):
    """Simulates a run that crashed right after mother-db row 5
    (cod_id=1010084) was checkpointed, leaving a stray duplicate match
    behind for it (standing in for a write that never made it cleanly to
    disk) and never reaching rows 6-10. Restart must clear row 5's
    matches before re-matching it, so the final db exactly reproduces an
    uninterrupted full run over tests/data/test-mother.db -- no missing
    or duplicate matches."""
    monkeypatch.chdir(tmp_path)

    match(output_db='golden.db', output_csv='golden.csv', restart=False, **MATCH_KWARGS)
    golden = match_rows('golden.db')

    # rows 1-4 survived the simulated crash intact; row 5 left a stray
    # duplicate behind; rows 6-10 were never reached.
    already_processed_cod_ids = {2300704, 9011664, 9008867, 9009122, 1010084}
    golden_db = connect('golden.db')
    interrupted_db = connect('interrupted.db')
    for row in golden_db.select():
        if row.cod_id in already_processed_cod_ids:
            interrupted_db.write(row.toatoms(), **row.key_value_pairs)
    row5 = next(golden_db.select(cod_id=1010084))
    interrupted_db.write(row5.toatoms(), **row5.key_value_pairs)

    write_checkpoint(checkpoint_path_for('interrupted.db'), 5)

    match(output_db='interrupted.db', output_csv='interrupted.csv', restart=True, **MATCH_KWARGS)

    assert match_rows('interrupted.db') == golden
