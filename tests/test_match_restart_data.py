import csv
from pathlib import Path

from db_ogre_match.match import checkpoint_path_for, match, read_csv_rows, write_checkpoint

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


def match_rows(csv_path):
    return sorted(
        tuple(row[k] for k in MATCH_COMPARISON_KEYS)
        for row in read_csv_rows(csv_path)
    )


def _write_partial_csv(csv_path, rows):
    """Write rows (as read_csv_rows would return them, including its own
    synthetic 'id' key) back out as a matches.csv-shaped file -- dropping
    'id', which is never one of match.py's own written columns."""
    fieldnames = [k for k in rows[0].keys() if k != 'id']
    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row[k] for k in fieldnames})


def test_restart_redoes_checkpointed_row_and_reproduces_a_full_run(tmp_path, monkeypatch):
    """Simulates a run that crashed right after mother-db row 5
    (cod_id=1010084) was checkpointed, leaving a stray duplicate match
    behind for it (standing in for a write that never made it cleanly to
    disk) and never reaching rows 6-10. Restart must clear row 5's
    matches before re-matching it, so the final CSV exactly reproduces an
    uninterrupted full run over tests/data/test-mother.db -- no missing
    or duplicate matches."""
    monkeypatch.chdir(tmp_path)

    match(output_csv='golden.csv', restart=False, **MATCH_KWARGS)
    golden = match_rows('golden.csv')

    # rows 1-4 survived the simulated crash intact; row 5 left a stray
    # duplicate behind; rows 6-10 were never reached.
    already_processed_cod_ids = {2300704, 9011664, 9008867, 9009122, 1010084}
    golden_rows = read_csv_rows('golden.csv')
    interrupted_rows = [row for row in golden_rows if row['cod_id'] in already_processed_cod_ids]
    row5 = next(row for row in golden_rows if row['cod_id'] == 1010084)
    interrupted_rows.append(row5)

    _write_partial_csv('interrupted.csv', interrupted_rows)
    write_checkpoint(checkpoint_path_for('interrupted.csv'), 5)

    match(output_csv='interrupted.csv', restart=True, **MATCH_KWARGS)

    assert match_rows('interrupted.csv') == golden
