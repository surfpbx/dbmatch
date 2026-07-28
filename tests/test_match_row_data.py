from pathlib import Path
import pytest
from ase.db import connect
from ase.io import read
from match import Matcher, match_row

DATA_DIR = Path(__file__).resolve().parent / "data"
MATCH_KEYS = [
    'area', 'strain', 'sub_hkl', 'flm_hkl',
    'sub_transform', 'flm_transform', 'n_sub_reps', 'n_flm_reps',
]


@pytest.fixture(scope="module")
def substrate():
    return read(str(DATA_DIR / "CdSe.cif"))


@pytest.fixture(scope="module")
def matcher(substrate):
    return Matcher(
        substrate=substrate,
        max_substrate_index=1,
        max_film_index=1,
        max_strain=0.05,
        max_area=500,
    )


@pytest.fixture(scope="module")
def mother_db():
    return connect(str(DATA_DIR / "test-mother.db"))


@pytest.fixture(scope="module")
def golden_db():
    return connect(str(DATA_DIR / "test-matches.db"))


def result_to_tuple(result):
    return tuple(result[k] for k in MATCH_KEYS)


def test_match_row_finds_matches_for_the_identical_substrate(mother_db, golden_db, matcher):
    row = next(mother_db.select(cod_id=9011664))  # CdSe itself

    _, _, results = match_row(row, matcher)

    actual = sorted(result_to_tuple(r) for r in results)
    expected = sorted(
        tuple(golden_row.key_value_pairs[k] for k in MATCH_KEYS)
        for golden_row in golden_db.select(cod_id=9011664)
    )

    assert actual == expected
