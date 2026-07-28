import os
from pathlib import Path
import pytest
from ase.db import connect
from ase.io import read
from match import Matcher, match_row

pytestmark = pytest.mark.skipif(
    os.environ.get("RUN_DATA_TESTS") != "1",
    reason="opt-in: runs real MillerSearch scans over all test-mother.db rows (~1 min); set RUN_DATA_TESTS=1 to run",
)

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


def golden_tuples(golden_db, cod_id):
    return sorted(
        tuple(row.key_value_pairs[k] for k in MATCH_KEYS)
        for row in golden_db.select(cod_id=cod_id)
    )


def test_match_row_matches_golden_fixture_for_every_row(mother_db, golden_db, matcher):
    for row in mother_db.select():
        _, _, results = match_row(row, matcher)

        actual = sorted(result_to_tuple(r) for r in results)
        expected = golden_tuples(golden_db, row.cod_id)

        assert actual == expected
