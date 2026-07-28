from unittest.mock import MagicMock
from match import match_row


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
