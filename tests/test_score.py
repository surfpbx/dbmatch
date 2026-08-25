import csv

import pytest
from ase import Atoms
from ase.db import connect
from db_ogre_match import config
from db_ogre_match.match import write_match_metadata
from db_ogre_match.score import build_tier_values, comp_score, geom_score_for_material, score, sg_score


def match_row(sub_hkl, flm_hkl, bravais='HEX', n_sub_reps=1, n_flm_reps=1):
    return {
        'sub_hkl': sub_hkl,
        'flm_hkl': flm_hkl,
        'bravais': bravais,
        'n_sub_reps': n_sub_reps,
        'n_flm_reps': n_flm_reps,
    }


def total_score(rows, sg_number, formula):
    tier_values = build_tier_values()
    geom = geom_score_for_material(rows, tier_values)['geom_score']
    sg = sg_score(sg_number)
    comp = comp_score(formula)
    return (config.score.w_geom * geom
            + config.score.w_sg * sg
            + config.score.w_comp * comp)


def test_total_score_is_1_when_all_major_facets_hit_and_sg_and_comp_match():
    rows = [match_row(hkl, hkl) for hkl in config.score.major_sub_facets]

    score = total_score(rows, sg_number=186, formula='CdSe')

    assert score == pytest.approx(1.0)


def test_total_score_is_0_5_when_only_two_of_three_major_facets_hit():
    hkls = sorted(config.score.major_sub_facets)[:2]
    rows = [match_row(hkl, hkl) for hkl in hkls]

    score = total_score(rows, sg_number=1, formula='SiO2')

    assert score == pytest.approx(0.5)


def test_total_score_is_0_when_no_matches_and_sg_and_comp_dont_match():
    score = total_score([], sg_number=1, formula='SiO2')

    assert score == pytest.approx(0.0)


def _make_matches_csv(matches_rows, cod_ids=(1,)):
    """Write matches.csv (matches_rows: list of dicts) plus a mother.db
    with one Atoms('H') row per cod_id, and matches.csv's metadata
    sidecar pointing at it -- score() needs mother_db to fetch each
    material's atoms/formula, matches.csv itself carries neither."""
    mother_db = connect('mother.db')
    for cod_id in cod_ids:
        mother_db.write(Atoms('H'), cod_id=cod_id)

    with open('matches.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(matches_rows[0].keys()))
        writer.writeheader()
        writer.writerows(matches_rows)

    write_match_metadata('matches.csv', mother_db='mother.db', substrate='substrate.cif')


def test_score_accepts_a_custom_scoring_fn_and_serializes_its_list_valued_extras(tmp_path, monkeypatch):
    """scoring_fn gets the raw cod_id matches group and can return whatever
    extra keys it wants (here a list) -- score() must carry them into both
    the output db and CSV, ';'-joining list values the same way
    geom_score_for_material's own selected_match_ids/major_facets_scored
    already rely on."""
    monkeypatch.chdir(tmp_path)
    _make_matches_csv([{'cod_id': 1, 'reduced_formula': 'H'}])

    def custom_scoring_fn(rows):
        assert rows[0]['cod_id'] == 1
        return {'total_score': 0.5, 'my_tag': ['a', 'b']}

    score('matches.csv', output_basename='scores', scoring_fn=custom_scoring_fn)

    scored = connect('scores.db').get(cod_id=1)
    assert scored.total_score == pytest.approx(0.5)
    assert scored.my_tag == 'a;b;'
    assert 'a;b;' in open('scores.csv').read()


def test_score_writes_a_single_element_numeric_id_list_without_raising(tmp_path, monkeypatch):
    """A material whose scoring_fn returns a one-element list of ids (e.g.
    geom_score_for_material's own selected_match_ids, when only one match
    row was kept) must not serialise to a bare digit string like '20007' --
    ase.db's key_value_pairs check rejects strings that look like ints."""
    monkeypatch.chdir(tmp_path)
    _make_matches_csv([{'cod_id': 1, 'reduced_formula': 'H'}])

    def custom_scoring_fn(rows):
        return {'total_score': 0.5, 'selected_match_ids': [20007]}

    score('matches.csv', output_basename='scores', scoring_fn=custom_scoring_fn)

    scored = connect('scores.db').get(cod_id=1)
    assert scored.selected_match_ids == '20007;'


def test_score_mother_db_override_takes_precedence_over_metadata_sidecar(tmp_path, monkeypatch):
    """score() must not require the metadata sidecar to still exist/be
    correct -- an explicit mother_db overrides whatever matches.csv's own
    sidecar says, both for the atoms/formula lookup and for what gets
    forwarded into the output db's own metadata."""
    monkeypatch.chdir(tmp_path)
    _make_matches_csv([{'cod_id': 1, 'reduced_formula': 'H'}])  # sidecar points at mother.db (H)

    override_db = connect('mother_override.db')
    override_db.write(Atoms('He'), cod_id=1)

    score(
        'matches.csv', output_basename='scores', mother_db='mother_override.db',
        scoring_fn=lambda rows: {'total_score': 0.5},
    )

    scored_db = connect('scores.db')
    scored = scored_db.get(cod_id=1)
    assert scored.formula == 'He'
    scored_db.count()  # metadata is only readable after some query
    assert scored_db.metadata['mother_db'] == str((tmp_path / 'mother_override.db').resolve())


def test_score_raises_if_scoring_fn_omits_total_score(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    _make_matches_csv([{'cod_id': 1, 'reduced_formula': 'H'}])

    with pytest.raises(ValueError):
        score(
            'matches.csv', output_basename='scores',
            scoring_fn=lambda rows: {'not_total_score': 1.0},
        )
