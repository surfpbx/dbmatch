import pytest
from ase import Atoms
from ase.db import connect
from db_ogre_match import config
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


def test_score_accepts_a_custom_scoring_fn_and_serializes_its_list_valued_extras(tmp_path, monkeypatch):
    """scoring_fn gets the raw cod_id matches group and can return whatever
    extra keys it wants (here a list) -- score() must carry them into both
    output_db and output_csv, ';'-joining list values the same way
    geom_score_for_material's own selected_match_ids/major_facets_scored
    already rely on."""
    monkeypatch.chdir(tmp_path)
    matches_db = connect('matches.db')
    matches_db.write(Atoms('H'), cod_id=1, reduced_formula='H')

    def custom_scoring_fn(rows):
        assert rows[0]['cod_id'] == 1
        return {'total_score': 0.5, 'my_tag': ['a', 'b']}

    score('matches.db', output_csv='scores.csv', output_db='scores.db', scoring_fn=custom_scoring_fn)

    scored = connect('scores.db').get(cod_id=1)
    assert scored.total_score == pytest.approx(0.5)
    assert scored.my_tag == 'a;b'
    assert 'a;b' in open('scores.csv').read()


def test_score_raises_if_scoring_fn_omits_total_score(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    matches_db = connect('matches.db')
    matches_db.write(Atoms('H'), cod_id=1, reduced_formula='H')

    with pytest.raises(ValueError):
        score(
            'matches.db', output_csv='scores.csv', output_db='scores.db',
            scoring_fn=lambda rows: {'not_total_score': 1.0},
        )
