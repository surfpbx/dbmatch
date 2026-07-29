import pytest
from db_ogre_match import config
from db_ogre_match.score import build_tier_values, comp_score, geom_score_for_material, sg_score


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
