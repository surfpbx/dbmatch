import pytest

from db_ogre_match import config
from db_ogre_match.config import apply_overrides, load_rc_file


@pytest.fixture
def restore_config():
    """apply_overrides mutates config.match/.score/.refine in place --
    shared module state other tests rely on -- so snapshot and restore
    every attribute touched by a test using this fixture."""
    snapshots = {
        name: dict(vars(namespace))
        for name, namespace in config._NAMESPACES.items()
    }
    yield
    for name, snapshot in snapshots.items():
        namespace = config._NAMESPACES[name]
        for key, value in snapshot.items():
            setattr(namespace, key, value)


def _write_rc(path, text):
    path.write_text(text)
    return str(path)


def test_load_rc_file_parses_scalars_as_python_literals(tmp_path):
    path = _write_rc(tmp_path / 'dbm.config', """
match.max_strain : 0.2
match.max_area : 250
match.output_csv : matches2.csv
score.w_geom : 0.5
""")

    overrides = load_rc_file(path)

    assert overrides == {
        'match.max_strain': 0.2,
        'match.max_area': 250,
        'match.output_csv': 'matches2.csv',
        'score.w_geom': 0.5,
    }


def test_load_rc_file_ignores_blank_lines_and_comments(tmp_path):
    path = _write_rc(tmp_path / 'dbm.config', """
# a full-line comment
match.max_strain : 0.2   # an inline comment

# score.w_geom : 0.9   -- this whole line is commented out
""")

    overrides = load_rc_file(path)

    assert overrides == {'match.max_strain': 0.2}


def test_load_rc_file_parses_a_set_and_a_dict_of_tuples(tmp_path):
    path = _write_rc(tmp_path / 'dbm.config', """
score.sub_compatible_sg : {186, 194}
score.major_facets_by_bravais : {'CUB': ('1,0,0', '1,1,1')}
""")

    overrides = load_rc_file(path)

    assert overrides == {
        'score.sub_compatible_sg': {186, 194},
        'score.major_facets_by_bravais': {'CUB': ('1,0,0', '1,1,1')},
    }


def test_apply_overrides_patches_the_named_config_attributes(restore_config):
    apply_overrides({'match.max_strain': 0.2, 'score.w_geom': 0.5})

    assert config.match.max_strain == 0.2
    assert config.score.w_geom == 0.5


def test_apply_overrides_raises_for_an_unknown_namespace(restore_config):
    with pytest.raises(ValueError):
        apply_overrides({'nonexistent.max_strain': 0.2})


def test_apply_overrides_raises_for_an_unknown_attribute(restore_config):
    with pytest.raises(ValueError):
        apply_overrides({'match.not_a_real_key': 0.2})


def test_load_rc_file_and_apply_overrides_round_trip_every_current_default(tmp_path, restore_config):
    """A dbm.config written by configure.write_rc_file (every line
    commented out, value = repr(current default)) must, once every line
    is uncommented, apply_overrides back to exactly the values config.py
    already held -- i.e. round-trips as a true no-op."""
    from db_ogre_match.configure import _rc_file_lines

    lines = _rc_file_lines()
    uncommented = '\n'.join(line[2:] if line.startswith('# ') else line for line in lines)
    path = tmp_path / 'dbm.config'
    path.write_text(uncommented)

    before = {
        f'{name}.{key}': value
        for name, namespace in config._NAMESPACES.items()
        for key, value in vars(namespace).items()
    }

    apply_overrides(load_rc_file(str(path)))

    after = {
        f'{name}.{key}': value
        for name, namespace in config._NAMESPACES.items()
        for key, value in vars(namespace).items()
    }
    assert after == before
