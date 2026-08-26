import pytest

from db_ogre_match import config
from db_ogre_match.configure import write_rc_file


def test_write_rc_file_writes_to_the_given_path(tmp_path):
    path = str(tmp_path / 'dbm.config')

    write_rc_file(path)

    assert (tmp_path / 'dbm.config').exists()


def test_write_rc_file_every_config_line_is_commented_out(tmp_path):
    path = tmp_path / 'dbm.config'
    write_rc_file(str(path))

    config_lines = [
        line for line in path.read_text().splitlines()
        if '.' in line and ':' in line and not line.startswith('##')
    ]
    assert config_lines  # sanity: found at least the match/score/refine lines
    assert all(line.startswith('# ') for line in config_lines)


def test_write_rc_file_includes_every_namespace_and_key(tmp_path):
    path = tmp_path / 'dbm.config'
    write_rc_file(str(path))

    text = path.read_text()
    for name, namespace in [('match', config.match), ('score', config.score), ('refine', config.refine)]:
        for key in vars(namespace):
            assert f'{name}.{key} :' in text


def test_write_rc_file_refuses_to_overwrite_an_existing_file(tmp_path):
    path = tmp_path / 'dbm.config'
    path.write_text('# a user has already customized this')

    with pytest.raises(FileExistsError):
        write_rc_file(str(path))

    assert path.read_text() == '# a user has already customized this'


def test_write_rc_file_overwrites_when_forced(tmp_path):
    path = tmp_path / 'dbm.config'
    path.write_text('# a user has already customized this')

    write_rc_file(str(path), force=True)

    assert '# a user has already customized this' not in path.read_text()
