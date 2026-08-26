"""Write a matplotlibrc-style project config template (config.RC_FILENAME,
'dbm.config' by default) to the current directory -- every
match/score/refine config key, commented out with its current default
already filled in. Uncommenting and editing a line overrides that
default for every `dbm` command run from this directory; see
config.load_rc_file/apply_overrides for how the file is parsed and
applied. The file is only ever read once, when db_ogre_match's config
module is first imported (like matplotlib's own matplotlibrc) -- editing
it, or changing directory, has no effect on a process already running.
"""

import argparse
import os

from db_ogre_match import config

_HEADER = f"""\
## {config.RC_FILENAME} -- db_ogre_match project configuration
##
## Every line below is commented out, with the package's own current
## default value already filled in. Uncomment (remove the leading '#')
## and edit a line to override that default for any `dbm` command run
## from this directory.
##
## This file is only ever read once, when db_ogre_match's config module
## is first imported -- like matplotlib's own matplotlibrc, editing this
## file (or changing directory) after that point has no effect until the
## next fresh run.
##
## Values are parsed as Python literals where possible (ast.literal_eval)
## -- numbers, True/False, quoted strings, and container literals like
## {{...}}/(...) all work directly; anything that doesn't parse as a
## literal (e.g. a bare, unquoted filename) is kept as a plain string.
"""

_SECTIONS = [
    ('match', config.match),
    ('score', config.score),
    ('refine', config.refine),
]


def _rc_file_lines():
    lines = [_HEADER]
    for section, namespace in _SECTIONS:
        lines.append(f'## {section} ##')
        for key, value in vars(namespace).items():
            lines.append(f'# {section}.{key} : {value!r}')
        lines.append('')
    return lines


def write_rc_file(path=config.RC_FILENAME, force=False):
    """
    Write the config template to path (config.RC_FILENAME, 'dbm.config',
    in the current directory by default).

    Refuses to overwrite an existing file unless force=True, so a
    already-customized config file is never silently clobbered.
    """
    if os.path.exists(path) and not force:
        raise FileExistsError(
            f'{path} already exists -- pass force=True (--force on the CLI) '
            f'to overwrite it, or remove/rename it first'
        )

    with open(path, 'w') as f:
        f.write('\n'.join(_rc_file_lines()))

    print(f'Wrote {path}.')


def add_arguments(parser):
    """Add config's CLI arguments to parser (shared by this file's own
    __main__ block and by cli.py's `dbm config` subcommand)."""
    parser.add_argument(
        '-o', '--output', default=config.RC_FILENAME,
        help='path to write the config template to (default: %(default)s)'
    )
    parser.add_argument(
        '--force', action='store_true',
        help='overwrite an existing config file (default: refuse if one already exists)'
    )
    return parser


def main(args):
    """Run write_rc_file from a parsed add_arguments() namespace -- shared
    by this file's own __main__ block and by cli.py's `dbm config`."""
    write_rc_file(args.output, force=args.force)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    main(parser.parse_args())
