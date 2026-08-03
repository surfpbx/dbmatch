"""dbm: single console-script entry point dispatching to db_ogre_match's
match / score / refine pipeline stages.

Two-phase argparse dispatch: which subcommand was requested is determined
first (parse_known_args against a bare, help-less subparser per stage), and
only *that* stage's module is imported afterwards. This matters because
refine.py imports OgreInterface.generate/OgreInterface.surface_matching at
module level (pulling in matscipy/scikit-opt, deliberately-optional extra
deps) -- dbm match/dbm score must never trigger that import.
"""

import argparse

_SUBCOMMAND_HELP = {
        'match': "match every structure in a mother database against a substrate. For options, check `dbm match --help`",
        'score': "score matched materials by geometric/space-group/chemistry compatibility. For options, check `dbm score --help`",
        'refine': "energetically refine a material's selected match. For options, check `dbm refine --help`",
}


def _load_command_module(command):
    """Import and return the db_ogre_match submodule for command, lazily --
    called only after the subcommand name is already known, so choosing
    'match'/'score' never imports 'refine' (see module docstring)."""
    if command == 'match':
        from db_ogre_match import match as mod
    elif command == 'score':
        from db_ogre_match import score as mod
    elif command == 'refine':
        from db_ogre_match import refine as mod
    return mod


def main(argv=None):
    """Entry point for the `dbm` console script.

    First pass: a bare top-level parser with help-less placeholder
    subparsers, just enough to find out which subcommand was requested
    without importing any stage module yet. Second pass: only the selected
    stage's module is imported, and its own add_arguments/main (the same
    functions its standalone `python <stage>.py` CLI uses) build the real
    parser and run it against the leftover argv.
    """
    parser = argparse.ArgumentParser(
        prog='dbm', description='db-ogre-match pipeline: match -> score -> refine',
    )
    subparsers = parser.add_subparsers(dest='command', required=True)
    for name, help_text in _SUBCOMMAND_HELP.items():
        subparsers.add_parser(name, add_help=False, help=help_text)

    args, remaining = parser.parse_known_args(argv)
    mod = _load_command_module(args.command)

    # description=mod.__doc__ surfaces each stage's own module docstring
    # under `dbm <command> --help`, same as e.g. refine.py's standalone
    # --help already does with its own folder-layout diagram.
    sub_parser = argparse.ArgumentParser(
        prog=f'dbm {args.command}', description=mod.__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    mod.add_arguments(sub_parser)
    mod.main(sub_parser.parse_args(remaining))


if __name__ == '__main__':
    main()
