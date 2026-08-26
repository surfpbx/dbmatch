"""Convert a CSV written by match.py/score.py/utils.py back to an ase db,
or an ase db to CSV, recognizing the direction from src/dst's own file
extensions: .csv -> .db uses utils.csv_to_db (an Atoms('H') placeholder
per row -- CSV carries no atomic structure to restore), .db -> .csv uses
utils.db_to_csv. Any other extension pairing raises ValueError. Each
metadata sidecar (see match.write_match_metadata) is carried across the
conversion in either direction.
"""

import argparse
import os

from db_ogre_match.utils import csv_to_db, db_to_csv

_CONVERTERS = {
    ('.csv', '.db'): csv_to_db,
    ('.db', '.csv'): db_to_csv,
}


def convert(src, dst):
    """Convert src to dst, dispatching on their own file extensions --
    .csv -> .db or .db -> .csv, nothing else. Extensions are matched
    exactly (not case-insensitively): ase.db.connect itself only
    recognizes a lowercase '.db' extension, so a differently-cased match
    here could still fail one step down. Writes only; returns nothing."""
    src_ext = os.path.splitext(src)[1]
    dst_ext = os.path.splitext(dst)[1]

    converter = _CONVERTERS.get((src_ext, dst_ext))
    if converter is None:
        raise ValueError(
            f"don't know how to convert {src_ext or '(no extension)'!r} to "
            f"{dst_ext or '(no extension)'!r} -- only .csv -> .db and "
            f".db -> .csv are supported"
        )

    converter(src, dst)
    print(f'Converted {src} -> {dst}.')


def add_arguments(parser):
    """Add convert's CLI arguments to parser (shared by this file's own
    __main__ block and by cli.py's `dbm convert` subcommand)."""
    parser.add_argument('src', help='path to convert from (.csv or .db)')
    parser.add_argument('dst', help='path to convert to (.db or .csv)')
    return parser


def main(args):
    """Run convert from a parsed add_arguments() namespace -- shared by
    this file's own __main__ block and by cli.py's `dbm convert`."""
    convert(args.src, args.dst)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    main(parser.parse_args())
