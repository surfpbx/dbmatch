"""Match every structure in a mother database against a substrate, writing
results to CSV."""

import argparse
import csv
import json
import os
from ase.db import connect
from ase.io import read
from ase.visualize import view
from db_ogre_match import config
from db_ogre_match.ogre_custom import MillerSearch
from tqdm import tqdm


class Matcher:
    """
    Runs a MillerSearch scan of a fixed substrate against a given film.
    """

    def __init__(
        self,
        substrate,
        max_substrate_index,
        max_film_index,
        max_strain,
        max_area
    ):
        self.substrate = substrate
        self.max_substrate_index = max_substrate_index
        self.max_film_index = max_film_index
        self.max_strain = max_strain
        self.max_area = max_area

    def __call__(self, film):
        ms = MillerSearch(
            substrate=self.substrate,
            film=film,
            refine_structure=False,
            max_film_index=self.max_film_index,
            max_substrate_index=self.max_substrate_index,
            max_strain=self.max_strain,
            max_area_mismatch=config.match.max_area_mismatch,
            max_area=self.max_area
        )
        return ms.run_scan()


def match_row(row, matcher):
    """
    Run matcher against a single mother-db row, returning its key-value
    pairs and match results.
    """
    atoms = row.toatoms()
    kwp = row.key_value_pairs
    results = matcher(atoms)
    return kwp, results


def _coerce(value):
    """Undo CSV's stringify-everything: try int, then float, then
    'True'/'False' -> bool, else leave as str."""
    for t in (int, float):
        try:
            return t(value)
        except ValueError:
            pass
    if value in ('True', 'False'):
        return value == 'True'
    return value


def read_csv_rows(csv_path):
    """
    Read csv_path back as a list of dicts, each cell coerced back to the
    most specific Python type it round-trips to (see _coerce) -- csv
    itself only ever writes/reads strings. Each dict also gets an 'id'
    key, its 1-based row position in the file -- a matches.csv row has no
    other stable identity, but this one is stable in practice: the only
    thing that ever rewrites matches.csv after the fact
    (clear_stale_matches, on restart) only ever touches the tail of the
    file, since mother_db rows are matched in strictly ascending id
    order, so it never shifts an earlier row's position.
    """
    with open(csv_path, newline='') as f:
        rows = [
            {k: _coerce(v) for k, v in row.items()}
            for row in csv.DictReader(f)
        ]
    for i, row in enumerate(rows, start=1):
        row['id'] = i
    return rows


class CsvWriter:
    """
    Appends dict rows to a CSV file, inferring the header from the first
    row. A later row with a key not in that header has the extra key
    silently dropped.
    """

    def __init__(self, path, append):
        self.file = open(path, 'a' if append else 'w', newline='')
        self.writer = None

    def writerow(self, row):
        """Write row, creating the DictWriter and header on the first call."""
        if self.writer is None:
            self.writer = csv.DictWriter(
                self.file, fieldnames=list(row.keys()), extrasaction='ignore'
            )
            if self.file.tell() == 0:
                self.writer.writeheader()
        self.writer.writerow(row)


def store_matches(csv_writer, kwp, results):
    """
    Write each match result for a row to csv_writer.
    """
    for result in results:
        kwp.update(result)
        csv_writer.writerow(kwp)


def _dotfile_path_for(output_csv, suffix):
    """Return a dotfile path alongside output_csv (e.g. its checkpoint or
    metadata sidecar), so it doesn't clutter directory listings."""
    directory, name = os.path.split(output_csv)
    return os.path.join(directory, f'.{name}{suffix}')


def checkpoint_path_for(output_csv):
    """
    Return the checkpoint file path that goes alongside output_csv, as a
    dotfile so it doesn't clutter directory listings.
    """
    return _dotfile_path_for(output_csv, '.checkpoint')


def metadata_path_for(output_csv):
    """
    Return the metadata sidecar path that goes alongside output_csv, as a
    dotfile -- CSV has no native metadata slot the way an ase db does, so
    mother_db's/substrate's paths (needed downstream by score.py/
    refine.py) are recorded here instead.
    """
    return _dotfile_path_for(output_csv, '.meta.json')


def write_match_metadata(output_csv, mother_db, substrate):
    """Record mother_db's/substrate's absolute paths in output_csv's
    metadata sidecar, so a matches.csv row can be traced back to its
    source mother-db row, and refine.py can find the substrate again,
    without separately tracking either path."""
    with open(metadata_path_for(output_csv), 'w') as f:
        json.dump(
            {'mother_db': os.path.abspath(mother_db), 'substrate': os.path.abspath(substrate)},
            f,
        )


def read_match_metadata(output_csv):
    """The dict written by write_match_metadata for output_csv, or {} if
    no metadata sidecar exists."""
    path = metadata_path_for(output_csv)
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return json.load(f)


def read_checkpoint(checkpoint_path):
    """
    Return the last completed row id, or None if no checkpoint exists.
    """
    if not os.path.exists(checkpoint_path):
        return None

    with open(checkpoint_path) as f:
        return int(f.read().strip())


def write_checkpoint(checkpoint_path, row_id):
    """
    Atomically record row_id as the last row fully processed.
    """
    tmp_path = f'{checkpoint_path}.tmp'
    with open(tmp_path, 'w') as f:
        f.write(str(row_id))
    os.replace(tmp_path, checkpoint_path)


def clear_stale_matches(csv_path, cod_id):
    """
    Delete any match rows already written for cod_id from csv_path.

    Used before re-matching the checkpointed row on restart: that row's
    matches may or may not have actually made it to disk before the
    previous run stopped. Restart always redoes that one row from
    scratch, so its own prior output (partial, complete, or none at all)
    must be cleared first to avoid duplicates.
    """
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0:
        return

    with open(csv_path, newline='') as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = [row for row in reader if row['cod_id'] != str(cod_id)]

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_matching(
    db,
    csv_writer,
    matcher,
    selection=None,
    checkpoint_path=None
):
    """
    Match every selected row in db, storing results and advancing the
    checkpoint after each row (even ones that raised and were skipped).

    Returns the number of match rows written.
    """
    n_matches = 0
    for row in tqdm(db.select(selection), total=db.count(selection)):
        try:
            kwp, results = match_row(row, matcher)
        except:
            print(row.cod_id)
        else:
            store_matches(csv_writer, kwp, results)
            n_matches += len(results)

        if checkpoint_path is not None:
            write_checkpoint(checkpoint_path, row.id)

    return n_matches


def match(
    substrate,
    mother_db,
    max_substrate_index,
    max_film_index,
    max_strain,
    max_area,
    output_csv,
    restart,
):
    """
    Match mother_db against substrate and write results to output_csv in
    the current directory, optionally resuming from output_csv's
    checkpoint file -- restart re-matches the checkpointed row itself too
    (not just the ones after it), first clearing any of its matches
    already in output_csv, since there's no guarantee its matches made it
    to disk before the previous run stopped.

    mother_db's and substrate's absolute paths are recorded in a metadata
    sidecar next to output_csv (see write_match_metadata), so that a
    matches row can be traced back to its source mother-db row, and
    refine.py can find the substrate again, without separately tracking
    either path.
    """
    print(f'\nMatching {mother_db} against {substrate}...')

    csv_path = output_csv

    checkpoint_path = checkpoint_path_for(csv_path)

    restart_id = None
    if restart:
        restart_id = read_checkpoint(checkpoint_path)
        if restart_id is None:
            raise ValueError(f'no checkpoint file found at {checkpoint_path}')

    do_restart = restart_id is not None
    selection = f'id>={restart_id}' if do_restart else None

    db = connect(mother_db)

    if do_restart:
        clear_stale_matches(csv_path, db.get(id=restart_id).cod_id)

    # opened only after clear_stale_matches's own rewrite (if any) is done,
    # so its append-mode file position reflects the file's true final size
    csv_writer = CsvWriter(csv_path, append=do_restart)

    write_match_metadata(csv_path, mother_db, substrate)

    substrate_atoms = read(substrate)

    matcher = Matcher(
        substrate_atoms,
        max_substrate_index,
        max_film_index,
        max_strain,
        max_area
    )

    n_matches = run_matching(
        db,
        csv_writer,
        matcher,
        selection,
        checkpoint_path
    )
    csv_writer.file.close()

    print(
        f'Matching finished. Found {n_matches} matches, written to {csv_path}.'
    )


def add_arguments(parser):
    """Add match's CLI arguments to parser (shared by this file's own
    __main__ block and by cli.py's `dbm match` subcommand)."""
    parser.add_argument(
        'substrate',
        help='path to the substrate structure file'
    )
    parser.add_argument(
        'mother_db',
        help='path to the mother database'
    )
    parser.add_argument(
        '--max-substrate-index', type=int, default=config.match.max_sub_index,
        help='max substrate miller index (default: %(default)s)'
    )
    parser.add_argument(
        '--max-film-index', type=int, default=config.match.max_flm_index,
        help='max film miller index (default: %(default)s)'
    )
    parser.add_argument(
        '--max-strain', type=float, default=config.match.max_strain,
        help='max strain (default: %(default)s)'
    )
    parser.add_argument(
        '--max-area', type=float, default=config.match.max_area,
        help='max area (default: %(default)s)'
    )
    parser.add_argument(
        '-o', '--output-csv', default=config.match.output_csv,
        help='name of the output CSV file (default: %(default)s)'
    )
    parser.add_argument(
        '-r', '--restart', action='store_true',
        help=(
            'resume from the checkpoint file next to the output CSV, appending to it '
            'and re-matching the checkpointed row itself '
            '(default: always start from the beginning)'
        )
    )
    return parser


def main(args):
    """Run match from a parsed add_arguments() namespace -- shared by this
    file's own __main__ block and by cli.py's `dbm match`."""
    match(
        args.substrate,
        args.mother_db,
        args.max_substrate_index,
        args.max_film_index,
        args.max_strain,
        args.max_area,
        args.output_csv,
        args.restart,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    main(parser.parse_args())
