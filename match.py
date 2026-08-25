"""Match every structure in a mother database against a substrate, writing
results to db and CSV."""

import argparse
import csv
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
    Run matcher against a single mother-db row, returning its atoms,
    key-value pairs, and match results.
    """
    atoms = row.toatoms()
    kwp = row.key_value_pairs
    results = matcher(atoms)
    return atoms, kwp, results


class CsvWriter:
    """
    Appends dict rows to a CSV file, inferring the header from the first
    row. A later row with a key not in that header has the extra key
    silently dropped from the CSV (it's still written in full to the db).
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


def store_matches(newdb, csv_writer, atoms, kwp, results):
    """
    Write each match result for a row to both newdb and csv_writer.
    """
    for result in results:
        kwp.update(result)
        newdb.write(atoms, **kwp)
        csv_writer.writerow(kwp)


def checkpoint_path_for(output_db):
    """
    Return the checkpoint file path that goes alongside output_db, as a
    dotfile so it doesn't clutter directory listings.
    """
    directory, name = os.path.split(output_db)
    return os.path.join(directory, f'.{name}.checkpoint')


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


def clear_stale_matches(newdb, cod_id):
    """
    Delete any match rows already written for cod_id from newdb.

    Used before re-matching the checkpointed row on restart: that row's
    matches may or may not have actually made it to disk before the
    previous run stopped. Restart always redoes that one row from
    scratch, so its own prior output (partial, complete, or none at all)
    must be cleared first to avoid duplicates.
    """
    stale_ids = [row.id for row in newdb.select(cod_id=cod_id)]
    if stale_ids:
        newdb.delete(stale_ids)


def run_matching(
    db,
    newdb,
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
            atoms, kwp, results = match_row(row, matcher)
        except:
            print(row.cod_id)
        else:
            store_matches(newdb, csv_writer, atoms, kwp, results)
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
    output_basename,
    restart,
):
    """
    Match mother_db against substrate and write results to
    <output_basename>.db/<output_basename>.csv in the current directory,
    optionally resuming from the output db's checkpoint file -- restart
    re-matches the checkpointed row itself too (not just the ones after
    it), first clearing any of its matches already in the output db,
    since there's no guarantee its matches made it to disk before the
    previous run stopped.

    mother_db's and substrate's absolute paths are recorded in the output
    db's metadata (as 'mother_db'/'substrate'), so that a matches row can
    be traced back to its source mother-db row, and refine.py can find the
    substrate again, without separately tracking either path.
    """
    print(f'\nMatching {mother_db} against {substrate}...')

    db_path = f'{output_basename}.db'
    csv_path = f'{output_basename}.csv'

    checkpoint_path = checkpoint_path_for(db_path)

    restart_id = None
    if restart:
        restart_id = read_checkpoint(checkpoint_path)
        if restart_id is None:
            raise ValueError(f'no checkpoint file found at {checkpoint_path}')

    do_restart = restart_id is not None
    selection = f'id>={restart_id}' if do_restart else None

    if do_restart and os.path.exists(f'{db_path}.lock'):
        print(
            f'WARNING: {db_path}.lock already exists -- if the previous run was '
            f'killed while writing (crash, OOM, kill -9) this is a stale lock '
            f'left behind, and ase.db will wait on it forever without any error. '
            f'Delete it yourself once you are sure no other process is writing '
            f'to {db_path} (rm {db_path}.lock), then rerun.'
        )

    db = connect(mother_db)
    newdb = connect(db_path, append=do_restart)
    csv_writer = CsvWriter(csv_path, append=do_restart)

    if do_restart:
        clear_stale_matches(newdb, db.get(id=restart_id).cod_id)

    newdb.metadata = {
        'mother_db': os.path.abspath(mother_db),
        'substrate': os.path.abspath(substrate),
    }

    substrate_atoms = read(substrate)

    matcher = Matcher(
        substrate_atoms,
        max_substrate_index,
        max_film_index,
        max_strain,
        max_area
    )

    n_matches = run_matching(
        db, newdb,
        csv_writer,
        matcher,
        selection,
        checkpoint_path
    )
    csv_writer.file.close()

    print(
        f'Matching finished. Found {n_matches} matches, '
        f'written to {db_path} and {csv_path}.'
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
        '-o', '--output', default=config.match.output_basename,
        help=(
            'basename for the output files -- results are written to '
            '<output>.db and <output>.csv (default: %(default)s)'
        )
    )
    parser.add_argument(
        '-r', '--restart', action='store_true',
        help=(
            'resume from the checkpoint file next to the output database, appending to it '
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
        args.output,
        args.restart,
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    add_arguments(parser)
    main(parser.parse_args())
