"""Match every structure in a mother database against a substrate, writing
results to db and CSV."""

import argparse
import csv
import os
from ase.db import connect
from ase.io import read
from ase.visualize import view
from db_ogre_match import config
from db_ogre_match.miller_custom import MillerSearch
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
    row.
    """

    def __init__(self, path, append):
        self.file = open(path, 'a' if append else 'w', newline='')
        self.writer = None

    def writerow(self, row):
        """Write row, creating the DictWriter and header on the first call."""
        if self.writer is None:
            self.writer = csv.DictWriter(self.file, fieldnames=list(row.keys()))
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
    Return the checkpoint file path that goes alongside output_db.
    """
    return f'{output_db}.checkpoint'


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


def match_database(
    mother_db,
    substrate,
    max_substrate_index,
    max_film_index,
    max_strain,
    max_area,
    output_db,
    output_csv,
    restart,
):
    """
    Match mother_db against substrate and write results under db/ and csv/,
    optionally resuming from the output db's checkpoint file.
    """
    print(f'\nMatching {mother_db} against {substrate}...')

    os.makedirs('db', exist_ok=True)
    os.makedirs('csv', exist_ok=True)

    db_path = os.path.join('db', output_db)
    csv_path = os.path.join('csv', output_csv)

    checkpoint_path = checkpoint_path_for(db_path)

    restart_id = None
    if restart:
        restart_id = read_checkpoint(checkpoint_path)
        if restart_id is None:
            raise ValueError(f'no checkpoint file found at {checkpoint_path}')

    do_restart = restart_id is not None
    selection = f'id>{restart_id}' if do_restart else None

    db = connect(mother_db)
    newdb = connect(db_path, append=do_restart)
    csv_writer = CsvWriter(csv_path, append=do_restart)
    # TODO: record mother_db's absolute path in newdb.metadata, the same way
    # score.py now records matches_db's path in the scores db's metadata --
    # so a matches row can be traced back to its source mother-db row without
    # separately tracking which mother_db a given matches_db came from.

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


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'mother_db',
        help='path to the mother database'
    )
    parser.add_argument(
        'substrate',
        help='path to the substrate structure file'
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
        '-o', '--output-db', default=config.match.output_db,
        help="name of the output database file, written under 'db/' (default: %(default)s)"
    )
    parser.add_argument(
        '--output-csv', default=config.match.output_csv,
        help="name of the output CSV file, written under 'csv/' (default: %(default)s)"
    )
    parser.add_argument(
        '-r', '--restart', action='store_true',
        help=(
            'resume from the checkpoint file next to the output database, appending to it '
            '(default: always start from the beginning)'
        )
    )
    args = parser.parse_args()

    match_database(
        args.mother_db,
        args.substrate,
        args.max_substrate_index,
        args.max_film_index,
        args.max_strain,
        args.max_area,
        args.output_db,
        args.output_csv,
        args.restart,
    )
