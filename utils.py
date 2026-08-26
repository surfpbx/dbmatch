"""Utility helpers for converting ASE databases to CSV, and back."""

import csv
import json
from ase import Atoms
from ase.db import connect
from db_ogre_match.match import metadata_path_for, read_csv_rows, read_match_metadata


def db_to_csv(db_path, csv_path):
    """
    Write every row's key_value_pairs from db_path to a fresh csv_path.
    If db_path has its own .metadata, it's carried over to csv_path's own
    metadata sidecar (see match.write_match_metadata/read_match_metadata)
    so a later csv_to_db(csv_path, ...) can restore it.
    """
    db = connect(db_path)
    rows = [row.key_value_pairs for row in db.select()]
    fieldnames = sorted({key for row in rows for key in row})

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    db.count()  # metadata is only readable after some query
    if db.metadata:
        with open(metadata_path_for(csv_path), 'w') as f:
            json.dump(db.metadata, f)


def csv_to_db(csv_path, db_path):
    """
    Write every row of csv_path to a fresh db_path, one row per line,
    using an empty Atoms('H') placeholder for each row -- CSV carries no
    atomic structure to restore, this is not real structural data.
    key_value_pairs come from read_csv_rows's own type-coerced dict,
    minus 'id' (its own synthetic row-position id, not a real CSV column)
    and 'formula' (an ase db reserved/computed name) -- both of which
    ase.db's key_value_pairs check would otherwise reject. If csv_path
    has its own metadata sidecar (see match.write_match_metadata), it's
    restored as db_path's own .metadata.
    """
    newdb = connect(db_path, append=False)
    for row in read_csv_rows(csv_path):
        kvp = {k: v for k, v in row.items() if k not in ('id', 'formula')}
        newdb.write(Atoms('H'), **kvp)

    metadata = read_match_metadata(csv_path)
    if metadata:
        newdb.metadata = metadata
