"""Utility helpers for converting ASE databases to CSV."""

import csv
from ase.db import connect


def db_to_csv(db_path, csv_path):
    """
    Write every row's key_value_pairs from db_path to a fresh csv_path.
    """
    rows = [row.key_value_pairs for row in connect(db_path).select()]
    fieldnames = sorted({key for row in rows for key in row})

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
