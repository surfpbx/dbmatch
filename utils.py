import csv
from ase.db import connect


def db_to_csv(db_path, csv_path):
    rows = [row.key_value_pairs for row in connect(db_path).select()]
    fieldnames = sorted({key for row in rows for key in row})

    with open(csv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
