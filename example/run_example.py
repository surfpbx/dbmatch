"""Run match_database against the mother-dataset/substrate symlinks kept
alongside this script."""

import os
from pathlib import Path
from db_ogre_match import match_database

#EXAMPLE_DIR = Path(__file__).resolve().parent
#
#os.chdir(EXAMPLE_DIR)

match_database(
    mother_db='mother-dataset.db',
    substrate='substrate.cif',
    max_substrate_index=1,
    max_film_index=1,
    max_strain=0.05,
    max_area=500,
    output_db='example.db',
    output_csv='example.csv',
    restart=False,
)
