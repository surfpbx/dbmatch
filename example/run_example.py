import os
import sys
from pathlib import Path

EXAMPLE_DIR = Path(__file__).resolve().parent

sys.path.insert(0, str(EXAMPLE_DIR.parent))
from match import main

os.chdir(EXAMPLE_DIR)

main(
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
