"""Run the full match-then-score workflow against the example mother-dataset/
substrate symlinks kept alongside this script."""

import os
from db_ogre_match import match_database
from db_ogre_match.score import score_materials

match_database(
    mother_db='mother-dataset.db',
    substrate='substrate.cif',
    max_substrate_index=1,
    max_film_index=1,
    max_strain=0.05,
    max_area=100,
    output_db='example_matches.db',
    output_csv='example_matches.csv',
    restart=False,
)

score_materials(
    matches_db='db/example_matches.db',
    output_csv='example_scores.csv',
    output_db='example_scores.db',
)
