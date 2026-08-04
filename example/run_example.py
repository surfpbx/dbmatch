"""Run the full match-then-score-then-refine workflow against the example
mother-dataset/substrate symlinks kept alongside this script.

The refine stage needs OgreInterface's IonicSurfaceMatcher, which pulls in
matscipy/scikit-opt -- see README.md's Installation section. It also takes
noticeably longer than match/score (real in-plane PES + z-shift energy
optimization per termination combo), rather than the ~1-2s match/score take.
"""

import os
from db_ogre_match.match import match_database
from db_ogre_match.score import score_materials
from db_ogre_match.refine import refine_material

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

# substrate/matches_db are both read from example_scores.db's own metadata
refine_material(
    cod_id=2300704,
    scores_db='db/example_scores.db',
)
