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
from ase.visualize import view

# Try matching the substrate (CdSe) with the 10 materials in the mother db
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

# Assign scores to the materials based on the found matches
score_materials(
    matches_db='db/example_matches.db',
    output_csv='example_scores.csv',
    output_db='example_scores.db',
)

# now inspect the results sorted by total_score in descending order
# by running "ase db db/example_scores.db -s total_score- -c +cod_id"

# entry 2300704 (MnTe) has a perfect score of 1.0. Let's refine it!
refine_material(
    selection='cod_id=2300704',
    scores_db='db/example_scores.db',
)

print('\nRefinement complete! To visualize the generated interfaces, run `ase gui TeMn-2300704/interfaces.db')
