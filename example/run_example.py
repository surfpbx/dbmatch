"""Run the full match-then-score-then-refine workflow against the example
mother-dataset/substrate symlinks kept alongside this script.

The refine stage needs OgreInterface's IonicSurfaceMatcher, which pulls in
matscipy/scikit-opt -- see README.md's Installation section. It also takes
noticeably longer than match/score (real in-plane PES + z-shift energy
optimization per termination combo), rather than the ~1-2s match/score take.
"""

import os
from db_ogre_match.match import match
from db_ogre_match.score import score
from db_ogre_match.refine import refine
from ase.visualize import view

# Try matching the substrate (CdSe) with the 10 materials in the mother db
match(
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
score(
    matches_db='example_matches.db',
    output_csv='example_scores.csv',
    output_db='example_scores.db',
)

# now inspect the results (already sorted by total_score descending) by
# running "ase db example_scores.db -c +cod_id"

# entry 2300704 (MnTe) has a perfect score of 1.0. Let's refine it!
refine(
    selection='cod_id=2300704',
    scores_db='example_scores.db',
)

print('\nRefinement complete! To visualize the generated interfaces, run `ase gui TeMn-2300704/interfaces.db')
