import warnings

from db_ogre_match import ogre_custom

# match_database/score_materials/refine_material are deliberately not
# re-exported here -- import each from its own module
# (db_ogre_match.match/.score/.refine). The ogre_custom import above is kept
# only for its side effect: it's what actually pulls in spglib (via
# OgreInterface's miller/surfaces/lattice_match submodules).
#
# noisy third-party deprecation warnings (spglib, OgreInterface) unrelated to
# this package's own code. Registered after importing ogre_custom.py (which
# pulls in spglib) since spglib.utils re-registers its own filter for the
# "dict interface" warning at import time, and warnings.filterwarnings()
# prepends to the filter list -- so ours must come last to take priority.
warnings.filterwarnings(
    'ignore',
    message='Set OLD_ERROR_HANDLING to false and catch the errors directly.',
    category=DeprecationWarning,
)
warnings.filterwarnings(
    'ignore',
    message='dict interface is deprecated. Use attribute interface instead',
    category=DeprecationWarning,
)
