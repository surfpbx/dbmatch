import warnings

from db_ogre_match.match import match_database

# noisy third-party deprecation warnings (spglib, OgreInterface) unrelated to
# this package's own code. Registered after importing match.py (which pulls
# in spglib) since spglib.utils re-registers its own filter for the "dict
# interface" warning at import time, and warnings.filterwarnings() prepends
# to the filter list -- so ours must come last to take priority.
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
