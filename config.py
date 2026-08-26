from types import SimpleNamespace

# ----------------------
# Matching stage options
# ----------------------
match = SimpleNamespace(
    max_sub_index=1,          # max substrate miller index scanned
    max_flm_index=1,          # max film miller index scanned
    max_strain=0.05,          # max strain allowed between substrate/film lattices
    max_area=500,             # max interface supercell area
    max_area_mismatch=0.08,   # max relative area mismatch between substrate/film supercells
    output_csv='matches.csv', # output CSV filename
)

# facets treated as "major" for each bravais class, keyed by bravais symbol
_major_facets_by_bravais = {
    'CUB': ('1,0,0', '1,1,1'),              # cubic
    'TET': ('1,0,0', '0,0,1', '1,1,0'),     # tetragonal
    'ORT': ('1,0,0', '0,1,0', '0,0,1'),     # orthorhombic
    'HEX': ('0,0,1', '0,1,0', '1,1,0'),     # hexagonal
    'TRG': ('0,0,1', '0,1,0', '1,1,0'),     # trigonal
}

# ---------------------
# Scoring stage options
# ---------------------
score = SimpleNamespace(
    output_basename='scores',      # output filename, minus extension -- results are written to
                                    # <output_basename>.db and <output_basename>.csv
    w_geom=0.75,                   # weight of geom_score in total_score
    w_sg=0.15,                     # weight of sg_score in total_score
    w_comp=0.10,                   # weight of comp_score in total_score
    base_major_major=1/3,          # tier value for a major-substrate/major-film facet pairing
    sub_major_minor_ratio=4.0,     # how much more a major substrate facet is worth than a minor one (film tier fixed)
    flm_major_minor_ratio=2.0,     # how much more a major film facet is worth than a minor one (substrate tier fixed)
    major_facets_by_bravais=_major_facets_by_bravais,   # major facets by bravais class, used to tier each match's film facet
    major_sub_facets={'0,0,1', '0,1,0', '1,1,0'},       # the fixed substrate's own major facets (independent of film bravais)
    sub_compatible_sg={186, 194, 216, 225},             # space groups compatible with the substrate (binary sg_score bonus)
    sub_compatible_elements={                           # elements compatible with the substrate (binary comp_score bonus)
        'Zn', 'Cd', 'Hg',
        'Al', 'Ga', 'In',
        'N', 'P', 'As', 'Sb', 'Bi',
        'S', 'Se', 'Te',
    },
)

# ------------------------
# Refinement stage options
# ------------------------
refine = SimpleNamespace(
    layers=10,                 # number of atomic layers in each generated surface slab
    vacuum=20.0,               # vacuum padding (Angstrom) around each generated surface slab
    z_shift_n_points=31,       # number of points sampled in the interfacial-distance scan (both the
                               # starting interfacial_distance and the scan bounds are derived per-combo
                               # from ionic radii -- see refine._contact_distance/_z_shift_range)
    contact_distance_factor=1.25,  # scales _contact_distance's raw ionic-radii sum -- pure Shannon/
                                    # ionic radii tend to sit a bit closer than a real relaxed contact
                                    # distance, so this nudges the starting point (and z-shift scan
                                    # range, since both derive from the same D) outward
    pes_colormap_bound=1.0,    # eV/A^2: fixed symmetric bound for the PES plot's colorbar, replacing
                               # OgreInterface's own data-min/max auto-scaling -- see
                               # refine._run_surface_matching
)

# ---------------------------------------------------------
# Project-local config file (see configure.py's `dbm config`)
# ---------------------------------------------------------
import ast
import os

RC_FILENAME = 'dbm.config'

_NAMESPACES = {'match': match, 'score': score, 'refine': refine}


def load_rc_file(path):
    """
    Parse an rc file into a {'namespace.key': value} dict.

    Blank lines and '#'-prefixed comments (including inline, matplotlibrc-
    style -- everything from the first '#' on a line onward is stripped)
    are skipped. Each value is parsed with ast.literal_eval (so ints,
    floats, strings, True/False, and container literals like {...}/(...)
    all round-trip as their real Python type); a value that isn't a valid
    literal (e.g. a bare, unquoted filename) is kept as the raw stripped
    string instead.
    """
    overrides = {}
    with open(path) as f:
        for line in f:
            line = line.split('#', 1)[0].strip()
            if not line:
                continue
            key, _, value = line.partition(':')
            key, value = key.strip(), value.strip()
            try:
                value = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                pass
            overrides[key] = value
    return overrides


def apply_overrides(overrides):
    """
    Patch match/score/refine in place from a {'namespace.key': value}
    dict (as returned by load_rc_file).

    Raises ValueError for a key whose namespace or attribute doesn't
    exist -- a typo'd or stale config key should fail loudly rather than
    be silently ignored.
    """
    for dotted_key, value in overrides.items():
        namespace_name, _, attr = dotted_key.partition('.')
        namespace = _NAMESPACES.get(namespace_name)
        if namespace is None or not hasattr(namespace, attr):
            raise ValueError(f'unknown config key {dotted_key!r} in {RC_FILENAME}')
        setattr(namespace, attr, value)


_rc_path = os.path.join(os.getcwd(), RC_FILENAME)
if os.path.exists(_rc_path):
    apply_overrides(load_rc_file(_rc_path))
