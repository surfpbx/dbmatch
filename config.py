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
    output_db='matches.db',   # output database filename, written under 'db/'
    output_csv='matches.csv', # output CSV filename, written under 'csv/'
)

# facets treated as "major" for each bravais class, keyed by bravais symbol
_major_facets_by_bravais = {
    'CUB': ('1,0,0', '1,1,1'),              # cubic
    'TET': ('1,0,0', '0,0,1', '1,1,0'),     # tetragonal
    'ORT': ('1,0,0', '0,1,0', '0,0,1'),     # orthorhombic
    'HEX': ('0,0,1', '0,1,0', '1,1,0'),     # hexagonal
    'TRG': ('0,0,1', '0,1,0', '1,1,0'),     # trigonal
}

# ----------------------
# Scoring stage options
# ----------------------
score = SimpleNamespace(
    output_db='scores.db',         # output database filename, written under 'db/'
    output_csv='scores.csv',       # output CSV filename, written under 'csv/'
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

# -----------------------
# Refinement stage options
# -----------------------
refine = SimpleNamespace(
    layers=5,                 # number of atomic layers in each generated surface slab
    vacuum=10.0,               # vacuum padding (Angstrom) around each generated surface slab
    interfacial_distance=3.0, # starting substrate-film separation (Angstrom) before z-shift optimization
    z_shift_min=1.0,           # lower bound (Angstrom) of the interfacial-distance scan
    z_shift_n_points=31,       # number of points sampled in the interfacial-distance scan
)
