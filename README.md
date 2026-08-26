# db-ogre-match

Match every structure in a "mother" database against a fixed substrate,
finding epitaxial interfaces via [OgreInterface](https://github.com/DerekDardzinski/OgreInterface),
score each material by how well and how compatibly it matches, then
energetically refine the best matches of a chosen material.

The pipeline has three stages. `match.py` writes plain CSV only (no atoms
duplicated per match -- see below); `score.py`/`refine.py` read and write
[ASE databases](https://docs.ase-lib.org/ase/db/db.html) and CSV files
(`refine.py` writes only a database, per material):

```
substrate + mother_db  --match.py -->  matches csv  --score.py -->  scores db/csv  --refine.py -->  <formula>-<cod_id>/interfaces.db
```

- **`match.py`** scans every structure in the mother database against the
  substrate with `OgreInterface`'s `MillerSearch`, writing one row per
  facet-pairing match found -- CSV only, no atoms: the film structure is
  unchanged from mother_db, so `score.py`/`refine.py` just re-read it from
  there by `cod_id` instead of it being duplicated once per match.
- **`score.py`** reads a matches CSV back in, grouped by material
  (`cod_id`), and assigns each material a `total_score` based on how many
  facets it hits, how compatible its space group is, and whether it
  contains substrate-compatible elements.
- **`refine.py`** takes a scores-database selection (e.g. `total_score>0.7`),
  and for every material it resolves to, rebuilds the matches `score.py`
  selected for it and runs `OgreInterface`'s ionic surface matching
  (in-plane + interfacial-distance energy optimization) over every
  substrate/film termination combination of each one.

## Installation

Below, the install instructions with both plain Python and Conda.
`OgreInterface` is **not** on PyPI and is not declared as a dependency in
`pyproject.toml` -- it has to be set up separately. 
To avoid dependency conflicts, **don't** `pip install` `OgreInterface` itself;
clone it, install only the needed dependencies directly, and put the clone on the
import path.

### Plain Python (`venv`)

```bash
# create a virtual environment 
python3 -m venv venv
source venv/bin/activate

# install db_ogre_match itself, plus OgreInterface's actual runtime deps
# (unpinned, and trimmed to just what ogre_custom.py's import path needs
# for match.py/score.py -- see the refine.py note below)
cd venv
git clone https://github.com/surfpbx/dbmatch.git
pip install -e ./dbmatch pymatgen matplotlib scipy pandas

# OgreInterface has no PyPI release -- clone it instead of pip-installing it
git clone https://github.com/DerekDardzinski/OgreInterface.git

# put it on the import path permanently, without touching PYTHONPATH:
# a .pth file in site-packages is read by Python on every startup in this
# environment, venv or conda alike
echo "$(pwd)/OgreInterface" > "$(python3 -c 'import site; print(site.getsitepackages()[0])')/ogreinterface.pth"

# only needed if you'll also run refine.py: its IonicSurfaceMatcher pulls in
# two more of OgreInterface's deps that match.py/score.py never import
pip install matscipy scikit-opt
```

### Conda

Same recipe -- conda just provides the Python interpreter/environment,
`pip` still does the actual installs:

```bash
conda create -n db-ogre-match python=3.10
conda activate db-ogre-match

mkdir db-ogre-match
cd db-ogre-match
git clone https://github.com/surfpbx/dbmatch.git
pip install -e ./dbmatch pymatgen matplotlib scipy pandas

git clone https://github.com/DerekDardzinski/OgreInterface.git

echo "$(pwd)/OgreInterface" > "$(python3 -c 'import site; print(site.getsitepackages()[0])')/ogreinterface.pth"

# only needed if you'll also run refine.py -- see the Plain Python note above
pip install matscipy scikit-opt
```

### Verify it worked

```bash
python3 -c "import db_ogre_match, OgreInterface; print('OK')"
dbm --help
pytest
```

`import OgreInterface` has to succeed before `db_ogre_match` will work --
`ogre_custom.py` imports directly from it.

## Input format

`mother_db` is an ASE database where every row is one candidate material,
with at least these key-value pairs:

| key              | used for                                          |
|------------------|----------------------------------------------------|
| `cod_id`         | grouping match rows back into materials when scoring |
| `reduced_formula`| the chemistry-compatibility score component        |
| `bravais`        | classifying that material's own facets as major/minor when scoring (one of `CUB`, `TET`, `ORT`, `HEX`, `TRG`) |
| `sg_number`      | the space-group-compatibility score component       |

Any other key-value pairs (lattice parameters, space group symbol,
stoichiometry, ...) are just carried through untouched into the matches
and scores output. `tests/data/test-mother.db` is a small worked example.

`substrate` is a structure file readable by `ase.io.read` (e.g. a `.cif`).

## Usage

### As a library

```python
from db_ogre_match.match import match
from db_ogre_match.score import score

match(
    substrate='substrate.cif',
    mother_db='mother.db',
    max_substrate_index=1,
    max_film_index=1,
    max_strain=0.05,
    max_area=500,
    output_csv='matches.csv',    # written to the current directory
    restart=False,
)

score(
    matches_db='matches.csv',
    output_basename='scores',    # writes scores.db and scores.csv to the current directory
)
```

`score` needs a mother database to re-read each material's atoms/formula
by `cod_id` (matches.csv itself carries no atoms). By default it reads
`mother_db`'s path from matches.csv's own metadata sidecar (recorded
there by `match`); pass `mother_db=` explicitly (`--mother-db` on the
CLI) to avoid relying on that sidecar -- e.g. if it could be lost or
moved -- or to point at a different mother_db than the one `match`
originally used:

```python
score(
    matches_db='matches.csv',
    output_basename='scores',
    mother_db='mother.db',       # overrides matches.csv's own metadata sidecar
)
```

`score`'s own geometric/space-group/chemistry scoring is just the built-in
default -- pass `scoring_fn` to score by any other criterion instead:

```python
def my_scoring_fn(rows):
    # rows: one cod_id's matches group (group_csv(matches_db)'s own value
    # type) -- one dict per facet-pairing match.py found for this material,
    # carrying both match.py's own per-match fields (area, strain, sub_hkl,
    # flm_hkl, sub_transform, flm_transform, n_sub_reps, n_flm_reps) and the
    # material's own mother-db fields (cod_id, reduced_formula, sg_number,
    # bravais, ...), the same on every row. Loop over rows and combine
    # whatever match.py fields matter for this scoring criterion, e.g.:
    best_area = min(row['area'] for row in rows)
    return {'total_score': ..., 'my_extra_field': [...]}   # total_score is mandatory

score(
    matches_db='matches.csv',
    output_basename='scores',
    scoring_fn=my_scoring_fn,
)
```

`scoring_fn` must return a dict with at least `total_score`; any other
keys are carried into the output db/CSV alongside the material's
own mother-db fields (list-valued extras are `;`-joined into a string, the
same way `geom_score_for_material`'s own `selected_match_ids`/
`major_facets_scored` already are). `w_geom`/`w_sg`/`w_comp` and the rest
of `score`'s tunable keyword arguments are only used to build the default
`scoring_fn` when one isn't given -- they're ignored otherwise. There's no
CLI flag for `scoring_fn` (an arbitrary Python callable has no CLI
encoding); `dbm score`/`python score.py` always use the built-in default.

After inspecting the output scores database, we can pick out materials for
refinement based on a selection string:

```python
from db_ogre_match.refine import refine

refine(
    selection='cod_id=2300704',   # or e.g. 'total_score>0.7'
    scores_db='scores.db',
)
```

`selection` is any query string `ase.db`'s own `select()` accepts (a single
`cod_id`, a threshold like `total_score>0.7`, comma-combined conditions, ...),
resolved against `scores_db`. Here's a few examples:

```python
'geom_score=1.0'                        # pick out anything with perfect geometric score
'S,total_score>0.6'                     # anything containing sulphur and scoring higher than 0.6
'material_class=halides,bravais=CUB'    # all cubic halides
'n_major_facets_scored=3'               # anything that matches 3 major facets of the substrate
```

`refine` writes
`<reduced_formula>-<cod_id>/interfaces.db` for each material it resolves to,
one row per substrate/film termination combination of each match `score.py`
selected for that material, plus one plot per match and one PES/z-shift plot
pair per combination -- see `refine.py`'s module docstring for the full
folder layout. It returns the list of root folders written, one per refined
`cod_id`.

See [ASE's own database documentation](https://docs.ase-lib.org/ase/db/db.html#querying)
for the full query syntax `selection` accepts, and more generally for how to
query/manipulate `db` files yourself -- `scores_db` is already written
sorted by `total_score` descending, best match first (e.g. `ase db
scores.db -c +cod_id` to list them).

`substrate`, `matches_db`, and `mother_db` all default to `None`, which
reads them from `scores_db`'s own metadata (forwarded there by `score`
from `match`'s own metadata sidecar next to matches.csv) instead of
requiring them to be passed/tracked separately -- `mother_db` is used to
re-read each selected match's film atoms by `cod_id`, since matches.csv
itself carries no atoms. Pass any of the three explicitly to override.

### From the command line

Installing the package (see Installation above) also installs a `dbm`
console app, with one subcommand per stage:

```bash
# match the substrate on the whole mother database, with custom options
dbm match substrate.cif mother.db --max-area 500 -o matches.csv
# score all materials based on match quality
dbm score matches.csv -o scores
# energetical refinement of the matches of COD entry 2300704
dbm refine cod_id=2300704
```

`match`'s `-o`/`--output-csv` is a plain CSV filename (its only output).
`score`'s `-o`/`--output` is a basename -- results go to `<output>.db` and
`<output>.csv` in the current directory. `score`'s `--mother-db` overrides
the mother database it reads from matches.csv's own metadata sidecar by
default (see "As a library" above).

Each stage also still runs standalone, with identical flags:

```bash
python match.py substrate.cif mother.db --max-area 500 -o matches.csv
python score.py matches.csv -o scores
python refine.py 'cod_id=2300704'   # we recommend using quotes, especially for multiple-condition selections
```

`refine`'s substrate/matches-csv/mother-db are read from `--scores-db`'s
own metadata by default (see "As a library" above); pass
`--substrate`/`--matches-db`/`--mother-db` to override any of them.

To inspect a specific result afterward instead of running refine again,
`dbm refine --view <id>` visualizes one `interfaces.db` row's optimized
structure (`ase.visualize.view`) and opens its `PES.png`/`z_shift.png`
(`xdg-open`) -- `<id>` is that row's own `ase db` id. Run it from inside the
material's own `<reduced_formula>-<cod_id>/` folder, or pass
`--interfaces-db` to point at a different one:

```bash
cd TeMn-2300704 && dbm refine --view 3
```

### Converting between CSV and db

`dbm convert` reads `src`'s own file extension and writes `dst`'s
opposite: `.csv` -> `.db` or `.db` -> `.csv`. Any other extension pairing
(same extension, unrecognized, missing) raises an error.

```bash
dbm convert matches.csv matches.db   # inspect with `ase db`/`ase gui` instead of a CSV
dbm convert matches.db matches.csv   # and back
```

CSV rows carry no atomic structure, so `.csv -> .db` writes a placeholder
`Atoms('H')` for every row -- not real structural data, just a valid atom
for the db format to store `key_value_pairs` against. Whichever side has
a metadata sidecar (see `match`'s own `.meta.json`, above) or `.metadata`
carries over to the other -- this works for `scores.db` too, not just
matches.csv.

Run `dbm --help` / `dbm <stage> --help` (or equivalently `python
<stage>.py --help`) for the full list of options (miller index cutoffs,
strain/area tolerances, score weights, slab layers/vacuum,
...) -- both forms share the exact same argument definitions.

### Resuming a long matching run

Lattice matching on large databases can take a few hours, so 
`match`/`match.py` checkpoint after every mother-db row processed.
Pass `restart=True` (`-r`/`--restart` on the CLI) to resume from the
checkpoint file next to the output CSV instead of starting over.
Restart re-matches the checkpointed row itself too, not just the ones
after it -- first clearing any of its matches already in the output
CSV -- since there's no guarantee they were actually written to
disk before the previous run stopped. `score` has no equivalent --
scoring is cheap enough to just rerun from scratch (a few ms per
material; see the module docstring in `score.py` for the model).

### Full example

`example/run_example.py` runs the whole pipeline end to end -- match, score,
*and* refine -- against the `tests/data/` fixtures (symlinked in as
`mother-dataset.db`/`substrate.cif`):

```bash
cd example
python run_example.py
```

The refine step needs `matscipy`/`scikit-opt` installed too (see
Installation above).

## Configuration

`config.py` centralizes the tunable defaults for all three stages, under
three namespaces:

- `config.match` -- miller index cutoffs, strain/area tolerances, and the
  default output CSV filename for `match.py`.
- `config.score` -- score weights, facet-tier values, the
  substrate-compatible space groups/elements/major-facets, and the default
  output basename for `score.py`. These only parameterize `score.py`'s
  own built-in `default_scoring_function` -- a custom `scoring_fn` (see
  "As a library" above) ignores them entirely.
- `config.refine` -- slab layers/vacuum, and the number of points sampled
  in the interfacial-distance scan for `refine.py` (the scan's bounds, and
  the starting interfacial distance itself, are derived per termination
  combo from the ionic radii of the two atoms facing each other across the
  interface -- see `refine._contact_distance`/`_z_shift_range`), plus the
  fixed +/-`pes_colormap_bound` eV/Å² range each combo's `PES.png` colorbar
  is capped to (see `refine._run_surface_matching`).

`match`/`match.py`, `score`/`score.py`, and
`refine`/`refine.py` all read their keyword/CLI defaults straight
from here. Edit the values in `config.py` to change the defaults everywhere
at once (e.g. to tune scoring for a different substrate), or override any
of them per-call as a keyword argument, or per-run via the CLI flags.

## Testing

```bash
pytest
```

Runs against the small fixtures in `tests/data/`, including golden-fixture
comparisons that re-run real code and diff the result against a checked-in
db: `test_run_matching_data.py` re-runs the `MillerSearch` scan against
`tests/data/test-matches.db`, `test_ogre_custom.py` checks that
`interface_from_row`'s closed-form reconstruction matches
`InterfaceGenerator`'s own search, and `test_workflow_data.py` runs the
full `match` -> `score` workflow and diffs the result
against `tests/data/test-scores.db`, material by material. `score.py`'s
scoring model itself is also checked directly against known total_score
values in `tests/test_score.py`, without touching disk.

## Project layout

| file | purpose |
|------|---------|
| `match.py` | matching stage: `match`, the CLI, checkpoint/restart |
| `score.py` | scoring stage: `score`, the CLI |
| `refine.py` | refinement stage: `refine`, the CLI |
| `convert.py` | `convert`, the CLI -- CSV <-> db by file extension |
| `cli.py` | the `dbm` console app -- `match`/`score`/`refine`/`convert` subcommands |
| `ogre_custom.py` | `OgreInterface.MillerSearch` subclass used by `match.py`, plus the closed-form `Interface` reconstruction used by `refine.py` |
| `config.py` | shared tunable defaults for all three pipeline stages |
| `utils.py` | `db_to_csv`/`csv_to_db`, the standalone conversion helpers `convert.py` dispatches to |
| `example/` | runnable example against the `tests/data/` fixtures |
| `tests/` | pytest suite, including golden-fixture data |
