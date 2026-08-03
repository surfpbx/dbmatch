# db-ogre-match

Match every structure in a "mother" database against a fixed substrate,
finding epitaxial interfaces via [OgreInterface](https://github.com/DerekDardzinski/OgreInterface),
score each material by how well and how compatibly it matches, then
energetically refine the best matches of a chosen material.

The pipeline has three stages, each reading and writing plain [ASE
databases](https://docs.ase-lib.org/ase/db/db.html) and CSV files (`refine.py`
writes only a database, per material):

```
substrate + mother_db  --match.py -->  matches db/csv  --score.py -->  scores db/csv  --refine.py -->  <formula>-<cod_id>/interfaces.db
```

- **`match.py`** scans every structure in the mother database against the
  substrate with `OgreInterface`'s `MillerSearch`, writing one row per
  facet-pairing match found.
- **`score.py`** reads a matches database back in, grouped by material
  (`cod_id`), and assigns each material a `total_score` based on how many
  facets it hits, how compatible its space group is, and whether it
  contains substrate-compatible elements.
- **`refine.py`** takes a single material's `cod_id`, rebuilds the matches
  `score.py` selected for it, and runs `OgreInterface`'s ionic surface
  matching (in-plane + interfacial-distance energy optimization) over every
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
from db_ogre_match import match_database
from db_ogre_match.score import score_materials

match_database(
    substrate='substrate.cif',
    mother_db='mother.db',
    max_substrate_index=1,
    max_film_index=1,
    max_strain=0.05,
    max_area=500,
    output_db='matches.db',      # written to db/matches.db
    output_csv='matches.csv',    # written to csv/matches.csv
    restart=False,
)

score_materials(
    matches_db='db/matches.db',
    output_csv='scores.csv',     # written to csv/scores.csv
    output_db='scores.db',       # written to db/scores.db
)
```

Both functions create `db/` and `csv/` under the current directory if
they don't already exist, and write output there. `score_materials` scores
and writes one material at a time as it goes, rather than buffering
everything in memory -- rows land in whatever `cod_id` order the source db
yields, not sorted by score.

Refining a single material's matches, once it's been scored:

```python
from db_ogre_match.refine import refine_material

refine_material(
    cod_id=2300704,
    scores_db='db/scores.db',
)
```

`refine_material` writes `<reduced_formula>-<cod_id>/interfaces.db`, one row
per substrate/film termination combination of each match `score.py` selected
for that material, plus one plot per match and one PES/z-shift plot pair per
combination -- see `refine.py`'s module docstring for the full folder layout.
Both `substrate` and `matches_db` default to `None`, which reads them from
`scores_db`'s own metadata (forwarded there by `score_materials` from
`match_database`) instead of requiring them to be passed/tracked separately;
pass either explicitly to override.

### From the command line

Installing the package (see Installation above) also installs a `dbm`
console app, with one subcommand per stage:

```bash
dbm match substrate.cif mother.db --max-area 500 -o matches.db --output-csv matches.csv
dbm score db/matches.db -o scores.db --output-csv scores.csv
dbm refine 2300704
```

Each stage also still runs standalone, with identical flags -- note
`match.py`'s positional order is `substrate` then `mother_db`, same as
`dbm match`'s:

```bash
python match.py substrate.cif mother.db --max-area 500 -o matches.db --output-csv matches.csv
python score.py db/matches.db -o scores.db --output-csv scores.csv
python refine.py 2300704
```

`refine`'s substrate/matches db are read from `--scores-db`'s own metadata
by default (see "As a library" above); pass `--substrate`/`--matches-db` to
override either.

Run `dbm --help` / `dbm <stage> --help` (or equivalently `python
<stage>.py --help`) for the full list of options (miller index cutoffs,
strain/area tolerances, score weights, slab layers/vacuum/interfacial-distance,
...) -- both forms share the exact same argument definitions.

### Resuming a long matching run

`match_database`/`match.py` checkpoint after every mother-db row processed.
Pass `restart=True` (`-r`/`--restart` on the CLI) to resume from the
checkpoint file next to the output database instead of starting over.
`score_materials` has no equivalent -- scoring is cheap enough to just
rerun from scratch (a few ms per material; see the module docstring in
`score.py` for the model).

### Full example

`example/run_example.py` runs the whole pipeline end to end against the
`tests/data/` fixtures (symlinked in as `mother-dataset.db`/
`substrate.cif`):

```bash
cd example
python run_example.py
```

## Configuration

`config.py` centralizes the tunable defaults for all three stages, under
three namespaces:

- `config.match` -- miller index cutoffs, strain/area tolerances, and
  default output filenames for `match.py`.
- `config.score` -- score weights, facet-tier values, the
  substrate-compatible space groups/elements/major-facets, and default
  output filenames for `score.py`.
- `config.refine` -- slab layers/vacuum, starting interfacial distance, and
  the interfacial-distance scan bounds for `refine.py`.

`match_database`/`match.py`, `score_materials`/`score.py`, and
`refine_material`/`refine.py` all read their keyword/CLI defaults straight
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
full `match_database` -> `score_materials` workflow and diffs the result
against `tests/data/test-scores.db`, material by material. `score.py`'s
scoring model itself is also checked directly against known total_score
values in `tests/test_score.py`, without touching disk.

## Project layout

| file | purpose |
|------|---------|
| `match.py` | matching stage: `match_database`, the CLI, checkpoint/restart |
| `score.py` | scoring stage: `score_materials`, the CLI |
| `refine.py` | refinement stage: `refine_material`, the CLI |
| `cli.py` | the `dbm` console app -- `match`/`score`/`refine` subcommands |
| `ogre_custom.py` | `OgreInterface.MillerSearch` subclass used by `match.py`, plus the closed-form `Interface` reconstruction used by `refine.py` |
| `config.py` | shared tunable defaults for all three stages |
| `utils.py` | `db_to_csv`, a standalone db -> CSV dump helper |
| `example/` | runnable example against the `tests/data/` fixtures |
| `tests/` | pytest suite, including golden-fixture data |
