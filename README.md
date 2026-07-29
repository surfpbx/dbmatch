# db-ogre-match

Match every structure in a "mother" database against a fixed substrate,
finding epitaxial interfaces via [OgreInterface](https://github.com/DerekDardzinski/OgreInterface),
then score each material by how well and how compatibly it matches.

The pipeline has two stages, each reading and writing plain [ASE
databases](https://wiki.fysik.dtu.dk/ase/ase/db/db.html) and CSV files:

```
mother_db + substrate  --match.py-->  matches db/csv  --score.py-->  scores db/csv
```

- **`match.py`** scans every structure in the mother database against the
  substrate with `OgreInterface`'s `MillerSearch`, writing one row per
  facet-pairing match found.
- **`score.py`** reads a matches database back in, grouped by material
  (`cod_id`), and assigns each material a `total_score` based on how many
  facets it hits, how compatible its space group is, and whether it
  contains substrate-compatible elements.

## Installation

```bash
pip install -e .
```

This installs `ase` and `tqdm`, the only two dependencies declared in
`pyproject.toml`.

### OgreInterface

`OgreInterface` is **not** on PyPI and is not declared as a dependency here
-- it has to be installed separately:

```bash
git clone https://github.com/DerekDardzinski/OgreInterface.git
pip install -e ./OgreInterface
```

That pulls in its own (heavier) dependencies: `pymatgen`, `scipy`,
`matscipy`, `matplotlib`, etc. Alternatively, if you already have a
checkout of `OgreInterface` somewhere, you can skip the install and just
put its parent directory on `PYTHONPATH`:

```bash
export PYTHONPATH="/path/to/OgreInterface:$PYTHONPATH"
```

Either way, `import OgreInterface` has to succeed before `db_ogre_match`
will work -- `miller_custom.py` imports directly from it.

**Known issue:** a fresh `pip install -e ./OgreInterface` can resolve a
`numpy`/`pandas` combination that's binary-incompatible (`OgreInterface`
pins `pandas<=2.0.0`, built pre-NumPy-2.0, but doesn't cap `numpy`, and
`pymatgen`/`matscipy` pull in `numpy>=2`), which fails on
`import OgreInterface` with `numpy.dtype size changed, may indicate
binary incompatibility`. This is unrelated to anything `db_ogre_match`
actually needs from `OgreInterface` -- pandas only gets touched because
`OgreInterface/data/ionic_radii.py` reads a CSV via pandas at import time,
in a module chain that `miller.py` pulls in regardless of whether the
function that uses it (`estimate_atomic_radius`, part of the unrelated
Lennard-Jones surface-matching code) is ever called. If you hit this,
pin `numpy<2` after installing `OgreInterface`; note that conflicts with
`matscipy`'s own `numpy>=2.0.0` requirement, so treat it as a workaround,
not a real fix -- this needs fixing upstream in `OgreInterface`.

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
    mother_db='mother.db',
    substrate='substrate.cif',
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

### From the command line

```bash
python match.py mother.db substrate.cif --max-area 500 -o matches.db --output-csv matches.csv
python score.py db/matches.db -o scores.db --output-csv scores.csv
```

Run `python match.py --help` / `python score.py --help` for the full list
of options (miller index cutoffs, strain/area tolerances, score weights,
...).

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

`config.py` centralizes the tunable defaults for both stages, under two
namespaces:

- `config.match` -- miller index cutoffs, strain/area tolerances, and
  default output filenames for `match.py`.
- `config.score` -- score weights, facet-tier values, the
  substrate-compatible space groups/elements/major-facets, and default
  output filenames for `score.py`.

Both `match_database`/`match.py` and `score_materials`/`score.py` read
their keyword/CLI defaults straight from here. Edit the values in
`config.py` to change the defaults everywhere at once (e.g. to tune
scoring for a different substrate), or override any of them per-call as a
keyword argument to `match_database`/`score_materials`, or per-run via the
CLI flags.

## Testing

```bash
pytest
```

Runs against the small fixtures in `tests/data/`, including a golden-
fixture comparison (`test_run_matching_data.py`) that re-runs the real
`MillerSearch` scan and checks it against `tests/data/test-matches.db`.
`score.py`'s scoring model itself is checked directly against known
total_score values in `tests/test_score.py`, without touching disk.

## Project layout

| file | purpose |
|------|---------|
| `match.py` | matching stage: `match_database`, the CLI, checkpoint/restart |
| `score.py` | scoring stage: `score_materials`, the CLI |
| `miller_custom.py` | `OgreInterface.MillerSearch` subclass used by `match.py` |
| `config.py` | shared tunable defaults for both stages |
| `utils.py` | `db_to_csv`, a standalone db -> CSV dump helper |
| `example/` | runnable example against the `tests/data/` fixtures |
| `tests/` | pytest suite, including golden-fixture data |
