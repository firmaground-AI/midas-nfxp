# MIDAS NFXP - API

Read MIDAS **GTS NX** result files (`.nfxp`, written by the MEC solver) in Python: displacements, reactions, stresses, strains, structural forces and factors of safety for **2D (plane strain), axisymmetric and 3D** models, every saved stage and increment.

```python
from midas_nfxp import NfxpFile

with NfxpFile('Excavation_CS.nfxp') as r:          # picks up Excavation_CS.mec next to it
    print(r)                                       # <NfxpFile ...: 2D, 7 stages, 7 steps>
    ids, uy = r.nodal('T2 TRANSLATION', step=-1)   # last saved step
    ids, sv = r.element('EFFECTIVE NORMAL-YY', step='exca 2')
    steps, u = r.history('TOTAL TRANSLATION', 1702)
    r.export_csv('results')                        # one CSV per step and element family
```

- **One file, one dependency.** `midas_nfxp.py` needs only NumPy. pandas and matplotlib are optional.
- **Fast on large models.** The file is memory-mapped and indexed when you open it; values are read only when you ask for them. An 8 GB, 250-step 3D result opens in about 2 seconds.
- **Coordinates and element types** come from the `.mec` solver deck that GTS NX writes next to the result file. They are used for centroids, the `x, y, z` columns, filtering by element type and detecting the model type.
- **Command line** for a quick look or a bulk CSV export, with no Python code needed.

> **Unofficial.** This project is not affiliated with or endorsed by MIDAS IT. The `.nfxp` format is undocumented; the reader is based on reverse engineering and was checked against the solver's own log output (see [Validation](docs/validation.md)). Verify critical numbers in GTS NX.

## Install

Requires Python 3.8+ and NumPy.

```bash
pip install git+https://github.com/firmaground-AI/midas-nfxp.git
```

Or copy `midas_nfxp.py` next to your script. Optional extras: `pip install pandas matplotlib`.

## Getting the files

Run the analysis in GTS NX as usual. The solver writes `<model>_<analysis>.nfxp` (results) and `<model>_<analysis>.mec` (the solver deck) into the model's folder. Keep them together with the same base name. Without the `.mec` everything still works except coordinates, element types and the model type.

## Quick start

```python
import numpy as np
from midas_nfxp import NfxpFile

r = NfxpFile(r'C:\models\Piled raft_analysis.nfxp')
print(r.model_type)                 # '3D'
for s in r.steps:                   # one Step per saved increment
    print(s.index, s.stage_label, s.increment, s.label, s.fos)

# Settlement at the end of the last stage
ids, uz = r.nodal('T3 TRANSLATION', step=-1)
k = np.argmin(uz)
print('max settlement %.1f mm at node %d' % (uz[k] * 1e3, ids[k]))

# Every element result of one step as columns, with centroids, per element family
t = r.table(step=-1, location='element')
piles = t['14_CBAR']                # keys like '31_CHEXA+CTETRA', '41_CQUAD4', '14_CBAR'
print(piles['z'][:5], piles['FORCE:AXIAL FORCE'][:5])
r.close()
```

**Steps** can be selected by index (`-1` is the last), stage label (its last increment), stage number, or a `(stage, increment)` tuple. Values are in the solver's units: **N, m and Pa** for MEC decks, even when the GTS NX interface shows kN. Compression is negative.

## Command line

```bash
midas-nfxp info   model.nfxp                       # model type, stages, steps, FoS, components
midas-nfxp export model.nfxp --steps last-of-stage # CSV per step and element family + summary.csv
midas-nfxp history model.nfxp "T2 TRANSLATION" 1702
midas-nfxp history model.nfxp "AXIAL FORCE" 5012 --element --category FORCE
```

If you didn't install the package, use `python midas_nfxp.py ...` instead of `midas-nfxp ...`.

## What it reads

| | |
|---|---|
| Model types | 2D plane strain, axisymmetric, 3D |
| Analyses | Linear static, nonlinear static, construction stage, SRM (strength reduction) and SAM (slope stability by stress analysis), with every saved increment. Seepage and consolidation stages use the same records and are read too, but haven't been checked separately against GTS NX. |
| Locations | Nodes, element centroids, element corner nodes (`corner=k`) |
| Elements | Solids, plane strain / axisymmetric, shells, beams, trusses, interfaces, pile and pile-tip elements, each kept as its own family |
| Extras | Active and deleted elements per stage, SRM FoS, SAM FoS per slip surface, run date |

It does **not** read Civil NX results, or eigenvalue, response-spectrum, moving-load, heat-of-hydration or load-combination result files. These use other record types, and opening one raises `NfxpFormatError` instead of returning wrong numbers.

## Documentation

- [User guide and API reference](docs/usage.md): every class and method, step selection, units, and recipes (pile axial force, wall bending moment, stress against depth, stage histories, pandas)
- [Result component names](docs/result-components.md): the component names you pass to `nodal()` and `element()`, by model type
- [File format notes](docs/file-format.md): the `.nfxp` record layout as far as it is known
- [Validation](docs/validation.md): how the reader was checked, and the known limitations

## Examples

```bash
python examples/quickstart.py        model.nfxp
python examples/pile_axial_force.py  model.nfxp --step "Surface load" --plot
python examples/factor_of_safety.py  model.nfxp
python examples/stress_vs_depth.py   model.nfxp 0 --plot
```

## Tests

```bash
pip install pytest
python -m pytest
```

The test suite builds small synthetic `.nfxp` and `.mec` files, so it needs no model files. To also check your own results, list them in `NFXP_TEST_FILES` (separated by `;` on Windows and `:` elsewhere):

```bash
set NFXP_TEST_FILES=C:\models\a.nfxp;C:\models\b.nfxp
python -m pytest tests/test_real_files.py -v
```

## Reporting a problem

If a file fails to open or gives numbers that disagree with GTS NX, please open an issue with:
- the output of `midas-nfxp info model.nfxp`
- the GTS NX version and analysis type
- the value GTS NX shows and the value the reader gives

Please don't attach client models.

## License

[MIT](LICENSE). GTS NX and MIDAS are trademarks of MIDAS Information Technology Co., Ltd.
