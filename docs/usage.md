# User guide and API reference

- [Concepts](#concepts)
- [Opening a file](#opening-a-file)
- [Steps and stages](#steps-and-stages)
- [Reading results](#reading-results)
- [Tables and CSV export](#tables-and-csv-export)
- [Mesh and model type](#mesh-and-model-type)
- [Command line](#command-line)
- [Recipes](#recipes)
- [Troubleshooting](#troubleshooting)

## Concepts

A GTS NX analysis writes its results as a sequence of **steps**. A step is one saved solution:
- one increment of a stage in nonlinear, construction stage and SRM analyses
- the single solution of a linear static analysis

Each step holds **result blocks**. A block is one component (for example `T2 TRANSLATION` or `EFFECTIVE NORMAL-YY`) at one **location**:

| location | meaning | how to read it |
|---|---|---|
| `node` | value per node | `nodal()` |
| `element` | value at the element centroid | `element()` |
| `corner` | value at element node *k* (labelled `NDk-` by GTS NX) | `element(..., corner=k)` |

Each block also has a **category**: `DISP`, `FORCE`, `STRESS`, `STRAIN` or `OTHER`. The same component name can appear in two categories, for example `NORMAL-XX` as a stress and as a strain. In that case pass `category=`.

Element results are written per **element family**. The family code comes from the result id:

| code | family | cards seen in `.mec` decks |
|---|---|---|
| 31 | 3D solid | CHEXA, CPENTA, CTETRA, CPYRAM |
| 43 | 2D plane strain | CQUAD4, CTRIA3, CQUAD8, CTRIA6 |
| 48 | axisymmetric | CQUADX, CTRIAX |
| 41 | shell / plate | CQUADR, CTRIAR |
| 42 | membrane / geogrid (3D) | CQUAD4 |
| 14 | beam | CBAR, CREBAR |
| 13 | truss / embedded truss | CROD, CREROD |
| 12 | spring | CBUSH |
| 21 | interface | CILINE2 (2D), CIQUAD4 (3D) |
| 23, 24 | pile interface, pile tip | CPILE2 |

These were collected from about 100 real models. A family that isn't listed is still read; it just keeps its numeric code. Always use the `type` column that `table()` adds from the `.mec` deck rather than assuming a card.

**Units** are the solver deck's: N, m, Pa (N/m² and N/m for shell/beam resultants) for MEC decks written by GTS NX, even if the GUI was set to kN. **Signs** follow GTS NX: compression is negative, and displacement components are in global X, Y, Z (T1, T2, T3). For axisymmetric models X is radial and Y axial. Reactions (`SPC_FORCE`) of axisymmetric models are for the **full 360°**, so no 2π factor is needed.

## Opening a file

```python
from midas_nfxp import NfxpFile, NfxpFormatError

r = NfxpFile(path, mec=None, find_mec=True)
```

| argument | |
|---|---|
| `path` | the `.nfxp` file |
| `mec` | path to the solver deck. The default looks for `<same name>.mec` next to the `.nfxp`. |
| `find_mec` | `False` skips the deck. Results still work; coordinates and element types don't. |

`NfxpFile` is a context manager (`with NfxpFile(p) as r:`). Call `r.close()` otherwise. On Windows an open file can't be deleted or overwritten by GTS NX.

It raises `NfxpFormatError` for an empty file, a file that isn't a MecSolver result, or an unsupported result type (Civil NX, eigenvalue, spectrum, moving load, heat of hydration, combinations).

It warns, without raising, in two cases:
- The file ends without the end-of-file record: the run was interrupted or is still writing. The last step may be incomplete.
- An unknown record follows valid results. Everything before it is indexed.

**Attributes**

| | |
|---|---|
| `r.steps` | list of `Step`, in file order |
| `r.stages` | `OrderedDict` stage number → stage label |
| `r.model_type` | `'2D'`, `'axisymmetric'` or `'3D'` (from the `.mec`). Without a deck: `'3D'` if any `T3 TRANSLATION` is non-zero, else `'2D or axisymmetric'`. |
| `r.complete` | `True` if the end-of-file record was found |
| `r.run_date` | date and time of the run as written by the solver |
| `r.mec_path` | the deck in use, or `None` |
| `r.mesh` | the parsed deck (see [Mesh and model type](#mesh-and-model-type)) |

## Steps and stages

```python
for s in r.steps:
    print(s)       # <Step 6: stage 7 'exca 3-SRM', INCR=9 (FOS=1.5156)>
```

| `Step` attribute | |
|---|---|
| `index` | position in `r.steps` (0-based) |
| `stage`, `stage_label` | stage number and name |
| `increment` | increment number within the stage, 1-based as in GTS NX |
| `step_no` | the solver's global step counter |
| `factor` | load factor of the increment (for SRM increments, the trial FoS) |
| `label` | GTS NX's increment label, e.g. `INCR=4 (LOAD=1.000)` |
| `fos` | SRM factor of safety: the solver's SRM record when present, else parsed from `FOS=` in the label. `None` for other steps. |
| `srm_fos` | the SRM record only |
| `sam` | `None`, or `dict(polyline_ids, fos)` for slope-analysis (SAM) results: one FoS per slip surface |
| `analysis` | solver's analysis code (undocumented) |

**Selecting a step.** Every method that takes `step=` accepts:

| key | selects |
|---|---|
| `-1`, `0`, `5` | by index; negative counts from the end |
| `'Excavation 2'` | the **last** increment of the stage with that label (case-insensitive) |
| `('Excavation 2', 3)` | increment 3 of that stage; `(stage, -1)` is the last one |
| `(4, 1)` | stage number 4, increment 1 |
| a `Step` object | itself |

```python
r.step('in situ')                    # Step object
r.stage_steps('exca 1-SRM')          # all steps of a stage (label or number)
last_of_each_stage = [r.stage_steps(st)[-1] for st in r.stages]
```

## Reading results

### `components(step=-1, location=None)`

Lists what a step contains as `(location, category, component)` tuples. Results can differ from stage to stage, because structural elements are only present once they are activated.

```python
for loc, cat, comp in r.components(-1, 'node'):
    print(loc, cat, comp)
```

### `nodal(component, step=-1, category=None)` → `Field(ids, values)`

```python
ids, uy = r.nodal('T2 TRANSLATION', step='exca 3')
ids, ry = r.nodal('T2 SPC_FORCE', step=0)          # reactions
```

`Field` is a named tuple, so `f.ids` and `f.values` work as well as unpacking. `ids` are GTS NX node ids (`int32`), and `values` are `float64` (the file stores `float32`). Component names are case-insensitive. A wrong name raises `KeyError` that lists similar names.

### `element(component, step=-1, category=None, corner=None, element_type=None)` → `Field(ids, values)`

```python
ids, s = r.element('EFFECTIVE NORMAL-ZZ', step=-1)                 # centroids
ids, s = r.element('EFFECTIVE NORMAL-ZZ', step=-1, corner=3)       # 3rd node of each element
ids, n = r.element('AXIAL FORCE', step=-1, category='FORCE', element_type='CBAR')
ids, m = r.element('BENDING MOMENT-XX', category='FORCE', element_type=['CQUAD4', 'CTRIA3'])
```

If several families carry the same component (for example pile and pile tip), their values are concatenated. `element_type` filters by the `.mec` card name and needs the deck.

### `history(component, id, location='node', category=None, corner=None, steps=None)` → `(step_indices, values)`

One node or element through the analysis. Steps where the entity is inactive (not yet constructed, or excavated) or the result is absent are skipped, so always use the returned step indices.

```python
idx, u = r.history('TOTAL TRANSLATION', 1702)
idx, s = r.history('EFFECTIVE NORMAL-YY', 845, location='element')
labels = [r.steps[i].stage_label for i in idx]
```

### Active and deleted elements

```python
r.active_elements('exca 2')     # element ids active in that step (or None if not written)
r.deleted_elements('exca 2')    # element ids removed so far
```

## Tables and CSV export

### `table(step=-1, location='node', corners=False)`

Returns every result at one location as columns: `{group: OrderedDict(column -> numpy array)}`.

- `location='node'` gives a single group, `'node'`.
- `location='element'` gives one group per element family. With the `.mec` deck the keys are `'<family>_<cards>'`, e.g. `'31_CHEXA+CPENTA+CTETRA'`, `'41_CQUAD4'`, `'14_CBAR'`. Without it, the keys are the integer family codes.
- Columns are `id`, then `type`, `x`, `y`, `z` (element families, with the deck) or `x`, `y`, `z` (nodes, with the deck), then one column per result named `CATEGORY:COMPONENT`, e.g. `STRESS:EFFECTIVE NORMAL-ZZ` or `FORCE:AXIAL FORCE`.
- For elements, `x, y, z` is the centroid: the mean of the element's node coordinates.
- A family is written in several chunks (per element shape and per result category). `table()` merges them, so a value that a chunk doesn't carry is `NaN`.
- `corners=True` adds the corner columns (`STRESS:ND1- EFFECTIVE NORMAL-ZZ`, …).

```python
t = r.table(-1, 'element')
soil = t[[k for k in t if str(k).startswith('31_')][0]]
deep = soil['z'] < -20
print(soil['STRESS:EFFECTIVE NORMAL-ZZ'][deep].mean())
```

### `dataframe(step=-1, location='node', corners=False)`

The same as pandas DataFrames: a single DataFrame for nodes, or a dict of DataFrames for elements. This needs pandas.

### `export_csv(out_dir, steps=None, corners=False)` → list of summary rows

For each step it writes:

| file | content |
|---|---|
| `<idx>_S<stage>_<stage label>_I<incr>_nodes.csv` | all nodal results |
| `<idx>_S<stage>_<stage label>_I<incr>_elements_<family>.csv` | all element results of a family |
| `<idx>_S<stage>_<stage label>_I<incr>_sam_fos.csv` | SAM FoS per slip surface, if any |
| `summary.csv` | one row per step: step, stage, stage_label, increment, label, factor, fos, sam_min_fos, max_total_translation, at_node |

`steps` is a list of step keys (see [Selecting a step](#steps-and-stages)). The default is all steps, which can be a lot of files for a long nonlinear run.

```python
r.export_csv('out', steps=[r.stage_steps(st)[-1].index for st in r.stages])
```

## Mesh and model type

`r.mesh` (parsed from the `.mec` deck on first use) is a dict:

| key | |
|---|---|
| `node_ids` | sorted node ids |
| `xyz` | `(n, 3)` coordinates, in the same order |
| `element_type` | `{element id: card name}`, e.g. `'CHEXA'` |
| `connectivity` | `{element id: [node ids]}` |
| `model_type` | `'3D'` if there are solid cards; `'axisymmetric'` if there are CQUADX/CTRIAX; otherwise `'2D'` (or `'3D'` if any node has z ≠ 0) |

`read_mec(path)` returns the same dict without opening a result file.

## Command line

```text
midas-nfxp info    NFXP [--mec MEC] [--step STEP]
midas-nfxp export  NFXP [--mec MEC] [--out DIR] [--steps all|last|last-of-stage|0,5,9] [--corners]
midas-nfxp history NFXP COMPONENT ID [--element] [--category CAT]
midas-nfxp --version
```

- `info` prints the model type, run date, every step (with FoS and SAM FoS) and the components of `--step` (default: last).
- `export` writes the CSVs described above. The default is `last-of-stage`, and the default folder is `<nfxp name>_results`.
- `history` prints CSV (step, stage, increment, label, value) to stdout. Redirect it to a file with `> node1702.csv`.

## Recipes

### Pile axial force along the length (3D)

```python
import numpy as np
t = r.table(step=-1, location='element')
piles = t['14_CBAR']
x, y, z = piles['x'], piles['y'], piles['z']
N = piles['FORCE:AXIAL FORCE'] / 1e3                       # kN
for px, py in sorted(set(zip(np.round(x, 2), np.round(y, 2)))):
    on = (np.round(x, 2) == px) & (np.round(y, 2) == py)
    o = np.argsort(-z[on])                                  # head to toe
    print((px, py), np.c_[z[on][o], N[on][o]].round(2))
```

`examples/pile_axial_force.py` does the same, writes a CSV and plots it. For the force at the element ends instead of the centroid, read `FORCE:ND1- AXIAL FORCE` and `FORCE:ND2- AXIAL FORCE` with `corners=True`. The node order is in `r.mesh['connectivity'][element_id]`.

### Wall or lining bending moment (2D)

```python
t = r.table(step='exca 3', location='element')
beams = t[[k for k in t if str(k).startswith('14_')][0]]
o = np.argsort(-beams['y'])
for yy, m in zip(beams['y'][o], beams['FORCE:BENDING MOMENT-Y'][o] / 1e3):
    print('%.2f  %.1f kN.m/m' % (yy, m))
```

In 2D models the in-plane beam results are `BENDING MOMENT-Y` and `SHEAR FORCE-Z` (element axes); `BENDING MOMENT-Z` and `SHEAR FORCE-Y` are zero.

### Raft or slab moments (3D shells)

```python
ids, mxx = r.element('BENDING MOMENT-XX', step=-1, category='FORCE')
print('max |Mxx| = %.0f kN.m/m' % (np.abs(mxx).max() / 1e3))
```

### Settlement trough along a line

```python
t = r.table(step=-1)['node']
on = np.isclose(t['z'], 0) & np.isclose(t['y'], 0)           # surface, along x at y = 0
o = np.argsort(t['x'][on])
trough = np.c_[t['x'][on][o], t['DISP:T3 TRANSLATION'][on][o] * 1e3]
```

### Factor of safety per stage

```python
for st, label in r.stages.items():
    s = r.stage_steps(st)[-1]
    if s.fos is not None:
        print(label, 'SRM', round(s.fos, 3))
    if s.sam is not None:
        print(label, 'SAM', round(s.sam['fos'].min(), 3))
```

The SRM FoS is the factor of the last converged increment of the SRM stage. GTS NX also writes it in the increment label (`INCR=19 (FOS=2.5520)`).

### Base reaction (axisymmetric)

```python
ids, fy = r.nodal('T2 SPC_FORCE', step=0)
y = r._node_xyz(ids)[:, 1]
print(fy[np.isclose(y, y.min())].sum())     # full 360 deg, N
```

### A node through all stages

```python
idx, u = r.history('T3 TRANSLATION', 1937)
for i, v in zip(idx, u):
    print(r.steps[i].stage_label, r.steps[i].increment, round(v * 1e3, 2))
```

### pandas

```python
nodes = r.dataframe(-1, 'node')
elems = r.dataframe(-1, 'element')          # dict of DataFrames per family
nodes.nsmallest(5, 'DISP:T3 TRANSLATION')
```

## Troubleshooting

| symptom | cause and fix |
|---|---|
| `KeyError: no node result 'UY'` | Use GTS NX's names: `T2 TRANSLATION`. Run `components()` or `midas-nfxp info` to list them. The error message also lists similar names. |
| `KeyError: ... exists as ['STRAIN', 'STRESS']; pass category=` | The name exists in two categories. Add `category='STRESS'`. |
| `RuntimeError: no .mec deck` | Coordinates, `element_type` and `type` columns need the deck. Put the `.mec` next to the `.nfxp` with the same base name, or pass `mec=`. |
| `x, y, z` are `NaN` for some rows | The deck isn't the one this result came from (the mesh changed). Re-run the analysis, or point `mec=` at the matching deck. |
| `NfxpFormatError: ... not supported` | Civil NX, eigenvalue, spectrum, moving-load, heat-of-hydration or combination results. These aren't covered. |
| warning: `file ends without the end-of-file record` | The run is still writing or crashed. Earlier steps are fine, but the last one may be partial. |
| Values look 1000× too big | MEC decks are in N, m, Pa. Divide by 1e3 for kN and kPa. |
| A structural element is missing in early steps | It isn't active yet. Check `r.active_elements(step)`. |
| `IndentationError` when pasting into the Python prompt | Lines at the `>>>` prompt must start in column 0. Paste multi-line recipes into a `.py` file and run that instead. |
