# Validation and limitations

## How the reader was checked

GTS NX doesn't document the `.nfxp` format, so every check compares what the reader returns with something the solver reports independently, or with physics.

| check | models | result |
|---|---|---|
| Maximum translation per saved step vs. the solver `.log` | 15 models (2D, axisymmetric, 3D; construction stage, nonlinear, SRM) | **177 of 177 steps match** |
| SRM factor of safety vs. the `.log` | staged 2D excavation with SRM after each stage | 2.552 / 2.291 / 1.516: match |
| SAM factor of safety vs. GTS NX | same model | 2.500 / 2.397 / 1.561: match |
| Reactions only at constrained nodes | all 15 | non-zero `SPC_FORCE` occurs only at `SPC1` nodes of the deck |
| Displacement continuity | all 15 | nodal displacements are continuous across element boundaries; there are no jumps from misassigned ids |
| In-situ vertical stress against depth | level-ground models | linear, R² ≥ 0.99, slope 16–25 kN/m³ |
| Axisymmetric base reaction | axisymmetric excavation | 2.827e10 N against a full-revolution weight of 2.819e10 N, so reactions are for 360° |
| Element tables vs. per-component reads | 3D model with 11 solid chunks | `table()` columns are identical to `element()` |
| Opening every result file on one engineering workstation | 240 files | 207 GTS NX files read (116 2D, 81 3D, 9 axisymmetric, 1 without a deck). The other 33 (Civil NX, eigenvalue, spectrum, moving load, heat of hydration, combination) raise `NfxpFormatError` as intended. |

The automated tests in `tests/` write synthetic files with known values and check that every API path returns them: indexing, step selection, chunk merging, FoS, corner values, interrupted files and the CLI. `tests/test_real_files.py` runs structural sanity checks on any real files you point it at.

## Known limitations

- **Format version.** The checks above used MEC solver V1.9.8. Older or newer solvers may add record types. The reader then stops at the unknown record with a warning, and keeps every step before it.
- **Unsupported analyses:** Civil NX, eigenvalue, response spectrum, moving load, heat of hydration and load combinations.
- **Seepage and consolidation** stages are read with the same record types, but their specific components (heads, flow, pore pressure dissipation) haven't been compared with GTS NX one by one.
- **Element axes** (ECSYS records) are skipped. Beam and shell forces are given in element axes as GTS NX writes them, and the reader doesn't rotate them to global axes.
- **Centroids** are the mean of the element's node coordinates from the deck. For higher-order elements this includes the mid-side nodes.
- **Precision:** values are stored as `float32` (about 7 significant figures) and returned as `float64`.
- **Units** are those of the solver deck (N, m, Pa for MEC decks), whatever the GUI's unit setting.
- **The `.log` file** isn't a reliable cross-reference for ids (see [file-format.md](file-format.md)).

If you find a file the reader gets wrong, please open an issue (see the README). Include the `midas-nfxp info` output, but not the model.
