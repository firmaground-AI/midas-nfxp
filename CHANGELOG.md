# Changelog

## 0.1.0 (2026-10-27)

First release.

- `NfxpFile`: memory-mapped reader for GTS NX / MEC `.nfxp` results (2D, axisymmetric, 3D).
- `nodal`, `element` (centroid or corner, filter by element type), `history`, `components`.
- `table` / `dataframe`: all results of a step as columns per element family, with coordinates from the `.mec` deck.
- `export_csv` and the `midas-nfxp` command line (`info`, `export`, `history`).
- SRM and SAM factors of safety, active and deleted elements per step.
- Synthetic-file test suite; optional checks on real files via `NFXP_TEST_FILES`.