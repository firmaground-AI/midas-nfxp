# Result component names

Pass these names to `nodal()` / `element()`. They are GTS NX's own labels, matched case-insensitively. Which ones exist depends on the model type, the element families present and the output options set in GTS NX. The lists below come from real 3D, 2D and axisymmetric models. To see what a particular file has, run:

```bash
midas-nfxp info model.nfxp --step -1
```

or call `r.components(step)` in Python.

## Nodal results (all model types)

| category | components |
|---|---|
| `DISP` | `TOTAL TRANSLATION`, `T1 TRANSLATION`, `T2 TRANSLATION`, `T3 TRANSLATION`, `T12 TRANSLATION`, `T23 TRANSLATION`, `T13 TRANSLATION`, `TOTAL ROTATION`, `R1 ROTATION`, `R2 ROTATION`, `R3 ROTATION` |
| `FORCE` | `TOTAL APPLIED_FORCE`, `T1/T2/T3 APPLIED_FORCE`, `TOTAL APPLIED_MOMENT`, `R1/R2/R3 APPLIED_MOMENT`, `TOTAL SPC_FORCE`, `T1/T2/T3 SPC_FORCE`, `TOTAL SPC_MOMENT`, `R1/R2/R3 SPC_MOMENT` |

`T1/T2/T3` are global X/Y/Z. In 2D and axisymmetric models, Y is vertical (axial) and X is horizontal (radial); in 3D, Z is usually vertical. `SPC_FORCE` are support reactions. They are non-zero only at constrained nodes and, for axisymmetric models, are for the full 360°.

## Soil and solid elements (centroid; also at corners as `NDk- <name>`)

**3D solids (family 31)**

| category | components |
|---|---|
| `STRESS` | `TOTAL NORMAL-XX/YY/ZZ`, `EFFECTIVE NORMAL-XX/YY/ZZ`, `SHEAR-XY`, `SHEAR-YZ`, `SHEAR-ZX`, `PORE STRESS`, `EXCESSIVE PORE STRESS`, `PRINCIPAL-A/B/C`, `PRINCIPAL-A/B/C COSINE-X/Y/Z`, `VON MISES`, `MAX SHEAR`, `MEAN EFFFECTIVE` (sic, three F's), `MEAN TOTAL`, `SAFETY FACTOR`, `EQUIVALENT`, `PLASTIC STATUS` |
| `STRAIN` | `NORMAL-XX/YY/ZZ`, `SHEAR-XY/YZ/ZX`, `PRINCIPAL-A/B/C`, `PRINCIPAL-A/B/C COSINE-X/Y/Z`, `VON MISES`, `MAX SHEAR`, `VOLUMETRIC`, `DEVIATORIC`, `EFFECTIVE PLASTIC`, `EQUIVALENT` |

**2D plane strain (family 43)**

| category | components |
|---|---|
| `STRESS` | `TOTAL NORMAL-XX/YY/ZZ`, `EFFECTIVE NORMAL-XX/YY/ZZ`, `SHEAR-XY`, `PORE STRESS`, `EXCESSIVE PORE STRESS`, `MAJOR PRINCIPAL`, `MINOR PRINCIPAL`, `VON MISES`, `MAX SHEAR`, `ZERO SHEAR ANGLE`, `MEAN EFFECTIVE`, `MEAN TOTAL`, `SAFETY FACTOR`, `EQUIVALENT`, `PLASTIC STATUS` |
| `STRAIN` | `NORMAL-XX`, `NORMAL-YY`, `SHEAR-XY`, `MAJOR PRINCIPAL`, `MINOR PRINCIPAL`, `VON MISES`, `MAX SHEAR`, `VOLUMETRIC`, `DEVIATORIC`, `EFFECTIVE PLASTIC`, `EQUIVALENT` |

Note that the 3D label is `MEAN EFFFECTIVE` and the 2D one is `MEAN EFFECTIVE`, exactly as GTS NX writes them.

**Axisymmetric (family 48)**

| category | components |
|---|---|
| `STRESS` | `NORMAL-XX` (radial), `NORMAL-YY` (axial), `SHEAR-XY`, `HOOP`, `MAJOR PRINCIPAL`, `MINOR PRINCIPAL`, `VON MISES`, `MAX SHEAR`, `SAFETY FACTOR`, `EQUIVALENT`, `PLASTIC STATUS` |
| `STRAIN` | `NORMAL-XX`, `NORMAL-YY`, `SHEAR-XY`, `HOOP`, `MAJOR PRINCIPAL`, `MINOR PRINCIPAL`, `VON MISES`, `MAX SHEAR`, `EFFECTIVE PLASTIC`, `EQUIVALENT`, `EFFECTIVE CREEP` |

In axisymmetric models the stress names have no `TOTAL`/`EFFECTIVE` prefix. `NORMAL-XX` exists as both a stress and a strain, so pass `category='STRESS'` or `'STRAIN'`.

## Structural elements

**Beams (family 14), `FORCE`:** `AXIAL FORCE`, `SHEAR FORCE-Y`, `SHEAR FORCE-Z`, `TORQUE`, `BENDING MOMENT-Y`, `BENDING MOMENT-Z`. These are in element axes. In 2D models the in-plane values are `BENDING MOMENT-Y` and `SHEAR FORCE-Z`. The end values are available as corners: `ND1- AXIAL FORCE` and `ND2- AXIAL FORCE`.

**Beam stresses (family 14), `OTHER`:** `AXIAL`, `POINT C/D/E/F`, `COMBINED MINIMUM`, `COMBINED MAXIMUM`, `SHEAR-XY`, `SHEAR-XZ`, `TORSIONAL`, `VON MISES`.

**Shells (family 41)**

| category | components |
|---|---|
| `FORCE` | `MEMBRANE FORCE-XX/YY/XY`, `BENDING MOMENT-XX/YY/XY`, `TRANSVERSE SHEAR FORCE-XZ/YZ` (per unit width) |
| `STRESS` | `NORMAL-XX/YY/ZZ`, `SHEAR-XY` and `MAJOR/MINOR PRINCIPAL`, `VON MISES`, `MAX SHEAR`, `ZERO SHEAR ANGLE`, `FIBER DISTANCE`, each at `BOTTOM`, `MIDDLE` and `TOP` (e.g. `NORMAL-XX TOP`), plus `MAJOR PRINCIPAL MAXIMUM`, `MINOR PRINCIPAL MINIMUM`, `VON MISES MAXIMUM`, `MAX SHEAR MAXIMUM` |
| `STRAIN` | the same names as the stresses, without `FIBER DISTANCE` |

**2D plates, geogrids and membranes:** `MEMBRANE FORCE-XX/YY/XY`, `TRANSVERSE SHEAR FORCE-XZ/YZ` (`FORCE`).

**Interfaces (family 21):**
- `FORCE`: `NORMAL-X`, `TANGENTIAL-Y`, `TANGENTIAL-Z`
- `STRAIN`: `NORMAL RELATIVE DISPLACEMENT-X`, `TANGENTIAL RELATIVE DISPLACEMENT-Y/Z`

**Springs (family 12):** `FORCE-X/Y/Z`, `MOMENT-X/Y/Z`, `TRANSLATION-…`, `ROTATION-X/Y…`, `MAX TRANSLATION`, `MAX ROTATION`.

## Names that appear in several places

A name can exist in more than one category or element family, for example `NORMAL-XX`, `SHEAR-XY`, `VON MISES` or `EQUIVALENT`.

- Several categories: `element()` raises `KeyError ... pass category=`.
- Several families in the same category: `element()` returns all of them concatenated. Filter them with `element_type=`, or use `table()`, which keeps the families apart.
