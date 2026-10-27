# `.nfxp` file format notes

The format is undocumented. This layout was worked out from files written by the MEC solver (V1.9.8) and checked against the solver's `.log` and `.out` summaries. If the reader meets a record type it doesn't know, it stops and warns rather than guessing.

All values are **little-endian**. Integers are `int32`, result values `float32`, FoS values `float64`. Labels are 80-byte, NUL-padded **UTF-8** (stage names can be non-ASCII).

## Overall layout

```text
header                      0x00 .. 0x51
STAGE 1                     10017 record
  INCREMENT 1               13 record          (absent in linear static)
    IDS kind 1              node ids
    RESULT ...              nodal results
    IDS kind 3              active elements
    DELETED                 elements removed so far
    ECSYS                   element axes        (structural families)
    IDS kind 2              element ids of family A, chunk 1
    RESULT ...              its results (centroid, then corner NDk-)
    IDS kind 2              family A, chunk 2 / family B ...
    RESULT ...
    SRM FOS / SAM FOS       (SRM / SAM stages)
  INCREMENT 2
    ...
STAGE 2
  ...
EOF                         24 B
```

A result block applies to the most recent id list of the matching kind: kind 1 for nodal blocks, kind 2 for element and corner blocks.

## Records

| type | name | size | layout |
|---|---|---|---|
| — | header | 0x52 B | `MecSolver` at 0x10; run date and time at 0x2C..0x3B |
| 10017 | STAGE | 108 B | `type, 0, next_stage_offset, ?, 12, stage_no, analysis_code, label[80]`. `next_stage_offset = -1` marks the **last stage**, not the end of the file. |
| 13 | INCREMENT | 100 B | `type, global_step_no, factor(f32), analysis_code, label[80], increment_index` (0-based). The label reads like `INCR=4 (LOAD=1.000)` or `INCR=19 (FOS=2.5520)`. |
| 1 / 2 / 3 | IDS | 8 + 4n | `kind, n, ids[n]`. 1 = node ids; 2 = element ids of the family chunk whose results follow; 3 = active elements in this step |
| 4 | DELETED | 8 + 4n | `4, n, ids[n]`: elements removed up to this step |
| 10003 | ECSYS | 8 + 24n | `10003, n, n × 6 f32`: element axis vectors |
| 14 | RESULT | 152 + 4n | 18-int header, `label[80]`, `values[n]` (f32). See below. |
| 10025 | SAM FOS | 24 + 12n | `10025, ?, post_id, ?, ?, n, polyline_ids[n] (i32), fos[n] (f64)` |
| 10032 | SRM FOS | 16 | `10032, pad, fos (f64)` |
| 0 | padding | 4 | skipped |
| 10017 | EOF | 24 B | `10017, 0, -1, -1, -1, 0` |

### RESULT header

| int | meaning |
|---|---|
| 0 | 14 (record type) |
| 1 | result id. `id % 10^8 // 10^4` is the element family code (31 solid, 43 plane strain, 48 axisymmetric, 41 shell, 42 membrane, 14 beam, 13 truss, 12 spring, 21 interface, 23/24 pile). |
| 2 | location: 0 node, 1 element centroid, 2 element node (label starts with `NDk- `) |
| 3 | category: 0 other, 1 displacement, 3 force, 4 stress, 5 strain |
| 4–17 | not decoded |

`n` is the length of the id list the block refers to. It isn't repeated in the header.

## Details that were easy to get wrong

- **A stage record is followed directly by an increment record.** The two must not be read as a single record.
- **`next_stage_offset = -1` doesn't mean end of file.** The last stage still has its increments. The end is the separate 24-byte EOF record, which an interrupted run doesn't write.
- **Each element family writes its own id list**, often in several chunks: one per element shape, and again per result category. A reader that keeps only the last id list attaches results to the wrong elements.
- **The SRM factor of safety is often only in the increment label** (`FOS=…`). The load factor field holds the trial factor, and the 10032 record isn't always written.
- **Linear static** results have a stage record but no increment record.
- **The `.log` file** sometimes lists node and element numbers as internal sequence numbers rather than ids, and its "MAXIMUM STRESS" isn't any of the centroid components. Don't use it to check ids.

## Record types not covered

Civil NX result files, eigenvalue, response-spectrum, moving-load, heat-of-hydration and load-combination results use other record types. `NfxpFile` raises `NfxpFormatError` for them when no usable steps are found.
