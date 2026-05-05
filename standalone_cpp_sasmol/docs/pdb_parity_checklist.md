# PDB Parity Checklist

Python `zazmol` remains the behavior oracle for PDB I/O. This checklist is
based on a read-only survey of:

- `src/python/pdb_io.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_read_pdb.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_write_pdb.py`
- `src/python/test_sasmol/test_file_io/test_unit_file_io_Files_helpers.py`

The fixed-column reference used for column ranges is summarized in
[`pdb_fixed_column_reference.md`](pdb_fixed_column_reference.md).

The C++ parser should not aim to be a narrow protein-only PDB reader. SASMOL
PDB behavior is intentionally tolerant because real SASSIE/ZAZZIE workflows
include proteins, RNA, mixed/odd atom names, nonstandard records, generated
model trajectories, and legacy files.

## Read Behavior To Preserve First

- Missing files return a clean file error.
- `ATOM` and `HETATM` records are coordinate records.
- Non-coordinate records before the first coordinate frame are preserved as
  header content where feasible.
- Single-frame files with no final `END` are accepted.
- Files ending with blank trailing lines are accepted.
- END-separated pseudo-trajectories are accepted when atom counts are
  consistent.
- MODEL/ENDMDL trajectories are accepted when structure is consistent.
- Malformed MODEL/ENDMDL/END combinations fail cleanly:
  - missing terminal `END` for MODEL trajectories
  - inconsistent atom counts per frame
  - mixed MODEL/END structures that imply ambiguous frame boundaries
- Descriptor parsing preserves Python defaults:
  - missing occupancy defaults
  - missing beta defaults
  - missing segname can fall back to chain
  - missing element/charge remain blank-style values
- Element resolution must preserve the Python conflict and CHARMM-name logic.
- `pdbscan` mode preserves alternate-location and optional fields more
  literally.
- `check_zero_coor` must apply the all-zero coordinate column guard.
- `CONECT` parsing is only populated in scan mode today and must be treated
  deliberately.

## Write Behavior To Preserve First

- Writes one selected coordinate frame.
- Supports normal `w`/`a` style behavior.
- Supports `MODEL` output and terminal `END`/`ENDMDL` behavior.
- Supports optional `conect` output through remapped original/current indices.
- Clamps overwide atom serials and residue ids to legacy-width strings.
- Coordinates are formatted as 8.3 columns.
- Optional PDB fields are not filled by default.
- `fill_missing_optional=True` fills optional fields but still raises for
  missing required fields.
- `check_for_all_zero_columns` nudges only the first atom of all-zero axes by
  `1e-10`.

## Minimum Fixture Parity

Read fixtures:

| Fixture | Behavior |
| --- | --- |
| `pdb_common/1ATM.pdb` | one atom, one frame |
| `pdb_common/1ATM-1to2.pdb` | one atom, two frames |
| `pdb_common/2AAD.pdb` | one-frame amino-acid fixture |
| `sasmol/file_io/2AAD-1to3-END.pdb` | END-separated three-frame pseudo-trajectory |
| `sasmol/file_io/2AAD-1to3-MODEL.pdb` | MODEL/ENDMDL three-frame trajectory |
| `sasmol/file_io/1AA-NoEND.pdb` | accepted single-frame no-END file |
| `pdb_common/dimcd_fixed_atoms.pdb` | accepted trailing blank lines |
| `pdb_common/rna-1to10.pdb` | RNA multi-frame fixture |

Failure fixtures:

| Fixture | Expected behavior |
| --- | --- |
| `sasmol/file_io/2AAD-1to3-END_wrong_number_atoms.pdb` | frame atom-count error |
| `sasmol/file_io/2AAD-1to3-MODEL_wrong_number_atoms.pdb` | frame atom-count error |
| `sasmol/file_io/2AAD-1to3-MODEL_wrongnumber_mix_END.pdb` | ambiguous MODEL/END error |
| `sasmol/file_io/2AAD-1to3-MODEL_mix_END_noterminating.pdb` | missing terminating END error |
| `pdb_common/1PSI.pdb` | invalid MODEL/ENDMDL structure |

Write fixtures:

| Fixture | Behavior |
| --- | --- |
| `pdb_common/1ATM.pdb` | one-frame write |
| `pdb_common/1ATM-1to2.pdb` | selected second-frame write |
| `pdb_common/2AAD.pdb` | amino-acid write |
| `pdb_common/rna-1to10.pdb` | selected RNA frame write |
| `pdb_common/1CRN.pdb` | protein write |
| optional-field synthetic tests | default blank behavior and fill behavior |

## Completed C++ Coverage

The first PDB parity surface is implemented and covered by fixture/generated
tests:

1. Add C++ helper functions for fixed-column slicing and numeric parsing that
   return `IoStatus`/small result types instead of throwing or exiting.
2. Added tests for `check_for_all_zero_columns` and `create_conect_pdb_lines`
   parity because those are isolated and already important.
3. Added PDB frame-boundary pre-scan tests for the fixture categories above.
4. Added fixed-column ATOM/HETATM record parsing helper tests using SASMOL
   field names/defaults.
5. Added one-frame `read_pdb` descriptor/coordinate population for `1ATM` and
   `2AAD`.
6. Added multi-frame coordinate loading while preserving first-frame
   descriptors for `1ATM-1to2`, END-separated `2AAD`, and MODEL/ENDMDL `2AAD`.
7. Added PDB read parity for moltype classification, all-zero coordinate guard,
   and `pdbscan` CONECT parsing.
8. Added single-frame `write_pdb` round-trip coverage for `1ATM` and `2AAD`.
9. Added selected-frame write support, MODEL/ENDMDL output, and CONECT output.
10. Added `write_all_frames` multi-model output with C++ readback coverage.
11. Added larger fixture coverage for `rna-1to10.pdb` multi-frame reads,
    selected-frame RNA writes, and `1CRN.pdb` protein read/write round trips.
12. Added direct fixture tests for no-terminal-`END` input and trailing-blank
    tolerance.
13. Added explicit missing-descriptor handling for PDB writes:
    default writes fail cleanly when optional descriptor vectors are missing,
    `fill_missing_optional` supplies Python-compatible optional defaults, and
    missing required descriptors remain errors.
14. Added core PDB element resolution for blank element fields: conflict atom
    names, isotope/salt aliases, known heavy elements, first-letter fallback,
    digit-prefixed fallback, and invalid-name errors.
15. Added table-driven element-resolution tests using the existing Python
    SASMOL property fixtures for H/C/N/O/S/P, other elements, miscellaneous
    aliases, and conflict atoms.

## Proposed Next Validation Slice

Completed Python cross-reader validation for generated C++ PDB outputs:

1. Generated C++ written PDB files for `1ATM`, `2AAD`, and `1ATM-1to2`.
2. Validated them with Python `zazmol` using
   `tools/validate_pdb_with_python.py`.
3. Recorded the result in `pdb_python_validation.md`.
4. Kept this separate from normal CTest so Python is not an unconditional C++
   build dependency.

## Current Hardening Inventory (May 5, 2026)

The PDB surface is in a useful state. It already covers the core fixed-column
reader/writer behavior, multi-frame frame-boundary handling, selected/all-frame
writes, generated malformed input cases, element resolution fixtures,
`pdbscan`-mode `CONECT` parsing, passive BIOMT metadata capture, and Python
cross-reader validation for generated C++ outputs.

The next PDB work should be hardening, not broadening. Do not add a new parser
mode, topology inference, BIOMT transform behavior, or stricter PDB schema.

Recommended next implementation slice:

1. Make `PdbReader::read_pdb` failure-atomic.
   Parse into a temporary `Molecule`, including coordinates, descriptors,
   `CONECT`, all-zero guard, and BIOMT metadata. Assign to the caller's molecule
   only after the whole read succeeds.
2. Add a generated malformed PDB test proving an existing destination molecule is
   unchanged after a mid-parse format error.
3. Add a generated missing-optional-field read test to lock ordinary mode
   defaults versus `pdbscan` literal preservation.
4. Run normal and ASAN C++ tests.

Suggested stop point:

- Stop after read-failure atomicity and the two generated read-side hardening
  tests. Leave new PDB syntax support, CHARMM/topology behavior, and Python
  cross-reader expansion for separate decisions.
