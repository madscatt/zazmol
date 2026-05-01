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

## Proposed First C++ Slice

Do not implement the whole PDB parser at once. Start with a strict contract
surface and helper tests:

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
