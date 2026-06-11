# PDBx/mmCIF Support Plan

This plan adds PDBx/mmCIF input support without altering the polished legacy
PDB reader/writer behavior in `pdb_io.py`.

## Scope

Base mmCIF support means reading coordinate and atom-descriptor data into the
same operational SASMOL molecule fields populated by `read_pdb()`. It does not
mean implementing a full mmCIF-aware `pdbscan`/`pdbrx` workflow in this phase.

Phase 1 is Python SASMOL support under `src/python`. Phase 2 is standalone C++
SASMOL support under `standalone_cpp_sasmol`, implemented beside the legacy PDB
reader so the existing `PdbReader`/`PdbWriter` behavior remains untouched.

## Dependencies And Fixtures

- Use the official RCSB/wwPDB `mmcif` Python package.
- A local `mmcif 1.1.1` source snapshot is stored under `third_party/rcsb/`.
- Representative mmCIF fixtures are stored under
  `src/python/test_sasmol/data/mmcif_common/`.
- The first fixture set covers small protein, DNA, RNA, heterogen-containing
  protein, multi-model/NMR, and larger multi-chain cases.

## Phase 1 Python Architecture

1. Add `src/python/mmcif_io.py` with an `MMCIF` mixin. Completed.
2. Keep `src/python/pdb_io.py` behavior unchanged. Completed.
3. Update `src/python/file_io.py` so `Files` inherits `PDB`, `MMCIF`, and
   `DCD`.
4. Add a top-level `read_structure(filename, format='auto', **kwargs)` helper
   that dispatches to `read_pdb()` or `read_mmcif()`.
5. Keep `write_pdb()` as the legacy PDB writer. Add mmCIF writing only after
   input behavior is stable and explicitly requested.

## Phase 2 C++ Architecture

The standalone C++ library follows the same base-I/O boundary as Python:

1. Add `MmcifReader` beside `PdbReader`.
2. Keep format dispatch above the PDB/mmCIF readers if/when a C++ dispatch API
   is added; do not route mmCIF through the legacy PDB parser.
3. Reuse the same `mmcif_common` fixtures for C++ parity tests.
4. Use a focused C++ `_atom_site` reader for base molecule I/O. The vendored
   official RCSB/wwPDB Python parser remains the Python implementation and a
   useful validation oracle, but the C++ base reader does not vendor or wrap a
   separate large C++ parser stack in this phase.
5. Preserve current C++ PDB parity tests unchanged.

## Field Mapping Policy

Base `read_mmcif()` should populate the same molecule fields as `read_pdb()`:

| SASMOL field | mmCIF source policy |
| --- | --- |
| `_atom` | `_atom_site.group_PDB` |
| `_name` | prefer `_atom_site.auth_atom_id`, fallback `_atom_site.label_atom_id` |
| `_loc` | `_atom_site.label_alt_id`, with `.` and `?` treated as blank |
| `_resname` | prefer `_atom_site.auth_comp_id`, fallback `_atom_site.label_comp_id` |
| `_chain` | prefer `_atom_site.auth_asym_id`, fallback `_atom_site.label_asym_id`; preserve full strings internally |
| `_segname` | derive from the chosen chain/asym id when no explicit segment id exists |
| `_resid`, `_original_resid` | prefer `_atom_site.auth_seq_id` only when compatible with current integer descriptor semantics |
| `_rescode` | `_atom_site.pdbx_PDB_ins_code`, with `.` and `?` treated as blank |
| `_coor` | `_atom_site.Cartn_x`, `_atom_site.Cartn_y`, `_atom_site.Cartn_z` |
| `_occupancy` | `_atom_site.occupancy`, preserving current PDB defaults when absent |
| `_beta` | `_atom_site.B_iso_or_equiv`, preserving current PDB defaults when absent |
| `_element` | `_atom_site.type_symbol`, then existing element-resolution logic if blank |
| `_charge` | `_atom_site.pdbx_formal_charge`, with `.` and `?` treated as blank |
| `_index` | SASMOL internal one-based read order |
| `_original_index` | `_atom_site.id` only when compatible with current integer descriptor semantics |
| `_moltype` | existing `read_pdb()` residue-name classification |

Group coordinate frames by `_atom_site.pdbx_PDB_model_num`. Require consistent
atom ordering and atom counts across models, matching current PDB trajectory
assumptions.

## Known Policy Risks

- mmCIF chain/asym ids can be longer than one character. Preserve them
  internally, but legacy `write_pdb()` cannot faithfully represent them in the
  single-character PDB chain field.
- mmCIF atom and residue identifiers can be non-numeric. Base input should fail
  clearly or preserve a side mapping when values cannot fit current integer
  `original_index`/`original_resid` descriptors.
- Legacy `segname` has no exact `_atom_site` correspondence in the official
  PDB-to-PDBx table, so base mmCIF input should follow SASMOL's current
  chain-as-segment fallback.
- Legacy `CONECT` has no direct PDB-to-PDBx correspondence. Semantic
  connections in `_struct_conn` belong to a later mmCIF-aware `pdbscan`/`pdbrx`
  plan, not base molecule loading.
- Raw PDB header lines should not be invented from mmCIF structured metadata in
  the base reader. Rich metadata parsing belongs to `pdbscan`/`pdbrx`.

## Acceptance Gate

The Python Phase 1 implementation is acceptable when:

1. `read_mmcif()` can parse the `mmcif_common` fixtures with the vendored RCSB
   parser.
2. For ordinary numeric-id entries, populated SASMOL descriptors match the
   fields expected from `read_pdb()` behavior.
3. Multi-model coordinates are grouped correctly for `2K39.cif`.
4. Long or non-numeric identifier cases fail with explicit diagnostics or are
   represented by a documented side-mapping policy.
5. Existing legacy PDB tests continue to pass unchanged.

## Phase 1 Status

Implemented in Python on 2026-06-11:

- `src/python/mmcif_io.py` adds base `read_mmcif()`.
- `src/python/file_io.py` adds `read_structure()` format dispatch.
- `src/python/system.py` routes constructor file input through
  `read_structure()`.
- `setup.py` packages the vendored official `mmcif` Python package.
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_read_mmcif.py`
  covers PDB parity for `1CRN`, dispatch, constructor loading, DNA/RNA/
  heterogen fixtures, and multi-model `2K39`.

## Phase 2 Status

Implemented in standalone C++ on 2026-06-11:

- `standalone_cpp_sasmol/include/sasmol/file_io.hpp` exposes `MmcifReader` and
  `MmcifReadOptions` beside the existing PDB API.
- `standalone_cpp_sasmol/src/file_io.cpp` adds implementation-local mmCIF
  tokenization, `_atom_site` loop extraction, model grouping, descriptor
  mapping, and failure-atomic molecule population.
- The C++ reader follows the Phase 1 base-field policy: `auth_*` identifiers
  are preferred with `label_*` fallbacks, full chain/asym ids are preserved
  internally, `segname` is derived from the chosen chain/asym id, `_struct_conn`
  is not interpreted, and `original_index`/`original_resid` require values that
  fit existing integer descriptor semantics.
- `standalone_cpp_sasmol/test_sasmol_cpp/test_mmcif_reader.cpp` covers `1CRN`
  parity with the legacy PDB fixture, DNA/RNA/heterogen fixtures, multi-model
  `2K39`, larger multi-chain `1KP8`, and failure-atomic handling for invalid
  atom identifiers.
- Validation gate passed on 2026-06-11:
  `ctest --test-dir standalone_cpp_sasmol/build --output-on-failure` reported
  17/17 standalone C++ tests passing.
