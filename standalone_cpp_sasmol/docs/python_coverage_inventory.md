# Python Coverage Inventory

Python `src/python` remains the read-only behavior oracle for the standalone C++
SasMol port. This inventory records the current C++ coverage status by Python
module surface so future work can move in small, reviewed slices.

Status labels:

- `implemented`: represented in the standalone C++ API and tests
- `partial`: core behavior exists, but the Python surface or legacy nuance is
  not fully represented
- `deferred`: intentionally postponed pending fixtures, caller need, or policy
- `intentional difference`: C++ deliberately keeps a safer or narrower contract
- `missing`: no direct C++ port yet

## Module Coverage

| Python surface | Python oracle tests | Current C++ home | Status | Risk | Recommended action |
| --- | --- | --- | --- | --- | --- |
| `file_io.Files` PDB read/write | `test_file_io/*read_pdb*`, `*write_pdb*` | `sasmol/file_io.hpp` | `implemented` | low | Continue hardening only with fixtures. |
| `file_io.Files` DCD open/read/step/close/write | `test_file_io/*dcd*` | `sasmol/file_io.hpp` | `implemented` | medium | Keep streaming APIs explicit for large trajectories. |
| DCD fixed/free atoms, unit-cell writing, true random seek | DCD fixture notes | `sasmol/file_io.hpp` | `deferred` | high | Require policy and fixtures before coding. |
| BIOMT metadata preservation | `test_file_io/*biomt_metadata*` | `sasmol/file_io.hpp`, `sasmol/molecule.hpp` | `implemented` | low | Keep `read_pdb` passive. |
| `calculate.Calculate` mass/formula/charge/COM/Rg/RMSD/minmax/PMI | `test_calculate/*` | `sasmol/calculate.hpp` | `implemented` | low | Add only fixture-backed edge hardening. |
| `operate.Move` translate/center/rotate/PMI alignment | `test_operate/*` | `sasmol/operate.hpp` | `implemented` | low | Existing coverage is strong. |
| `operate.Move.align(...)` legacy mode surface | `test_operate/*align*` | `sasmol/operate.hpp` | `partial` | high | Port only after mode-specific oracle review. |
| `operate.set_average_vdw` | `test_operate/test_unit_operate_set_average_vdw.py` | descriptor/property helpers | `missing` | low-medium | Port as isolated helper after FASTA/renumber slices. |
| `subset.Mask` named basis, masks, indices, copy, duplicate, merge, descriptors | `test_subset/*` | `sasmol/selection.hpp`, `sasmol/subset.hpp` | `implemented` | medium | Keep expanding only through explicit grammar decisions. |
| `subset.Mask.get_subset_mask` open Python eval behavior | `test_subset/*get_subset_mask*` | `sasmol/selection.hpp` | `intentional difference` | high | Survey real expressions before any grammar expansion. |
| `subset.Mask.apply_biomt`, `copy_apply_biomt` | future dedicated fixtures | `sasmol/subset.hpp` | `deferred` | high | Keep separate from metadata and assembly helpers. |
| C++ BIOMT assembly helpers | C++ BIOMT tests | `sasmol/subset.hpp` | `implemented` | medium | Document as assembly helpers, not Python selected-transform parity. |
| `charmm_topology.CharmmTopology` parser records | topology oracle fixtures | `sasmol/topology.hpp` | `implemented` | high | Validate new records against Python fixtures first. |
| CHARMM residue atom lists, atom-only patches, order choice, reorder planning/copy/in-place | topology C++ tests | `sasmol/topology.hpp` | `implemented` | high | Keep mutation explicit and no-partial-mutation. |
| CHARMM36 topology upgrade | planned GitHub issue | `sasmol/topology.hpp` | `deferred` | high | Handle as joint Python/C++ tracked work. |
| `topology.Topology.create_fasta` | utilities/system FASTA tests | `sasmol/molecule.hpp` plus future helper | `partial` | low | Next implementation slice. |
| `topology.Topology.renumber` | `test_system/*index*`, `*resid*` | future topology/system helper | `missing` | low-medium | Port after `create_fasta`. |
| `topology.Topology.make_constraint_pdb` | dedicated oracle fixture needed | future topology/file-IO helper | `missing` | medium | Require tiny Python output fixture before coding. |
| `topology.Topology.make_backbone_pdb_from_fasta` | dedicated oracle fixture needed | future topology builder | `missing` | medium-high | Pause before coding; behavior is structural synthesis. |
| `multiprocessing_sasmol.Multiprocessing_SasMol` | no core C++ parity tests | optional adapter layer | `missing` | high | Defer unless standalone C++ needs orchestration parity. |
| Python property tables: AMU, VDW, scattering lengths, CHARMM names, amino-acid SLD | `test_properties/*` | `sasmol/properties.hpp` | `partial` | medium | Expand table-by-table only when a caller or helper needs it. |
| VMD/view helpers | `test_system/*send_coordinates_to_vmd*` | optional adapter | `deferred` | medium | Keep outside portable core. |

## Recommended Port Order

1. `topology.Topology.create_fasta`

   Low risk. The C++ molecule already owns FASTA storage, and the behavior can
   be tested with small residue-name fixtures.

2. `topology.Topology.renumber`

   Low to medium risk. This is deterministic descriptor mutation, but it must
   preserve failure-before-mutation behavior where C++ can validate inputs.

3. `operate.set_average_vdw`

   Low to medium risk. This is isolated, but may require a focused property-table
   expansion for VDW values.

4. `topology.Topology.make_constraint_pdb`

   Medium risk. It combines selection, descriptor changes, and PDB writing, so a
   tiny Python oracle output should come first.

5. `topology.Topology.make_backbone_pdb_from_fasta`

   Medium to high risk. This builds molecular structure from FASTA and should
   not be ported without dedicated fixtures for supported molecule types.

6. Legacy `operate.Move.align(..., mode=...)` details

   High risk. Review mode-specific Python behavior before adding compatibility
   APIs.

7. Python-style selected BIOMT `apply_biomt` / `copy_apply_biomt`

   High risk. Keep separate from existing BIOMT metadata preservation and C++
   assembly helpers.

## Policy-Gated Work

The following should remain explicit design decisions, not incidental parity
cleanup:

- broad Python eval-style selection compatibility
- contextual selection aliases such as generic `backbone`
- DCD fixed/free atom variants
- DCD unit-cell writing
- true random-access DCD seeking
- CHARMM36 topology upgrade
- multiprocessing and VMD/view adapter behavior

## Slice Rhythm

For each future port slice:

- inspect the Python oracle behavior first
- add or identify a tiny Python fixture when output is structural or formatted
- implement only the named surface
- preserve clear status/error returns and no partial mutation
- run normal and sanitizer C++ tests
- update this inventory and the relevant parity note
- commit and push before starting the next slice
