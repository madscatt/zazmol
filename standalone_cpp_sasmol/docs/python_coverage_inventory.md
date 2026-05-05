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
| DCD fixed/free atoms, unit-cell writing, true random seek | DCD fixture notes, C++ guard tests | `sasmol/file_io.hpp` | `intentional difference` | high | Fixed/free variants and unit-cell writes fail with explicit `unsupported` status; true random seek is not exposed. Add fixtures only if policy changes. |
| BIOMT metadata preservation | `test_file_io/*biomt_metadata*` | `sasmol/file_io.hpp`, `sasmol/molecule.hpp` | `implemented` | low | Keep `read_pdb` passive. |
| `calculate.Calculate` mass/formula/charge/COM/Rg/RMSD/minmax/PMI | `test_calculate/*` | `sasmol/calculate.hpp` | `implemented` | low | Add only fixture-backed edge hardening. |
| `linear_algebra` helper functions | `test_linear_algebra/*` | `sasmol/linear_algebra.hpp` | `implemented` | medium | Includes legacy `find_u` and `cmp` helpers for alignment parity. |
| `operate.Move` translate/center/rotate/PMI alignment | `test_operate/*` | `sasmol/operate.hpp` | `implemented` | low | Existing coverage is strong. |
| `operate.Move.align(...)` legacy mode surface | `test_operate/*align*` | `sasmol/operate.hpp` | `implemented` | high | C++ exposes Python initialization/production semantics as explicit plan helpers rather than string mode dispatch. |
| `operate.set_average_vdw` | `test_operate/test_unit_operate_set_average_vdw.py` | `sasmol/operate.hpp` | `implemented` | low-medium | C++ stores the usable legacy radius column as scalar `atom_vdw()`. |
| `subset.Mask` named basis, masks, indices, copy, duplicate, merge, descriptors | `test_subset/*` | `sasmol/selection.hpp`, `sasmol/subset.hpp` | `implemented` | medium | Keep expanding only through explicit grammar decisions. |
| `subset.Mask.get_subset_mask` open Python eval behavior | `test_subset/*get_subset_mask*` | `sasmol/selection.hpp` | `intentional difference` | high | Survey real expressions before any grammar expansion. |
| `subset.Mask.apply_biomt`, `copy_apply_biomt` | C++ selected-transform tests | `sasmol/subset.hpp` | `implemented` | high | C++ keeps these as one-transform selected-coordinate helpers, separate from metadata and assembly helpers. |
| C++ BIOMT assembly helpers | C++ BIOMT tests | `sasmol/subset.hpp` | `implemented` | medium | Document as assembly helpers, not Python selected-transform parity. |
| `charmm_topology.CharmmTopology` parser records | topology oracle fixtures | `sasmol/topology.hpp` | `implemented` | high | Validate new records against Python fixtures first. |
| CHARMM residue atom lists, atom-only patches, order choice, reorder planning/copy/in-place | topology C++ tests | `sasmol/topology.hpp` | `implemented` | high | Keep mutation explicit and no-partial-mutation. |
| CHARMM36 topology upgrade | planned GitHub issue | external tracking issue | `deferred` | high | Out of scope for this parity pass; handle as joint Python/C++ tracked work. |
| `topology.Topology.create_fasta` | utilities/system FASTA tests | `sasmol/topology.hpp`, `sasmol/molecule.hpp` | `implemented` | low | Keep C++ result typed: sequence vector plus formatted FASTA string. |
| `topology.Topology.renumber` | `test_topology/test_intg_topology_Topology.py`, `test_system/*index*`, `*resid*` | `sasmol/topology.hpp` | `implemented` | low-medium | Preserve explicit options for index/resid starts. |
| `topology.Topology.make_constraint_pdb` | `test_topology/test_intg_topology_Topology.py` | `sasmol/topology.hpp`, `sasmol/file_io.hpp` | `implemented` | medium | C++ splits descriptor application from the PDB-writing wrapper. |
| `topology.Topology.make_backbone_pdb_from_fasta` | `test_topology/test_intg_topology_Topology.py` | `sasmol/topology.hpp`, `sasmol/file_io.hpp` | `implemented` | medium-high | C++ exposes both molecule builder and PDB-writing wrapper. |
| `multiprocessing_sasmol.Multiprocessing_SasMol` | no core C++ parity tests | none | `intentional difference` | low | Experimental Python orchestration helper; no C++ port planned. |
| Python property tables: AMU, VDW, scattering lengths, amino-acid SLD | `test_properties/*` | `sasmol/properties.hpp` | `implemented` | medium | C++ table accessors are fixture-checked against Python oracle data. |
| `charmm_topology.CharmmTopology.charmm_names` | `test_properties/test_unit_properties_Atomic_charmm_names.py` | `sasmol/topology.hpp` | `implemented` | medium | Pure CHARMM atom-name classification table fixture-checked against Python data. |
| VMD/view helpers | `test_system/*send_coordinates_to_vmd*` | optional adapter | `deferred` | medium | Keep outside portable core. |

## Recommended Port Order

The remaining in-scope work needs explicit guardrails. Recommended order:

1. Bounded selection grammar expansion

   Survey real selection expressions first, then add only named grammar features
   with Python/C++ parity tests. Do not recreate Python `eval`.

2. VMD/view optional adapter boundary

   Keep portable core independent of VMD. If added, expose a small adapter API
   around coordinate buffers and clear validation errors.

## Policy-Gated Work

The following should remain explicit design decisions, not incidental parity
cleanup:

- broad Python eval-style selection compatibility
- contextual selection aliases such as generic `backbone`
- DCD fixed/free atom variants
- DCD unit-cell writing
- true random-access DCD seeking
- VMD/view adapter behavior

## Slice Rhythm

For each future port slice:

- inspect the Python oracle behavior first
- add or identify a tiny Python fixture when output is structural or formatted
- implement only the named surface
- preserve clear status/error returns and no partial mutation
- run normal and sanitizer C++ tests
- update this inventory and the relevant parity note
- commit and push before starting the next slice
