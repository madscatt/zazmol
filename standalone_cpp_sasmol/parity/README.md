# Parity Ledger

Python `zazmol` is the read-only behavior standard for this standalone C++
library. The parity ledger tracks which Python-facing concepts have been
reviewed, implemented, deferred, or intentionally expressed differently in C++.

## Status Labels

- `required`: needed for v1 parity
- `implemented`: present in the standalone C++ library
- `deferred`: real behavior, but outside the current milestone
- `different`: intentionally different in C++, with the reason documented
- `question`: needs design or behavior review before porting

## Current Decisions

| Python area | C++ target | Status | Notes |
| --- | --- | --- | --- |
| `system` molecule state | `sasmol::Molecule` | implemented | Initial descriptor, coordinate, and integrity model only. |
| `config.COORD_DTYPE` | `sasmol::coord_type` | implemented | `float`, matching storage-friendly coordinate intent. |
| `config.CALC_DTYPE` | `sasmol::calc_type` | implemented | `double`, matching calculation-friendly intent. |
| PDB I/O | `PdbReader` / `PdbWriter` | implemented | Implemented for tolerant fixed-column reads, END and MODEL trajectories, selected/all-frame writes, CONECT output, optional-field fill, element resolution, and larger fixture validation through Python `zazmol`. |
| DCD I/O | `DcdReader` / `DcdWriter` | implemented | Sequential normal full-coordinate reader/writer implemented with reopen-and-scan single-frame reads, whole-trajectory convenience reads, malformed-input handling, lifecycle checks, C++ round trips, and Python validation tooling. |
| calculations | free functions / module namespaces | implemented | Coordinate bounds, DCD streaming min/max, mass, center-of-mass, radius-of-gyration, RMSD, formula, residue charge, and PMI slices implemented. |
| operate | free functions / `AlignmentPlan` | implemented | Translate, center, rotations, PMI alignment, and explicit-index basis alignment implemented with in-place and pure variants. Python keyword-mode `align` wrapper remains deferred. |
| selection/subset | safe parser plus subset helpers | implemented | Bounded expression parser, `all`/`heavy` bases, mask bridge helpers, coordinate subset operations, descriptor get/set, copy, duplicate, and merge slices implemented. BIOMT remains deferred. |
| topology | explicit CHARMM helpers | deferred | Explicit CHARMM type assignment helpers implemented; full CHARMM topology parser/patch/reorder behavior remains deferred. |
| GPU/MPI backends | optional backends | deferred | Data views are being shaped now; dependencies wait. |

Run `../tools/generate_parity_ledger.py` from the repository root to produce a
fresh method survey from Python source.
