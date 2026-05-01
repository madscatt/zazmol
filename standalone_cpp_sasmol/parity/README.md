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
| PDB I/O | `PdbReader` / `PdbWriter` | required | Contract API added; parser deferred until tolerant fixture parity is reviewed. |
| DCD I/O | `DcdReader` / `DcdWriter` | required | Reader implemented for normal full-coordinate small fixtures; writer round-trips through C++ and Python readers for `1ATM`, `2AAD`, and `rna-1to10`. |
| calculations | free functions / module namespaces | deferred | Port only after file and descriptor parity are stable. |
| GPU/MPI backends | optional backends | deferred | Data views are being shaped now; dependencies wait. |

Run `../tools/generate_parity_ledger.py` from the repository root to produce a
fresh method survey from Python source.
