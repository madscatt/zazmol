# Descriptor Policy

Python `zazmol` remains the behavioral inventory for atom-aligned molecular
descriptors. The descriptors already present there are treated as hard-earned
library surface, not incidental implementation detail.

## Current Built-In Atom Descriptors

The standalone C++ `Molecule` should support the stable descriptor families
used by Python `zazmol`:

- PDB record field: `record()`; this corresponds to Python `_atom`
- atom indices: `index()`, `original_index()`
- residue indices: `resid()`, `original_resid()`
- PDB strings: `name()`, `loc()`, `resname()`, `chain()`, `rescode()`,
  `occupancy()`, `beta()`, `segname()`, `element()`, `charge()`
- molecule/type strings: `moltype()`, `charmm_type()`
- numeric atom data: `atom_charge()`, `atom_vdw()`, `mass()`
- residue annotations: `residue_flag()`, `residue_charge()`
- coordinate and connectivity data: `coor()`, `conect()`

`charmm_type()` is optional because Python treats `_charmm_type` as an optional
force-field descriptor. `residue_flag()` is part of the stable atom-aligned core
state and is sized with `natoms`.

## Extension Descriptors

Developer or one-off descriptors can live in typed extension maps:

- `extra_string_descriptors()`
- `extra_int_descriptors()`
- `extra_calc_descriptors()`

Each entry must be atom-aligned when present. Copy, subset, merge, and integrity
checks should preserve and validate these maps. Stable descriptors that become
regularly used should graduate from extension maps into explicit `Molecule`
fields after review.

## Maintenance Rule

When Python `zazmol` gains or reveals a stable atom-aligned descriptor, update
this policy and choose one of:

- add it as an explicit C++ `Molecule` field
- keep it in extension maps temporarily with a stated reason
- document why it is intentionally not part of standalone C++ parity
