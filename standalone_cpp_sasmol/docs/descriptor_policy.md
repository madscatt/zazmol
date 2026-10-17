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

## PDB I/O Scope

PDB read/write populates and emits the descriptors actually represented by PDB
ATOM/HETATM records. On read, `residue_flag()` is initialized to zero for each
atom, matching Python's `False` default. `charmm_type()` is intentionally left
empty by PDB reading because CHARMM atom type data is not a standard PDB column;
topology/CHARMM workflows should populate it explicitly when that information is
available.

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

## Unavailable PDB atom serials

An exact `*****` ATOM/HETATM serial uses the one-based atom position in
`original_index`; `index` remains sequential. This substitutes an unavailable
source serial rather than recovering it. Preserved source serials must not
collide with substituted positions. CONECT endpoints must not depend on those
substitutions. The complete reader checks the first-frame map and returns
`IoCode::format_error` on ambiguity, leaving the destination molecule unchanged.
The standalone record helper requires a positive `atom_ordinal_context` to
accept stars; it cannot validate whole-file collisions by itself. Other invalid
fields and VMD-hex rules are unchanged. No new descriptor or writer behavior
is introduced. Python uses its existing `ValueError` mechanism for these cases.
