# C++ `merge_two_molecules` Design

Python `zazmol.subset.Mask.merge_two_molecules` is the behavioral reference,
but the standalone C++ implementation should use the C++ data model directly
instead of recreating Python's open-ended attribute mechanics.

## Current Python Behavior

- `mol1` must contain at least one atom.
- `mol2` may contain zero atoms.
- Coordinates are copied only from frame 0 of each input molecule.
- The merged molecule has one frame.
- `mol1.index()` is preserved.
- `mol2.index()` is regenerated sequentially from the last `mol1.index()`.
- Core atom descriptors are concatenated in atom order.
- Optional descriptors are merged only when both molecules provide them.
- Missing optional descriptors may be reported when requested.
- `CONECT` data is copied without aliasing the first molecule.

## C++ Design Choice

Python can decide whether a descriptor exists with `hasattr`. C++ `Molecule`
does not have open-ended attributes; it owns typed vectors for descriptors,
coordinates, formula data, unit cell data, and connectivity. Therefore the C++
merge contract should validate vector lengths and preserve typed state instead
of emulating Python missing-attribute behavior.

This keeps the API readable, memory-safe, and backend-neutral while preserving
the behavior that matters to users: atom order, coordinates, descriptor values,
index regeneration for the second molecule, and non-aliased connectivity.

## Implemented API

```cpp
struct MergeOptions {
  bool report_skipped_descriptors{false};
};

[[nodiscard]] SubsetResult merge_two_molecules(
    const Molecule& mol1, const Molecule& mol2, Molecule& merged,
    MergeOptions options = {});

[[nodiscard]] Molecule merged_two_molecules(
    const Molecule& mol1, const Molecule& mol2,
    MergeOptions options = {});
```

`merge_two_molecules` is the non-throwing worker. `merged_two_molecules` is the
convenience value-returning wrapper and should throw `std::invalid_argument`
when the worker reports errors.

## Validation Rules

- Return an error if `mol1.natoms() == 0`.
- Return an error if either molecule has zero frames but nonzero atoms.
- Return an error if coordinate storage does not match
  `natoms * number_of_frames * 3`.
- Return an error if `mol1.index()` does not cover all `mol1` atoms.
- Validate core descriptors before mutating `merged`.
- Do not partially mutate `merged` on failure.
- Treat empty optional vectors as absent only when both inputs are empty.
- Merge optional numeric vectors only when both vectors have atom-length data.
- If exactly one optional vector is present, skip it and optionally report that.

## Merge Rules

- Resize a temporary output molecule to `mol1.natoms() + mol2.natoms()` and one
  frame.
- Copy frame 0 coordinates from `mol1` followed by frame 0 coordinates from
  `mol2`.
- Concatenate string, integer, and calculation descriptors that are part of the
  C++ core `Molecule` model.
- Preserve `mol1.index()` exactly.
- Assign `mol2` indices as `last(mol1.index()) + 1 ... + mol2.natoms()`.
- Concatenate `original_index()` and `original_resid()` unchanged.
- Copy `conect()` by value. The stored linkage values remain original PDB atom
  indices, consistent with current C++ PDB writer remapping behavior.
- Copy molecule-level scalar/metadata fields conservatively:
  - formula is cleared unless recalculated later
  - total mass is the sum of input total masses only when both are nonzero
  - unit cell is copied from `mol1`
  - FASTA is concatenated only when either input has FASTA text

## Implemented Slice

The core merge is implemented without BIOMT or force-field-specific inference:

- coordinates
- core descriptors
- atom charges, VDW, mass
- residue charge
- connectivity
- index regeneration
- one-frame merged output
- `charmm_type()`
- typed extension descriptor maps

Focused tests cover:

- 1-atom + 2AAD fixture merge
- mol2 index regeneration
- frame 0 only when mol2 has multiple frames
- mol1 empty error
- coordinate shape/integrity error without output mutation
- optional numeric descriptor skip/report behavior
- non-aliased `CONECT`
- typed extension descriptor merging and skipped-descriptor reporting
- value-returning wrapper behavior

## Deferred

- caller-defined required descriptor sets
- BIOMT-aware merge behavior
- formula/residue-charge recalculation policies beyond direct concatenation

`charmm_type()` is no longer deferred; it is a typed C++ descriptor. New
one-off or experimental atom-aligned data should use typed extension descriptor
maps first, then graduate into explicit molecule fields only after review.
