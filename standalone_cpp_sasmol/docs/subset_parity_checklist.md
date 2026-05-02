# Subset Parity Checklist

Python `zazmol.subset` combines mask construction, coordinate extraction,
coordinate replacement, molecule copying, duplication, BIOMT, and other
higher-level operations. The standalone C++ port should build this area from
explicit index operations first, then add mask-shaped compatibility wrappers
where useful.

## Current Parity Checkpoint

Subset has a useful first-pass C++ surface for atom selection workflows:

- index and 0/1 mask forms for coordinate extraction, coordinate replacement,
  molecule copy, and value-returning copy helpers
- descriptor get/set helpers for built-in string, integer, and calculation
  descriptors
- extension descriptor map preservation during copy and merge
- duplicate and merge helpers with structured errors
- no-mutation failure behavior for bad masks, bad indices, bad frames, shape
  mismatches, and invalid merge sources

The remaining subset work is not a blocker for ordinary selection/copy/merge
use. BIOMT remains the main deferred behavior because it is more than a wrapper:
it changes coordinates through supplied transforms and should be ported with
dedicated fixtures.

At this checkpoint, subset is complete for v1 ordinary selection-driven
workflows: index and 0/1 mask conversion, coordinate extraction/replacement,
copy/value-copy, descriptor get/set, duplication, merge, and structured
failure handling.

## First Implemented Slice

Implemented explicit-index operations:

- `get_coordinates_using_indices(molecule, frame, indices)`
- `set_coordinates_using_indices(molecule, source, frame, indices)`
- `copy_molecule_using_indices(source, destination, indices, frame)`
- `copied_molecule_using_indices(source, indices, frame)`

Implemented 0/1 mask compatibility wrappers:

- `get_indices_from_mask(molecule, mask)`
- `get_coordinates_using_mask(molecule, frame, mask)`
- `set_coordinates_using_mask(molecule, source, frame, mask)`
- `with_coordinates_using_indices(target, source, frame, indices)`
- `with_coordinates_using_mask(target, source, frame, mask)`
- `copy_molecule_using_mask(source, destination, mask, frame)`
- `copied_molecule_using_mask(source, mask, frame)`

Selection masks can be produced by `selection` helpers:

- `select_mask(molecule, expression)`
- `select_named_basis_mask(molecule, basis_name)`

Implemented descriptor get/set helpers:

- `get_string_descriptor_using_indices(...)`
- `get_string_descriptor_using_mask(...)`
- `set_string_descriptor_using_indices(...)`
- `set_string_descriptor_using_mask(...)`
- `get_int_descriptor_using_indices(...)`
- `get_int_descriptor_using_mask(...)`
- `set_int_descriptor_using_indices(...)`
- `set_int_descriptor_using_mask(...)`
- `get_calc_descriptor_using_indices(...)`
- `get_calc_descriptor_using_mask(...)`
- `set_calc_descriptor_using_indices(...)`
- `set_calc_descriptor_using_mask(...)`

Implemented molecule duplication:

- `duplicate_molecule(molecule, number_of_duplicates)`

Implemented first molecule merge slice:

- `merge_two_molecules(mol1, mol2, merged, options)`
- `merged_two_molecules(mol1, mol2, options)`

Behavior notes:

- failures return structured errors for non-throwing APIs
- empty selections, bad frames, bad atom indices, bad mask lengths, and
  non-0/1 mask values fail without mutation
- coordinate replacement requires the source atom count to match the selected
  destination atom count
- copy creates a one-frame molecule using the requested source frame
- copied descriptors follow the explicit selected atom order
- `CONECT` entries are filtered to selected atoms only
- descriptor APIs use typed descriptor enums instead of arbitrary descriptor
  object mutation
- built-in and extension descriptor policy is documented in
  `docs/descriptor_policy.md`
- duplication returns independent value copies; zero requested duplicates returns
  an empty vector
- merge copies frame 0 coordinates into a one-frame molecule
- merge preserves `mol1.index()` and regenerates `mol2.index()` sequentially
- merge validates input shape before mutating the output molecule
- merge copies connectivity by value

## Deferred

- BIOMT operations

The current merge behavior is documented in `docs/merge_two_molecules_design.md`.
Future merge expansion should be reviewed separately rather than folded into
ordinary mask/copy cleanup.
