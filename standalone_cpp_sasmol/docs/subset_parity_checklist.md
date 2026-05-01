# Subset Parity Checklist

Python `zazmol.subset` combines mask construction, coordinate extraction,
coordinate replacement, molecule copying, duplication, BIOMT, and other
higher-level operations. The standalone C++ port should build this area from
explicit index operations first, then add mask-shaped compatibility wrappers
where useful.

## First Implemented Slice

Implemented explicit-index operations:

- `get_coordinates_using_indices(molecule, frame, indices)`
- `set_coordinates_using_indices(molecule, source, frame, indices)`
- `copy_molecule_using_indices(source, destination, indices, frame)`
- `copied_molecule_using_indices(source, indices, frame)`

Behavior notes:

- failures return structured errors for non-throwing APIs
- empty selections, bad frames, and bad atom indices fail without mutation
- coordinate replacement requires the source atom count to match the selected
  destination atom count
- copy creates a one-frame molecule using the requested source frame
- copied descriptors follow the explicit selected atom order
- `CONECT` entries are filtered to selected atoms only

## Deferred

- 0/1 mask wrappers matching Python names
- descriptor set/get using mask or indices
- `duplicate_molecule`
- `merge_two_molecules`
- BIOMT operations
