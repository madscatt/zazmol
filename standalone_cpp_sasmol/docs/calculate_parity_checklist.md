# Calculate Parity Checklist

Python `zazmol.calculate` remains the behavior oracle for scientific
calculations. C++ calculation parity should grow in small, fixture-backed
batches so numerical choices stay easy to review.

## First Implemented Slice

Implemented coordinate bounds calculations:

- `calculate_minimum_and_maximum(molecule)`
- `calculate_minimum_and_maximum(molecule, frames)`
- `calculate_minimum_and_maximum_all_steps(molecule)`
- `calc_minmax_all_steps(molecule)`

The first C++ slice covers:

- all loaded frames by default
- selected frame lists
- empty molecule and bad-frame error handling
- Python fixture parity for `2AAD.pdb`, `1CRN.pdb`, and multi-frame
  `1ATM-1to2.pdb`
- alias parity for `calc_minmax_all_steps`

## Deferred

- DCD streaming overloads for `calculate_minimum_and_maximum_all_steps`
- mass-dependent calculations: mass, center of mass, radius of gyration,
  principal moments of inertia
- RMSD/alignment-dependent calculations
- molecular formula and residue charge
