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

## Mass Slice

Implemented mass/property groundwork:

- `amu()` standard atomic weight table using Python SASMOL values
- `calculate_mass(molecule)`
- updates per-atom `mass` and `total_mass`
- reports unknown elements explicitly while preserving zero mass for those
  entries
- rejects element descriptor length mismatches

The C++ tests cover:

- `2AAD.pdb` per-atom masses and total mass
- `rna.pdb` and `1CRN.pdb` total masses
- unknown element reporting
- descriptor mismatch handling

## Deferred

- DCD streaming overloads for `calculate_minimum_and_maximum_all_steps`
- center of mass, radius of gyration, and principal moments of inertia
- RMSD/alignment-dependent calculations
- molecular formula and residue charge
