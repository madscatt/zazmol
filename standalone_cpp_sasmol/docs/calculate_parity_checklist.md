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

## Center Of Mass Slice

Implemented `calculate_center_of_mass(molecule, frame)`:

- auto-calculates mass when total mass is not populated
- returns calculation-precision coordinates
- rejects unknown masses and out-of-range frames
- fixture parity for `2AAD.pdb`, `rna.pdb`, and `1CRN.pdb`

## Radius Of Gyration Slice

Implemented `calculate_radius_of_gyration(molecule, frame)`:

- follows Python SASMOL convention: mass-weighted center of mass, then
  unweighted mean squared coordinate distance over atom count
- fixture parity for `1ATM.pdb`, `2AAD.pdb`, `rna.pdb`, and `1CRN.pdb`
- reuses center-of-mass validation for unknown masses and frame errors

## RMSD Slice

Implemented `calculate_root_mean_square_deviation(first, second)`:

- follows Python SASMOL convention: sum squared differences over all loaded
  coordinate values divided by atom count
- synthetic parity for one-atom and two-atom examples
- fixture parity for identical `1ATM.pdb`, `1CRN-rot.pdb`, and
  `1CRN-rot-shift.pdb`
- shape mismatches are explicit errors instead of printed warnings

## Formula And Residue Charge Slice

Implemented descriptor/topology calculations:

- `calculate_molecular_formula(molecule)`
- `calculate_residue_charge(molecule)`
- formula is stored as `std::map<std::string, std::size_t>`, matching Python's
  dictionary-style result
- residue charge is stored per atom as `calc_type`
- fixture parity for `1ATM.pdb` and `2AAD.pdb` formula counts
- Python unit-test parity for two-residue charge sums
- descriptor mismatches are explicit errors

## Principal Moments Of Inertia Slice

Implemented `calculate_principal_moments_of_inertia(molecule, frame)`:

- follows Python SASMOL's center-of-mass inertia tensor convention
- returns inertia tensor, sorted eigenvalues, eigenvectors, and an explicit
  `singular` flag for rank-deficient tensors
- keeps eigen decomposition local to `calculate.cpp`; no Eigen dependency is
  introduced in this slice
- fixture parity for `2AAD.pdb`
- synthetic coverage for a simple non-linear three-atom system
- one-atom tensors are reported as singular rather than assigned unstable
  principal axes

## DCD Min/Max Streaming Slice

Implemented DCD trajectory overloads:

- `calculate_minimum_and_maximum_all_steps(trajectory_filename)`
- `calc_minmax_all_steps(trajectory_filename)`
- streams frame coordinates through `DcdReader` without allocating a full
  trajectory molecule
- keeps the existing loaded-molecule overloads unchanged
- fixture parity for `1ATM.dcd`
- loaded-versus-streamed parity for `2AAD.dcd`

## Deferred

- alignment-dependent calculations
