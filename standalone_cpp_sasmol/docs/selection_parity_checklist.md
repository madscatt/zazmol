# Selection Parity Checklist

Python `zazmol.subset.get_subset_mask` uses `eval` over per-atom descriptor
arrays. The standalone C++ library must not reproduce arbitrary evaluation as
the default selection mechanism. Selection support should grow through an
explicit, safe grammar and typed helper APIs.

## First Implemented Slice

Implemented explicit index helpers:

- `indices_all(molecule)`
- `indices_by_name(molecule, name)`
- `indices_by_resname(molecule, resname)`
- `indices_by_resid_range(molecule, first_resid, last_resid)`

Implemented a bounded expression parser:

- `select_indices(molecule, expression)`

Supported grammar:

- `all`
- descriptor comparisons of the form `descriptor[i] OP literal`
- `and`, `or`, and parentheses
- string descriptors: `name`, `resname`, `chain`, `segname`, `element`,
  `moltype`
- integer descriptors: `resid`, `index`, `original_index`, `original_resid`
- string operators: `==`, `!=`
- integer operators: `==`, `!=`, `<`, `<=`, `>`, `>=`

Failure policy:

- unsupported descriptors fail with an error
- unsupported tokens or grammar fail with an error
- descriptor length mismatches fail with an error
- zero selected atoms fail with an error
- failures return no indices

This is intentionally narrower than Python `eval`; it is meant to fail gently
and stop rather than guessing.

## Deferred

- mask objects and 0/1 mask arrays
- copy/get/set coordinate using mask or indices
- distance, angle, and hydrogen-bond selections
- named selection aliases such as backbone or calpha
- full Python expression compatibility
