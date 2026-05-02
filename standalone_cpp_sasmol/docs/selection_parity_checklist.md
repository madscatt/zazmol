# Selection Parity Checklist

Python `zazmol.subset.get_subset_mask` uses `eval` over per-atom descriptor
arrays. The standalone C++ library must not reproduce arbitrary evaluation as
the default selection mechanism. Selection support should grow through an
explicit, safe grammar and typed helper APIs.

## Current Parity Checkpoint

Selection has a safe first-pass C++ surface:

- explicit index helpers for common descriptor selections
- bounded expression parsing for the Python-style basis expressions we have
  intentionally admitted
- named `all` and `heavy` basis helpers
- 0/1 mask bridge helpers for compatibility with subset APIs
- structured failure returns with no partial indices on failure

The C++ selection layer is intentionally narrower than Python `eval`. That is a
feature, not an accidental gap: unsupported grammar should fail clearly rather
than executing arbitrary Python-like expressions.

## First Implemented Slice

Implemented explicit index helpers:

- `indices_all(molecule)`
- `indices_by_name(molecule, name)`
- `indices_by_resname(molecule, resname)`
- `indices_by_resid_range(molecule, first_resid, last_resid)`

Implemented named basis helpers:

- `basis_expression("all")` -> `not name[i] == None`
- `basis_expression("heavy")` -> `not name[i][0] == "H"`
- `select_named_basis(molecule, basis_name)`

Implemented selection-to-mask bridge helpers:

- `mask_from_indices(molecule, indices)`
- `select_mask(molecule, expression)`
- `select_named_basis_mask(molecule, basis_name)`

These return 0/1 atom masks compatible with the subset mask APIs while keeping
the explicit-index selection API available for callers that do not need Python
mask-shaped compatibility.

`backbone` is intentionally not supported as a generic named basis. Current
Python/SASSIE practice uses different backbone definitions for protein overlap,
protein constraints, and nucleic acid workflows, so a generic C++ alias would
hide real caller intent.

Implemented a bounded expression parser:

- `select_indices(molecule, expression)`

Supported grammar:

- `all`
- descriptor comparisons of the form `descriptor[i] OP literal`
- unary `not` before a comparison or parenthesized expression
- first-character string checks of the form `descriptor[i][0] OP literal`
- `and`, `or`, and parentheses
- string descriptors: `record`, `name`, `loc`, `resname`, `chain`, `rescode`,
  `occupancy`, `beta`, `segname`, `element`, `charge`, `moltype`,
  `charmm_type`
- integer descriptors: `resid`, `index`, `original_index`, `original_resid`
- boolean/integer descriptor: `residue_flag`
- string operators: `==`, `!=`
- integer operators: `==`, `!=`, `<`, `<=`, `>`, `>=`
- `None` comparisons for string descriptors, primarily to support historical
  expressions such as `not name[i] == None`

Failure policy:

- unsupported descriptors fail with an error
- unsupported tokens or grammar fail with an error
- descriptor length mismatches fail with an error
- zero selected atoms fail with an error
- failures return no indices

This is intentionally narrower than Python `eval`; it is meant to fail gently
and stop rather than guessing.

## Deferred

- distance, angle, and hydrogen-bond selections
- contextual named selections such as `backbone` or `calpha`
- arbitrary Python calls, slicing, arithmetic, and list membership
- full Python expression compatibility
