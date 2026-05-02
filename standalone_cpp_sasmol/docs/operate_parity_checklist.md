# Operate Parity Checklist

Python `zazmol.operate` remains the behavior oracle for coordinate mutation,
rotation, and alignment behavior. This area is intentionally split into small
steps because sign conventions, row/column multiplication order, and in-place
mutation are user-visible.

## API Policy

Public operate functions should offer both:

- in-place mutation for Python SasMol parity and performance
- pure copy-returning forms for safer composition and testing

Internal implementation should be shared on explicit coordinate views and
owned SasMol types. Public APIs must not expose Eigen or backend-specific
storage.

## First Implemented Slice

Implemented primitive frame operations:

- `translate(molecule, frame, value, options)`
- `translated(molecule, frame, value, options)`
- `center(molecule, frame)`
- `centered(molecule, frame)`
- `rotate(molecule, frame, axis, theta)`
- `rotated(molecule, frame, axis, theta)`
- `rotate_general_axis(molecule, frame, theta, unit_axis)`
- `rotated_general_axis(molecule, frame, theta, unit_axis)`
- `rotate_euler(molecule, frame, phi, theta, psi)`
- `rotated_euler(molecule, frame, phi, theta, psi)`

Behavior notes:

- `translate(..., {.point = true})` moves the frame center of mass to the
  requested point, matching Python's `point=True` behavior.
- `rotate(axis)` preserves Python's column-vector multiplication convention.
- `rotate_general_axis` and `rotate_euler` preserve Python's row-vector
  multiplication convention.
- `rotate_general_axis` intentionally accepts non-unit axes because Python
  tests and legacy behavior use non-unit values despite the method name and
  docstring.

## Deferred

## PMI Alignment Slice

Implemented PMI alignment operations:

- `align_pmi_on_axis(molecule, frame, pmi_eigenvector, alignment_axis)`
- `pmi_aligned_on_axis(molecule, frame, pmi_eigenvector, alignment_axis)`
- `align_pmi_on_cardinal_axes(molecule, frame)`
- `pmi_aligned_on_cardinal_axes(molecule, frame)`

Behavior notes:

- the method centers the selected frame before alignment, matching Python
  behavior
- the rotation axis uses Python's `cross(target_axis, pmi_axis)` convention
  because this is coupled to the row-vector general-axis rotation path
- tests assert absolute axis alignment, not eigenvector sign identity
- tests assert eigenvalue preservation across PMI alignment rather than
  arbitrary eigenvector signs
- PMI alignment leaves the selected frame centered at the origin
- singular PMI tensors are rejected explicitly

## Deferred

## Basis Alignment Slice

Implemented explicit-index basis alignment:

- `initialize_alignment(moving, reference, moving_basis_indices,
  reference_basis_indices, frame)`
- `align(moving, plan, frame)`
- `aligned(moving, plan, frame)`

Behavior notes:

- this ports the rigid-body fit core used by Python `align`
- C++ callers provide explicit atom indices rather than Python basis strings
- the transform is applied to the whole moving molecule frame, not just the
  basis atoms
- alignment centers each selected basis on its own mass-weighted basis center
  of mass, matching Python's subset-copy alignment path
- pure `aligned(...)` returns a copy and leaves the source molecule unchanged

Deferred from Python `align`:

- initialization/production mode keyword interface
