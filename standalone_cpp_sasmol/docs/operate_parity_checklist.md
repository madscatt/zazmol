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

- `align_pmi_on_axis`
- `align_pmi_on_cardinal_axes`
- full basis/Kabsch-style `align`
- sign-convention tests for PMI alignment
