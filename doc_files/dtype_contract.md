### SASMOL Dtype Contract

This project uses the repository name `zazmol`, but the code package remains
`sasmol`.

The dtype policy is intentionally narrow and explicit:

- molecular coordinate storage uses `config.COORD_DTYPE`
- derived numerical calculations use `config.CALC_DTYPE`

Current definitions:

```python
COORD_DTYPE = numpy.float32
CALC_DTYPE = numpy.float64
```

#### Coordinate Storage: `COORD_DTYPE`

Use `config.COORD_DTYPE` for coordinate storage and coordinate-transfer
buffers.

This includes:

- `Molecule._coor`
- coordinates read from PDB/DCD files
- coordinates written to PDB/DCD/VMD interfaces
- coordinates produced by subset, copy, merge, and duplicate operations
- trajectory buffers and GPU-facing coordinate transfer arrays

The reason is practical and historical: coordinate-heavy workflows dominate
memory use, especially when SASSIE calculations operate on large trajectories
or pass coordinate data toward CUDA-backed code. Doubling coordinate storage to
`float64` is not acceptable without a separate, explicit performance and memory
review.

#### Derived Calculations: `CALC_DTYPE`

Use `config.CALC_DTYPE` for derived numerical calculations, reductions, and
scalar/vector/tensor outputs that accumulate over coordinates.

This includes:

- coordinate sums used by tests
- center of mass
- radius of gyration
- principal moments of inertia
- matrix and vector helper calculations
- residue charges and similar derived physical properties

The reason is numerical stability: reductions over many `float32` coordinates
can accumulate meaningful error if the accumulation itself is also `float32`.

#### Review Rules

New numerical code should not use bare `dtype=float`, `numpy.float32`,
`numpy.float64`, `np.float32`, or `np.float64` directly unless there is a local,
documented reason.

Use:

- `config.COORD_DTYPE` for coordinate storage and coordinate buffers
- `config.CALC_DTYPE` for derived calculations and reductions

Do not change the dtype contract, expected numerical outputs, or test
tolerances just to make tests pass. Any dtype-policy change must document:

- old behavior
- new behavior
- reason for the change
- affected code paths
- affected tests or reference data
- comparison to legacy behavior where possible

#### GitHub Issue Text

Suggested issue summary:

`sasmol` historically mixed `float32` coordinate storage with Python/NumPy
default floating-point calculations. The port now makes that implicit behavior
explicit through `config.COORD_DTYPE` and `config.CALC_DTYPE`: coordinates and
coordinate-transfer buffers remain `float32`, while derived calculations and
reductions use `float64`. This preserves memory behavior for large trajectories
and CUDA-facing workflows while avoiding avoidable accumulation error in
calculated properties.

Acceptance criteria:

- coordinate storage and coordinate-transfer code uses `config.COORD_DTYPE`
- derived calculation and reduction code uses `config.CALC_DTYPE`
- tests assert the intended dtype boundaries
- direct `numpy.float32` / `numpy.float64` usage is limited to configuration and
  tests that verify configuration
- full normal, large, and huge test suites pass
