# Eigen Usage Policy

Eigen may be used as a CPU-side helper library, but it must not define the core
SasMol data model or GPU interface.

## Role Of Eigen

Allowed:

- small vector/matrix math
- local numerical helpers
- CPU-side convenience operations

Not allowed:

- primary storage for simulation data
- data structures passed to Metal or CUDA
- replacing existing buffer-based representations

Eigen is a helper layer, not the core representation.

## GPU Interface Rule

Before any GPU dispatch:

- data must be in plain contiguous buffers
- no Eigen objects cross the CPU to GPU boundary
- no Eigen-dependent layouts are assumed

## Precision Model

- storage coordinates follow `COORD_DTYPE` semantics, typically `float32`
- calculations and reductions follow `CALC_DTYPE` semantics, typically
  `float64`

Backend kernels may use their native precision policy, documented per backend.

Conversions:

- occur only at explicit backend boundaries
- do not happen implicitly within core logic

Validation:

- uses absolute and relative error thresholds
- does not require bitwise equality

Future backend note: CUDA `float32` kernels may be worth evaluating for
scattering, GV, and Debye-style workloads where experimental uncertainty and
throughput tradeoffs differ from coordinate/topology parity calculations. No
CUDA precision policy is being set here; core SasMol parity remains
`COORD_DTYPE` for coordinate storage and `CALC_DTYPE` for calculations and
reductions.

## Performance Guideline

Eigen is suitable for:

- small/local computations
- clarity and correctness

Eigen is not intended for:

- large-N loops such as GV, Debye, or PBC kernels
- performance-critical inner loops
- GPU-bound execution paths

Core kernels must use explicit loops and buffers.

## Licensing And Redistribution

Eigen is licensed under MPL2, the Mozilla Public License 2.0.

Policy requirement:

- define `EIGEN_MPL2_ONLY` to disable non-MPL2/LGPL features

Implications:

- Eigen is safe to include in public repositories under this policy
- MPL2 file-level copyleft does not extend to the project's own source files
- include Eigen's license file if Eigen is vendored

Eigen is header-only and can be vendored into the repository or used as an
external dependency.

## Dependency Policy

- keep Eigen usage minimal
- avoid deep coupling to Eigen types
- prefer explicit, backend-neutral data structures

Eigen types:

- may be used in private `.cpp` implementation details
- must not appear in public headers unless explicitly approved

Rationale:

- preserves ABI stability
- keeps compile times reasonable
- maintains clean GPU and API boundaries

## Final Rule

Eigen is a CPU-side helper only. All core data structures and
performance-critical paths must remain explicit, backend-neutral buffers.
