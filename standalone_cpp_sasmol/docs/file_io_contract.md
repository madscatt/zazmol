# File I/O Contract

This milestone defines the public C++ file I/O surface before porting parser
internals. Python `zazmol` and the DCD C-extension behavior remain the parity
standards.

## PDB

PDB parsing should be tolerant by default. The C++ implementation must preserve
the hard-earned Python behavior before it becomes the default path:

- accept common non-ideal PDB variants rather than enforcing a protein-only
  schema
- preserve multi-model semantics, including `MODEL`, `ENDMDL`, `END`, and mixed
  legacy file shapes
- keep `CONECT` records and remapping behavior
- preserve descriptor defaults where fields are absent or partially populated
- keep element-resolution behavior explicit and test-backed
- preserve the all-zero coordinate column guard used before writing

The C++ API exposes `PdbReadOptions` and `PdbWriteOptions`. The current
implementation covers the first tolerant read/write parity surface with fixture
tests and Python cross-reader validation tooling; additional unusual PDB shapes
should still be added only with fixtures and reviewed behavior notes.

The detailed PDB implementation gate is tracked in
[`pdb_parity_checklist.md`](pdb_parity_checklist.md).

## DCD

DCD is treated as a binary compatibility problem first and a cleanup problem
second. The primary reader API is sequential:

```cpp
sasmol::DcdReader reader;
reader.open_dcd_read("trajectory.dcd");
reader.read_header(header);
reader.read_next_frame_coordinates(coordinates);
reader.close_dcd_read();
```

For large trajectories, callers should process each frame and then reuse or
discard the frame buffer. `read_next_frame(molecule)` is also sequential, but it
writes into a molecule trajectory container. `read_dcd(filename, molecule)` is a
whole-trajectory convenience wrapper and may allocate storage proportional to
`natoms * nframes * 3`; it should be reserved for fixtures, tests, and
moderate-size data.

Writing supports normal full-coordinate DCD files through explicit
open/header/step/close calls, all-frame `write_dcd(...)`, and bounded
`write_dcd_frames(filename, molecule, start, end)` using zero-based
`[start, end)` frame ranges.

Random-frame behavior should not be hidden as normal seekable access. If
supported, it should be explicit reopen/scan behavior through
`read_single_dcd_step(...)` or a similarly named API. This matches the observed
practical behavior of the current stack and avoids promising random access
until the legacy C-extension contract proves it is safe.

The DCD parity pass must cover:

- header fields and frame counts
- atom counts
- unit-cell/no-unit-cell variants
- CHARMM/NAMD-style flags where present
- sequential frame reading
- explicit close/reopen behavior for cycling or random-frame consumers
- write/read round trips
- selected-frame range writes
- failure modes that return library status instead of terminating callers

The detailed DCD implementation gate is tracked in
[`dcd_parity_checklist.md`](dcd_parity_checklist.md).

## Error Handling

The standalone C++ library returns `IoStatus` for file I/O operations. It should
not call `exit`, print directly to stdout/stderr, or impose application logging
policy. Applications can translate `IoStatus` into their own logging,
exceptions, GUI messages, or workflow errors.
