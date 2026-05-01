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

The initial C++ API therefore exposes `PdbReadOptions` and `PdbWriteOptions`,
but parsing and writing return `IoCode::not_implemented` until fixture parity
tests are in place.

## DCD

DCD is treated as a binary compatibility problem first and a cleanup problem
second. The primary reader API is sequential:

```cpp
sasmol::DcdReader reader;
reader.open_dcd_read("trajectory.dcd");
reader.read_header(header);
reader.read_next_frame(molecule);
reader.close_dcd_read();
```

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
- failure modes that return library status instead of terminating callers

## Error Handling

The standalone C++ library returns `IoStatus` for file I/O operations. It should
not call `exit`, print directly to stdout/stderr, or impose application logging
policy. Applications can translate `IoStatus` into their own logging,
exceptions, GUI messages, or workflow errors.
