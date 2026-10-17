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

## Unavailable PDB atom serials

An exact `*****` ATOM/HETATM serial uses the one-based atom position in
`original_index`; `index` remains sequential. This substitutes an unavailable
source serial rather than recovering it. Preserved source serials must not
collide with substituted positions. CONECT endpoints must not depend on those
substitutions. The complete reader checks the first-frame map and returns
`IoCode::format_error` on ambiguity, leaving the destination molecule unchanged.
The standalone record helper requires a positive `atom_ordinal_context` to
accept stars; it cannot validate whole-file collisions by itself. Other invalid
fields and VMD-hex rules are unchanged. No new descriptor or writer behavior
is introduced. Python uses its existing `ValueError` mechanism for these cases.


### Star-serial validation checkpoint (2026-10-17)

Validated on the Mac mini against baseline `989a3dd` in separate private
Python 3.9 environments, without changing production installations:

- Focused Python reader checks: 36 passed.
- Complete SasMol Python suite: baseline 913 passed / 1 skipped; candidate
  919 passed / 1 skipped, including six new reader regressions. The optional
  archive-layout metadata check skipped in both isolated layouts also passed
  in a paired recheck against the actual neighboring application checkout.
  Large DCD tests, including the 6.4 GB fixture, were enabled.
- Standalone C++ Debug CTest: 17/17 passed in each build, including the new
  starred-serial, collision, connectivity and unchanged-destination checks.
- SASSIE PDBScan/PDBRx checks: 1,531 passed with each SasMol version.
- A 100,001-atom PDB/PSF identity probe accepted unavailable serials while
  still rejecting a changed atom name.
- Sphinx documentation built with two existing docstring diagnostics.
- After correcting only the integration harness and private ProDy setup,
  17 affected SASSIE integration checks produced identical outcomes:
  15 passed and 2 unrelated failures in each environment, with no collection
  errors or skips. This does not claim an entirely passing SASSIE suite.
  No SASSIE test fixes are included in this patch.

The paired rerun used original test/fixture paths and unchanged expectations;
only the test-helper namespace was registered for collection. Scientific
packages came from installed distributions. Hash checks found no changes to
12,338 protected SASSIE repository files, production SASSIE/SasMol packages,
production ProDy files or the production version-stack binding during that
rerun. No tolerances, input contracts or scenario catalogs were changed.

Machine-local evidence is retained under
`/Volumes/ten_tee_bee/pdbrx_validation/sasmol_star_serial_20261016`, including
`candidate_source_hashes.json`, the per-suite commands and results, and
`affected_integration_recheck_20261017/{selected_nodes.json,audit.json,report.txt}`.
These paths describe validation evidence, not portable runtime dependencies.
