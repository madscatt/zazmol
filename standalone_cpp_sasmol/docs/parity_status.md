# Standalone C++ SasMol Parity Status

This checkpoint summarizes the standalone C++ library state after the first
round of module-by-module parity work. Python `zazmol` remains the read-only
behavior oracle.

## Current Strengths

- Core `sasmol::Molecule` data model with explicit descriptor vectors,
  contiguous frame-major coordinates, connectivity, formula, FASTA, unit cell,
  moltype reporting, and typed extension descriptor maps.
- Tolerant PDB reader/writer coverage for fixed-column records, one-frame and
  multi-frame inputs, `END` and `MODEL` trajectories, selected/all-frame writes,
  optional-field handling, `CONECT` output, all-zero coordinate guarding, element
  resolution, and Python cross-reader validation of generated C++ files.
- Sequential DCD reader/writer coverage for normal full-coordinate files,
  explicit reopen-and-scan single-frame reads, whole-trajectory convenience
  reads, caller-owned frame-buffer streaming, generated multi-frame DCD
  hardening, normal writes, C++ round trips, Python cross-reader tooling,
  malformed input/status handling, and lifecycle misuse checks.
- Calculation coverage for mass, molecular formula, residue charge, center of
  mass, radius of gyration, RMSD, min/max, DCD streaming min/max, and principal
  moments of inertia.
- Operate coverage for translate, center, axis/general/Euler rotations, PMI
  alignment, explicit-index basis alignment, pure copy-returning variants, row
  and column convention preservation, and failure-before-mutation behavior.
- Selection/subset coverage for a bounded safe expression grammar, named `all`
  and `heavy` bases, 0/1 mask bridge helpers, coordinate get/set/copy, molecule
  copy/duplicate/merge, descriptor get/set, extension descriptors, and structured
  failure returns.
- Topology support for explicit CHARMM type assignment from already-trusted
  atom-aligned type vectors or atom-name/type tables, with no PDB-name
  inference and no partial mutation on validation failure.
- Normal and sanitizer C++ test runs are part of the checkpoint rhythm.

## Intentional Boundaries

- This is not a drop-in Python binding yet; it is a standalone C++ core.
- Python `eval` selection compatibility is intentionally not reproduced. The C++
  parser is bounded and fails rather than guessing.
- Generic `backbone` selection remains deferred because Python/SASSIE usage is
  context-dependent across protein and nucleic-acid workflows.
- `rotate_general_axis` preserves legacy Python behavior for non-unit axes; only
  unit-axis use is treated as a true rigid rotation.
- PDB reading does not infer CHARMM atom types.
- CHARMM topology parsing, patching, atom completeness checks, and atom
  reordering remain separate design work.
- BIOMT subset operations remain deferred.
- Fixed/free atom DCD variants remain unsupported unless fixtures and policy are
  added.
- Large DCD workflows should use explicit streaming APIs; whole-trajectory DCD
  reads remain convenience behavior for bounded data.
- GPU/MPI/Python bindings remain architectural considerations, not v1
  dependencies.

## Recommended Next Slice

Do a small selection/subset documentation cleanup before new feature work:

- update stale deferred notes that still imply mask/subset bridge operations are
  missing
- add a concise selection/subset checkpoint like the operate and topology
  checkpoints
- then decide whether to tackle BIOMT, opt-in large DCD streaming stress tests,
  or broader Python expression coverage as a reviewed design choice

The next implementation-heavy step should not be a CHARMM parser or arbitrary
selection evaluator without a separate plan.
