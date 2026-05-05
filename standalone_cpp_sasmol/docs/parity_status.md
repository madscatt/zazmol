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
  hardening, normal writes, selected-frame range writes, C++ round trips,
  Python cross-reader tooling, malformed input/status handling, and lifecycle
  misuse checks.
- Calculation coverage for mass, molecular formula, residue charge, center of
  mass, radius of gyration, RMSD, min/max, DCD streaming min/max, and principal
  moments of inertia.
- Linear algebra coverage for cross products, matrix multiplication, legacy
  alignment matrix construction, vector helpers, Python-style comparison,
  signed angles, dihedral angles, and ordinary three-point angles.
- Operate coverage for translate, center, axis/general/Euler rotations, PMI
  alignment, explicit-index and basis-expression alignment initialization and
  production helpers, pure copy-returning variants, row and column convention
  preservation, legacy average vdW radius assignment, and
  failure-before-mutation behavior.
- Selection/subset coverage for a bounded safe expression grammar, surveyed
  SASSIE/VMD-like basis translation, named `all` and `heavy` bases, contextual
  `backbone`/`calpha` compatibility helpers, 0/1 mask bridge helpers,
  coordinate get/set/copy, molecule copy/duplicate/merge, dihedral subset masks,
  descriptor get/set, extension descriptors, Python-style selected BIOMT
  apply/copy-apply helpers, and structured failure returns.
- Overlap coverage for coordinate-vector and molecule-frame overlap checks.
- Property table coverage for atomic masses, amino-acid SLD values, element,
  nucleotide, DNA, RNA, and protein scattering lengths, and van der Waals radii
  including Python `None` entries represented as empty optionals.
- Topology support for explicit CHARMM type/charge assignment from
  already-trusted atom-aligned vectors or atom-name tables; Python-parity
  CHARMM atom-name classification lists; Python-parity CHARMM topology parsing
  for reviewed record types; pure residue atom-list, patch, atom-order choice,
  residue reorder, whole-molecule reorder planning, reordered-copy, explicit
  in-place reorder helpers, and FASTA creation from residue descriptors;
  Python-style index/resid renumbering; constraint descriptor/PDB generation;
  FASTA backbone molecule/PDB generation; no PDB-name inference and no partial
  mutation on validation failure.
- Normal and sanitizer C++ test runs are part of the checkpoint rhythm.

## Intentional Boundaries

- This is not a drop-in Python binding yet; it is a standalone C++ core.
- Python `eval` selection compatibility is intentionally not reproduced. The C++
  parser is bounded and fails rather than guessing.
- Generic `backbone` selection remains rejected without context because
  Python/SASSIE usage is context-dependent across protein and nucleic-acid
  workflows. C++ exposes explicit protein, nucleic, and nucleic-overlap contexts.
- `rotate_general_axis` preserves legacy Python behavior for non-unit axes; only
  unit-axis use is treated as a true rigid rotation.
- PDB reading does not infer CHARMM atom types.
- CHARMM topology support is still a minimal parity port, not a full topology
  engine. It preserves Python's current behavior and known limits, including
  atom-only patching and no automatic PDB-read typing or reordering.
- BIOMT support is deliberately split into passive metadata preservation,
  Python-style selected one-transform helpers, and optional coordinate-only
  assembly helpers. PDB reading still does not apply transforms automatically.
- VMD/view support is an optional adapter: the portable C++ API mirrors
  Python's frame-to-X/Y/Z coordinate extraction, and the promoted legacy IMD/VMD
  C sender is compiled only when requested.
- Fixed/free atom DCD variants and DCD unit-cell writing are deliberate
  unsupported-status paths with tests. True random-access DCD seeking is not
  exposed; single-frame access remains explicit reopen-and-scan behavior.
- Python multiprocessing support is experimental orchestration code and is not
  planned for the standalone C++ core.
- CHARMM36 support is tracked as separate joint Python/C++ design work, not as
  part of this parity pass.
- Large DCD workflows should use explicit streaming APIs; whole-trajectory DCD
  reads remain convenience behavior for bounded data.
- GPU/MPI/Python bindings remain architectural considerations, not v1
  dependencies.

## Recommended Next Slice

No broad in-scope core module remains as an obvious next port target. The next
slice should be chosen from a concrete caller need or fixture-backed parity gap:
selection grammar expansion from surveyed real expressions, real-VMD runtime
smoke testing if that workflow matters, or unusual PDB/DCD/topology fixture
hardening.
