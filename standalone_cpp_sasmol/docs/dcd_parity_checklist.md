# DCD Parity Checklist

This checklist is the implementation gate for the standalone C++ `DcdReader`
and `DcdWriter`. It is based on a read-only survey of:

- `src/python/dcd_io.py`
- `src/python/extensions/dcdio/dcdio.c`
- `src/python/extensions/dcdio/dcdio.h`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_open_dcd_read.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_read_dcd.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_read_dcd_step.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_read_single_dcd_step.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_write_dcd.py`
- `src/python/test_sasmol/test_file_io/test_intg_file_io_Files_write_dcd_frames.py`

Python `zazmol` and the existing `_dcdio` extension remain the behavior oracle.

## Current Behavior To Preserve First

- `open_dcd_read(filename)` returns a stream descriptor equivalent to:
  file handle, atom count, frame count, endian flag, and CHARMM flags.
- Missing files raise cleanly before header parsing.
- Truncated or malformed headers raise a Python exception rather than crashing
  the caller.
- `read_dcd(filename)` reads the whole trajectory into one molecule coordinate
  array.
- `read_dcd_step(open_descriptor, frame)` reads the next frame from the current
  stream position. The `frame` argument is passed through to the C function as
  the `first` flag, but practical use is sequential.
- `read_single_dcd_step(filename, frame)` opens the file, reads the header,
  scans forward, reads the selected frame, closes the file, and stores one
  frame. This is reopen-and-scan behavior, not true random access.
- `write_dcd(filename)` writes every frame.
- `write_dcd_frames(filename, start, end)` writes `[start, end)` frames.
- DCD coordinate arrays are written as separate X, Y, Z float arrays.

## Binary Details To Test Before Implementation

- Initial Fortran record marker must be `84`, with endian reversal detection.
- Header ID must be `CORD`.
- Header fields to preserve:
  - `NSET`
  - `ISTART`
  - `NSAVC`
  - `DELTA`
  - `NAMNF`
  - atom count
  - reverse-endian flag
  - CHARMM flags
- CHARMM detection depends on nonzero data in the last integer of the 84-byte
  header block.
- CHARMM extra-block and 4D flags must be recognized and skipped.
- X, Y, and Z frame blocks each use record markers of `4 * natoms` for normal
  full-frame reads.
- The old code has machinery for fixed/free atoms through `NAMNF`; even if v1
  does not implement fixed atoms, it must fail with a precise status rather than
  silently corrupting coordinates.
- Large-file behavior matters. The old tests exercise generated 1.0 GB through
  6.4 GB DCDs behind `SASMOL_HUGETEST`.

## Minimum Fixture Parity

Small required fixtures:

| Fixture | Expected atoms | Expected frames | Expected flags |
| --- | ---: | ---: | --- |
| `src/python/test_sasmol/data/dcd_common/1ATM.dcd` | 1 | 2 | `reverseEndian=0`, `charmm=5` |
| `src/python/test_sasmol/data/dcd_common/2AAD.dcd` | 15 | 3 | `reverseEndian=0`, `charmm=5` |
| `src/python/test_sasmol/data/dcd_common/rna-1to10.dcd` | 10632 | 10 | `reverseEndian=0`, `charmm=5` |

Coordinate parity must include:

- full-trajectory read for all three fixtures
- sequential `read_next_frame` for frame 1, middle frame where available, and
  final frame
- explicit reopen-and-scan single-frame read for the same frames
- sums using calculation precision and point samples matching Python tests
  within the existing DCD tolerances

Write parity must include:

- PDB to DCD to readback for `1ATM`, `2AAD`, `rna`, `rna-1to10`, and `1CRN`
- `write_dcd_frames` for first frame, last frame, and full range where fixtures
  exist
- generated files must round-trip through the C++ reader and Python `zazmol`
  before replacing any legacy path

## Implementation Guardrails

- Keep the public C++ reader sequential by default.
- Do not expose hidden random access unless the implementation truly supports
  it. Reopen-and-scan should remain explicit.
- Return `IoStatus`; do not `exit`, print progress, or impose logging policy.
- Treat any unsupported DCD variant as `IoCode::unsupported` or
  `IoCode::format_error`, with enough detail for callers to log.
- Keep coordinate storage conversion explicit: DCD blocks are X/Y/Z arrays;
  `Molecule` storage is frame-major atom `xyz` triplets.
- Do not require GPU, MPI, Eigen, or Python bindings for DCD v1.

## Proposed Next Implementation Slice

Completed first slice: **header-only C++ DCD parsing tests and implementation**.

1. Added binary little/big-endian record-marker helpers.
2. Implemented `DcdReader::open_dcd_read` and `read_header`.
3. Populated `DcdHeader` with atom count, frame count, unit-cell flag, CHARMM
   flag, endian flag, `istart`, `nsavc`, `delta`, and `namnf`.
4. Added C++ tests for the three small DCD fixtures above.
5. Stopped before frame-coordinate parsing.

## Proposed Next Implementation Slice

Add sequential frame-coordinate parsing for normal full-coordinate DCD files:

1. Implement CHARMm unit-cell/extra-block skipping before each frame.
2. Read X, Y, and Z float blocks with record-marker validation.
3. Convert DCD X/Y/Z arrays into `Molecule` frame-major atom `xyz` triplets.
4. Add tests for first/middle/final frame samples and coordinate sums for the
   three small fixtures.
5. Keep fixed/free atom DCDs as explicit `IoCode::unsupported`.

This is the next risky binary step and should stay separate from DCD writing.
