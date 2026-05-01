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

1. Implemented CHARMm unit-cell/extra-block skipping before each frame.
2. Read X, Y, and Z float blocks with record-marker validation.
3. Converted DCD X/Y/Z arrays into `Molecule` frame-major atom `xyz` triplets.
4. Added tests for first/middle/final frame samples for the three small
   fixtures.
5. Keep fixed/free atom DCDs as explicit `IoCode::unsupported`.

This slice remains separate from DCD writing.

## Proposed Next Implementation Slice

Completed reader hardening slice before writer behavior:

1. Added coordinate sum parity for the same frames already sampled.
2. Added `read_single_dcd_step` reopen-and-scan behavior using one-based frame
   numbering to match Python `zazmol`.
3. Added malformed/truncated-frame tests that return `IoStatus` errors.
4. Kept fixed/free atom DCDs as explicit `IoCode::unsupported`.

## Proposed Next Implementation Slice

Completed whole-trajectory convenience reading:

1. Added `DcdReader::read_dcd(filename, molecule)` as a convenience wrapper over
   open/read_header/repeated `read_next_frame`/close.
2. Tested whole-file read parity for `1ATM.dcd`, `2AAD.dcd`, and
   `rna-1to10.dcd`.
3. Preserved the same sequential internals rather than adding hidden random
   access.
4. Stopped before DCD writing.

## Proposed Next Implementation Slice

Completed first DCD writing slice with the same narrowness:

1. Added `DcdWriter::open_dcd_write`, `write_dcd_header`, `write_dcd_step`, and
   `close_dcd_write` behavior for normal full-coordinate files only.
2. Wrote non-unit-cell headers matching the existing Python writer path for
   normal full-coordinate files.
3. Round-tripped C++ written DCDs through the C++ reader for `1ATM` and `2AAD`
   first.
4. Deferred huge-file and cross-reader Python validation to a separate slice.

## Proposed Next Implementation Slice

Completed cross-reader preparation and convenience write wrapper:

1. Added `DcdWriter::write_dcd(filename, molecule)` as
   open/header/all-steps/close.
2. Tested C++ writer output read by C++ reader for `rna-1to10`.
3. Added a Python-side validation script to verify C++
   generated DCDs with Python `zazmol` before claiming full parity.
4. Kept huge-file tests separate.

Python-side validation command:

```bash
python3 standalone_cpp_sasmol/tools/validate_dcd_with_python.py path/to/file.dcd
```

## Proposed Next Implementation Slice

Completed deliberate Python cross-reader parity fixture tooling:

1. Added a C++ `dcd_roundtrip_writer` tool to generate C++ written DCD files
   from existing DCD fixtures.
2. Validate generated files with Python `zazmol` using
   `tools/validate_dcd_with_python.py`.
3. Record results in docs without making Python an unconditional CTest
   dependency.
4. Then decide whether to start PDB I/O contracts or expand DCD writer options.

## Memory-Safety Status

The DCD reader/writer uses `std::ifstream`, `std::ofstream`, `std::vector`, and
`std::array` for ownership. It does not use raw owning pointers. The binary
parsing path includes explicit checks for:

- negative header counts
- record-marker size overflow
- unsupported fixed/free atom tables
- unsupported unit-cell writes
- malformed/truncated input returning `IoStatus`

The current C++ tests pass under AddressSanitizer and UndefinedBehaviorSanitizer
using the documented `SASMOL_ENABLE_SANITIZERS=ON` build.

## Additional DCD Hardening

Completed malformed-input and misuse coverage:

- truncated headers return `IoStatus`
- bad `CORD` magic returns `IoStatus`
- truncated frame reads return `IoStatus`
- `read_dcd` closes the reader on failure
- `read_single_dcd_step` closes its internal reader when reading past EOF
- unit-cell DCD writing is explicit `IoCode::unsupported`
- out-of-range DCD write frames return `IoStatus`

These paths pass in both the normal build and the sanitizer build.

## Proposed Next Implementation Slice

Completed targeted binary-header edge tests:

1. Invalid title block size/count.
2. Negative header counts.
3. Oversized/unsupported record marker paths where feasible without allocating
   huge memory.
4. Kept these as status-return tests rather than exceptions or process exits.

These paths pass in both the normal build and the sanitizer build.

## Proposed Next Implementation Slice

Add explicit lifecycle and state-machine hardening:

1. Repeated open/close calls are safe.
2. Re-reading the header resets frame position deliberately.
3. Read/write calls after close return `IoCode::not_open`.
4. Temporary files are cleaned by tests.
