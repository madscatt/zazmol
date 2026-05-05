# BIOMT Plan (Metadata + Transform, Fixture-First)

Python BIOMT behavior now has two explicit layers that must stay separate in
standalone C++ parity work:

1. `read_pdb()` records BIOMT metadata from `REMARK 350` header lines.
2. coordinate transforms are explicit helper calls and are not auto-applied.

This document defines the guardrails and fixture-backed contracts for both
layers so c5.3 can proceed without guessing.

## Scope For The Next Implementation Slices

- Add BIOMT metadata parity before BIOMT transform implementation parity.
- Keep metadata parsing passive: record only, no coordinate mutation.
- Keep transform APIs explicit and manual-only.
- Do not combine BIOMT work with selection grammar expansion, topology parser
  work, DCD policy changes, or merge policy expansion.

## Guardrails

- Python `zazmol` remains the read-only behavior oracle.
- Public APIs must use SasMol-owned types and remain backend-neutral.
- Failures must be explicit and non-partial.
- No hidden global state or logging side effects.
- No mutation on failed validation.
- No large binary fixtures in git; use generated test inputs.
- No automatic transform application during PDB read.

## Python Delta (May 5, 2026)

Behavior reference changed in Python commit `da1d399`:

- Added passive BIOMT metadata parsing from `REMARK 350` headers during
  `read_pdb()`.
- Added BIOMT metadata storage accessors on `Molecule`.
- Clarified `apply_biomt` and `copy_apply_biomt` as explicit transform helpers
  that do not parse headers or assemble full biological units.

## Phase 3 Design Decision

BIOMT metadata should be molecule-owned, not parser-owned. This mirrors Python's
`Molecule.biomt()` storage and keeps metadata available after `read_pdb()`
returns.

The PDB reader should be responsible for passive `REMARK 350` parsing because
the metadata is part of file input behavior. The reader must populate the
molecule's BIOMT map after a successful PDB read, just as Python does after
header capture. Reading must not apply BIOMT transforms or resize/duplicate
coordinate storage.

Transform helpers remain subset/operate-style coordinate APIs and should be
implemented later, after metadata parity is complete.

## Approved C++ Metadata API Shape

Metadata record:

```cpp
struct BiomtRecord {
  std::vector<std::string> subdivs;
  std::string auth_bio_unit;
  std::string soft_bio_unit;
  std::vector<std::array<std::array<calc_type, 3>, 3>> rot;
  std::vector<std::array<calc_type, 3>> trans;
};
```

Metadata container:

```cpp
using BiomtMap = std::map<int, BiomtRecord>;
```

Molecule accessors:

```cpp
[[nodiscard]] BiomtMap& biomt() noexcept;
[[nodiscard]] const BiomtMap& biomt() const noexcept;
void set_biomt(BiomtMap value);
```

PDB parse helper:

```cpp
[[nodiscard]] BiomtMap parse_biomt_header_records(
    const std::vector<std::string>& header_lines);
```

Parser boundary:

- expose the helper from `sasmol/file_io.hpp` so parser behavior can be unit
  tested without constructing a PDB file
- call the same helper from `PdbReader::read_pdb`
- store the result on `Molecule::biomt()`
- clear BIOMT metadata on `Molecule::resize`, matching the rest of molecule
  state reset behavior

## Transform API Shape (Deferred Review Target)

Transform helper candidates:

In-place worker:

```cpp
[[nodiscard]] SubsetResult apply_biomt_transforms(
    Molecule& molecule, std::size_t frame,
    const std::vector<BiomtTransform>& transforms);
```

Pure value-returning helper:

```cpp
[[nodiscard]] Molecule biomt_transformed(
    const Molecule& molecule, std::size_t frame,
    const std::vector<BiomtTransform>& transforms);
```

Candidate transform type:

```cpp
struct BiomtTransform {
  std::array<std::array<calc_type, 3>, 3> rotation;
  std::array<calc_type, 3> translation;
};
```

Notes:

- transform API names remain review targets, not finalized ABI.
- Frame argument is explicit; no implicit current-frame behavior.
- Transform sequence order is caller-defined and preserved.

## Metadata Semantics To Lock Before Implementation

1. Parsing trigger:
   Parse only `REMARK 350` header lines.
2. Passive behavior:
   Store metadata only; do not alter coordinates.
3. BIOMOLECULE grouping:
   Group by integer biomolecule id.
4. Chain list support:
   Handle both `APPLY THE FOLLOWING TO CHAINS:` and `AND CHAINS:` continuations.
5. BIOMT row assembly:
   Accept only complete row triplets (`BIOMT1/2/3`) per transform id.
6. Robustness:
   Ignore malformed lines rather than guessing.
7. Read behavior:
   Empty or absent BIOMT headers produce an empty map.
8. Numeric policy:
   Rotation and translation values are parsed into `calc_type`.

## Transform Semantics To Lock Before Implementation

1. Transform order:
   Apply transforms in input order. Output assembly order must match transform
   order.
2. Frame scope:
   BIOMT should consume one explicit source frame per call.
3. Assembly expansion:
   Each transform contributes one transformed copy of the source coordinates.
4. Coordinate precision:
   Calculations use `calc_type`, storage follows `coord_type` policy.
5. Validation:
   Reject bad frame, empty transform set, non-finite matrix/translation entries,
   and non-rigid rotation matrices if policy requires rigidity checks.
6. Failure mode:
   Non-throwing worker returns structured errors with no partial mutation.

## Test Scaffolding In This Step

`tests/test_biomt_fixture_semantics.cpp` adds generated fixture-style checks for:

- deterministic transform-order behavior
- explicit frame targeting
- source-coordinate immutability while assembling transformed outputs

These tests do not call BIOMT production APIs yet. They lock fixture semantics
and expected transformed coordinates for the implementation step.

`tests/test_biomt_metadata_fixture_semantics.cpp` adds generated fixture-style
checks for:

- empty metadata when no BIOMT header lines are present
- one complete transform record parsing
- multi-chain continuation parsing (`APPLY...` + `AND CHAINS...`)
- incomplete BIOMT triplets being ignored

These tests also avoid production BIOMT APIs and lock parsing expectations
before implementation.

## Next Implementation Slice

Implement metadata parity only:

1. Add `BiomtRecord` and `BiomtMap` to `sasmol/molecule.hpp`.
2. Add `Molecule::biomt()` accessors and reset behavior.
3. Add `parse_biomt_header_records(...)` in file I/O.
4. Wire `PdbReader::read_pdb` to store parsed metadata after successful header
   capture.
5. Convert metadata fixture semantics tests to call production APIs.
6. Add a PDB-read test proving identity BIOMT metadata does not mutate
   coordinates.

Stop after this slice. Do not implement coordinate transform helpers in the same
commit.

## Explicit Non-Goals In This Step

- no BIOMT implementation in `src/subset.cpp`
- no BIOMT metadata implementation in `src/file_io.cpp` or parser internals yet
- no parser/selection alias expansion
- no topology parser or PDB-driven force-field inference
