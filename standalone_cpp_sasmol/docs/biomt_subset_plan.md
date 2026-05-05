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

## Proposed C++ API Shape (Review Target)

Metadata record candidate:

```cpp
struct BiomtRecord {
  std::vector<std::string> subdivs;
  std::string auth_bio_unit;
  std::string soft_bio_unit;
  std::vector<std::array<std::array<calc_type, 3>, 3>> rot;
  std::vector<std::array<calc_type, 3>> trans;
};
```

Metadata container candidate:

```cpp
using BiomtMap = std::map<int, BiomtRecord>;
```

Metadata parse helper candidate:

```cpp
[[nodiscard]] BiomtMap parse_biomt_header_records(
    const std::vector<std::string>& header_lines);
```

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

- API names are review targets, not finalized ABI.
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

## Explicit Non-Goals In This Step

- no BIOMT implementation in `src/subset.cpp`
- no BIOMT metadata implementation in `src/file_io.cpp` or parser internals yet
- no parser/selection alias expansion
- no topology parser or PDB-driven force-field inference
