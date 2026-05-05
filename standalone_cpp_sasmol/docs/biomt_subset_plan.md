# BIOMT Subset Plan (Fixture-First)

Python `zazmol.subset` BIOMT behavior is a deferred feature in standalone C++.
This document defines the implementation guardrails and the first fixture-backed
contract so BIOMT can be added without guessing.

## Scope For The Next Implementation Slice

- Add BIOMT subset APIs only after fixture semantics are locked.
- Keep this as coordinate-transform behavior only.
- Do not combine BIOMT work with selection grammar expansion, topology parser
  work, or merge policy expansion.

## Guardrails

- Python `zazmol` remains the read-only behavior oracle.
- Public APIs must use SasMol-owned types and remain backend-neutral.
- Failures must be explicit and non-partial.
- No hidden global state or logging side effects.
- No mutation on failed validation.
- No large binary fixtures in git; use generated test inputs.

## Proposed API Shape (For Review Before Implementation)

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

## Semantics To Lock Before Implementation

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

## Explicit Non-Goals In This Step

- no BIOMT implementation in `src/subset.cpp`
- no parser/selection alias expansion
- no topology parser or PDB-driven force-field inference
