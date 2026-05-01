# Standalone C++ SasMol

This directory starts a new standalone C++20 SasMol library. It is intentionally
separate from the existing Python `zazmol` package and from the legacy
`src/cpp/sasmol` tree. Python `zazmol` remains the read-only behavioral
standard, fixture source, and naming reference.

## Design Goals

- Preserve familiar SasMol concepts and names where they make user code easier
  to recognize: `Molecule`, `read_pdb`, `write_pdb`, `read_dcd`,
  `get_subset_mask`, `calculate_center_of_mass`, and related snake_case APIs.
- Use modern, readable C++20 without forcing conventional C++ naming where it
  would make parity harder to see.
- Keep public APIs owned by this library. Optional implementation dependencies
  such as Eigen can be evaluated later without leaking into user-facing types.
- Keep LASSIE `SasCalc`, Debye, GV, CUDA, and scattering code separate. Those
  algorithms may later consume this library but are not part of this core.
- Keep coordinate storage explicit and contiguous so future CPU, GPU, MPI, or
  Python-binding layers can consume predictable views without changing
  `sasmol::Molecule`.

## Current Scope

The first milestone is the core data model and parity machinery:

- `sasmol::Molecule`
- Python-like descriptor vectors exposed through named accessors
- contiguous frame-major coordinate storage
- storage-friendly `coord_type` and calculation-friendly `calc_type`
- integrity reporting that names descriptor length mismatches
- file I/O contract types for PDB/DCD readers and writers
- a read-only parity ledger generator for surveying Python `zazmol`

This is not yet a feature-complete port of Python SasMol.

The PDB/DCD APIs are defined in
[`docs/file_io_contract.md`](docs/file_io_contract.md), but parser internals are
intentionally deferred until fixture parity is reviewed.

## Build

```bash
cmake -S standalone_cpp_sasmol -B standalone_cpp_sasmol/build
cmake --build standalone_cpp_sasmol/build
ctest --test-dir standalone_cpp_sasmol/build --output-on-failure
```

## Parity Survey

From the repository root:

```bash
python3 standalone_cpp_sasmol/tools/generate_parity_ledger.py --python-src src/python
```

The script reads Python source and writes a Markdown survey to stdout. It does
not modify Python `zazmol`.
