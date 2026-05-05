# Python Extension Migration README

This standalone C++ library should promote the useful behavior of the old
Python extension modules into named C++ API homes. It should not expose the old
Python C-extension ABI, NumPy object handling, or hidden implementation names as
the public C++ contract.

## Mapping

| Old Python-side extension area | Standalone C++ home | Status |
| --- | --- | --- |
| `_dcdio` / DCD extension code | `sasmol/file_io.hpp`, `DcdReader`, `DcdWriter` | Promoted as explicit reader/writer APIs with status returns and streaming helpers. |
| `matrix_math` | `sasmol/linear_algebra.hpp` | Promoted as small vector, angle, dihedral, and matrix helper APIs. |
| `mask` dihedral helpers | `sasmol/subset.hpp`, `get_dihedral_subset_mask` | Promoted as a subset helper because it returns atom masks. |
| `overlap` | `sasmol/overlap.hpp`, `has_overlap` | Promoted as a named geometry helper over raw coordinates or molecule frames. |
| VMD/view helpers | `sasmol/view.hpp` | Optional adapter. The C++ surface mirrors Python's coordinate extraction and can compile the promoted legacy IMD/VMD C sender when `SASMOL_ENABLE_VMD_ADAPTER=ON`. |

## Rules for Future Ports

- Promote behavior to the module where a Python developer would look for it:
  file formats in `file_io`, math in `linear_algebra`, masks/copies in
  `subset`, and geometry collision checks in `overlap`.
- Keep old extension names visible in this README and the parity ledger, not in
  public C++ type names.
- Public APIs should use SasMol-owned types such as `Molecule`, `Vec3`, result
  structs, and standard containers.
- Do not leak Python, NumPy, or extension-module implementation details into the
  standalone C++ public headers.
- Kernel-like code can be optimized internally, but callers should find it
  through the same conceptual module names they know from Python SasMol.

## Why This Shape

The old extension modules were hidden partly because they were implementation
details of the Python package. In standalone C++, they are ordinary C++ code and
should be discoverable. The promotion boundary is the behavior, not the old file
or module name: a Python developer should be able to recognize where the
functionality went without carrying forward the legacy extension layout.
