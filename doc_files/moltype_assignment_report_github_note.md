### Moltype Assignment Report Note

#### Prior Analysis

I took the deep look read-only. The key correction is: `read_pdb()` does **not** currently assign generic `nucleic` during normal PDB read. It assigns atom-level `moltype` from residue-name tables in [pdb_io.py](/Users/curtisj/git_working_copies/zazmol/src/python/pdb_io.py:704), using this order: `protein`, then `rna`, then `dna`, then `water`, else `other`.

That order matters. Shared residue names like `ADE`, `GUA`, and `CYT` are classified as `rna` because RNA is checked before DNA. DNA-specific names like `DA`, `DG`, `DC`, `DT`, and `THY` become `dna`. The `nucleic_resnames` list exists in [pdb_io.py](/Users/curtisj/git_working_copies/zazmol/src/python/pdb_io.py:327), but it is not used for the actual `read_pdb()` moltype assignment.

The SASSIE implication is real: several workflows branch on exact strings `protein`, `rna`, `dna`, and sometimes `nucleic`. Some already compensate locally. For example, `contrast_calculator_filter.py` has explicit logic for ambiguous `ADE/CYT/GUA` cases and comments that DNA can be classified as RNA. `torsion_angle_md.py` also has a `fix_moltypes()` path where user-provided DNA segnames override the read-in moltype. So changing `read_pdb()` classification itself would be risky.

My recommendation still holds: add a safe Option B report later, without changing `read_pdb()` behavior or return values. Something like `moltype_by_segname_report()` / `validate_moltype_by_segname()` could return a plain dict saying, per segment: assigned moltypes, residue-name evidence, ambiguous DNA/RNA evidence, mixed/other status, and human-readable messages. No logging, no mutation, no altered caller behavior. That gives SASSIE and users better visibility without disturbing decades of hard-earned PDB behavior.

For C++ parity, I’d mirror the current Python behavior first, then add the same non-mutating report API. The report can be smarter than assignment; assignment stays conservative and compatible.

#### Implemented

Added the safe Option B report in both Python SasMol and standalone C++ SasMol.

Python:

- Added `moltype_by_segname_report()` to the existing PDB/file-I/O mixin.
- The method returns a plain dictionary with `overall_status`, `segments`, and `messages`.
- It does not modify `moltype`, `resname`, `segname`, or any other molecule descriptor.
- It does not change `read_pdb()` behavior or return values.
- It reports segment-level assigned moltypes, residue names, residue counts, atom counts, ambiguous DNA/RNA residue names, DNA-specific residue evidence, RNA-specific residue evidence, RNA atom-name evidence, and explanatory messages.

C++:

- Added `sasmol::MoltypeSegmentReport` and `sasmol::MoltypeReport`.
- Added `sasmol::Molecule::moltype_by_segname_report() const`.
- The C++ report mirrors the Python status vocabulary and evidence fields.
- The C++ method is const and does not mutate molecule state.

Status vocabulary:

- Segment status: `clean`, `mixed`, `ambiguous_nucleic`, `all_other`
- Overall status: `clean`, `mixed_by_segname`, `ambiguous_nucleic`, `unknown`

The report intentionally detects ambiguity without reassigning moltypes. For example, a segment containing only `ADE`, `CYT`, and `GUA` can be reported as `ambiguous_nucleic` even though the existing `read_pdb()` assignment remains `rna` because of the historical residue-table order.

#### Regression Coverage

Added Python file-I/O unit tests covering:

- ambiguous DNA/RNA overlap residue names
- no mutation of assigned `moltype`
- DNA-specific residue evidence
- RNA atom-name evidence such as `O2'`
- mixed moltype segments

Added C++ molecule-core tests covering the same report behavior.

#### Verification

Python:

```text
python development_tools/test_runner.py --suite file_io
Ran 182 tests in 30.043s
OK (skipped=34)

python development_tools/test_runner.py
Ran 848 tests in 48.105s
OK (skipped=77)
```

C++:

```text
cmake -S standalone_cpp_sasmol -B standalone_cpp_sasmol/build
cmake --build standalone_cpp_sasmol/build
ctest --test-dir standalone_cpp_sasmol/build --output-on-failure
100% tests passed, 0 tests failed out of 9
```
