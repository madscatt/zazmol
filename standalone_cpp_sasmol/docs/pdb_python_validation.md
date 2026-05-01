# PDB Python Cross-Reader Validation

This records manual cross-reader validation for C++ generated PDB files. Python
`zazmol` is used as the external reader so this remains separate from normal
C++ CTest runs.

## Commands

```bash
cmake --build standalone_cpp_sasmol/build
standalone_cpp_sasmol/build/pdb_roundtrip_writer \
  src/python/test_sasmol/data/pdb_common/1ATM.pdb \
  /private/tmp/sasmol_cpp_1ATM.pdb
python3 standalone_cpp_sasmol/tools/validate_pdb_with_python.py \
  /private/tmp/sasmol_cpp_1ATM.pdb
```

Repeat for `2AAD.pdb` and `1ATM-1to2.pdb`.

## Results

Validated on the local development environment:

| Generated file | Python frames | Python atoms | Python coordinate sum |
| --- | ---: | ---: | ---: |
| `/private/tmp/sasmol_cpp_1ATM.pdb` | 1 | 1 | `157.395000` |
| `/private/tmp/sasmol_cpp_2AAD.pdb` | 1 | 15 | `2407.675995` |
| `/private/tmp/sasmol_cpp_1ATM-1to2.pdb` | 2 | 1 | `314.790001` |
