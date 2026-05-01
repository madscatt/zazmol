# DCD Python Cross-Reader Validation

This records manual cross-reader validation for C++ generated DCD files. Python
`zazmol` is used as the external reader so this remains separate from normal
C++ CTest runs.

## Commands

```bash
cmake --build standalone_cpp_sasmol/build
standalone_cpp_sasmol/build/dcd_roundtrip_writer \
  src/python/test_sasmol/data/dcd_common/1ATM.dcd \
  /private/tmp/sasmol_cpp_1ATM.dcd
python3 standalone_cpp_sasmol/tools/validate_dcd_with_python.py \
  /private/tmp/sasmol_cpp_1ATM.dcd
```

Repeat for `2AAD.dcd` and `rna-1to10.dcd`.

## Results

Validated on the local development environment:

| Generated file | Python frames | Python atoms | Python coordinate sum |
| --- | ---: | ---: | ---: |
| `/private/tmp/sasmol_cpp_1ATM.dcd` | 2 | 1 | `314.790001` |
| `/private/tmp/sasmol_cpp_2AAD.dcd` | 3 | 15 | `3644.293983` |
| `/private/tmp/sasmol_cpp_rna-1to10.dcd` | 10 | 10632 | `-430804.377611` |
