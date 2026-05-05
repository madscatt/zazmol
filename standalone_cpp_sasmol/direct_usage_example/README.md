# Direct Usage Example

This directory shows the smallest standalone use of SasMol: read a PDB file
and write a DCD file.

Two examples are included:

- `simple_pdb_in_dcd_out.cpp` uses the public C++ API directly.
- `simple_pdb_in_dcd_out.c` is a real C program. It calls the small local
  `sasmol_direct_usage.cpp` bridge because the standalone SasMol library does
  not currently publish a broad C ABI.

Build:

```bash
make
```

Run with the default small fixture:

```bash
make run
```

Run just one example:

```bash
make run_cpp
make run_c
```

Run with your own input file:

```bash
make run INPUT=/path/to/input.pdb
```
