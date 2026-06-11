# mmCIF Fixture Set

These fixtures were downloaded from the official RCSB coordinate endpoint on
2026-06-11 for planned PDBx/mmCIF input support.

## Fixtures

| File | Role | atom_site rows | Models | Auth chains | SHA256 |
| --- | --- | ---: | ---: | ---: | --- |
| `1CRN.cif` | Small protein baseline | 327 | 1 | 1 | `23787562c427d7c1abe5420e86d5f1d0a6c7007dec1e8ce85645a6d69c32e8ba` |
| `1BNA.cif` | DNA/nucleic-acid baseline | 566 | 1 | 2 | `76d550f32963de46f4f3c63caacb530f91a24f333fbc4312766d240b14e8cf30` |
| `1EHZ.cif` | RNA/nucleic-acid baseline | 1821 | 1 | 1 | `3021dd2b6461f850bb66d6748d3c22a5e1bb6cc891c625c15d793f38215b28ab` |
| `4HHB.cif` | Protein with heme/heterogens | 4779 | 1 | 4 | `6c977e3c48fcae60ef116dffff90ce2ae2dbc987f39c3f930d325361a1916812` |
| `2K39.cif` | Multi-model/NMR stress case | 142796 | 116 | 1 | `a6bf4eed6d432c37089a81aba2ef1ef1bf6500e819ac628c5a5258460e1bb824` |
| `1KP8.cif` | Larger multi-chain legacy-stress counterpart | 57085 | 1 | 14 | `1710ea732055b53e7a8b70d384345382ac5a643232341a11a7a2f66f23e722f7` |

## Source URLs

Each file uses the standard RCSB download URL pattern:

```text
https://files.rcsb.org/download/<PDB_ID>.cif
```

These files are intended for parser and descriptor-parity tests. They are not a
`pdbscan`/`pdbrx` metadata test set; connection, assembly, and dictionary-rich
metadata handling should be covered by a separate mmCIF-aware `pdbscan`/`pdbrx`
plan.
