# CHARMM Topology Python Oracle

This records the validation harness for future standalone C++ CHARMM topology
parser work.

Python `zazmol.charmm_topology.CharmmTopology` is the behavior oracle. Before
adding any C++ parser, tiny topology fixtures are parsed by Python and normalized
to JSON with:

```bash
/Users/curtisj/anaconda3/bin/python \
  standalone_cpp_sasmol/tools/validate_charmm_topology_with_python.py \
  standalone_cpp_sasmol/tests/data/topology/minimal_resi_atoms.rtf
```

The tool deliberately preserves Python's parsed shape:

- charges remain strings
- residue and patch records are dictionary entries keyed by residue or patch name
- `DELE` is nested by delete property such as `ATOM` or `ANGL`
- atom order and bond pair order are preserved
- dictionary keys are sorted only for stable JSON output

## Characterization Fixtures

- `minimal_mass_only.rtf`
- `minimal_resi_atoms.rtf`
- `minimal_resi_atoms_bonds.rtf`
- `minimal_pres_atoms_dele.rtf`
- `minimal_comments_blank_lines.rtf`
- `minimal_multiple_residues.rtf`

All six fixtures currently parse through Python with `errors == []`.

## Observed Python Shapes

Global records:

- `MASS` is stored as a list of three-token lists: index, atom type, mass.
- `DECL` is stored as a list of strings.
- `DEFA` is stored as a list of strings.
- `AUTO` is stored as a list of strings.

Residue and patch records:

- `RESI GLY 0.00` creates a top-level `GLY` dictionary.
- `PRES NTER 1.00` creates a top-level `NTER` dictionary.
- `TOTAL_CHARGE` is stored as a string.
- `ATOM` records are stored as ordered three-token lists: atom name, CHARMM type,
  charge string.
- `BOND` records are stored as ordered two-token lists.
- `DOUBLE` records are stored under the key `DOUB`.
- `DELE ATOM HN` is stored under `DELE.ATOM` as the string `HN`.
- `DELE ANGL HT1 N CA` is stored under `DELE.ANGL` as a token list.

Comment behavior:

- blank lines are ignored.
- full-line comments beginning with `!` are ignored.
- inline comments in pair-token records such as `BOND` stop pair parsing when the
  parser reaches the `!` token.
- inline comments after `ATOM` records do not affect the stored atom record
  because Python keeps only `words[1:4]`.

## Parser Guardrail

Do not implement the C++ parser from the production topology file first. Each
new C++ parser slice must compare against this Python oracle on small fixtures
before production topology summaries are considered.
