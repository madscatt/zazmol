# CHARMM Topology Plan

Python `zazmol.charmm_topology` is the behavioral reference for future C++
topology work. It is not just a source of atom type names; it parses topology
records, builds residue atom lists, applies patches, checks atom completeness,
and can reorganize atoms to match CHARMM topology order.

## Current C++ Boundary

The standalone C++ `Molecule` already has:

- explicit `charmm_type()` storage
- explicit `atom_charge()` storage
- explicit `residue_flag()` storage
- typed extension descriptor maps for one-off data

PDB reading intentionally does not populate `charmm_type()`. Ordinary
user-supplied PDB files often do not follow CHARMM atom naming closely enough
for safe automatic assignment.

## Current Parity Checkpoint

Topology support is intentionally narrow and safe at this stage. The C++ core
can store CHARMM-related descriptors, can accept explicit caller-provided
CHARMM type assignments, and can parse the first reviewed global-record subset
from CHARMM topology files. It can also parse `RESI`/`PRES` headers and their
ordered `ATOM`, `BOND`, `DOUBLE`, `ANGL`, and `THET` rows as data-only records.
It does not infer types, apply patches, parse dihedral/improper records, or
mutate molecules from topology files yet.

Implemented:

- `Molecule::charmm_type()` as an optional atom-aligned descriptor
- `assign_charmm_types(molecule, types)` for already-known atom-aligned type
  lists
- `assign_atom_charges(molecule, charges)` for already-known atom-aligned charge
  values
- `assign_charmm_types_from_atom_table(molecule, assignments)` for explicit
  atom-ordered `(atom name, CHARMM type)` rows, such as data parsed from a PSF
  or caller table
- `assign_charmm_types_and_atom_charges_from_atom_table(molecule, assignments)`
  for explicit atom-ordered `(atom name, CHARMM type, atom charge)` rows
- hand-built `CharmmResidueDefinition` records plus
  `validate_charmm_residue_atoms(...)` for exact-match, missing atom, extra atom,
  and duplicate atom reports
- `parse_charmm_topology_globals(...)` for Python-matched `MASS`, `DECL`,
  `DEFA`, and `AUTO` records, preserving values as strings
- `parse_charmm_topology(...)` for those global records plus Python-matched
  `RESI`/`PRES` headers and ordered `ATOM`, `BOND`, `DOUBLE`, `ANGL`, and
  `THET` rows, preserving total charge and atom charges as strings
- no-partial-mutation failure behavior for length mismatches, atom-name
  mismatches, and molecule name-vector mismatches

This is enough for workflows that already have trustworthy CHARMM type data.
It is not a topology engine. The parser slices are data-only and do not assign
descriptors to a molecule. Dihedral/improper parsing, patch application,
completeness checks, and possible atom reordering remain separate reviewed
steps.

Recommended next step: validate the next parser slice against tiny Python-oracle
fixtures before any production topology summary work. Do not parse dihedral,
improper, patch-delete, or reorder behavior without a separate fixture-backed
slice.

The Python oracle harness for future parser work is recorded in
`docs/charmm_topology_python_oracle.md`.

## Safety Rule

Never infer CHARMM atom types from PDB atom names alone.

Typing must come from an explicit topology source, caller-provided mapping, or
future topology subsystem. Ambiguous or unmatched atoms must be reported rather
than guessed.

## Recommended Implementation Path

1. **Explicit Assignment Helper** implemented

   Add a small helper that takes caller-provided atom-aligned values and assigns
   `charmm_type()` or `atom_charge()` after validating length. This supports
   workflows that already know the force-field descriptors without adding
   topology parsing yet.

   Current helpers:

   - `assign_charmm_types(molecule, types)`
   - `assign_atom_charges(molecule, charges)`

2. **PSF/Caller Atom Table Helper**

   Implemented helpers that mirror the SASSIE TAMC usage: consume atom-ordered
   table values from an explicit force-field source such as a PSF, validate them
   against `Molecule::name()`, and assign descriptors only if every atom
   matches.

   Current helpers:

   - `assign_charmm_types_from_atom_table(molecule, assignments)`
   - `assign_charmm_types_and_atom_charges_from_atom_table(molecule, assignments)`

   This is stricter than the historical Python caller: no partial mutation
   after a mismatch.

3. **Residue Definition Model** implemented

   Added plain C++ data structures for hand-built residue definitions:

   - `CharmmAtomDefinition`
   - `CharmmResidueDefinition`
   - `CharmmResidueValidation`
   - `validate_charmm_residue_atoms(molecule_atom_names, residue)`

   This validates residue atom-name sets without parsing topology files, applying
   patches, mutating molecules, or reordering atoms.

4. **Topology Parser Subsystem**

   First slice implemented:

   - `parse_charmm_topology_globals(...)`
   - parses only `MASS`, `DECL`, `DEFA`, and `AUTO`
   - stores Python-equivalent string tokens instead of numeric-coercing masses
   - reports malformed global records through `errors`

   Second slice implemented:

   - `parse_charmm_topology(...)`
   - parses the same global records plus `RESI`, `PRES`, and `ATOM`
   - stores residue and patch total charges as strings
   - stores `ATOM` records as ordered `(atom name, CHARMM type, charge)` string
     triples
   - deliberately ignored `BOND`, `DELE`, and other section records in that
     slice

   Third slice implemented:

   - parses `BOND` as ordered two-token string pairs
   - parses `DOUBLE` under Python's `DOUB` behavior as ordered two-token string
     pairs
   - stops pair parsing at inline comments beginning with `!`

   Fourth slice implemented:

   - parses `ANGL` as ordered three-token string triples
   - parses `THET` as ordered three-token string triples
   - stops triple parsing at inline comments beginning with `!`

   Future slices should port Python `CharmmTopology` behavior as its own module:

   - parse `DIHE`, `IMPR`, `CMAP`, `DONO`, `ACCE`, `IC`, and `DELE`
   - build residue atom lists
   - support reviewed residue patches such as `NTER`, `CTER`, `GLYP`, `PROP`,
     and disulfide/HIS variants
   - validate missing, extra, duplicate, and ambiguous atoms
   - optionally reorder molecules only through an explicit API

## Non-Goals For Now

- no hidden CHARMM typing during PDB read
- no simple `(resname, atom name)` type inference as if it were Python parity
- no best-effort guessing for unknown residues
- no automatic atom reordering as a side effect of reading PDB
- no partial mutation when topology validation fails

## Test Requirements

Future topology work should use small fixtures first:

- one residue with exact atom-name match
- one residue with missing atom
- one residue with extra atom
- terminal patch cases
- HIS/CYS special cases only after the simple cases pass
- atom reordering only through an explicitly named API
