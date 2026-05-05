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
can store CHARMM-related descriptors and can accept explicit caller-provided
CHARMM type assignments, but it does not infer types and does not parse CHARMM
topology files yet.

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
- no-partial-mutation failure behavior for length mismatches, atom-name
  mismatches, and molecule name-vector mismatches

This is enough for workflows that already have trustworthy CHARMM type data.
It is not a topology engine. A full CHARMM topology parser would be a separate
design step because Python `CharmmTopology` includes residue parsing, patches,
completeness checks, and possible atom reordering.

Recommended next step: pause topology feature expansion unless a caller needs a
minimal parser for a reviewed subset of CHARMM topology records. Do not start
parser work without a separate reviewed plan and fixtures.

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

   Port Python `CharmmTopology` behavior as its own module:

   - parse `MASS`, `DECL`, `DEFA`, `AUTO`, `RESI`, `PRES`, `ATOM`, `BOND`,
     `DOUBLE`, `IMPR`, `CMAP`, `DONO`, `ACCE`, `IC`, and `DELE`
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
