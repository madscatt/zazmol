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
ordered `ATOM`, `BOND`, `DOUBLE`, `ANGL`, `THET`, `DIHE`, `IMPR`, and `CMAP`
rows, plus `DONO`, `ACCE`, `IC`, and `DELE` rows, as data-only records. It does
not infer types, apply patches, or mutate molecules from topology files yet.

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
- `compare_list_ignore_order(...)` matching Python's helper behavior used by
  `check_charmm_atomic_order_reorganize`
- `setup_charmm_residue_atoms(...)` as a pure data helper that builds
  residue/patch names to ordered atom-name lists from parsed topology entries
- `setup_cys_patch_atoms_simple(...)` as a pure helper for the Python `DISU`
  convention, removing `HG1` from the CYS atom list
- `patch_charmm_residue_atoms(...)` as a pure helper for Python's atom-only
  residue patch behavior, returning a patched topology entry and atom-name list
- `choose_charmm_residue_atom_order(...)` as a pure helper that extracts the
  per-residue topology-order decision from Python's mutating
  `check_charmm_atomic_order_reorganize`
- `parse_charmm_topology_globals(...)` for Python-matched `MASS`, `DECL`,
  `DEFA`, and `AUTO` records, preserving values as strings
- `parse_charmm_topology(...)` for those global records plus Python-matched
  `RESI`/`PRES` headers and ordered `ATOM`, `BOND`, `DOUBLE`, `ANGL`, `THET`,
  `DIHE`, `IMPR`, `CMAP`, `DONO`, `ACCE`, `IC`, and `DELE` rows, preserving
  total charge and atom charges as strings
- no-partial-mutation failure behavior for length mismatches, atom-name
  mismatches, and molecule name-vector mismatches

This is enough for workflows that already have trustworthy CHARMM type data.
It is not a topology engine. The parser slices are data-only and do not assign
descriptors to a molecule. Patch application, completeness checks, and possible
atom reordering remain separate reviewed steps.

Recommended next step: pause before moving from passive parsing to any topology
operation. Do not apply patches, build residue atom lists, or reorder atoms
without a separate fixture-backed design slice.

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

   Fifth slice implemented:

   - parses `DIHE` as ordered four-token string records
   - parses `IMPR` as ordered four-token string records
   - parses `CMAP` as ordered four-token string records, matching Python's
     current four-token chunk behavior
   - stops four-token parsing at inline comments beginning with `!`

   Sixth slice implemented:

   - parses `DONO` as ordered single-token string records
   - parses `ACCE` as ordered single-token string records
   - stops single-token parsing at inline comments beginning with `!`

   Seventh slice implemented:

   - parses `IC` as ordered string vectors preserving Python's `words[1:10]`
     shape
   - does not validate internal-coordinate geometry or require a full nine-token
     record, matching Python's passive storage behavior

   Eighth slice implemented:

   - parses `DELE ATOM` as ordered atom-name strings, preserving Python's
     first-token-only behavior after `DELE ATOM`
   - parses `DELE ANGL` as ordered string vectors, preserving Python's
     `words[2:]` shape
   - ignores unknown `DELE` types, matching Python's current parser behavior

   Ninth slice implemented:

   - ports Python's `compare_list_ignore_order` helper as
     `compare_list_ignore_order(...)`
   - preserves the Python helper behavior exactly; this is less strict than
     `validate_charmm_residue_atoms(...)` and exists for future reorder parity

   Tenth slice implemented:

   - ports Python's `setup_charmm_residue_atoms` behavior as a pure
     `setup_charmm_residue_atoms(...)` data helper
   - ports Python's `setup_cys_patch_atoms_simple` behavior as
     `setup_cys_patch_atoms_simple(...)`
   - adds `DISU` from the CYS atom list when CYS is present
   - avoids Python's hard failure on tiny topology fixtures without CYS by
     omitting `DISU` when no CYS atom list exists

   Eleventh slice implemented:

   - ports Python's `patch_charmm_residue_atoms` atom behavior as pure
     `patch_charmm_residue_atoms(...)`
   - applies `DELE ATOM` records from the patch to the residue atom list
   - replaces same-named residue atoms before adding patch atoms
   - preserves Python's patch placement policy: `NTER`, `GLYP`, and `PROP`
     insert patch atoms at the front; `CTER` appends patch atoms
   - still does not patch bonds, angles, or any other topology records, matching
     Python's documented limitation

   Twelfth slice implemented:

   - extracts Python's per-residue order decision into pure
     `choose_charmm_residue_atom_order(...)`
   - uses residue atom lists plus pure patch application to choose the topology
     atom-name order
   - preserves Python's CYS-to-DISU and HIS-to-HSE/HSD/HSP fallback behavior
   - preserves Python's terminal patch selection for `NTER`, `GLYP`, `PROP`,
     and `CTER`
   - still does not mutate molecule descriptors or coordinates

   Future slices should port Python `CharmmTopology` behavior as its own module:

   - integrate this order choice into a no-partial-mutation molecule reorder
     workflow
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
