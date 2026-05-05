# Development Progress Notes

This file records design follow-ups that should not be lost while the
standalone C++ SasMol parity work proceeds in small implementation slices.

## Selection / Eval Scope

The first C++ selection parser intentionally covers only a bounded subset of
Python `get_subset_mask` expressions. A future pass should survey real Python
SasMol/SASSIE selection usage and determine the maximum practical scope that
can be supported safely.

That future expansion must keep the current failure contract:

- unsupported or ambiguous expressions fail with a clear error
- failures return no partial indices
- the parser must stop rather than guessing
- arbitrary Python-style `eval` behavior must not become the default C++ path

## CHARMM Topology And Type Assignment

`Molecule::charmm_type()` is a first-class descriptor, but it must not be
auto-guessed from ordinary PDB input. Many user-supplied PDB files do not follow
CHARMM atom naming closely enough for automatic assignment to be safe. A wrong
force-field type can be worse than a missing one because downstream calculations
may treat it as authoritative.

CHARMM typing should therefore stay explicit:

- PDB read leaves `charmm_type()` empty
- topology/CHARMM workflows may populate it only after they have a real topology
  source or caller-provided atom table
- SASSIE currently assigns CHARMM types from a PSF atom table, then uses those
  atom-aligned values for torsion parameter matching
- the C++ helper mirrors that PSF/caller-table workflow by validating atom names
  and refusing partial mutation on mismatch
- helper APIs should report unmatched or ambiguous atoms instead of guessing

The standalone C++ topology port now includes the minimal Python-parity path:

- reviewed CHARMM topology record parsing with Python oracle fixtures
- residue atom-list construction
- atom-only residue patch application
- residue atom order choice, including Python's terminal patch, CYS/DISU, and
  HIS variant behavior
- pure residue and whole-molecule reorder planning
- pure reordered molecule copy
- explicit in-place reorder wrapper that swaps only after planning and copying
  succeed

This is still not a general topology engine. Python currently has important
limits, and the C++ port should not silently improve them in ways that split the
two trees. CHARMM36 support needs a joint Python/C++ design pass.

The staged topology direction is captured in `docs/charmm_topology_plan.md`.
The follow-up issue draft is captured in `docs/charmm36_issue_draft.md`.

## Moltype Assignment Reporting

Python `read_pdb()` historically assigns atom-level `moltype` from residue-name
tables in the order `protein`, `rna`, `dna`, `water`, then `other`. Shared
nucleic residue names such as `ADE`, `GUA`, and `CYT` therefore remain assigned
as `rna` by default. That behavior is intentionally unchanged because SASSIE
and long-lived workflows branch on the exact existing strings and sometimes
apply their own higher-level corrections.

The safer parity direction is a non-mutating report API. Python now exposes
`moltype_by_segname_report()` on the PDB/file-I/O mixin, and standalone C++
exposes `sasmol::Molecule::moltype_by_segname_report() const`. Both report
mixed moltype segments, DNA/RNA-overlap residue names, and DNA/RNA-specific
evidence without logging, relabeling, or changing `read_pdb()` behavior.

## Python Extension Promotion

The old Python extension modules should be treated as behavior sources, not as
public C++ naming or ABI templates. DCD behavior belongs in `file_io`,
matrix/vector helpers belong in `linear_algebra`, dihedral mask behavior belongs
in `subset`, and overlap checks belong in a named `overlap` module. VMD/view
integration remains an optional adapter outside the portable core.

The current mapping for Python developers is documented in
`docs/python_extension_migration_readme.md`.
