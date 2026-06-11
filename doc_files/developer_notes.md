### Developer Notes

This page collects project-maintenance notes for the Python 3 `sasmol` port in
the `zazmol` repository.

- [`testing and parity`](testing.md)
- [`python 3 parity ledger`](parity_ledger.md)
- [`dtype contract`](dtype_contract.md)
- [`moltype assignment report note`](moltype_assignment_report_github_note.md)
- [`PDBx/mmCIF support plan`](mmcif_support_plan.md)
- [`original developer notes`](../development_tools/notes.md)

## PDBx/mmCIF Numbering Policy Note

Legacy PDB input deliberately separates depositor numbering from SASMOL's
internal atom numbering. `read_pdb()` stores parsed file atom serials in
`_original_index`, stores parsed residue ids in `_original_resid`, and assigns
`_index` from the read order so files with atom-serial overflow or unusual
numbering remain usable internally. `CONECT` output depends on this distinction:
stored connectivity is keyed by original atom serials and is remapped through
`original_index -> index` when writing.

PDBx/mmCIF support must preserve the same external/internal distinction, but it
cannot assume a single fixed-column PDB number exists. The `_atom_site` category
usually carries both depositor/auth identifiers and dictionary-normalized label
identifiers. The mmCIF reader should define an explicit policy for which fields
populate `original_index`, `original_resid`, `index`, and `resid`, and should
record enough mapping information for later `pdbscan`/`pdbrx` work to relate
auth atom/residue ids, label atom/residue ids, chain/asym ids, insertion codes,
model numbers, and connection records without silently changing SASMOL's
current numbering contract.

Two additional PDBx/mmCIF field-policy caveats should remain visible during
implementation. First, legacy PDB `segname` has no exact `_atom_site`
correspondence in the official PDB-to-PDBx table. Base mmCIF I/O should preserve
SASMOL's current fallback behavior by deriving `segname` from the chosen
chain/asym identifier when no explicit segment identifier is available, rather
than introducing a new segment concept in the base reader. Second, the official
PDB-to-PDBx table marks legacy `CONECT` fields as not applicable; related
semantic connection records live in `_struct_conn`, but interpreting those
records belongs to a later mmCIF-aware `pdbscan`/`pdbrx` plan, not the base
molecule reader.

Standalone C++ Phase 2 follows the same policy. `MmcifReader` populates
`original_index` from `_atom_site.id`, assigns `index` from read order, uses
depositor/auth residue numbering when it fits the current integer descriptors,
and rejects incompatible identifiers with a format error rather than silently
renumbering depositor fields. It preserves full chain/asym id strings inside
the molecule, derives `segname` from that chosen chain/asym id, and leaves
`_struct_conn`/CONECT-equivalent semantics for the separate pdbscan/pdbrx plan.
