# CHARMM36 Issue Draft

Title:

Update Python and C++ SasMol topology handling for CHARMM36 compatibility

Body:

The current SasMol topology behavior is based on the historical Python
`zazmol.charmm_topology.CharmmTopology` implementation and the CHARMM27-style
topology files used by older SASSIE workflows. The standalone C++ port now
preserves the reviewed Python behavior for the minimal topology path:

- global and residue topology record parsing
- residue atom-list construction
- atom-only patch application
- residue atom-order choice
- residue and molecule reorder planning
- reordered molecule copy
- explicit in-place reorder after successful validation

This issue tracks a separate upgrade path for CHARMM36. It should be handled in
both the Python and C++ trees together so the C++ port does not drift away from
the Python behavior oracle.

Guardrails:

- Do not infer CHARMM atom types from ordinary PDB atom names.
- Do not make `read_pdb()` assign force-field types or reorder atoms as a side
  effect.
- Do not silently guess unknown residues, atom aliases, or patch choices.
- Preserve current Python behavior unless a test-backed Python and C++ change is
  made together.
- Keep mutation behind explicit APIs and preserve no-partial-mutation failure
  behavior in C++.
- Validate parser and patch changes against small fixtures before trying full
  production topology files.

Open design questions:

- Which CHARMM36 topology and parameter files should be the reference fixtures?
- Which existing Python `CharmmTopology` limitations are intentional compatibility
  behavior, and which should be fixed before or during CHARMM36 support?
- Are new patch records, residue aliases, terminal patches, or nucleic-acid
  conventions needed beyond the CHARMM27 behavior?
- Should bond/angle/dihedral/IC patching remain passive data storage, or should
  Python and C++ gain a real topology graph update path?
- What caller workflows need atom reordering, and should they use Python's
  historical in-place path or the safer C++ plan/copy/in-place split?

Suggested validation plan:

1. Add tiny Python oracle fixtures for one CHARMM36 residue with no patch.
2. Add fixtures for N-terminal and C-terminal patch behavior.
3. Add fixtures for one nucleic-acid residue if CHARMM36 naming differs from
   current assumptions.
4. Add fixtures for any changed HIS/CYS/disulfide behavior.
5. Compare Python and C++ normalized parse output before adding mutation.
6. Compare Python and C++ chosen atom order before adding molecule reorder.
7. Only after small fixtures pass, summarize full CHARMM36 topology files and
   compare high-level counts and representative residue records.

Expected outcome:

Python and C++ SasMol should share a documented CHARMM36-compatible topology
contract, with fixtures proving the parser, patch, order-choice, and explicit
reorder behavior. The upgrade should not narrow legacy CHARMM27 behavior unless
the change is deliberate, documented, and covered by tests.
