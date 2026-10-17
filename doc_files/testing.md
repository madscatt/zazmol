### Testing And Parity

This project uses the repository name `zazmol`, but the code package remains `sasmol`.
The active Python regression suite lives under `src/python/test_sasmol`.

#### Standard Test Workflow

Run the full Python suite from the repository root:

```bash
python development_tools/test_runner.py
```

Run one supported slice:

```bash
python development_tools/test_runner.py --suite subset
python development_tools/test_runner.py --suite properties
python development_tools/test_runner.py --suite system
```

List the supported suite names:

```bash
python development_tools/test_runner.py --list
```

The test runner builds `sasmol` first and runs tests from `build/lib.*` so the results reflect the current workspace code instead of a previously installed copy in `site-packages`.

#### What Counts As Parity

For this Python 3 port, parity means:

1. Public, caller-relied behavior from legacy `sasmol` is available in `sasmol`, either under the same name or under an intentional moved API that is documented and tested.
2. Scientific behavior is preserved unless a change is explicitly approved and documented.
3. Existing validated sample data, reference outputs, and test tolerances are preserved.
4. The active Python regression suite passes from the build tree.
5. Intentional exclusions and differences are written down instead of being rediscovered later.

#### Hard Compatibility Rules

Do not change casually:

- numerical precision
- numerical tolerances
- algorithms
- dtype behavior
- array shape or ordering
- atom ordering
- residue numbering
- segment names
- chain names
- frame indexing
- coordinate units or layout
- file-format numeric formatting
- write-read round-trip behavior
- public return types
- caller-relied error behavior

When a numerical or scientific test fails, compare to legacy behavior where possible and report the cause before changing code or tests.

#### PDB Reader Compatibility

The PDB reader is intentionally tolerant because it supports long-lived SASSIE
workflows beyond canonical protein-only PDB files, including PDB scan, PDB
repair/build workflows, simulation setup, trajectory extraction, and merge
utilities.

Do not make the reader strictly PDB-spec compliant as a cleanup change. Preserve
the tested behavior for blank trailing lines, single-frame files without
``END``, multi-frame files separated by ``END`` or ``MODEL``/``ENDMDL``,
non-CHARMM atom names, HETATM records, and non-protein systems unless a fixture
and legacy comparison show that a behavior change is intentional.

#### Dtype Contract

Coordinate storage and coordinate-transfer buffers use `config.COORD_DTYPE`
(`numpy.float32`). Derived calculations and reductions use `config.CALC_DTYPE`
(`numpy.float64`).

See the project dtype contract for the full policy:

- [`dtype contract`](dtype_contract.md)

#### Recent Coverage Additions

Focused Python 3 parity coverage was added for stable public behavior in:

- `calculate.calculate_molecular_formula`
- `calculate.calc_minmax_all_steps` alias parity against
  `calculate_minimum_and_maximum_all_steps`
- `subset.set_coor_using_mask`
- `pdb_io` helper behavior:
  `check_for_all_zero_columns`, `create_conect_pdb_lines`
- `utilities.parse_fasta`
- `utilities.check_integrity`

These additions intentionally avoid questionable legacy aliases and
experimental internals so tests improve confidence without freezing uncertain
behavior.

#### Runner Note

`development_tools/test_runner.py` does not currently expose a `utilities`
suite selector. Run utility tests directly with `unittest` when needed.

#### Current Intentional Exclusions

These legacy areas are not part of the active parity target right now:

- `SasSys`
- `SasAss`
- `SasSol`
- `SasHybrid`
- removed `manual_tests`
- removed historical `src/python/broken` files

#### Extension Tests

Extension tests are a separate workstream and are not part of the main Python suite yet.
Treat them per extension:

1. decide whether the extension is still supported
2. classify its tests as active, port-needed, or historical/demo
3. only then bring supported extension tests into a repeatable workflow

#### Current Ledger

See the current parity status ledger here:

- [`python 3 parity ledger`](parity_ledger.md)


#### Prepared selection coverage

`test_subset/test_unit_subset_Mask_get_subset_mask_prepared.py` is a `unittest`
module discovered automatically by both `--suite subset` and the default full
suite. No separate runner or production installation is needed:

```bash
~/anaconda3/bin/python -B development_tools/test_runner.py --suite subset
~/anaconda3/bin/python -B development_tools/test_runner.py
```

It compares the separate prepared method against unchanged `get_subset_mask()`:
all selection descriptors, atom order and index conversion, mask types/dtypes,
current values after mutation/atom-count changes, reader-produced molecules,
existing subclass dispatch and a small set of existing return/failure boundaries.
These are focused compatibility checks, not an expanded error-reporting contract.
Developers remain responsible for expression suitability and handling results.
Timing measurements belong outside the tests; they must verify equivalent results
and distinguish selection-call improvements from whole-workflow improvements.
