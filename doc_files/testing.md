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

#### Dtype Contract

Coordinate storage and coordinate-transfer buffers use `config.COORD_DTYPE`
(`numpy.float32`). Derived calculations and reductions use `config.CALC_DTYPE`
(`numpy.float64`).

See the project dtype contract for the full policy:

- [`dtype contract`](dtype_contract.md)

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
