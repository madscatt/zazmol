### Python 3 Parity Ledger

This project uses the repository name `zazmol`, but the code package remains
`sasmol`.

This ledger records the current status of Python 3 parity against frozen
read-only `legacy_sasmol`.

#### Scope

This ledger covers the active Python package and its active Python regression
suite.

It does not claim:

- external SASSIE caller audit is complete
- historical/manual scripts are preserved in the active repo
- standalone extension-side script coverage is part of the active suite

#### Current Status

The active Python suite passes from the build tree:

- `830` tests run
- `77` skipped

The skipped tests are size-gated `large` and `huge` cases controlled through
`SASMOL_LARGETEST` and `SASMOL_HUGETEST`.

Current expanded validation status:

- normal suite passes
- `large` suite passes
- `huge` suite passes

Recent targeted parity/coverage additions completed:

- direct tests for `calculate_molecular_formula`
- alias parity test for `calc_minmax_all_steps` vs
  `calculate_minimum_and_maximum_all_steps`
- focused behavior tests for `subset.set_coor_using_mask`
- focused helper tests for `pdb_io`:
  `check_for_all_zero_columns`, `create_conect_pdb_lines`
- utility tests for `parse_fasta` and `check_integrity`

#### Ported And Passing

The following active Python areas have been reviewed against legacy behavior,
repaired where needed, and validated through the current suite:

- `calculate`
- `charmm_topology`
- `dcd_io`
- `file_io`
- `linear_algebra`
- `operate`
- `pdb_io`
- `properties`
- `subset`
- `system`
- `topology`
- `utilities`
- `view`

Notable restored or repaired parity items include:

- import/build compatibility shims
- PMI alignment behavior and tests
- `properties` physical-data methods restored from legacy:
  - `element_scattering_lengths`
  - `nucleotide_scattering_lengths`
  - `dna_scattering_lengths`
  - `rna_scattering_lengths`
  - `protein_scattering_lengths`
  - `van_der_waals_radii`
- compatibility wrappers kept for legacy shorthand names:
  - `element_sl`
  - `nucleotide_sl`
  - `dna_sl`
  - `rna_sl`
  - `protein_sl`
  - `van_der_Waals_radii`
- `subset` helper-level unit coverage restored for:
  - `get_indices_from_mask`
  - `get_dihedral_subset_mask`
  - `get_subset_mask`
- `linear_algebra` helper parity restored for:
  - `vec_sub`
  - `vec_scale`
- legacy-style `filename()` / `setFilename()` restored on `system`
- legacy streaming min/max capability restored in `calculate`:
  - `calculate_minimum_and_maximum_all_steps(...)`
  - `calc_minmax_all_steps(...)`
- large/huge test harness repaired for Python 3 workflow
- unconditional DCD success-path stdout noise removed
- historical DCD progress dots moved behind the Python-side `DEBUG` setting
- Python-side logging defaults changed to quiet mode
- explicit dtype contract added:
  - coordinate storage and transfer buffers use `config.COORD_DTYPE`
  - derived calculations and reductions use `config.CALC_DTYPE`

#### Moved Or Renamed By Design

These legacy behaviors are intentionally preserved under clearer or moved APIs:

- legacy `sasproperties.charmm_names` is now covered through
  `charmm_topology.py`
- legacy `create_fasta` behavior is now in `topology.py`
- current readable `properties` method names expand legacy shorthand
  `_sl` into `*_scattering_lengths`

#### Intentionally Excluded

These legacy areas are not part of the active parity target:

- `SasSys`
- `SasAss`
- `SasSol`
- `SasHybrid`
- removed `manual_tests`
- removed historical `src/python/broken` files
- historical extension-side script tests removed from this repo:
  - `extensions/dcdio/test.py`
  - `extensions/dcdio/test2.py`
  - `extensions/mask/test_map_mask.py`
  - `extensions/old_dcdio/test_dcdio.py`

These remain available in git history and, where relevant, in frozen
`legacy_sasmol`.

#### Reviewed And Not Restored

The following legacy items were reviewed and are not currently judged important
active parity gaps:

- `type()` / `setType()`
- `molcharge()`
- `corr()`

Rationale:

- current code uses `moltype()` / `setMoltype()`
- no meaningful implementation was present for `molcharge()` / `corr()`
- no strong evidence of active caller reliance was found inside this repo

#### Remaining Important Unknowns

No known active, in-repo, non-extension core-package parity gap remains open
after the current normal, large, and huge suite passes.

The largest remaining parity risk is outside this repo:

1. External SASSIE caller audit

We have not yet proven that downstream SASSIE code never depends on:

- legacy-only method names
- legacy logging/output side effects
- legacy return/error behavior not covered by this repo's tests

This is the main remaining item before claiming total ecosystem parity with high
confidence.

2. Standalone extension API hardening

Compiled extensions are exercised indirectly by the active Python suite, which
is the primary validation target.

What is not yet claimed:

- a separate modern automated regression layer for direct extension imports and
  direct low-level extension calls

At present, this is optional hardening, not a blocker for the base Python 3
port.

3. Historical low-value or unused legacy APIs

The old `type()` / `setType()`, `molcharge()`, and `corr()` APIs remain reviewed
and not restored. If a downstream SASSIE audit finds active caller reliance,
that should become a targeted compatibility issue rather than a broad
port-completion blocker.

4. Test-runner suite naming gap

`development_tools/test_runner.py` does not currently provide a `utilities`
suite selector. Utility coverage is validated via direct `unittest` module
execution in the current workflow.

#### Practical Definition Of "Base Port Complete"

For this repo, the base Python 3 port can be considered complete when all of
the following are true:

1. the active Python regression suite passes from the build tree
2. meaningful legacy core-package gaps have been restored or intentionally
   excluded
3. scientific and numerical behavior has not been casually changed
4. intentional differences are written down
5. remaining uncertainty is explicitly limited to external caller audit and
   optional extension hardening


#### October 16, 2026 — issue #106: separate prepared selections

Codex; `dev`, baseline `9c76e3b4064e03a9561d54b625e9578bf6db8473`;
uncommitted implementation and validation on the Mac mini.

Added only `Mask.get_subset_mask_prepared(self, basis_filter)` to the library
source. The old method's signature/body and all other `subset.py` bytes remain
unchanged (verified by removing the additive method and comparing to baseline).
No existing caller, SASSIE file, production installation or default was changed.
Preparation occurs once per invocation inside the first atom's existing exception
boundary. No cross-call cache, molecular result cache, new dependency or error
communication/validation layer was added. Developers retain responsibility for
expression suitability and result handling. Arbitrary Python frame/local
introspection and translation audit events are not an equivalence claim; a custom
old-method override is not automatically used by the new separate method.

Added 11 `unittest` methods under `test_subset`; the existing subset/full runner
automatically discovers them. Tests cover all descriptor selections, exact mask
values/types/dtypes/order and indices, fresh data after descriptor/atom-count
changes, reader fixtures (1CRN, 1KP8, RNA), existing override dispatch and a small
set of existing return/failure boundaries. The initial run exposed a test-fixture
mistake (a nonexistent Molecule_Maker residue-flag setter); corrected only the new
test to mutate its existing flag list. Final checks passed without changing old
tests, fixtures, expected outputs or tolerances.

Validation used `/Users/curtisj/anaconda3/bin/python` and the documented build tree,
not the installed production package and not a source overlay. Commands/results:

- `python -B setup.py build`: succeeded, including the existing native build.
- `SASMOL_LARGETEST=y python -B development_tools/test_runner.py --no-build --suite subset`:
  248 passed, no skips/failures, 377.856 seconds.
- `SASMOL_LARGETEST=n SASMOL_HUGETEST=n python -B development_tools/test_runner.py --no-build`:
  907 run, 830 passed, 77 existing size-gated skips, no failures, 81.842 seconds.
- The 17 remaining existing large cases outside subset were selected from the
  recorded normal-suite skip inventory and run serially with `SASMOL_LARGETEST=y`:
  17 passed, no skips/failures, 171.166 seconds. Thus all 43 large cases were
  exercised, with 873 unique Python cases passed across runs. The 34 huge-trajectory
  cases remain gated and were not run; do not label them passes.
- Existing standalone C++ suite: normal CMake test build and `ctest --test-dir
  build/zazmol106_ctest --output-on-failure`, 17/17 passed.
- Four examples in the new method's docstring passed with `doctest`.
- Usage/testing docs and Sphinx subset source updated; Sphinx HTML built successfully
  with the same two baseline warnings. Only six generated artifacts that differed
  from the baseline generation were updated. API anchor, source link and inventory
  entry verified. Full generated-HTML diff whitespace warnings come from existing
  Sphinx template output; handwritten source/docs pass `git diff --check`.

Unprofiled selector timing: five measured paired repetitions after warm-up,
alternating old/new order, identical inputs and outside-timing mask checks. Budget
was ten minutes; completed in 31.14 seconds. These are selection timings, not an
end-to-end SASSIE speedup or a production acceptance claim.

| Workload | Atoms | Calls/sample | Old median s | Prepared median s | Ratio |
| --- | ---: | ---: | ---: | ---: | ---: |
| synthetic_100 | 100 | 1 | 0.000712 | 0.000096 | 7.44 |
| synthetic_100_repeated | 100 | 100 | 0.079140 | 0.009978 | 7.93 |
| synthetic_10000 | 10000 | 1 | 0.101996 | 0.009556 | 10.67 |
| synthetic_10000_repeated | 10000 | 10 | 0.990663 | 0.103460 | 9.58 |
| reader_1KP8 | 57085 | 1 | 0.480726 | 0.034542 | 13.92 |
| reader_1KP8_repeated | 57085 | 3 | 1.380733 | 0.101980 | 13.54 |

Process peak RSS for the whole timing experiment was 108,920,832 bytes; this is
not a per-method memory comparison. Microbenchmarks are machine/workload specific.
SASSIE-wide profiling and any caller adoption remain separate work.

Read-only final fingerprints matched all 684 recorded production package files.
Machine-local raw evidence (commands, logs, skip/test inventory, timing samples,
source baseline and audit JSON) is recorded at:
`/var/folders/85/03z1k3gj79972t52h_ft4rlc0000gn/T/zazmol106_prepared_vovlt7p0`.
No installation, Git commit/push, release or deployment was performed.


October 16 follow-up before requested GACP: The earlier Release C++ run had NDEBUG enabled, so its assertion-based checks were inactive. Rebuilt the existing C++ suite with `cmake -S standalone_cpp_sasmol -B build/zazmol106_debug -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Debug`, built that tree and ran `ctest --test-dir build/zazmol106_debug --output-on-failure`: all 17 tests passed with assertions enabled (22.44 seconds). Confirmed test compile flags omit NDEBUG. No C++ source/API changes were made. This Debug result supplies the assertion-based regression evidence missing from the earlier Release run.
