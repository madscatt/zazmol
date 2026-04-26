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

- `736` tests run
- `77` skipped

The skipped tests are primarily size-gated `large` and `huge` cases controlled
through `SASMOL_LARGETEST` and `SASMOL_HUGETEST`.

#### Ported And Passing

The following active Python areas have been reviewed against legacy behavior,
repaired where needed, and validated through the current suite:

- `system`
- `file_io`
- `calculate`
- `operate`
- `subset`
- `properties`
- `linear_algebra`

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
- historical files under `src/python/broken`
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
