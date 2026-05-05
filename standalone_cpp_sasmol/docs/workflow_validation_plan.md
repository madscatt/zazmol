# Workflow Validation Plan

The active next work is validation, not broad new porting.

## Active

1. Selection compatibility validation

   Keep checking surveyed SASSIE/ZAZZIE basis strings against the safe C++
   translator. Add only fixture-backed grammar features. Do not recreate Python
   `eval`.

2. End-to-end molecule workflow smoke tests

   Validate real operation chains such as PDB read -> SASSIE basis selection ->
   subset copy -> DCD write/read, and DCD streaming calculations against loaded
   trajectory calculations.

## Backlog

- Real VMD runtime smoke test

  Build with `SASMOL_ENABLE_VMD_ADAPTER=ON`, open a real VMD IMD listener, and
  send a tiny molecule frame. Normal tests should continue to use the mock
  sender.

- Python binding planning

  Decide which standalone C++ APIs to expose first, how to preserve Python
  names, how NumPy ownership/conversion should work, and how `IoStatus` /
  structured errors should map into Python.

- Policy-gated feature work

  CHARMM36, fixed/free DCD variants, DCD unit-cell writing, true random DCD
  seeking, and `rg/x2` frame-filter helpers should remain separate decisions.
