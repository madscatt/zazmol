## Scientific Computing Guidance

- Do not introduce shortcuts, parameter reductions, approximations, or workflow simplifications solely to save time. Preserve the requested physics, parameter space, and analysis unless explicitly directed otherwise.

- When a shortcut or optimization is possible, propose it separately and explain the expected impact. Do not apply it without approval.

- Prefer correctness, reproducibility, and scientific validity over runtime reduction.

- For long-running jobs, implement a periodic heartbeat ("pulse") indicating the process is still active. Include elapsed time and useful progress metrics when available.

- Clearly identify all assumptions. Never silently replace requested behavior with a simpler alternative.

- Algorithmic and computational improvements are encouraged, but should be presented as recommendations for review before adoption.

## Installed Version Tracking

Installed versions come from build receipts carried with verified artifacts,
not from the current checkout HEAD. Preserve independent release identities;
unknown installations stay unrecorded and partial differing builds stay mixed.
Git hooks must not publish installed revisions. Read `README.md` and the
paired application/framework policies before changing this boundary. Deployment
and live legacy-hook cleanup require their own task scope; local build checks
do not establish deployed acceptance.
