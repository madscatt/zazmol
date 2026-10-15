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

Normal supported installs update a configured instance automatically after
owner verification. The one-time binding belongs to the selected Python
environment at `etc/sassie/version-stack.json`; standalone environments need
no website. Preserve the existing verified web-deployment publisher. An update
failure must fail the command and identify the completed install separately.
The binding schema and operator workflow are owned by
`genapp_zazzie/docs/source_revision_metadata.md` in the paired checkout.
Live binding setup and deployment require their own approved scope.
