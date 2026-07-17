## Scientific Computing Guidance

- Do not introduce shortcuts, parameter reductions, approximations, or workflow simplifications solely to save time. Preserve the requested physics, parameter space, and analysis unless explicitly directed otherwise.

- When a shortcut or optimization is possible, propose it separately and explain the expected impact. Do not apply it without approval.

- Prefer correctness, reproducibility, and scientific validity over runtime reduction.

- For long-running jobs, implement a periodic heartbeat ("pulse") indicating the process is still active. Include elapsed time and useful progress metrics when available.

- Clearly identify all assumptions. Never silently replace requested behavior with a simpler alternative.

- Algorithmic and computational improvements are encouraged, but should be presented as recommendations for review before adoption.
