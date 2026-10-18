# Kakeya TTC sampling experiments

This directory has exactly two experiment categories:

- [`samp_only/`](samp_only/) uses TTC's opt-in `--exact-box-volumes` oracle.
  Every body is an axis-aligned box with an analytic product-of-sides volume;
  this isolates error in TTC's union sampling stage. It has two sampler
  directories, [`ball_walk/`](samp_only/ball_walk/) and
  [`hit_and_run/`](samp_only/hit_and_run/), each with selected failing SMT2
  inputs, a combined result CSV, cluster runner, and PDF/PGF/PNG figures.
  Shared catalog, accumulator, and projection scripts are directly under
  `samp_only/`.
- [`vol_and_samp/`](vol_and_samp/) uses TTC's practical volume estimator. The
  body-volume and union-sampling stages can both affect the result. It has the
  same benchmark, figure, and script layout.

The previous Kakeya directory is a local archive and is not included with
these benchmark sources.
All experiment failures compare TTC's returned union volume against the
analytically known union volume at `epsilon = 0.2`. A failed TTC run is an
observed seed-42 result, not a measured failure probability.

The current campaign under
[`samp_only/hit_and_run/`](samp_only/hit_and_run/README.md) uses TTC's
random-direction hit-and-run sampler and exact box volumes. Its 20 target
walks span 5 through 1,000,000. The previous ball-walk jobs were stopped at
the user's request; their completed measurements are selected in
`samp_only/ball_walk/`, and their full history is preserved in
the local archive. The hit-and-run result table distinguishes completed TTC runs
from target walks that are still running or pending.
