# Ball-walk sampler-only failures

`benchmarks/` contains 137 distinct box unions with at least one validated
ball-walk TTC union-volume error above 20%. `catalog.csv` records body count,
dimension, SHA-256, and the original archived path. `cases.csv` lists the
164 observed benchmark/walk combinations. All runs use exact box volumes,
seed 42, epsilon 0.2, and delta 0.2.

Run `creator.sh` under MPI or Slurm to produce only `.out.xz` and
`.timeout.xz` TTC artifacts. From the repository root, run
`python3 example/kakeya/samp_only/accumulate.py --sampler ball_walk` to
validate compressed logs, write `combined_results.csv`, and create the
combined walk-length versus relative-error figure in PDF, PGF, and PNG.
When given `out_ball_walk_<jobid>/`, the accumulator writes
`combined_results_<jobid>.csv` and `figures_<jobid>/` instead.
The legend labels body count and dimension. See the
[`sampler-only README`](../README.md) for cluster commands and the
sunflower projection command.
