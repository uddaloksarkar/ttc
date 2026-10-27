# Exact-volume sampler-only Kakeya benchmarks

Benchmark inputs and generation code are versioned here. Result CSVs, TTC
logs, and figures are generated artifacts and are not part of the source commit.

This directory has two sampler campaigns: [`ball_walk/`](ball_walk/) and
[`hit_and_run/`](hit_and_run/). A benchmark enters the cluster catalog only
after at least one validated TTC run on that box union exceeded the 20%
relative union-volume error threshold. The same input may
also have passing walks. All runs use TTC's exact axis-aligned box volume
oracle (`--exact-box-volumes`), seed 42, epsilon 0.2, and delta 0.2.
The sampler is `ball` or random-direction hit-and-run (`rdhr`).

Each campaign contains:

- `benchmarks/*.smt2`: source input geometries. The hit-and-run catalog uses
  SHA-256 names for selected inputs, but the runner stages those names in
  scratch rather than keeping duplicate files here.
- `catalog.csv`: benchmark filename, body count, dimension, full SHA-256,
  and preserved original location.
- `cases.csv`: benchmark/walk pairs for the cluster runner, including the
  analytic union volume and expected compressed-log filename.
- `measurements.csv`: validated historical measurements for selected inputs.
- `combined_results_<jobid>.csv` for cluster outputs (or `combined_results.csv`
  for a campaign run): the accumulator's unified result table. `run_status=0`
  means a validated result; pending or malformed cluster logs do not become
  measured failures. `evidence` identifies a prior validated CSV or a newly
  checked compressed TTC log.
- `hit_and_run/figures_<jobid>/` contains separate
  `boxes_walk_vs_relative_error_<jobid>` and
  `radial_walk_vs_relative_error_<jobid>` figures in PDF, PGF, and PNG.
  Box color darkens as width decreases; radial color darkens as body count
  increases. Both use the same walk-length axis and lightly shade the region
  below the 20% error threshold. The box plot excludes `d >= 4,000,000`.
  Ball-walk retains its combined walk-length plot.

The hit-and-run local pilot runners and their logs have been retired. The
existing catalog metadata and compressed `out_hit_and_run*` files are the
inputs for its current figures. See the [hit-and-run README](hit_and_run/README.md)
for generation, accumulation, and projections.

The older ball-walk experiment sources and logs are preserved at
[`example/kakeya_previous/samp_only_legacy/`](../../kakeya_previous/samp_only_legacy/).

## Run the selected benchmarks on a cluster

Build TTC so `build_static/ttc` exists, then stage the repository and benchmark
files on the cluster. From a Slurm allocation at the repository root, run
one or both campaigns with the desired number of ranks, for example:

```bash
srun -n 8 bash script_ttc.sh hit_and_run
srun -n 8 bash script_ttc.sh ball_walk
```

`script_ttc.sh` uses the same 20 walk lengths from 5 through 1,000,000 for
every selected input. The hit-and-run campaign has 40 catalog inputs and 20
radial inputs (1,200 commands); the ball-walk campaign has 137 (2,740
commands). It **only runs TTC**. Each
rank takes a disjoint block, copies SMT2 files to scratch, and returns
`<benchmark>_walk<length>.out.xz` and matching `.timeout.xz` files to its
output directory. Existing `.out.xz` files are skipped. `TTC_BIN`,
`SCRATCH`, and `TLIMIT_SECONDS` may be set in the environment. To inspect
the generated first command without running TTC, set `DRY_RUN=1`.

The 20 radial hit-and-run inputs run at all 20 walk lengths using
`--exact-parallelotope-volumes`. `script_ttc.sh hit_and_run` schedules the full
1,200-case sweep. To run only the 400 radial cases, use
`srun -n 4 bash script_ttc.sh hit_and_run_radial` with a TTC binary built
from the current source. The other 53 hit-and-run source inputs are in
[`hit_and_run/unselected_benchmarks/`](hit_and_run/unselected_benchmarks/),
listed in [`unselected_benchmarks.csv`](hit_and_run/unselected_benchmarks.csv).
They are outside the selected catalog and are not scheduled by
`script_ttc.sh hit_and_run`.

After the cluster outputs are copied back, run the accumulator from the
repository root:

```bash
python3 example/kakeya/samp_only/accumulate.py
```

To process one campaign directory, pass its path instead:

```bash
python3 example/kakeya/samp_only/accumulate.py example/kakeya/samp_only/hit_and_run
```

The sampler is inferred from a `ball_walk` or `hit_and_run` directory name.
For a differently named directory, specify `--sampler ball_walk` or
`--sampler hit_and_run`.

`script_ttc.sh` writes to `out_hit_and_run_<jobid>/` or
`out_ball_walk_<jobid>/` (or appends `_<jobid>` to `OUT_DIR`). Keep downloaded
cluster outputs beside their campaign directories. The accumulator reads the sibling
campaign's `catalog.csv`, `cases.csv`, `measurements.csv`, and `benchmarks/`,
then writes `combined_results_<jobid>.csv` and `figures_<jobid>/` into that
campaign directory.
For example:

```bash
python3 example/kakeya/samp_only/accumulate.py example/kakeya/samp_only/out_hit_and_run_123456
```

When given an output directory, it also includes compressed logs whose
walk lengths are absent from `cases.csv`.

When supplied an output directory, the accumulator covers the historical
`cases.csv` rows, the 400 planned radial cases, and any extra `.out.xz` files.
It validates box cases against their SMT2 bounds and radial cases against the
exact volumes in `radial_geometry.csv`. It writes one
`combined_results_<jobid>.csv` and refreshes one walk-length versus relative-error
figure in PDF, PGF, and PNG formats. When run on an output directory,
the plots use validated logs from that directory. A previously
validated result remains available if a later cluster rerun is incomplete;
`cluster_log_status` records the rerun's status separately.

## Project a benchmark represented in an output directory

For hit-and-run box and radial projections, use `hit_and_run/plot_projections.py`
with an `out_hit_and_run_<jobid>` directory and one or more `--benchmark`
filenames. See the [hit-and-run instructions](hit_and_run/README.md).
The figures show the SMT2 input geometry, not sampled TTC points.
