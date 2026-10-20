# Exact-volume sampler-only Kakeya benchmarks

Benchmark inputs and generation code are versioned here. Result CSVs, TTC
logs, and figures are generated artifacts and are not part of the source commit.

This directory has two sampler campaigns: [`ball_walk/`](ball_walk/) and
[`hit_and_run/`](hit_and_run/). A benchmark appears in a campaign's
`benchmarks/` only after at least one validated TTC run on that box union
exceeded the 20% relative union-volume error threshold. The same input may
also have passing walks. All runs use TTC's exact axis-aligned box volume
oracle (`--exact-box-volumes`), seed 42, epsilon 0.2, and delta 0.2.
The sampler is `ball` or random-direction hit-and-run (`rdhr`).

Each campaign contains:

- `benchmarks/*.smt2`: selected failing input geometries. A SHA-256 suffix
  prevents different historical inputs with the same name from colliding.
- `catalog.csv`: benchmark filename, body count, dimension, full SHA-256,
  and preserved original location.
- `cases.csv`: benchmark/walk pairs for the cluster runner, including the
  analytic union volume and expected compressed-log filename.
- `measurements.csv`: validated historical measurements for selected inputs.
- `combined_results.csv`: the accumulator's unified result table. `run_status=0`
  means a validated result; pending or malformed cluster logs do not become
  measured failures. `evidence` identifies a prior validated CSV or a newly
  checked compressed TTC log.
- `figures/walk_length_vs_relative_error.{pdf,pgf,png}`: all validated walks
  of all selected benchmark inputs on one graph. Colors and legend labels
  indicate **number of bodies and dimension**. The dashed line is the 20%
  failure threshold. Each input is connected only to its own measurements.

The current catalogs contain 137 ball-walk geometries (164 validated
measurements) and 40 RDHR geometries (151 validated measurements at catalog
creation). The RDHR `cases.csv` also includes longer walks that are still
running locally. Regenerate catalogs and figures when their logs validate:

```bash
python3 example/kakeya/samp_only/hit_and_run/sweep.py analyze
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py analyze
python3 example/kakeya/samp_only/build_catalog.py
python3 example/kakeya/samp_only/accumulate.py
```

The older ball-walk experiment sources and logs are preserved at
[`example/kakeya_previous/samp_only_legacy/`](../../kakeya_previous/samp_only_legacy/).
The older RDHR logs and scripts remain in `hit_and_run/` while their
long-running TTC processes are active; they are not cluster outputs.

## Run the selected benchmarks on a cluster

Build TTC so `build_static/ttc` exists, then stage the repository and benchmark
files on the cluster. From a Slurm allocation at the repository root, run
one or both campaigns with the desired number of ranks, for example:

```bash
srun -n 8 bash script_ttc.sh hit_and_run
srun -n 8 bash script_ttc.sh ball_walk
```

`script_ttc.sh` uses the same 20 walk lengths from 5 through 1,000,000 for
every selected input. The hit-and-run campaign has 40 inputs (800 commands);
the ball-walk campaign has 137 (2,740 commands). It **only runs TTC**. Each
rank takes a disjoint block, copies SMT2 files to scratch, and returns
`<benchmark>_walk<length>.out.xz` and matching `.timeout.xz` files to its
output directory. Existing `.out.xz` files are skipped. `TTC_BIN`,
`SCRATCH`, and `TLIMIT_SECONDS` may be set in the environment. To inspect
the generated first command without running TTC, set `DRY_RUN=1`.

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

Keep downloaded cluster outputs in `out_hit_and_run/` or `out_ball_walk/`
beside their campaign directories. The accumulator reads the sibling
campaign's `catalog.csv`, `cases.csv`, `measurements.csv`, and `benchmarks/`,
then writes `combined_results.csv` and figures into that campaign directory.
For example:

```bash
python3 example/kakeya/samp_only/accumulate.py example/kakeya/samp_only/out_hit_and_run
```

When given an output directory, it also includes compressed logs whose
walk lengths are absent from `cases.csv`.

The accumulator covers the historical `cases.csv` rows and any extra `.out.xz`
files in a supplied output directory. It checks the xz
stream, TTC exit footer and status, sampler and walk,
exact-box oracle, body count, each body volume against the SMT2 box bounds,
and the returned union estimate. It writes one `combined_results.csv` and refreshes
the combined PDF, PGF, and PNG graph in **each** campaign. A previously
validated result remains available if a later cluster rerun is incomplete;
`cluster_log_status` records the rerun's status separately.

## Project one chosen box union

Choose a filename from the relevant `catalog.csv`. For example:

```bash
python3 example/kakeya/samp_only/sunflower.py --sampler hit_and_run \
  --benchmark RDHR_dimeqbody_m0008_d05000000_3b26cafdb64e.smt2
```

The script parses the chosen SMT2 box union, projects each body onto the
`x0`–`x1` plane, draws its rectangle and center-to-corner rays in a
sunflower-like view, and writes PDF, PGF, and PNG under that campaign's
`figures/`. The figure is a projection of the input geometry, not a
projection of sampled TTC points.
