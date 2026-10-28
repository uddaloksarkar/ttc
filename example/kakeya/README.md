# Kakeya sampler benchmarks

This directory contains the sampler-only Kakeya benchmarks used with TTC's exact
per-body volume oracles. The selected suites are 26 axis-aligned box unions and
20 radial sea urchins for random-direction hit-and-run, plus 136 box unions for
ball walk. A benchmark can fail the union-volume tolerance even though each
body volume is supplied exactly.
The two box catalogs share only three SMT2 formulas by SHA-256, so their
benchmark sets remain separate.

## Files

| Path | Purpose |
| --- | --- |
| `hit_and_run/benchmarks/` | 26 box and 20 radial SMT2 inputs |
| `hit_and_run/geometry.csv` | Exact box body and union volumes; maximum per-body circumradius/inradius ratio |
| `hit_and_run/radial_geometry.csv` | Exact radial determinant, core, body, and union volumes |
| `ball_walk/benchmarks/` | 136 selected box SMT2 inputs |
| `ball_walk/geometry.csv` | Exact rational volume of every body and union |
| `*/catalog.csv` | Selected input names, dimensions, body counts, and SHA-256 hashes |
| `*/cases.csv` | Historical benchmark/walk combinations and exact union volumes |
| `*/measurements.csv` | Validated historical relative errors |
| `*/combined_results*.csv` | Validated TTC result tables |
| `hit_and_run/figures_2625803/` | Job-2625803 walk/error and geometry figures in PDF, PGF, and PNG |
| `ball_walk/figures/` | Ball-walk figures in PDF, PGF, and PNG |

The shared generator checks written box inequalities against independent
exact-volume calculations. It constructs hit-and-run 3D staircase boxes,
boxes with dimension equal to body count, and slanted radial sea urchins.
The ball-walk inputs are a curated historical set; the same command verifies
their SMT2 hashes and rebuilds their exact geometry metadata. With no options
it prepares both suites:

```bash
python3 example/kakeya/generate_benchmarks.py
```

Use `--sampler hit_and_run` or `--sampler ball_walk` to prepare one suite.
For hit-and-run, the optional modes are `box`, `staircase`, and `radial`. For example,
`staircase --counts 8 --denominator 500` selects the eight-dimensional box
case. Denominators of 4,000,000 or greater are excluded.

## Run and analyze

From the repository root, build TTC and run a sampler campaign under Slurm:

```bash
srun -n 8 bash script_ttc.sh hit_and_run
srun -n 8 bash script_ttc.sh ball_walk
```

The hit-and-run campaign has 46 inputs and 20 walk lengths, for 920 jobs. The
ball-walk campaign has 136 inputs and 20 walk lengths, for 2,720 jobs. Set
`TTC_BIN`, `SCRATCH`, or `TLIMIT_SECONDS` if needed. The runner writes compressed
logs to `<sampler>/out_<sampler>_<jobid>/` by default. To validate the
job-2625803 output and regenerate its
result table and figures:

```bash
python3 example/kakeya/accumulate.py example/kakeya/hit_and_run/out_hit_and_run_2625803
```

The accumulator verifies each available compressed log, checks the exact body
volume oracle, and writes `hit_and_run/combined_results_<jobid>.csv` and
`<sampler>/figures_<jobid>/`. Ball-walk local figures go to
`ball_walk/figures/`. Figures include the box and radial walk-length
curves and one two-line benchmark-index comparison. For that comparison the
target is `1000 * dimension`; when it exceeds the last validated walk, the
plotted error uses that last walk. `fixed_walk_index_<jobid>.csv` records both
walk lengths for every plotted benchmark.

To produce a geometry projection from a benchmark represented in an output
directory:

```bash
python3 example/kakeya/plot_projections.py \
  example/kakeya/hit_and_run/out_hit_and_run_2625803 \
  --benchmark RDHR_radial_urchin_m0016_n0016.smt2
```

The projections show the SMT2 geometry; they do not plot TTC samples.
For a ball-walk input without a compressed output directory, use
`python3 example/kakeya/plot_projections.py --sampler ball_walk --benchmark NAME.smt2`.
