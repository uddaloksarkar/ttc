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

## Run one benchmark directly

Build TTC first. The commands below use `./ttc`, as in the top-level README, and
are run from the repository root. They use seed 42, tolerance 0.2, and a
five-step random-direction hit-and-run walk.

**Box example** (two 3D bodies):

```bash
./ttc \
  example/kakeya/hit_and_run/benchmarks/RDHR_box_m0002_d00000010_s0001.smt2 \
  --seed 42 --epsilon 0.2 --delta 0.2 \
  --sampler rdhr --walklen-samp 5 --exact-box-volumes
```

Its exact union volume is `796/5`, recorded in
`hit_and_run/geometry.csv`. TTC prints its sampled union estimate on the
`s volume` line. The exact-box flag supplies exact **body** volumes; TTC still
estimates the volume of their union.

**Sea urchin example** (four slanted bodies in four dimensions):

```bash
./ttc \
  example/kakeya/hit_and_run/benchmarks/RDHR_radial_urchin_m0004_n0004.smt2 \
  --seed 42 --epsilon 0.2 --delta 0.2 \
  --sampler rdhr --walklen-samp 5 --exact-parallelotope-volumes
```

Its exact union volume is `849983697/3125000`, recorded in
`hit_and_run/radial_geometry.csv`. For a box from `ball_walk/benchmarks/`, use
`--sampler ball` with `--exact-box-volumes`. Other walk lengths can be passed
with `--walklen-samp`.

These one-off commands print to the terminal. `accumulate.py` expects the
wrapped, compressed output format from the recorded campaign, so bare TTC
stdout is not an input to it. To revalidate the existing job-2625803 logs and
regenerate its result table and figures, run:

```bash
python3 example/kakeya/accumulate.py \
  example/kakeya/hit_and_run/out_hit_and_run_2625803
```

The accumulator writes `hit_and_run/combined_results_2625803.csv` and
`hit_and_run/figures_2625803/`. The benchmark-index figure uses the target
walk `1000 * dimension`, capped at the last validated walk;
`hit_and_run/fixed_walk_index_2625803.csv` records the plotted walk for each
benchmark.

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
