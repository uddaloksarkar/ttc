# Exact-volume hit-and-run TTC benchmarks

The 60 selected SMT2 inputs (40 box benchmarks and 20 radial sea urchins) are
in `benchmarks/`; the other 53 inputs are in `unselected_benchmarks/`.
`catalog.csv`, `cases.csv`, `measurements.csv`, `unselected_benchmarks.csv`, and
`radial_geometry.csv` provide the metadata needed to validate cluster output.
The local `logs/` directory and its retired runner scripts have been removed.

## Generate inputs

From the repository root, regenerate and verify the complete inventory:

```bash
python3 example/kakeya/samp_only/hit_and_run/generate_benchmarks.py
```

The modes `box`, `staircase`, `factored`, and `radial` select one family. For
example:

```bash
python3 example/kakeya/samp_only/hit_and_run/generate_benchmarks.py staircase --counts 4,8 --denominator 500
```

The generator checks exact rational body and union volumes. Its fixed 3D box
check independently parses the written SMT2 bounds; the radial geometry and
exact volumes are recorded in `radial_geometry.csv`.
For all 93 box benchmarks (fixed 3D, staircase, and factored), `geometry.csv`
records the exact volume of one body and `max_body_R_over_r`, with a decimal
approximation. The ratio is the maximum over individual boxes; each box's
Euclidean circumradius reaches a corner, and its inradius is its smallest
half-width. `utility.py` contains the shared exact per-body volume calculations
for boxes and radial sea urchins. Regenerate this metadata with the no-argument
generator command above.

## Figures from a cluster output directory

Keep `out_hit_and_run_<jobid>/` beside `hit_and_run/`. To validate its
compressed TTC logs and update the two walk length versus relative error plots:

```bash
python3 example/kakeya/samp_only/accumulate.py example/kakeya/samp_only/out_hit_and_run_<jobid>
```

The results table is `hit_and_run/combined_results_<jobid>.csv`. The box and
radial plots are written to `hit_and_run/figures_<jobid>/` as PDF, PGF, and
PNG. Their log-scale walk axes have identical limits; the pale region below
`epsilon = 0.2` marks acceptable error. Darker radial curves have more
bodies, while darker box curves have smaller width. The box plot omits
denominators `d >= 4,000,000`. To plot
an input geometry represented in that output directory, choose a benchmark
filename from `catalog.csv` or `radial_geometry.csv` and run:

```bash
python3 example/kakeya/samp_only/hit_and_run/plot_projections.py \
  example/kakeya/samp_only/out_hit_and_run_<jobid> \
  --benchmark RDHR_radial_urchin_m0004_n0004.smt2
```

Repeat `--benchmark` to plot more inputs. The script accepts box and radial
benchmarks, checks that the output directory contains a matching `.out.xz`,
and writes the projections into `hit_and_run/figures_<jobid>/` with the job
ID in each filename. The projections show SMT2 geometry, not sampled
TTC points. Historical local result CSVs are retained only as inputs to
`build_catalog.py`; they are not read to make these figures.

For 3D boxes that differ only in their thin `x2` width, use `--plane 3d`
for a spatial view, or `--plane x1x2` for a flat projection. Pass the same
`--x2-limit` for a directly comparable set of plots:

```bash
python3 example/kakeya/samp_only/hit_and_run/plot_projections.py \
  example/kakeya/samp_only/out_hit_and_run_<jobid> \
  --plane 3d --x2-limit 0.11 \
  --benchmark RDHR_box_m0008_d00000010_s0001_b6176ccad001.smt2
```
