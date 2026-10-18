# Exact-volume hit-and-run TTC benchmarks

The consolidated failing-input catalog for this sampler is now
[`catalog.csv`](catalog.csv), with selected inputs in [`benchmarks/`](benchmarks/),
the unified [`combined_results.csv`](combined_results.csv), and the combined
[`walk-length error figure`](figures/walk_length_vs_relative_error.png).
Run [`script_ttc.sh`](../../../../script_ttc.sh) from the repository root on
the cluster to produce `.out.xz` files, then
run [`../accumulate.py`](../accumulate.py) to validate them and refresh the
PDF, PGF, and PNG graph. This README also documents the preserved local
fixed-3D campaign, whose original scripts and logs remain in use while
long-running TTC processes finish. Its `results.csv` is separate from
`combined_results.csv`.

This campaign uses TTC's random-direction hit-and-run sampler (`--sampler
rdhr`). Each benchmark is a union of axis-aligned boxes, and TTC receives
their analytic volumes with `--exact-box-volumes`. Runs use seed 42,
`--epsilon 0.2`, and `--delta 0.2`; `--walklen-samp` is varied on identical
SMT2 files. The volume oracle is therefore isolated from the sampling error.

The primary target is **100 TTC cases**: five independently verified box
unions with **2, 4, 8, 16, and 32 bodies**, each tested at the same 20 walk
lengths from 5 through 1,000,000. Every box has z half-width 1/5,000,000;
the x-y staircase geometry and exact union change with body count. For each
fixed input, only `--walklen-samp` changes between its 20 TTC runs.
`benchmark_matrix.csv` records every planned case,
its input, its log, and its measured status. The original 1/5,000,000 input
has exact body volume 0.000160 and exact union volume 0.0003184;
`target_status.csv` retains its detailed 20-walk status. A `pending`, `running`, or
`interrupted_or_unreported` entry is not a TTC result. Only rows marked
`run_status=0` in `results.csv` have a normal TTC exit, the requested RDHR
walk, all exact-volume body rows, a positive union estimate, and the expected
benchmark geometry.

The separate body-count sweep uses the same staircase construction at 2, 4,
8, and 16 bodies. At a fixed denominator and walk length, only the body count
changes. Each body's volume is the same within that sweep, while the exact
union is integrated from the box bounds by an independent verifier.
The [`dimension_equals_bodies/`](dimension_equals_bodies/) candidate keeps
dimension equal to body count. It is stored separately, preserving all
fixed-3D inputs and logs. At eight bodies and walk 1000, the dimension-8
input at half-width 1/50 fails with 36.32% error, while the corresponding
dimension-3 input passes with 19.00%; dimension 8 also fails at half-width
1/10 with 22.18% error. The dimension-8 thickness search found an
adjacent walk-1000 boundary: half-width 1/7 passes (19.80% error) and
1/8 fails (20.90%). At walk 10,000 and half-width 1/500, the matched
dimension-equals-body-count controls now show four bodies passing at 18.34%
error and eight bodies failing at 35.92% error.
See that directory's README and PDF/PGF/PNG comparison figures.

The two families serve different current targets. Among validated runs at
this snapshot, the preserved fixed-3D, two-body input at denominator
5,000,000 fails at the longest completed walk, 76,559 (49.74% error). The
dimension-equals-body-count matrix has a four-body failure through walk
11,143 (64.17% error) and validated 32-body failures through walk 34
(78.35% error); its million-step cases are still active. At walk
1,000 with eight bodies, dimension 8 permits the fattest measured failing
half-width, 1/8, with an adjacent 1/7 pass; the preserved dimension-3
boundary is 1/55 fail and 1/54 pass. `results.csv` in each directory is the
authoritative source as more TTC runs complete.
`figures/` contains PDF, PGF, and PNG exports. Data plots include completed
TTC runs only, and the projection shows the benchmark geometry. Curve labels
identify body counts. The two-body overlap plots compare TTC's
second-body deletion count divided by the first-body sample count with the
exact 1% intersection fraction. This is a logged deletion ratio, not an
independent estimate of overlap probability.
`target_run_coverage` separately shows completed and running target walks.
`benchmark_matrix_union_error` and `benchmark_matrix_coverage` cover the
same 100-case set and are exported as PDF, PGF, and PNG. The error figure
shows only completed TTC runs, with curves labeled by body count. The
coverage figure distinguishes completed, running, and
pending cases. Thus a generated input is not counted as an observed failure.

## Why the thin case fails

TTC starts every RDHR sample at the box's computed inner-ball center.
VolEsti then chooses a random direction and a uniform point on the chord
through the current point. The very narrow z coordinate truncates most such
chords before they can move far along x or y. This provides a mechanism for
samples remaining near the common center even though the true intersection
is only 1% of either box. The logged deletion ratio is near 100% on the
thin-box runs, and the union estimate stays close to one exact body volume.
This explanation is an inference from the sampler code and measured TTC
logs; the exact union and intersection come from the box geometry.

A same-input coordinate-direction hit-and-run (`--sampler cdhr`) control at
walk 5 returned 0.000304243 against the exact 0.0003184 union, or 4.45%
error, and logged 39,179 deletions out of 405,375 first-body samples.
Random-direction hit-and-run at the same walk returned 0.000159494, or
49.91% error, and logged 405,375 deletions. Both runs used the exact box
oracle and exited normally. The [CDHR control log](controls/CDHR_control_m0002_d05000000_walk0000005.log)
is kept separate from the RDHR target result table.
The [direction-control figure](figures/rdhr_vs_cdhr_control.pdf) also has
PGF and PNG exports.

## Reproduce

Start at the repository root with its `build/ttc` executable,
Python 3 with Matplotlib, and `pdflatex` available. `prepare` generates the
SMT2 boxes and independently checks the exact volumes in `geometry.csv`.
This fixed-3D campaign and the dimension-equals-body-count campaign have
separate inputs, logs, result tables, and figures. To regenerate the fixed-3D
100-case target matrix, use the five counts **2, 4, 8, 16, 32**, denominator
**5,000,000**, and the same 20 walk lengths for every input:

```bash
cd example/kakeya/samp_only/hit_and_run
python3 sweep.py prepare --counts 2,4,8,16,32 --denominators 5000000
python3 sweep.py run --counts 2,4,8,16,32 --denominators 5000000 --walks 5,10,18,34,65,124,236,449,853,1622,3083,5861,11143,21183,40271,76559,145544,276692,526015,1000000 --jobs 1 --timeout-seconds 1209600
python3 sweep.py analyze
python3 plot.py
```

The `cd` line assumes the shell starts at the repository root; omit it if
already in this directory. For a quick check, replace the long `--walks`
list with `5`. The wrapper invokes `build/ttc` with `--seed 42 --epsilon
0.2 --delta 0.2 --sampler rdhr --walklen-samp W --exact-box-volumes`.
For example, its walk-10000 eight-body command is:

```bash
../../../../build/ttc benchmarks/RDHR_box_m0008_d05000000_s0001.smt2 --seed 42 --epsilon 0.2 --delta 0.2 --sampler rdhr --walklen-samp 10000 --exact-box-volumes
```

The wrapper writes `logs/RDHR_box_mMMMM_dDDDDDDDD_s0001_walkWWWWWWW.log`.
`sweep.py analyze` rebuilds `results.csv`, `benchmark_matrix.csv`, and
`target_status.csv` from the logs. Interpret a case as measured only when
its result has `run_status=0`; `fails_20pct=1` denotes absolute relative
union error above 20%. `plot.py` reads validated results and writes PDF,
PGF, and PNG to `figures/`, including `benchmark_matrix_union_error` and
`benchmark_matrix_coverage` for the same 100 planned cases. It also makes
the body-count and projection figures from completed data. Long walks can
take days; use a smaller `--jobs` value when memory is limited. A completed
case is reused on rerun. Before replacing an incomplete log, the runner
copies it into `logs/attempts/` with a timestamp.

For the matched fixed-3D body-count example at walk 10000 and denominator
1000, and for the thicker walk-1000 comparison, run:

```bash
python3 sweep.py prepare --counts 2,4,8,16 --denominators 1000
python3 sweep.py run --counts 2,4,8,16 --denominators 1000 --walks 10000 --jobs 1 --timeout-seconds 1209600
python3 sweep.py analyze
python3 plot.py
```

To reproduce the other thickness controls, call `sweep.py prepare` and
`sweep.py run` with the desired positive comma-separated denominators and
walk lengths, then rerun `analyze` and `plot.py`. For the paired
dimension-equals-body-count runs and their PDF/PGF/PNG figures, follow
[`dimension_equals_bodies/README.md`](dimension_equals_bodies/README.md).

The initial walk-5 TTC result failed on all
five earlier two-body thickness controls with about 49.9% union error and
validated exact body volumes. Those controls are separate from the primary
body-count matrix. Thickness refinement remains secondary to finding
body-count failures at higher walks.

The million-step run is launched separately so it can continue while shorter
walks finish. Its live PID is recorded in `target_walk1000000_runner.pid`,
and its TTC output is in `logs/RDHR_box_m0002_d05000000_s0001_walk1000000.log`.
The PID file alone does not establish that the process is live; check the
process table and completed log when interpreting its status.
`monitor_target.py` refreshes `target_status.csv`, `results.csv`, and the
PDF/PGF/PNG figures whenever another one of these 20 target walks validates.
The current monitor's PID and output are in `target_monitor.pid` and
`target_monitor.log`; it exits after all 20 walks validate.

## Current measured findings

In the initial 2026-10-18 measurements, the thin two-box case fails at walks
5 through 124 with about 49.91% union error. A later completed exact-volume
RDHR run at walk **76,559** also fails, with **49.74%** union error. See
`target_status.csv` for all completed target walks; the million-step result
is still running.
At z denominator 1 and walk 5, two bodies pass the 20% limit with 19.91% error,
while four bodies fail with 28.49% error. At walk 10, the corresponding
errors are 13.23% and 22.02%. At denominator 10 and walk 100, four bodies
pass with 19.61% error, while eight fail with 24.40%. These matched
comparisons show that adding bodies can change a passing RDHR TTC estimate
into a failing one. At denominator 100 and walk 1000, two bodies pass with
10.46% error and four pass with 19.50%, while eight fail with 24.20%.
Sixteen bodies also fail there, with 26.65% error.
On that fixed-denominator family, the two- and four-body cases fail at walk
100 but pass at walk 1000; the eight-body case still fails at walk 1000.
The longer target walks remain active. At denominator 1000 and walk 10000,
the two-body case passes with 10.48% error and the four-body case passes
with 19.13%, while eight and sixteen bodies fail with 24.06% and 26.03%
error, respectively. This is a matched
same-dimension, same-walk, same-thickness result in which increasing bodies
changes a passing TTC estimate into a failing one.

The eight-body walk-1000 thickness sweep has an adjacent measured
passing/failing integer boundary: z half-width `1/54` passes with 19.916%
error, while the thinner `1/55` fails with 20.014% error. Thus `1/55` is
the fattest measured failing input among the integer denominators tested in
this fixed-3D family. This is a seed-42 result for the stated TTC settings;
it does not establish a monotone probability boundary. See the
[thickness figure](figures/fatness_boundary_m8_walk1000.pdf), also exported
as PGF and PNG.
