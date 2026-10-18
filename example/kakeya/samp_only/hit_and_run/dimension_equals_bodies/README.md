# Dimension equals body count: RDHR pilot

This is a separate candidate family. The existing fixed-3D TTC inputs,
completed logs, and figures in the parent directory are preserved.
Across body counts in this family, dimension changes too. Therefore those
rows establish failures for dimension-equals-body-count configurations,
while the parent's fixed-3D body-count comparison isolates the effect of
adding bodies.

There are five exact-volume box unions with 2, 4, 8, 16, and 32 bodies,
respectively. Each input has that same number of real variables. The first
two coordinates form the previous x-y staircase. For dimension at least 3,
coordinate 2 has half-width 1/5,000,000 and every further coordinate has
half-width 1. Thus body volume and the exact union are simple products of
the integrated x-y staircase with the common-coordinate widths. The
generator parses the written SMT2 file again and integrates its parsed
bounds before reporting the exact values.

## Reproduce the dimension-equals-body-count matrix

Start at the repository root; the first command enters this experiment's
directory. You need the built executable
`build/ttc`, Python 3 with Matplotlib, and `pdflatex` for PGF output. The
generator writes five SMT2 files in `benchmarks/` and independently checks
their exact body and union volumes in `geometry.csv`:

```bash
cd example/kakeya/samp_only/hit_and_run/dimension_equals_bodies
python3 probe.py prepare
```

The body counts are **2, 4, 8, 16, 32**. Each input has exactly that many
real dimensions. All have the same thin-coordinate denominator 5,000,000.
The 20 target walk lengths are **5, 10, 18, 34, 65, 124, 236, 449, 853,
1622, 3083, 5861, 11143, 21183, 40271, 76559, 145544, 276692,
526015, 1000000**. Thus the matrix contains 5 × 20 = 100 TTC cases. Use
one job at first; increase `--jobs` only if the machine has enough memory:

```bash
python3 probe.py run --counts 2,4,8,16,32 --walks 5,10,18,34,65,124,236,449,853,1622,3083,5861,11143,21183,40271,76559,145544,276692,526015,1000000 --jobs 1 --timeout-seconds 1209600
python3 probe.py analyze
python3 plot.py
```

For a quick end-to-end check, replace the long `--walks` value with `5`.
For a single case, use `--counts 8 --walks 1000` only in the separate
thickness workflow below: walk 1000 is **not** one of the matrix's 20 walks.
For a matrix case, for example, use `--counts 8 --walks 853`.

The wrapper runs exactly this TTC configuration for each requested case
(substitute a matrix SMT2 filename and target walk):

```bash
../../../../../build/ttc benchmarks/RDHR_dimeqbody_m0008_d05000000.smt2 --seed 42 --epsilon 0.2 --delta 0.2 --sampler rdhr --walklen-samp 853 --exact-box-volumes
```

`probe.py run` writes one log per case in `logs/`. It skips completed cases;
`probe.py analyze` rebuilds `results.csv` from the logs, including pending
rows. `run_status=0` means the validator found TTC's normal exit, the exact
box oracle, the requested sampler and walk, all body-volume rows, and a
positive union estimate. `fails_20pct=1` means absolute relative union
error exceeds 0.2. A full TTC `c exit 0` footer can validate a run even if
its wrapper stopped before appending `run_status=0`. Pending or partial logs
are not measured outcomes. Long walks, especially at high body counts, can
take days. Before replacing an incomplete log, the runner copies it into
`logs/attempts/` with a timestamp; completed logs are reused.

`plot.py` reads the validated rows and writes **PDF, PGF, and PNG** under
`figures/` for the matrix error, coverage, body-count/dimension comparison,
projection, and completed thickness comparisons. The central outputs are
`figures/matrix_error_vs_walk.{pdf,pgf,png}` and
`figures/matrix_coverage.{pdf,pgf,png}`. Both represent the same 100 planned
cases; the error plot shows only completed runs. The fixed-3D comparison
requires the parent directory's `results.csv`, regenerated with the parent
README's commands.

In the validated matrix, the four-body/four-dimensional case still fails
at walk 11143 with 64.17% relative union-volume error. The
sixteen-body/sixteen-dimensional case fails through walk 236 with 75.67% error.
The eight-body/eight-dimensional case fails through walk 5861 with 71.40% error.
The thirty-two-body/thirty-two-dimensional case fails through walk 34,
with 78.35% error at each completed walk.
These statements are supported by completed TTC logs and exact-volume
validation rows; longer walks remain pending. The
two-, four-, eight-, sixteen-, and thirty-two-body dimension-matched walk-1,000,000
cases are active and are not yet measured outcomes.

`fatness_compare.py 50` tests the eight-dimensional, eight-body input at
walk 1000 with half-width 1/50. The existing three-dimensional eight-body
input at that same thickness and walk passed with 19.00% union-volume error.
This comparison tests whether dimension matching permits a fatter failure.
Neither pending nor active runs are measured outcomes.
`figures/projection_m8_xy` shows the positive-quadrant x-y projection of
the eight-body staircase in PDF, PGF, and PNG. All thickness variants of
this eight-body family have this same projection; only their common third
coordinate width changes.

The dimension-eight denominator-50 TTC run completed its own normal-exit
footer, all eight exact-volume body rows, and a positive union estimate.
The wrapper ended before appending its usual status line, so the validator
records `completion_evidence=ttc_exit_footer` and checks the full footer and
body rows. The measured error is **36.32%**, compared with **19.00%** in
dimension 3. Thus dimension matching improves this benchmark at the fixed
walk-1000, denominator-50 configuration. The previous 3D files remain as
controls. The `fatness_comparison_m8_walk1000_d50` figure is available in
PDF, PGF, and PNG. The denominator-25 and denominator-10 runs also
completed and failed at 29.79% and 22.18% error, respectively. Half-width
1/10 is 5.5 times the fattest measured failing half-width in the preserved
3D family (1/55). Denominators 1 through 7 all completed and passed;
D7 has 19.803% error. Denominators 8 through 10 fail; D8 has 20.900%
error. **D8 is the fattest measured failing integer denominator** in this
dimension-8, eight-body, walk-1000 family, with an adjacent D7 pass.
Its half-width 1/8 is 6.875 times the preserved 3D family's fattest
measured failing half-width 1/55. These are seed-42 TTC outcomes, not a
claim that failure probability changes monotonically with thickness.
`fatness_boundary_dimension8_m8_walk1000` plots all completed dimension-8
thickness cases; `fatness_comparison_m8_walk1000` plots matched 3D and 8D
results where both TTC runs are complete.

Because this matched comparison improves the measured failure, the next
dimension-matched matrix runs have begun. The 2-body and 4/8-body walkers
cover the shorter target walks in `results.csv`; their live runner PIDs are
recorded in `matrix_2_runner.pid` and `matrix_4_8_runner.pid`. A PID file is
only a pointer: check `/proc` and the TTC logs for current state.
The 16-body walk-5, walk-10, walk-18, walk-34, walk-65, walk-124, and
walk-236 pilots completed and failed. The 32-body walk-5, walk-10, walk-18,
and walk-34 cases also completed and failed, each at 78.35% relative error.
The first 32-body walk-5 attempt had an insufficient one-hour wrapper
timeout; its partial log remains preserved under `logs/attempts/`. The
successful rerun used
the 14-day timeout under `matrix_32_walk5_runner.pid`.
Existing 3D results and figures have not been removed.

A matched walk-10000, half-width-1/500 comparison has completed. The
four-body/four-dimensional input passes with 18.34% relative union error;
the eight-body/eight-dimensional input fails with 35.92% error. Both use
exact body volumes, seed 42, epsilon and delta 0.2, and RDHR at walk 10000.
Their validated rows and TTC logs are in
`walk0010000_m0004_d00000500/` and `walk0010000_d00000500/`.
`figures/body_count_effect_walk10000_d500.{pdf,pgf,png}` shows both cases
against the 20% limit. This is a matched high-walk example in which adding
bodies changes a pass to a failure while keeping the thin-coordinate width
and sampling walk fixed; dimension changes with body count as requested.

`queue_remaining.py` covers the 46 matrix cases not assigned to the
short/medium-walk runners: the seven higher target walks for counts 2, 4,
8, and 16, plus count-32 walks from 18 through 1,000,000. It limits the
combined live RDHR TTC process count to 18 and uses at most three workers.
It archives any interrupted log before retrying it. The detached queue's
PID is in `remaining_queue.pid`, with decisions in `remaining_queue.log`.
It started the dimension-8 walk-1,000,000 case early; this is an active
attempt, not a measured failure.
The dimension-2 walk-1,000,000 control is active under user service
`ttc-rdhr-m2-w1m.service`, with TTC output in
`logs/RDHR_dimeqbody_m0002_walk1000000.log`. It will provide a matched
body-count baseline if it validates.
The dimension-4 walk-1,000,000 case was started separately under user
service `ttc-rdhr-m4-w1m.service`, with TTC output in
`logs/RDHR_dimeqbody_m0004_walk1000000.log`. It uses the same matrix
settings and is unmeasured until TTC exits normally and validates.
The dimension-16 walk-1,000,000 case was also started early as the user
service `ttc-rdhr-m16-w1m.service` so body-count scaling can be tested at
the million-step target. Its TTC output is
`logs/RDHR_dimeqbody_m0016_walk1000000.log`; inspect the live service with
`systemctl --user status ttc-rdhr-m16-w1m.service`. Its entry in
`results.csv` remains unvalidated until `probe.py analyze` checks a normal
TTC exit. This run uses the same generated benchmark and seed, epsilon,
delta, sampler, and exact body-volume settings as the matrix.
The dimension-32 walk-1,000,000 case was started under user service
`ttc-rdhr-m32-w1m.service`, with TTC output in
`logs/RDHR_dimeqbody_m0032_walk1000000.log`. It also remains unmeasured
until the full TTC result validates.

`large_counts/` contains separately staged 64-body/64-dimensional and
100-body/100-dimensional SMT2 unions. `prepare_large.py` parsed their
written bounds and independently integrated the exact body and union
volumes; both volumes fit finite double precision. No TTC failure is
claimed for either larger input until a full TTC run validates. The
original-form 64D walk-5 attempt is retained under
`large_counts/m64_walk5_runner.pid` and its partial log. The original-form
100D walk-5 pilot was started below the
18-process RDHR limit by `launch_100_when_slot.py`; its PID and launcher log
are retained in `large_counts/m100_walk5_runner.pid` and
`large_counts/m100_walk5_runner.log`.
Both original-form high-count attempts were stopped during DNFization to
release resources; neither high-count input has a validated TTC failure yet.

The high-count SMT2 file originally repeated the same 98 common-coordinate
bounds inside every one of the 100 bodies. `prepare_factored.py` writes an
equivalent formula with those bounds once, outside the x-y body disjunction.
It parses the written bounds and checks the same exact body and union
volumes. An eight-body, eight-dimensional walk-5 TTC control using this
factored form completed with all eight exact-volume rows and exactly the
same union estimate, 0.0051125, as the original form. The stopped original
64- and 100-body attempts remain in `large_counts/logs/` with
`run_status=-15`; neither finished DNFization. The equivalent factored
walk-5 pilots are active under `factored/m64_walk5_runner.pid` and
`factored/m100_walk5_runner.pid`. Their completion and errors must be
validated separately; the eight-body control does not prove either
high-count case has finished. Both factored inputs have completed
DNFization; TTC reports 64 real variables/64 polytopes and 100 real
variables/100 polytopes, respectively. Their body sampling and final union
estimates are still pending. The 64-body run has emitted its first two
exact-volume body rows; body 2 logged 217,171 deletions against 217,171
first-body samples. The exact intersection of the first two boxes occupies
7/18 (38.89%) of the first box, so the logged deletion ratio is much larger
at this point. It is an intermediate TTC count, not an independent overlap
probability estimate or a validated union-volume failure claim.

## Reproduce the thickness and larger-count pilots

From this directory, `fatness_compare.py` changes the common thin-coordinate
half-width to `1/D`, where `D` is its positional denominator. Each run
generates its own SMT2 file, verifies its exact volumes, runs TTC with the
same seed/oracle/sampler settings, and writes `result.csv` alongside the log.
The measured dimension-eight, eight-body walk-1000 boundary can be
regenerated with:

```bash
for D in 1 2 3 4 5 6 7 8 9 10 25 50; do
  python3 fatness_compare.py "$D" --count 8 --walk 1000
done
python3 plot.py
```

For the matched walk-10000 candidates, run these separately. The `--count`
value also sets the dimension:

```bash
python3 fatness_compare.py 500 --count 4 --walk 10000
python3 fatness_compare.py 500 --count 8 --walk 10000
python3 plot.py
```

The 64- and 100-body inputs are **outside the 100-case matrix**. Generate
and verify both, then run one pilot at a time:

```bash
python3 prepare_large.py
python3 run_large.py 64 --walk 5 --timeout-seconds 1209600
python3 run_large.py 100 --walk 5 --timeout-seconds 1209600
```

Their SMT2 files and exact volumes appear in `large_counts/benchmarks/` and
`large_counts/geometry.csv`; logs and one-row validation files appear in
`large_counts/logs/` and `large_counts/result_mNNNN_walkWWWWWWW.csv`.
`plot.py` currently plots the five-count matrix and thickness comparisons;
the 64/100-body pilots are reported in their validation files, not silently
added to the 100-case figure. Repeating a completed pilot reuses its log.

To regenerate the equivalent factored high-count inputs and run one TTC
case, use:

```bash
python3 prepare_factored.py
python3 run_factored.py 8 --walk 5 --timeout-seconds 1209600
python3 run_factored.py 64 --walk 5 --timeout-seconds 1209600
python3 run_factored.py 100 --walk 5 --timeout-seconds 1209600
python3 plot_factored.py
```

The validator writes `factored/result_m0100_walk0000005.csv`. To validate
an existing log without running TTC again, add `--analyze-only`. The
eight-body equivalence control is `factored/result_m0008_walk0000005.csv`.
`plot_factored.py` reads and validates the logs directly and writes PDF,
PGF, and PNG for `error_vs_body_count_walk5`, `coverage_walk5`, and
`xy_projection` under `factored/figures/`. The error graph omits incomplete
cases; the coverage graph names them.
`monitor_factored.py` can keep those figures current as the long runs
finish. The current monitor's PID and output are in `factored/monitor.pid`
and `factored/monitor.log`; it stops after all three walk-5 logs validate.

An optional preprocessing comparison uses the same factored SMT2 inputs
with `--no-cdd-simp`. This skips cddlib constraint canonicalization; every
axis-aligned box in these inputs already has exactly its necessary lower
and upper bounds. The completed eight-body walk-5 control has the same
eight exact-volume rows and union estimate, 0.0051125, with and without
this option. Its separate log and validation row are in `factored/nocdd/`.
A 64-body no-cdd walk-5 runtime trial was stopped after about 10 minutes
without a body row; its partial log and `run_status=-15` validation row are
preserved. It showed lower memory use, but no completed speed or union-error
result. The default 64-body run continues. Reproduce the no-cdd comparison
with:

```bash
python3 run_factored.py 8 --walk 5 --no-cdd-simp
python3 run_factored.py 64 --walk 5 --no-cdd-simp --timeout-seconds 1209600
```

The no-cdd logs never replace the default factored logs.

Initial validated walk-5 outcomes: the two-dimensional two-body case passes
at 8.42% error, whereas its fixed-3D control fails at 49.91%. The
four-dimensional four-body case fails at 64.196%, almost identical to its
fixed-3D control. The eight-dimensional eight-body case fails at 71.402%,
also nearly identical to its fixed-3D control. The 16-dimensional walk-5
case fails at 75.667%; the 32-dimensional walk-5 TTC run remains active.
See `figures/dimension_comparison_walk5` in
PDF, PGF, and PNG; it includes completed pairs only.
