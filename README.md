# TTC: Toolbox for Theory Counting

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**ttc**, the Toolbox for Theory Counting, is a tool for measuring the solution space of Satisfiability Modulo Theories (SMT) formulas. It supports counting over a number of theories and their combinations, including bit-vectors and uninterpreted functions, as well as volume computation for the theory of linear real arithmetic.

## Documentation

A documentation on usage and building is available at https://arijitsh.github.io/docs-ttc/

## Reproduce the Kakeya benchmarks and figures

The current sampler-only layout and cluster workflow are described in the
[`samp_only` README](example/kakeya/samp_only/README.md): separate
`ball_walk/` and `hit_and_run/` catalogs of validated failing inputs,
parallel `creator.sh` runners that return xz logs, one accumulator per
sampler, combined PDF/PGF/PNG error figures, and a selectable projected
sunflower view. The commands below document the preserved local RDHR sweeps
that still feed that catalog while long TTC jobs finish.

Run the commands in this section from the **repository root**. They use the
built `build/ttc` executable, Python 3 with Matplotlib, and `pdflatex` for
Matplotlib's PGF export. The current experiment uses random-direction
hit-and-run (`--sampler rdhr`) and TTC's exact axis-aligned box volume oracle
(`--exact-box-volumes`). Every run fixes seed 42, `epsilon = delta = 0.2`;
for a given SMT2 input, only `--walklen-samp` changes. A body-count comparison
also changes the input geometry, as recorded in the corresponding manifest.

The main set has five inputs with **2, 4, 8, 16, and 32 bodies**, each with
dimension equal to its body count. Each input is tested at the same **20 walk
lengths**, from 5 to 1,000,000, for 100 planned TTC runs. First generate and
check the exact body and union volumes by parsing the generated SMT2 files:

```bash
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py prepare
```

The inputs appear in
`example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/benchmarks/`;
their exact volumes are in `geometry.csv` in the same directory. For a quick
TTC smoke run on one of these inputs, run:

```bash
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py run --counts 8 --walks 5 --jobs 1
```

To run the **entire 100-case matrix**, use the following command. These long
walks and higher-dimensional inputs can require substantial CPU time and
memory; `--jobs 1` runs one TTC process at a time. Increase it only when the
machine has capacity. The per-case timeout below is 14 days, not a promise
that every case will finish within that time.

```bash
WALKS=5,10,18,34,65,124,236,449,853,1622,3083,5861,11143,21183,40271,76559,145544,276692,526015,1000000
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py run \
  --counts 2,4,8,16,32 --walks "$WALKS" --jobs 1 --timeout-seconds 1209600
```

The wrapper writes one TTC log per case in
`example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/logs/`. It
reuses validated completions and archives a partial log under that
directory's `attempts/` subdirectory before a rerun.
The equivalent manual command for the eight-body, walk-5 case is:

```bash
build/ttc example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/benchmarks/RDHR_dimeqbody_m0008_d05000000.smt2 \
  --seed 42 --epsilon 0.2 --delta 0.2 --sampler rdhr \
  --walklen-samp 5 --exact-box-volumes
```

After runs finish, validate the main matrix logs:

```bash
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py analyze
```

The analyzer writes `dimension_equals_bodies/results.csv` beneath the
hit-and-run directory. Only `run_status=0` rows are completed, validated TTC
measurements; `fails_20pct=1` means the absolute relative error in TTC's
union-volume estimate exceeds 20% of the analytic union volume. Running,
pending, interrupted, and timed-out cases are **not** measured failures. The
error graph plots completed cases only; its companion coverage graph shows
which of the same 100 planned cases have completed. Both, the eight-body
projection, and comparison figures are written as **PDF, PGF, and PNG** under
`example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/figures/`.
The main files are `matrix_error_vs_walk.{pdf,pgf,png}` and
`matrix_coverage.{pdf,pgf,png}`.

The dimension-matched plot script also reads the preserved **fixed-3D**
`results.csv` for its comparison figures. If that CSV is already present,
generate both sets of figures now with:

```bash
python3 example/kakeya/samp_only/hit_and_run/plot.py
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot.py
```

To reproduce both experiment sets from scratch, first generate and run the
fixed-3D inputs, then analyze both sets and plot them:

```bash
python3 example/kakeya/samp_only/hit_and_run/sweep.py prepare \
  --counts 2,4,8,16,32 --denominators 5000000
python3 example/kakeya/samp_only/hit_and_run/sweep.py run \
  --counts 2,4,8,16,32 --denominators 5000000 \
  --walks "$WALKS" --jobs 1 --timeout-seconds 1209600
python3 example/kakeya/samp_only/hit_and_run/sweep.py analyze
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/probe.py analyze
python3 example/kakeya/samp_only/hit_and_run/plot.py
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot.py
```

For the fixed-3D set, the matching error and coverage figures are
`benchmark_matrix_union_error.{pdf,pgf,png}` and
`benchmark_matrix_coverage.{pdf,pgf,png}` in `hit_and_run/figures/`.
The plot scripts can also be rerun on the existing validated CSV files
without rerunning TTC.

### Larger-body and thickness pilots

The **64- and 100-body** inputs are separate pilots, outside the 100-case
matrix. Their common-coordinate bounds are factored out of the body
disjunction to reduce SMT2 preprocessing. Regenerate and verify the inputs,
then run each pilot separately because they require substantial memory:

```bash
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/prepare_factored.py
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/run_factored.py 8 --walk 5
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/run_factored.py 64 --walk 5 --timeout-seconds 1209600
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/run_factored.py 100 --walk 5 --timeout-seconds 1209600
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot_factored.py
```

The eight-body run checks the factored formula against the original input.
Each run writes a log and a `result_mNNNN_walkWWWWWWW.csv` validation row
under `dimension_equals_bodies/factored/`; `run_status=0` is required before
calling it a measured failure. The factored error, coverage, and projection
figures are exported as PDF, PGF, and PNG under `factored/figures/`. To
revalidate an existing factored log without launching TTC, use, for example,
`run_factored.py 64 --walk 5 --analyze-only` with the full script path above.

The measured **eight-body, dimension-eight, walk-1000 thickness search**
varies the common thin-coordinate half-width `1/D`. Its scripts live beside
the main matrix; this command regenerates the tested integer denominators,
their TTC logs, validation rows, and PDF/PGF/PNG comparison figures:

```bash
for D in 1 2 3 4 5 6 7 8 9 10 25 50; do
  python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/fatness_compare.py "$D" --count 8 --walk 1000
done
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot.py
```

The adjacent measured boundary in this family is a **passing** half-width
`1/7` and a **failing** half-width `1/8`, with the same TTC settings. For its
fixed-3D control, regenerate the adjacent denominator-54 passing and
denominator-55 failing cases and their comparison figures:

```bash
python3 example/kakeya/samp_only/hit_and_run/sweep.py prepare --counts 8 --denominators 54,55
python3 example/kakeya/samp_only/hit_and_run/sweep.py run --counts 8 --denominators 54,55 --walks 1000 --jobs 1
python3 example/kakeya/samp_only/hit_and_run/sweep.py analyze
python3 example/kakeya/samp_only/hit_and_run/plot.py
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot.py
```

These thickness boundaries are for walk 1000; a million-step thickness
boundary has **not** been established.

For the completed **walk-10,000 body-count comparison** at half-width
`1/500`, regenerate the four- and eight-body inputs, TTC results, and
`body_count_effect_walk10000_d500.{pdf,pgf,png}`:

```bash
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/fatness_compare.py 500 --count 4 --walk 10000
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/fatness_compare.py 500 --count 8 --walk 10000
python3 example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/plot.py
```

The matched four-body case passes with 18.34% union error; the eight-body
case fails with 35.92%. Both use the same sampling walk and thin-coordinate
width, with dimension equal to body count.

See the [fixed-3D instructions](example/kakeya/samp_only/hit_and_run/README.md)
and [dimension-equals-body-count instructions](example/kakeya/samp_only/hit_and_run/dimension_equals_bodies/README.md)
for further controls, validation details, and all figure names. The
[Kakeya index](example/kakeya/README.md) links the preserved sampler-only
and combined volume-and-sampling experiments.

## Different Usage

### Usage: Bitvector Counting
To counts the models of formulas over the theory of fixed-size bit-vectors (logics BV and UFBV), use this mode. The formula is eagerly bit-blasted to a model-preserving CNF in memory and counted with Arjun + ApproxMC in-process, so no intermediate file is written. This engine is auto-selected for bit-vector logics.


### Usage: Volume Computation
Computes the volume of the solution space of linear real arithmetic (QF_LRA) formulas — that is, the measure of the region of real assignments that satisfy the formula. The polytopes are estimated with the volesti random-walk sampler. This engine is auto-selected for QF_LRA inputs with no projection variables.

### Usage: Projected Counting
Hybrid formulas mix discrete and continuous variables. Engine 2 handles them with projection-based approximate counting (PACT): it counts the satisfying assignments to the discrete (bit-vector / Boolean) projection variables, while the continuous part is discharged by the SMT solver.

This engine is auto-selected for non-BV, non-LRA logics that have bit-vector or Boolean projection variables. 

### Usage: Skolem Function Counting

To count the interpretations (Skolem functions) of an uninterpreted function that satisfy a formula, use the SkolemFC mode. It is intended for QF_UFBV formulas that constrain the behaviour of an uninterpreted function.

### Example
```
$ ./ttc example/test.smt2
...
s mc 1
...
```



## Issues, questions, bugs, etc.
Please [open an issue](https://github.com/meelgroup/ttc/issues/new) for
support or bug reports. Please get in touch with Arijit Shaw or Kuldeep Meel in case you want a new feature / theory to be covered.

## How to Cite
If you use TTC in your research, please cite the following:


- **Efficient Volume Computation for SMT Formulas**

  Arijit Shaw, Uddalok Sarkar, and Kuldeep S. Meel

  *In Proceedings of Knowledge Representation and Reasoning (KR), Nov 2025*

- **Approximate SMT Counting Beyond Discrete Domains**

  Arijit Shaw, and Kuldeep S. Meel

  *In Proceedings of Design Automation Conference (DAC), Jun 2025*

