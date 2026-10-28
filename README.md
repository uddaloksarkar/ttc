# TTC: Toolbox for Theory Counting

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**ttc**, the Toolbox for Theory Counting, is a tool for measuring the solution space of Satisfiability Modulo Theories (SMT) formulas. It supports counting over a number of theories and their combinations, including bit-vectors and uninterpreted functions, as well as volume computation for the theory of linear real arithmetic.

## Documentation

A documentation on usage and building is available at https://arijitsh.github.io/docs-ttc/

## Reproduce the Kakeya benchmarks and figures

See the [Kakeya README](example/kakeya/README.md) for the sampler-only
benchmark inventory, exact volumes, generation, and direct TTC examples.

From the repository root, regenerate and verify the hit-and-run inputs with:

```bash
python3 example/kakeya/generate_benchmarks.py
```

Current inputs live in `example/kakeya/hit_and_run/benchmarks/` and
`example/kakeya/ball_walk/benchmarks/`.

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

