#!/usr/bin/env python3
"""Run a reproducible dimension-matched exact-volume TTC comparison."""

import argparse
import csv
from pathlib import Path

import probe


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("denominator", type=int)
    parser.add_argument("--walk", type=int, default=1000)
    parser.add_argument("--count", type=int, default=8)
    args = parser.parse_args()
    if args.denominator < 1:
        raise ValueError("denominator must be positive")
    count, walk = args.count, args.walk
    if walk < 1 or count not in probe.COUNTS:
        raise ValueError("unsupported walk or body count")
    dirname = (f"fatness_d{args.denominator:08d}" if walk == 1000 and count == 8
               else f"walk{walk:07d}_d{args.denominator:08d}" if count == 8
               else f"walk{walk:07d}_m{count:04d}_d{args.denominator:08d}")
    home = Path(__file__).resolve().parent / dirname
    probe.DENOMINATOR = args.denominator
    probe.BENCHMARKS = home / "benchmarks"
    probe.LOGS = home / "logs"
    probe.BENCHMARKS.mkdir(parents=True, exist_ok=True)
    probe.LOGS.mkdir(parents=True, exist_ok=True)
    benchmark = probe.path_for(count)
    formula = probe.make_formula(count)
    if benchmark.exists() and benchmark.read_text() != formula:
        raise ValueError(f"existing benchmark changed: {benchmark}")
    benchmark.write_text(formula)
    probe.verify_written(benchmark, count)
    result = probe.execute(count, walk, 1209600)
    row = probe.parse_log(count, walk)
    row["thin_denominator"] = args.denominator
    with (home / "result.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    print(result, row["run_status"], row["relative_error"], flush=True)


if __name__ == "__main__":
    main()
