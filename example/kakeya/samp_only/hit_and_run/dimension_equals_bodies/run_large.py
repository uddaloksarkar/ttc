#!/usr/bin/env python3
"""Run a staged 64D or 100D exact-volume RDHR TTC benchmark."""

import argparse
import csv
from pathlib import Path

import probe

HERE = Path(__file__).resolve().parent
LARGE = HERE / "large_counts"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("count", type=int, choices=(64, 100))
    parser.add_argument("--walk", type=int, default=5)
    parser.add_argument("--timeout-seconds", type=int, default=1209600)
    args = parser.parse_args()
    if args.walk < 1:
        raise ValueError("walk must be positive")
    probe.BENCHMARKS = LARGE / "benchmarks"
    probe.LOGS = LARGE / "logs"
    probe.LOGS.mkdir(exist_ok=True)
    path = probe.path_for(args.count)
    if not path.exists():
        raise FileNotFoundError(f"run prepare_large.py first: {path}")
    probe.verify_written(path, args.count)
    result = probe.execute(args.count, args.walk, args.timeout_seconds)
    row = probe.parse_log(args.count, args.walk)
    target = LARGE / f"result_m{args.count:04d}_walk{args.walk:07d}.csv"
    with target.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    print(result, row["run_status"], row["relative_error"], flush=True)


if __name__ == "__main__":
    main()
