#!/usr/bin/env python3
"""Reproduce saved non-exact-volume TTC runs or extend their walk sweep."""

import argparse
import csv
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
OBSERVED = {(8,1): [10,20,50,100,500],
            (8,50): [10,20,50,100], (10,50): [10,20,50,100],
            (12,50): [10,20,50], (16,50): [10,20,50], (20,50): [10]}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--walks", nargs="+", type=int,
                        help="Override recorded walk lengths for selected cases")
    parser.add_argument("--names", nargs="*", help="Run only these benchmark names")
    parser.add_argument("--missing", action="store_true")
    parser.add_argument("--timeout", type=int, default=1800)
    args = parser.parse_args()
    if args.walks is not None and any(w < 1 for w in args.walks):
        raise SystemExit("Sampling walks must be positive")
    rows = list(csv.DictReader((HERE / "manifest.csv").open(newline="")))
    if args.names:
        rows = [r for r in rows if r["name"] in args.names]
    for row in rows:
        walks = args.walks or OBSERVED[(int(row["bodies"]), int(row["volume_walk_length"]))]
        for walk in walks:
            log = HERE / "benchmarks" / f"{row['name']}_samp{walk}.log"
            if args.missing and log.exists() and "run_status=0" in log.read_text():
                continue
            cmd = [str(PROJECT / "build" / "ttc"), str(HERE / row["path"]),
                   "--seed", row["seed"], "--epsilon", "0.2", "--delta", "0.2",
                   "--walklen-vol", row["volume_walk_length"],
                   "--walklen-samp", str(walk)]
            start = time.monotonic()
            try:
                result = subprocess.run(cmd, capture_output=True, text=True,
                                        timeout=args.timeout)
                log.write_text(result.stdout + result.stderr
                               + f"\nrun_status={result.returncode}\n"
                               + f"wall_seconds={time.monotonic()-start:.2f}\n")
                print(f"{row['name']} walk={walk}: exit={result.returncode}", flush=True)
            except subprocess.TimeoutExpired as exc:
                log.write_text((exc.stdout or b"").decode(errors="replace")
                               + (exc.stderr or b"").decode(errors="replace")
                               + f"\nrun_status=timeout\nwall_seconds={time.monotonic()-start:.2f}\n")
                print(f"{row['name']} walk={walk}: timeout", flush=True)


if __name__ == "__main__":
    main()
