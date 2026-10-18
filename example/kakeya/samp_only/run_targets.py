#!/usr/bin/env python3
"""Run selected members of the target suite and preserve complete TTC logs."""

import argparse
import csv
import math
import re
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]


def valid_completed_log(raw: str, expected_bodies: int) -> bool:
    if not re.search(r"^run_status=0$", raw, re.M):
        return False
    result = re.search(r"^s volume ([\d.eE+-]+)$", raw, re.M)
    if result is None:
        return False
    try:
        returned = float(result.group(1))
    except ValueError:
        return False
    if not math.isfinite(returned) or returned <= 0:
        return False
    samples = re.findall(
        r"^c\s+\d+\s+[\d.eE+-]+\s+\d+\s+(\d+)$", raw, re.M
    )
    return len(samples) == expected_bodies and all(int(n) > 0 for n in samples)


def run_category(root: Path) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--targets", nargs="+", type=int)
    parser.add_argument("--replicas", nargs="+", type=int)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--timeout", type=int, default=1800)
    parser.add_argument("--missing", action="store_true")
    args = parser.parse_args()
    rows = list(csv.DictReader((root / "target_manifest.csv").open(newline="")))
    if args.targets:
        rows = [r for r in rows if int(r["target_walk_length"]) in args.targets]
    if args.replicas:
        rows = [r for r in rows if int(r["replica"]) in args.replicas]
    if args.limit is not None:
        rows = rows[:args.limit]
    for row in rows:
        log = root / "benchmarks" / f"{row['name']}.log"
        if args.missing and log.exists() and valid_completed_log(
            log.read_text(), int(row["bodies"])
        ):
            continue
        command = [str(PROJECT / "build" / "ttc"), str(root / row["path"]),
                   "--seed", row["seed"], "--epsilon", row["epsilon"],
                   "--delta", row["delta"], "--sampler", "ball",
                   "--walklen-samp", row["target_walk_length"]]
        if row["volume_walk_length"] == "exact":
            command.append("--exact-box-volumes")
        else:
            command.extend(["--walklen-vol", row["volume_walk_length"]])
        started = time.monotonic()
        try:
            completed = subprocess.run(command, text=True, capture_output=True,
                                       timeout=args.timeout)
            body = completed.stdout + completed.stderr
            status = str(completed.returncode)
        except subprocess.TimeoutExpired as exc:
            body = (exc.stdout or b"").decode(errors="replace")
            body += (exc.stderr or b"").decode(errors="replace")
            status = "timeout"
        elapsed = time.monotonic() - started
        log.write_text(body + f"\nrun_status={status}\nwall_seconds={elapsed:.2f}\n")
        print(f"{row['name']}: {status}, {elapsed:.1f}s", flush=True)


if __name__ == "__main__":
    run_category(HERE)
