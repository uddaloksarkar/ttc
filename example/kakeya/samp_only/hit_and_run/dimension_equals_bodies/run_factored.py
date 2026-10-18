#!/usr/bin/env python3
"""Run and validate an equivalent factored exact-volume RDHR TTC input."""

import argparse
import csv
import re
import subprocess
import time
from pathlib import Path

import probe
from prepare_factored import HERE, COUNTS, verify


def paths(count: int, walk: int, no_cdd: bool = False) -> tuple[Path, Path, Path]:
    stem = f"RDHR_factored_m{count:04d}_d{probe.DENOMINATOR:08d}"
    output = HERE / "nocdd" if no_cdd else HERE
    return (HERE / "benchmarks" / f"{stem}.smt2",
            output / "logs" / f"RDHR_factored_m{count:04d}_walk{walk:07d}.log",
            output / f"result_m{count:04d}_walk{walk:07d}.csv")


def parse(count: int, walk: int, no_cdd: bool = False) -> dict[str, object]:
    benchmark, log, _ = paths(count, walk, no_cdd)
    body, union = verify(benchmark, count)
    raw = log.read_text(errors="replace") if log.exists() else ""
    status_match = re.search(r"^run_status=(\S+)\s*$", raw, re.M)
    footer = bool(re.search(r"^c exit 0\s*$", raw, re.M))
    status = (status_match[1] if status_match else
              "0" if footer else "running_or_interrupted" if raw else "pending")
    volumes = re.findall(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$", raw, re.M)
    estimates = re.findall(r"^s volume ([0-9.eE+-]+)\s*$", raw, re.M)
    if status == "0":
        headers = (f"c real variables: {count}", f"c polytopes: {count}",
                   "c sampler walk: rdhr", f"c sample walk length: {walk}",
                   f"c cdd simplification: {'no' if no_cdd else 'yes'}",
                   "c body volume oracle: exact axis-aligned box", "c exit 0")
        tolerance = max(5.1e-7, 1e-6 * float(body))
        if (any(header not in raw for header in headers)
                or [int(row[0]) for row in volumes] != list(range(1, count + 1))
                or any(abs(float(row[1]) - float(body)) > tolerance
                       for row in volumes)
                or len(estimates) != 1 or float(estimates[0]) <= 0):
            status = "invalid_output"
    estimate = float(estimates[0]) if status == "0" else None
    error = abs(estimate / float(union) - 1) if estimate is not None else None
    return dict(bodies=count, dimension=count, walk_length=walk,
                exact_body_volume=str(body), exact_union_volume=str(union),
                estimated_union_volume="" if estimate is None else f"{estimate:.17g}",
                relative_error="" if error is None else f"{error:.17g}",
                fails_20pct="" if error is None else int(error > 0.2),
                run_status=status, body_rows=len(volumes),
                completion_evidence="wrapper_status" if status_match else
                "ttc_exit_footer" if footer else "",
                benchmark=str(benchmark.relative_to(HERE)),
                log=str(log.relative_to(HERE)))


def live_case(benchmark: Path, walk: int, no_cdd: bool = False) -> bool:
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            args = [part.decode(errors="replace") for part in
                    (process / "cmdline").read_bytes().split(b"\0") if part]
        except (OSError, PermissionError):
            continue
        if (args and Path(args[0]).name == "ttc" and str(benchmark) in args
                and "--walklen-samp" in args
                and args[args.index("--walklen-samp") + 1] == str(walk)
                and ("--no-cdd-simp" in args) == no_cdd):
            return True
    return False


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("count", type=int, choices=COUNTS)
    parser.add_argument("--walk", type=int, default=5)
    parser.add_argument("--timeout-seconds", type=int, default=1209600)
    parser.add_argument("--analyze-only", action="store_true")
    parser.add_argument("--no-cdd-simp", action="store_true")
    args = parser.parse_args()
    if args.walk < 1:
        raise ValueError("walk must be positive")
    benchmark, log, result = paths(args.count, args.walk, args.no_cdd_simp)
    verify(benchmark, args.count)
    log.parent.mkdir(parents=True, exist_ok=True)
    row = parse(args.count, args.walk, args.no_cdd_simp)
    if not args.analyze_only and row["run_status"] != "0":
        if live_case(benchmark, args.walk, args.no_cdd_simp):
            print("already running", flush=True)
            return
        if log.exists():
            attempts = log.parent / "attempts"
            attempts.mkdir(exist_ok=True)
            archive = attempts / f"{log.stem}_before_resume_{int(time.time())}.log"
            archive.write_bytes(log.read_bytes())
        start = time.monotonic()
        with log.open("w") as stream:
            try:
                command = [str(probe.TTC), str(benchmark),
                           "--seed", "42", "--epsilon", "0.2",
                           "--delta", "0.2", "--sampler", "rdhr",
                           "--walklen-samp", str(args.walk),
                           "--exact-box-volumes"]
                if args.no_cdd_simp:
                    command.append("--no-cdd-simp")
                complete = subprocess.run(command,
                                          stdout=stream, stderr=subprocess.STDOUT,
                                          timeout=args.timeout_seconds)
                status = str(complete.returncode)
            except subprocess.TimeoutExpired:
                status = "timeout"
            stream.write(f"\nrun_status={status}\n"
                         f"wall_seconds={time.monotonic() - start:.2f}\n")
        row = parse(args.count, args.walk, args.no_cdd_simp)
    with result.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=row)
        writer.writeheader()
        writer.writerow(row)
    print(row["run_status"], row["relative_error"], flush=True)


if __name__ == "__main__":
    main()
