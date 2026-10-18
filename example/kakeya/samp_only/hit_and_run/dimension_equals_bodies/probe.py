#!/usr/bin/env python3
"""Pilot exact-volume RDHR benchmarks with dimension equal to body count."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARENT = HERE.parent
PROJECT = HERE.parents[4]
sys.path.insert(0, str(PARENT))
from generate import exact_union_quadrant, smt_negative, smt_positive  # noqa: E402

TTC = PROJECT / "build" / "ttc"
COUNTS = (2, 4, 8, 16, 32)
DENOMINATOR = 5000000
TARGET_WALKS = (5, 10, 18, 34, 65, 124, 236, 449, 853, 1622,
                3083, 5861, 11143, 21183, 40271, 76559, 145544,
                276692, 526015, 1000000)
BENCHMARKS = HERE / "benchmarks"
LOGS = HERE / "logs"


def path_for(count: int) -> Path:
    return BENCHMARKS / f"RDHR_dimeqbody_m{count:04d}_d{DENOMINATOR:08d}.smt2"


def log_for(count: int, walk: int) -> Path:
    return LOGS / f"RDHR_dimeqbody_m{count:04d}_walk{walk:07d}.log"


def common_factor(count: int) -> Fraction:
    return (Fraction(1) if count == 2 else
            Fraction(2, DENOMINATOR) * 2 ** (count - 3))


def exact_volumes(count: int) -> tuple[Fraction, Fraction]:
    factor = common_factor(count)
    return (400 * factor, 4 * exact_union_quadrant(count) * factor)


def make_formula(count: int) -> str:
    lines = ["(set-logic QF_LRA)"]
    lines += [f"(declare-fun x{j} () Real)" for j in range(count)]
    lines.append("(assert (or")
    for i in range(count):
        x = Fraction(count - 1 + 99 * i, count - 1)
        y = Fraction(100) / x
        radii = [x, y] + ([Fraction(1, DENOMINATOR)] +
                           [Fraction(1)] * (count - 3) if count >= 3 else [])
        lines.append("  (and")
        for j, radius in enumerate(radii):
            lines.append(f"    (<= {smt_negative(radius)} x{j})")
            lines.append(f"    (<= x{j} {smt_positive(radius)})")
        lines.append("  )")
    return "\n".join(lines + ["))", "(check-sat)"]) + "\n"


def parse_fraction(token: str) -> Fraction:
    match = re.fullmatch(r"\(/\s+(\d+)\s+(\d+)\)", token.strip())
    return Fraction(int(match[1]), int(match[2])) if match else Fraction(token)


def verify_written(path: Path, count: int) -> tuple[Fraction, Fraction]:
    raw = path.read_text()
    declared = re.findall(r"^\(declare-fun x(\d+) \(\) Real\)$", raw, re.M)
    if declared != [str(j) for j in range(count)]:
        raise ValueError(f"variable count differs: {path}")
    bounds = []
    for j in range(count):
        upper = [parse_fraction(x) for x in re.findall(
            rf"^\s*\(<= x{j} (.+)\)$", raw, re.M)]
        lower = [parse_fraction(x) for x in re.findall(
            rf"^\s*\(<= \(- (.+)\) x{j}\)$", raw, re.M)]
        if len(upper) != count or upper != lower:
            raise ValueError(f"body bounds differ: {path}, coordinate {j}")
        bounds.append(upper)
    xs, ys = bounds[:2]
    if xs != [Fraction(count - 1 + 99 * i, count - 1)
              for i in range(count)] or ys != [Fraction(100) / x for x in xs]:
        raise ValueError(f"staircase geometry differs: {path}")
    for j in range(2, count):
        expected = Fraction(1, DENOMINATOR) if j == 2 else Fraction(1)
        if bounds[j] != [expected] * count:
            raise ValueError(f"common coordinates differ: {path}")
    factor = Fraction(1)
    for j in range(2, count):
        factor *= 2 * bounds[j][0]
    body = 4 * xs[0] * ys[0] * factor
    endpoints = [Fraction(0)] + sorted(set(xs))
    quadrant = sum((right - left) * max(y for x, y in zip(xs, ys) if x >= right)
                   for left, right in zip(endpoints, endpoints[1:]))
    union = 4 * quadrant * factor
    if (body, union) != exact_volumes(count):
        raise ValueError(f"exact integral differs: {path}")
    return body, union


def prepare() -> None:
    BENCHMARKS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    rows = []
    for count in COUNTS:
        path = path_for(count)
        formula = make_formula(count)
        if path.exists() and path.read_text() != formula:
            raise ValueError(f"existing input differs: {path}")
        path.write_text(formula)
        body, union = verify_written(path, count)
        rows.append(dict(benchmark=path.name, bodies=count, dimension=count,
                         exact_body_volume=str(body), exact_union_volume=str(union)))
    with (HERE / "geometry.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Verified {len(rows)} dimension=body-count SMT2 geometries", flush=True)


def parse_log(count: int, walk: int) -> dict[str, object]:
    path = log_for(count, walk)
    raw = path.read_text(errors="replace") if path.exists() else ""
    marker = re.search(r"^run_status=(\S+)\s*$", raw, re.M)
    # A wrapper can be terminated after TTC has printed its complete normal
    # exit footer but before it appends run_status. The TTC footer plus all
    # other checks below is sufficient evidence of a completed TTC run.
    footer_complete = bool(re.search(r"^c exit 0\s*$", raw, re.M))
    status = (marker[1] if marker else
              ("0" if footer_complete else "running_or_interrupted" if raw else "pending"))
    body, union = exact_volumes(count)
    estimates = re.findall(r"^s volume ([0-9.eE+-]+)\s*$", raw, re.M)
    records = re.findall(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$", raw, re.M)
    if status == "0":
        headers = (f"c real variables: {count}", f"c polytopes: {count}",
                   "c sampler walk: rdhr", f"c sample walk length: {walk}",
                   "c body volume oracle: exact axis-aligned box", "c exit 0")
        tolerance = max(5.1e-7, 1e-6 * float(body))
        if (any(header not in raw for header in headers)
                or [int(record[0]) for record in records] != list(range(1, count + 1))
                or any(abs(float(record[1]) - float(body)) > tolerance
                       for record in records)
                or len(estimates) != 1 or float(estimates[0]) <= 0
                or path_for(count).read_text() != make_formula(count)):
            status = "invalid_output"
    estimate = float(estimates[0]) if status == "0" else None
    error = abs(estimate / float(union) - 1) if estimate is not None else None
    return dict(bodies=count, dimension=count, walk_length=walk,
                exact_body_volume=str(body), exact_union_volume=str(union),
                estimated_union_volume="" if estimate is None else f"{estimate:.17g}",
                relative_error="" if error is None else f"{error:.17g}",
                fails_20pct="" if error is None else int(error > 0.2),
                run_status=status, body_rows=len(records),
                completion_evidence=("wrapper_status" if marker else
                                     "ttc_exit_footer" if footer_complete else ""),
                benchmark=str(path_for(count).relative_to(HERE)),
                log=str(path.relative_to(HERE)))


def analyze() -> None:
    rows = [parse_log(count, walk) for count in COUNTS for walk in TARGET_WALKS]
    target = HERE / "results.csv"
    temporary = HERE / f"results.csv.tmp.{os.getpid()}.{threading.get_ident()}"
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(target)
    print(f"Validated {sum(row['run_status']=='0' for row in rows)}/100 TTC runs", flush=True)


def execute(count: int, walk: int, timeout: int) -> tuple[int, int, str, float]:
    existing = parse_log(count, walk)
    if existing["run_status"] == "0":
        return count, walk, "existing", 0.0
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            argv = [part.decode(errors="replace") for part in
                    (process / "cmdline").read_bytes().split(b"\0") if part]
        except (OSError, PermissionError):
            continue
        if (argv and Path(argv[0]).name == "ttc"
                and str(path_for(count)) in argv
                and "--walklen-samp" in argv
                and argv[argv.index("--walklen-samp") + 1] == str(walk)):
            return count, walk, "already_running", 0.0
    log = log_for(count, walk)
    if log.exists() and log.stat().st_size:
        attempts = LOGS / "attempts"
        attempts.mkdir(exist_ok=True)
        shutil.copy2(log, attempts / f"{log.stem}_before_resume_{time.time_ns()}.log")
    start = time.monotonic()
    with log.open("w") as stream:
        try:
            completed = subprocess.run([str(TTC), str(path_for(count)), "--seed", "42",
                                        "--epsilon", "0.2", "--delta", "0.2",
                                        "--sampler", "rdhr", "--walklen-samp", str(walk),
                                        "--exact-box-volumes"],
                                       stdout=stream, stderr=subprocess.STDOUT,
                                       timeout=timeout)
            status = str(completed.returncode)
        except subprocess.TimeoutExpired:
            status = "timeout"
        elapsed = time.monotonic() - start
        stream.write(f"\nrun_status={status}\nwall_seconds={elapsed:.2f}\n")
    return count, walk, status, elapsed


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "analyze"))
    parser.add_argument("--counts", default="2,4,8,16,32")
    parser.add_argument("--walks", default="5")
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--timeout-seconds", type=int, default=1209600)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "analyze":
        analyze()
    else:
        counts = [int(x) for x in args.counts.split(",")]
        walks = [int(x) for x in args.walks.split(",")]
        if not set(counts) <= set(COUNTS) or not set(walks) <= set(TARGET_WALKS):
            raise ValueError("unsupported count or walk")
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            futures = [pool.submit(execute, count, walk, args.timeout_seconds)
                       for walk in walks for count in counts]
            for future in as_completed(futures):
                count, walk, status, elapsed = future.result()
                print(f"m=d={count} W={walk} status={status} seconds={elapsed:.2f}",
                      flush=True)
                analyze()


if __name__ == "__main__":
    main()
