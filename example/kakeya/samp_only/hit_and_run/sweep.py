#!/usr/bin/env python3
"""Generate, run, and validate exact-box random-direction hit-and-run TTC cases."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import subprocess
import sys
import tempfile
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[3]
sys.path.insert(0, str(HERE))
from generate import exact_union_quadrant, make_formula  # noqa: E402
from verify_thickness_geometry import verify  # noqa: E402

BENCHMARKS = HERE / "benchmarks"
LOGS = HERE / "logs"
TTC = PROJECT / "build" / "ttc"
SPREAD = 100
SHORT = 1
EPSILON = 0.2
DELTA = 0.2
SEED = 42
TARGET_WALKS = (5, 10, 18, 34, 65, 124, 236, 449, 853, 1622,
                3083, 5861, 11143, 21183, 40271, 76559, 145544,
                276692, 526015, 1000000)
MATRIX_COUNTS = (2, 4, 8, 16, 32)
MATRIX_DENOMINATOR = 5000000


def numbers(value: str) -> list[int]:
    result = [int(item) for item in value.split(",")]
    if not result or any(item < 1 for item in result):
        raise argparse.ArgumentTypeError("expected positive comma-separated integers")
    return result


def stem(count: int, denominator: int) -> str:
    return f"RDHR_box_m{count:04d}_d{denominator:08d}_s{SHORT:04d}"


def input_path(count: int, denominator: int) -> Path:
    return BENCHMARKS / f"{stem(count, denominator)}.smt2"


def log_path(count: int, denominator: int, walk: int) -> Path:
    return LOGS / f"{stem(count, denominator)}_walk{walk:07d}.log"


def geometry(count: int, denominator: int) -> tuple[Fraction, Fraction]:
    return (Fraction(8 * SPREAD * SHORT, denominator),
            Fraction(8, denominator) * exact_union_quadrant(count, SPREAD, SHORT))


def prepare(counts: list[int], denominators: list[int]) -> None:
    BENCHMARKS.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    manifest = HERE / "geometry.csv"
    existing = ({row["name"]: row for row in csv.DictReader(manifest.open())}
                if manifest.exists() else {})
    rows = []
    for denominator in sorted(set(denominators)):
        for count in sorted(set(counts)):
            if count < 2:
                raise ValueError("at least two bodies are needed")
            path = input_path(count, denominator)
            formula = make_formula(count, 5, denominator, SPREAD, SHORT)
            if path.exists():
                if path.read_text() != formula:
                    raise ValueError(f"existing benchmark differs: {path}")
            else:
                path.write_text(formula)
            # The independent verifier checks a historical filename grammar.
            # Feed it the exact bytes written for this RDHR case under that
            # grammar, without changing or duplicating the benchmark itself.
            with tempfile.TemporaryDirectory() as temporary:
                verification_path = (Path(temporary) /
                                     f"SAMPONLY_thickness_m{count:04d}_walk0000005_"
                                     f"d{denominator:08d}.smt2")
                verification_path.write_bytes(path.read_bytes())
                parsed = verify(verification_path)
            body, union = geometry(count, denominator)
            if (int(parsed["bodies"]) != count
                    or Fraction(parsed["exact_body_volume"]) != body
                    or Fraction(parsed["exact_union_volume"]) != union):
                raise ValueError(f"independent geometry check failed: {path}")
            rows.append(dict(name=stem(count, denominator), bodies=count,
                             dimension=3, thin_denominator=denominator,
                             exact_body_volume=str(body), exact_union_volume=str(union),
                             benchmark=str(path.relative_to(HERE))))
    existing.update({row["name"]: row for row in rows})
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))
    print(f"Prepared and independently verified {len(rows)} box unions", flush=True)


def command(path: Path, walk: int) -> list[str]:
    return [str(TTC), str(path), "--seed", str(SEED),
            "--epsilon", str(EPSILON), "--delta", str(DELTA),
            "--sampler", "rdhr", "--walklen-samp", str(walk),
            "--exact-box-volumes"]


def is_live_case(path: Path, walk: int) -> bool:
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            argv = [part.decode(errors="replace") for part in
                    (process / "cmdline").read_bytes().split(b"\0") if part]
        except (OSError, PermissionError):
            continue
        if (argv and Path(argv[0]).name == "ttc" and str(path) in argv
                and "--sampler" in argv and "--walklen-samp" in argv
                and argv[argv.index("--sampler") + 1] == "rdhr"
                and argv[argv.index("--walklen-samp") + 1] == str(walk)):
            return True
    return False


def execute(case: tuple[int, int, int], timeout_seconds: float) -> tuple[tuple[int, int, int], str, float]:
    count, denominator, walk = case
    path = input_path(count, denominator)
    log = log_path(count, denominator, walk)
    if not path.exists():
        raise FileNotFoundError(path)
    if log.exists():
        if parse_log(log, count, denominator, walk)["run_status"] == "0":
            return case, "existing", 0.0
        if is_live_case(path, walk):
            raise RuntimeError(f"TTC is already running this case: {path}, walk {walk}")
        if log.stat().st_size:
            attempts = LOGS / "attempts"
            attempts.mkdir(exist_ok=True)
            shutil.copy2(log, attempts / f"{log.stem}_before_resume_{time.time_ns()}.log")
    start = time.monotonic()
    with log.open("w") as output:
        try:
            process = subprocess.run(command(path, walk), stdout=output,
                                     stderr=subprocess.STDOUT, timeout=timeout_seconds)
            status = str(process.returncode)
        except subprocess.TimeoutExpired:
            status = "timeout"
        elapsed = time.monotonic() - start
        output.write(f"\nrun_status={status}\nwall_seconds={elapsed:.2f}\n")
    return case, status, elapsed


def run(counts: list[int], denominators: list[int], walks: list[int], jobs: int,
        timeout_seconds: float) -> None:
    cases = [(count, denominator, walk)
             for walk in sorted(set(walks))
             for count in sorted(set(counts))
             for denominator in sorted(set(denominators))]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(execute, case, timeout_seconds) for case in cases]
        for future in as_completed(futures):
            case, status, elapsed = future.result()
            print(f"m={case[0]} D={case[1]} W={case[2]} status={status} "
                  f"seconds={elapsed:.2f}", flush=True)
            analyze()
    analyze()


LOG_NAME = re.compile(r"RDHR_box_m(\d+)_d(\d+)_s0001_walk(\d+)\.log$")
BODY_ROW = re.compile(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$", re.M)


def parse_log(path: Path, count: int, denominator: int, walk: int) -> dict[str, object]:
    raw = path.read_text(errors="replace")
    marker = re.search(r"^run_status=(\S+)\s*$", raw, re.M)
    status = marker.group(1) if marker else "running_or_interrupted"
    body, union = geometry(count, denominator)
    estimates = re.findall(r"^s volume ([0-9.eE+-]+)\s*$", raw, re.M)
    records = BODY_ROW.findall(raw)
    if status == "0":
        headers = (f"c sample walk length: {walk}", "c sampler walk: rdhr",
                   "c body volume oracle: exact axis-aligned box",
                   f"c polytopes: {count}", "c real variables: 3", "c exit 0")
        tolerance = max(5.1e-7, 1e-5 * float(body))
        if (input_path(count, denominator).read_text()
                != make_formula(count, 5, denominator, SPREAD, SHORT)
                or any(header not in raw for header in headers) or len(estimates) != 1
                or [int(item[0]) for item in records] != list(range(1, count + 1))
                or any(float(item[1]) <= 0 or
                       abs(float(item[1]) - float(body)) > tolerance
                       for item in records)):
            status = "invalid_output"
    estimate = float(estimates[0]) if status == "0" else None
    error = abs(estimate / float(union) - 1) if estimate is not None else None
    first_samples = int(records[0][3]) if records else None
    deleted_second = int(records[1][2]) if count == 2 and len(records) == 2 else None
    deletion_ratio = (deleted_second / first_samples
                      if status == "0" and first_samples and deleted_second is not None
                      else None)
    wall = re.search(r"^wall_seconds=([0-9.]+)", raw, re.M)
    return dict(name=path.stem, bodies=count, dimension=3,
                thin_denominator=denominator, walk_length=walk,
                exact_body_volume=str(body), exact_union_volume=str(union),
                estimated_union_volume="" if estimate is None else f"{estimate:.17g}",
                relative_error="" if error is None else f"{error:.17g}",
                fails_20pct="" if error is None else int(error > EPSILON),
                run_status=status, body_rows=len(records),
                first_body_samples="" if first_samples is None else first_samples,
                second_body_deleted="" if deleted_second is None else deleted_second,
                deletion_ratio="" if deletion_ratio is None else f"{deletion_ratio:.17g}",
                true_overlap_fraction="" if count != 2 else f"{SHORT / SPREAD:.17g}",
                wall_seconds=wall.group(1) if wall else "",
                benchmark=str(input_path(count, denominator).relative_to(HERE)),
                log=str(path.relative_to(HERE)))


def analyze() -> None:
    rows = []
    for path in sorted(LOGS.glob("*.log")):
        match = LOG_NAME.fullmatch(path.name)
        if not match:
            continue
        count, denominator, walk = map(int, match.groups())
        rows.append(parse_log(path, count, denominator, walk))
    if not rows:
        return
    target = HERE / "results.csv"
    temporary = HERE / f"results.csv.tmp.{os.getpid()}.{threading.get_ident()}"
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(target)
    complete = [r for r in rows if r["run_status"] == "0"]
    indexed = {(int(r["bodies"]), int(r["thin_denominator"]),
                int(r["walk_length"])): r for r in rows}
    live_walks = set()
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            argv = [part.decode(errors="replace") for part in
                    (process / "cmdline").read_bytes().split(b"\0") if part]
        except (OSError, PermissionError):
            continue
        if (argv and Path(argv[0]).name == "ttc"
                and str(input_path(2, 5000000)) in argv
                and "--sampler" in argv and argv[argv.index("--sampler") + 1] == "rdhr"
                and "--walklen-samp" in argv):
            live_walks.add(int(argv[argv.index("--walklen-samp") + 1]))
    targets = []
    for walk in TARGET_WALKS:
        result = indexed.get((2, 5000000, walk))
        status = result["run_status"] if result else "pending"
        if status != "0" and walk in live_walks:
            status = "running"
        elif status == "running_or_interrupted":
            status = "interrupted_or_unreported"
        targets.append(dict(walk_length=walk, bodies=2,
                            thin_denominator=5000000,
                            run_status=status,
                            relative_error=result["relative_error"] if result else "",
                            fails_20pct=result["fails_20pct"] if result else "",
                            log=str(log_path(2, 5000000, walk).relative_to(HERE))))
    target_temporary = HERE / (f"target_status.csv.tmp.{os.getpid()}."
                               f"{threading.get_ident()}")
    with target_temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=targets[0])
        writer.writeheader()
        writer.writerows(targets)
    target_temporary.replace(HERE / "target_status.csv")
    matrix = []
    for count in MATRIX_COUNTS:
        for walk in TARGET_WALKS:
            result = indexed.get((count, MATRIX_DENOMINATOR, walk))
            path = input_path(count, MATRIX_DENOMINATOR)
            status = result["run_status"] if result else "pending"
            if status != "0" and is_live_case(path, walk):
                status = "running"
            elif status == "running_or_interrupted":
                status = "interrupted_or_unreported"
            matrix.append(dict(benchmark=stem(count, MATRIX_DENOMINATOR), bodies=count,
                               thin_denominator=MATRIX_DENOMINATOR, walk_length=walk,
                               run_status=status,
                               relative_error=result["relative_error"] if result else "",
                               fails_20pct=result["fails_20pct"] if result else "",
                               input=str(path.relative_to(HERE)),
                               log=str(log_path(count, MATRIX_DENOMINATOR, walk).relative_to(HERE))))
    matrix_target = HERE / "benchmark_matrix.csv"
    matrix_temporary = HERE / f"benchmark_matrix.csv.tmp.{os.getpid()}.{threading.get_ident()}"
    with matrix_temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=matrix[0])
        writer.writeheader()
        writer.writerows(matrix)
    matrix_temporary.replace(matrix_target)
    print(f"Validated {len(complete)} of {len(rows)} RDHR runs", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run", "analyze"))
    parser.add_argument("--counts", type=numbers, default=[2, 4, 8, 16])
    parser.add_argument("--denominators", type=numbers, default=[1, 10, 100, 5000000])
    parser.add_argument("--walks", type=numbers, default=[5, 10, 100, 1000])
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--timeout-seconds", type=float, default=172800)
    args = parser.parse_args()
    if args.action == "prepare":
        prepare(args.counts, args.denominators)
    elif args.action == "run":
        run(args.counts, args.denominators, args.walks, args.jobs,
            args.timeout_seconds)
    else:
        analyze()


if __name__ == "__main__":
    main()
