#!/usr/bin/env python3
"""Queue matrix cases not assigned to the existing shorter-walk runners."""

from __future__ import annotations

import argparse
import csv
import shutil
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import probe

HERE = Path(__file__).resolve().parent


def planned_cases() -> list[tuple[int, int]]:
    high = probe.TARGET_WALKS[13:]
    cases = [(count, walk) for count in (2, 4, 8, 16) for walk in high]
    cases += [(32, walk) for walk in probe.TARGET_WALKS[2:]]
    cases.sort(key=lambda case: case[1] * case[0] ** 2)
    # Start one million-step dimension-matched candidate early.
    cases.remove((8, 1000000))
    return [(8, 1000000)] + cases


def active_rdhr_ttc() -> int:
    total = 0
    for process in Path("/proc").iterdir():
        if not process.name.isdigit():
            continue
        try:
            argv = [part.decode(errors="replace") for part in
                    (process / "cmdline").read_bytes().split(b"\0") if part]
        except (OSError, PermissionError):
            continue
        if argv and Path(argv[0]).name == "ttc" and "--sampler" in argv:
            total += argv[argv.index("--sampler") + 1] == "rdhr"
    return total


def archive_partial(count: int, walk: int) -> None:
    source = probe.log_for(count, walk)
    if not source.exists() or probe.parse_log(count, walk)["run_status"] == "0":
        return
    archive = probe.LOGS / "attempts"
    archive.mkdir(exist_ok=True)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    target = archive / f"{source.stem}_before_resume_{stamp}.log"
    shutil.copy2(source, target)
    print(f"archived {source.name} -> {target.name}", flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-active", type=int, default=18)
    parser.add_argument("--workers", type=int, default=3)
    parser.add_argument("--timeout-seconds", type=int, default=1209600)
    args = parser.parse_args()
    if args.max_active < 1 or args.workers < 1:
        raise ValueError("active and worker limits must be positive")
    pending = planned_cases()
    running = {}
    failures = []
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        while pending or running:
            for future, case in list(running.items()):
                if not future.done():
                    continue
                del running[future]
                try:
                    _, _, status, elapsed = future.result()
                    validated = probe.parse_log(*case)["run_status"]
                    print(f"finished m={case[0]} W={case[1]} wrapper={status} "
                          f"validated={validated} seconds={elapsed:.1f}", flush=True)
                    probe.analyze()
                    if validated != "0":
                        failures.append((case, status, validated))
                except Exception as error:
                    failures.append((case, "exception", str(error)))
                    print(f"error m={case[0]} W={case[1]} {error}", flush=True)
            if pending and len(running) < args.workers and active_rdhr_ttc() < args.max_active:
                case = pending.pop(0)
                if probe.parse_log(*case)["run_status"] == "0":
                    print(f"existing m={case[0]} W={case[1]}", flush=True)
                else:
                    archive_partial(*case)
                    future = pool.submit(probe.execute, *case, args.timeout_seconds)
                    running[future] = case
                    print(f"started m={case[0]} W={case[1]}", flush=True)
                    time.sleep(1)  # Let the TTC child appear before another capacity check.
            else:
                time.sleep(30)
    target = HERE / "queue_failures.csv"
    with target.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(("bodies", "walk_length", "wrapper_status", "validated_status"))
        writer.writerows((case[0], case[1], wrapper, validated)
                         for case, wrapper, validated in failures)
    print(f"queue finished; unresolved={len(failures)}", flush=True)


if __name__ == "__main__":
    main()
