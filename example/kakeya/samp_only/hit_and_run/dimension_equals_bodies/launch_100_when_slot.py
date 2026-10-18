#!/usr/bin/env python3
"""Start the staged 100-body pilot when fewer than 18 RDHR TTC runs are live."""

import os
import sys
import time
import csv
from pathlib import Path

from queue_remaining import active_rdhr_ttc

HERE = Path(__file__).resolve().parent
RESULT = HERE / "large_counts" / "result_m0100_walk0000005.csv"
RUNNER = HERE / "run_large.py"


def main() -> None:
    if RESULT.exists():
        with RESULT.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        if rows and rows[0]["run_status"] == "0":
            print("100-body pilot already has a completed result", flush=True)
            return
    while active_rdhr_ttc() >= 18:
        time.sleep(10)
    print("starting 100-body, 100-dimensional walk-5 pilot", flush=True)
    os.execv(sys.executable, [sys.executable, "-u", str(RUNNER), "100",
                              "--walk", "5", "--timeout-seconds", "1209600"])


if __name__ == "__main__":
    main()
