#!/usr/bin/env python3
"""Refresh factored pilot figures when a new TTC log validates."""

import subprocess
import sys
import time
from pathlib import Path

from prepare_factored import COUNTS
from run_factored import parse

HERE = Path(__file__).resolve().parent


def completed() -> tuple[int, ...]:
    return tuple(count for count in COUNTS if parse(count, 5)["run_status"] == "0")


def main() -> None:
    previous = None
    while True:
        current = completed()
        if current != previous:
            subprocess.run([sys.executable, str(HERE / "plot_factored.py")], check=True)
            print(f"refreshed factored figures: validated {current}", flush=True)
            previous = current
        if len(current) == len(COUNTS):
            print("all factored walk-5 pilots validated", flush=True)
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
