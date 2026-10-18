#!/usr/bin/env python3
"""Refresh validated fixed-3D target tables and figures as TTC runs finish."""

import csv
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def validated_walks() -> tuple[int, ...]:
    with (HERE / "target_status.csv").open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    return tuple(int(row["walk_length"]) for row in rows
                 if row["run_status"] == "0")


def main() -> None:
    previous = None
    while True:
        subprocess.run([sys.executable, str(HERE / "sweep.py"), "analyze"],
                       check=True, stdout=subprocess.DEVNULL)
        current = validated_walks()
        if current != previous:
            subprocess.run([sys.executable, str(HERE / "plot.py")], check=True)
            print(f"refreshed fixed-3D figures: {len(current)}/20 walks", flush=True)
            previous = current
        if len(current) == 20:
            print("all target walks validated", flush=True)
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
