#!/usr/bin/env python3
"""Refresh validated dimension-matched TTC tables and figures."""

import csv
import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main():
    previous = None
    while True:
        subprocess.run(["python3", str(HERE / "probe.py"), "analyze"], check=True)
        completed = tuple((row["bodies"], row["walk_length"], row["relative_error"])
                          for row in csv.DictReader((HERE / "results.csv").open())
                          if row["run_status"] == "0")
        fatness = tuple((p.name, p.stat().st_mtime_ns)
                        for p in sorted(HERE.glob("fatness_d*/result.csv")))
        state = (completed, fatness)
        if state != previous:
            subprocess.run(["python3", str(HERE / "plot.py")], check=True)
            previous = state
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} refreshed figures",
                  flush=True)
        time.sleep(60)


if __name__ == "__main__":
    main()
