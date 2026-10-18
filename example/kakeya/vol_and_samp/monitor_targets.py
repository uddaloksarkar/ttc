#!/usr/bin/env python3
"""Refresh the measured target figures as TTC result logs arrive."""

import csv
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "target_results.csv"
JOURNAL = HERE / "target_monitor.log"


def refresh():
    subprocess.run(["python3", str(HERE / "analyze_targets.py")], check=True,
                   stdout=subprocess.DEVNULL)
    with RESULTS.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    state = tuple(tuple(row.items()) for row in rows)
    complete = sum(row["run_status"] == "0" for row in rows)
    return state, complete, len(rows)


def main():
    previous = None
    while True:
        state, complete, total = refresh()
        if state != previous:
            subprocess.run(["python3", str(HERE / "plot_targets.py")], check=True)
            stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")
            with JOURNAL.open("a") as stream:
                stream.write(f"{stamp} {complete}/{total} valid TTC results; "
                             "PDF/PGF/PNG figures refreshed\n")
            previous = state
        if complete == total:
            return
        time.sleep(60)


if __name__ == "__main__":
    main()
