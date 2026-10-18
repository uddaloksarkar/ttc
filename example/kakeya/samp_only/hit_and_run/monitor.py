#!/usr/bin/env python3
"""Refresh RDHR result tables and figures as background TTC runs finish."""

import subprocess
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def main() -> None:
    previous = None
    while True:
        subprocess.run(["python3", str(HERE / "sweep.py"), "analyze"], check=True)
        state = (HERE / "results.csv").read_bytes()
        if state != previous:
            subprocess.run(["python3", str(HERE / "plot.py")], check=True)
            previous = state
            print(f"{time.strftime('%Y-%m-%d %H:%M:%S')} refreshed RDHR figures",
                  flush=True)
        time.sleep(60)


if __name__ == "__main__":
    main()
