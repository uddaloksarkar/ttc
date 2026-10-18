#!/usr/bin/env python3
"""Run the matching target suite with volesti body-volume estimates."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "samp_only"))
from run_targets import run_category  # noqa: E402

if __name__ == "__main__":
    run_category(HERE)
