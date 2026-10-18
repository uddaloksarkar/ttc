#!/usr/bin/env python3
"""Plot measured ordinary-volume TTC target failures."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "samp_only"))
from plot_targets import plot_category  # noqa: E402

if __name__ == "__main__":
    plot_category(HERE)
