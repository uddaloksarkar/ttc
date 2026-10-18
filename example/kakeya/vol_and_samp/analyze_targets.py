#!/usr/bin/env python3
"""Analyze completed target runs using ordinary body-volume estimates."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "samp_only"))
from analyze_targets import analyze_category  # noqa: E402

if __name__ == "__main__":
    analyze_category(HERE)
