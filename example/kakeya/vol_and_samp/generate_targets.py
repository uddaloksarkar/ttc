#!/usr/bin/env python3
"""Create matching target geometries using the ordinary body-volume oracle."""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "samp_only"))
from generate_targets import generate_category  # noqa: E402

if __name__ == "__main__":
    generate_category(HERE, "VOLANDSAMP")
