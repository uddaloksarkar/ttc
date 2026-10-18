#!/usr/bin/env python3
"""Regenerate the saved two-sided box unions with analytic true volumes."""

import csv
import sys
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from sea_urchin_smt_generator import bounds_for_body, write_union  # noqa: E402


def main() -> None:
    rows = list(csv.DictReader((HERE / "manifest.csv").open(newline="")))
    for row in rows:
        n = int(row["dimension"])
        if int(row["bodies"]) != 2*n or int(row["width"]) != 20:
            raise ValueError(f"Unsupported geometry: {row['name']}")
        forms = [f"x{i}" for i in range(n)]
        bodies = [bounds_for_body(n, axis, sign, Fraction(19), forms)
                  for axis in range(n) for sign in ("+", "-")]
        write_union(HERE / row["path"], n, bodies)
    print(f"Regenerated {len(rows)} benchmark inputs")


if __name__ == "__main__":
    main()
