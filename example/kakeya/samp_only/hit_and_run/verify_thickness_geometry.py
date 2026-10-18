#!/usr/bin/env python3
"""Independently verify every written fixed-3D thickness SMT benchmark."""

import csv
import re
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAME = re.compile(r"SAMPONLY_thickness_m(\d+)_walk(\d+)_d(\d+)"
                  r"(?:_r(\d+))?(?:_s(\d+))?\.smt2$")
SIDE_NAME = re.compile(r"SAMPONLY_side_m(\d+)_walk(\d+)_d(\d+)_r(\d+)_"
                       r"s(\d+)over(\d+)\.smt2$")


def radius(token):
    token = token.strip()
    fraction = re.fullmatch(r"\(/\s+(\d+)\s+(\d+)\)", token)
    if fraction:
        return Fraction(int(fraction.group(1)), int(fraction.group(2)))
    return Fraction(int(token))


def bounds(raw, coordinate):
    upper = re.findall(rf"^\s*\(<= x{coordinate} (.+)\)\s*$", raw, re.M)
    lower = re.findall(rf"^\s*\(<= \(- (.+)\) x{coordinate}\)\s*$", raw, re.M)
    upper, lower = [radius(x) for x in upper], [radius(x) for x in lower]
    if lower != upper:
        raise ValueError(f"asymmetric x{coordinate} bounds")
    return upper


def verify(path):
    match = NAME.fullmatch(path.name)
    side_match = SIDE_NAME.fullmatch(path.name) if not match else None
    if not match and not side_match:
        raise ValueError(f"unexpected benchmark name: {path.name}")
    if match:
        count, walk, denominator, spread, short = match.groups()
        count, walk, denominator = map(int, (count, walk, denominator))
        spread, short = int(spread or 100), Fraction(int(short or 1))
    else:
        count, walk, denominator, spread, numerator, side_denominator = map(
            int, side_match.groups())
        spread, short = int(spread), Fraction(numerator, side_denominator)
    raw = path.read_text()
    xs, ys, zs = (bounds(raw, coordinate) for coordinate in range(3))
    if not (len(xs) == len(ys) == len(zs) == count):
        raise ValueError(f"wrong body count: {path.name}")
    expected_x = [(short * (count - 1) + (spread - short) * i) / (count - 1)
                  for i in range(count)]
    if xs != expected_x:
        raise ValueError(f"wrong x staircase: {path.name}")
    if ys != [spread * short / x for x in xs]:
        raise ValueError(f"wrong y staircase: {path.name}")
    if zs != [Fraction(1, denominator)] * count:
        raise ValueError(f"wrong z thickness: {path.name}")
    volumes = [8 * x * y * z for x, y, z in zip(xs, ys, zs)]
    expected_body = Fraction(8 * spread, denominator) * short
    if volumes != [expected_body] * count:
        raise ValueError(f"wrong body volume: {path.name}")
    # Integrate the parsed boxes over x slabs, independent of the generator's
    # exact_union_quadrant implementation.
    endpoints = [Fraction(0)] + sorted(set(xs))
    quadrant = sum((right - left) * max(y for x, y in zip(xs, ys) if x >= right)
                   for left, right in zip(endpoints, endpoints[1:]))
    union = 8 * zs[0] * quadrant
    return {"name": path.stem, "bodies": count, "walk_length": walk,
            "denominator": denominator, "spread": spread,
            "short_radius": str(short),
            "exact_body_volume": str(expected_body),
            "exact_union_volume": str(union), "verified": 1}


def main():
    paths = sorted((HERE / "benchmarks").glob("SAMPONLY_thickness_*.smt2"))
    rows = [verify(path) for path in paths]
    if not rows:
        raise RuntimeError("no thickness SMT files found")
    output = HERE / "thickness_geometry_verification.csv"
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Verified {len(rows)} written fixed-3D thickness SMT inputs -> {output}")


if __name__ == "__main__":
    main()
