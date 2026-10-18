#!/usr/bin/env python3
"""Generate logically equivalent box unions with shared bounds factored out."""

import csv
import re
from fractions import Fraction
from pathlib import Path

import probe

HERE = Path(__file__).resolve().parent / "factored"
COUNTS = (8, 64, 100)


def formula(count: int) -> str:
    lines = ["(set-logic QF_LRA)"]
    lines += [f"(declare-fun x{j} () Real)" for j in range(count)]
    lines.append("(assert (and")
    for j in range(2, count):
        radius = Fraction(1, probe.DENOMINATOR) if j == 2 else Fraction(1)
        lines.extend((f"  (<= {probe.smt_negative(radius)} x{j})",
                      f"  (<= x{j} {probe.smt_positive(radius)})"))
    lines.append("  (or")
    for i in range(count):
        x = Fraction(count - 1 + 99 * i, count - 1)
        y = Fraction(100) / x
        lines.extend(("    (and",
                      f"      (<= {probe.smt_negative(x)} x0)",
                      f"      (<= x0 {probe.smt_positive(x)})",
                      f"      (<= {probe.smt_negative(y)} x1)",
                      f"      (<= x1 {probe.smt_positive(y)})",
                      "    )"))
    return "\n".join(lines + ["  )", "))", "(check-sat)"]) + "\n"


def verify(path: Path, count: int) -> tuple[Fraction, Fraction]:
    raw = path.read_text()
    if raw != formula(count):
        raise ValueError(f"factored formula differs: {path}")
    declared = re.findall(r"^\(declare-fun x(\d+) \(\) Real\)$", raw, re.M)
    if declared != [str(j) for j in range(count)]:
        raise ValueError(f"incorrect variable declarations: {path}")
    bounds = []
    for j in range(count):
        upper = [probe.parse_fraction(token) for token in re.findall(
            rf"^\s*\(<= x{j} (.+)\)$", raw, re.M)]
        lower = [probe.parse_fraction(token) for token in re.findall(
            rf"^\s*\(<= \(- (.+)\) x{j}\)$", raw, re.M)]
        expected_rows = count if j < 2 else 1
        if len(upper) != expected_rows or upper != lower:
            raise ValueError(f"incorrect bounds for x{j}: {path}")
        bounds.append(upper)
    xs, ys = bounds[:2]
    if xs != [Fraction(count - 1 + 99 * i, count - 1)
              for i in range(count)] or ys != [Fraction(100) / x for x in xs]:
        raise ValueError(f"incorrect x-y staircase: {path}")
    for j in range(2, count):
        expected = Fraction(1, probe.DENOMINATOR) if j == 2 else Fraction(1)
        if bounds[j] != [expected]:
            raise ValueError(f"incorrect shared bounds for x{j}: {path}")
    factor = Fraction(1)
    for j in range(2, count):
        factor *= 2 * bounds[j][0]
    body = 4 * xs[0] * ys[0] * factor
    endpoints = [Fraction(0)] + sorted(xs)
    quadrant = sum((right - left) * max(y for x, y in zip(xs, ys) if x >= right)
                   for left, right in zip(endpoints, endpoints[1:]))
    union = 4 * quadrant * factor
    if (body, union) != probe.exact_volumes(count):
        raise ValueError(f"factored geometry differs from original: {path}")
    return body, union


def main() -> None:
    benchmarks = HERE / "benchmarks"
    benchmarks.mkdir(parents=True, exist_ok=True)
    rows = []
    for count in COUNTS:
        path = benchmarks / f"RDHR_factored_m{count:04d}_d{probe.DENOMINATOR:08d}.smt2"
        expected = formula(count)
        if path.exists() and path.read_text() != expected:
            raise ValueError(f"existing input differs: {path}")
        if not path.exists():
            path.write_text(expected)
        body, union = verify(path, count)
        rows.append(dict(benchmark=path.name, bodies=count, dimension=count,
                         thin_denominator=probe.DENOMINATOR,
                         exact_body_volume=str(body), exact_union_volume=str(union)))
        print(f"verified equivalent factored input with {count} bodies", flush=True)
    with (HERE / "geometry.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
