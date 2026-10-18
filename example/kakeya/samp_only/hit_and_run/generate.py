#!/usr/bin/env python3
"""Generate fixed-dimension, exact-volume box unions with varying body counts.

For body i, x has radius L_i=1+99*i/(m-1), y has radius 100/L_i,
and z has radius 1/(20*walk). Thus every body has volume 800/(20*walk).
The boxes all contain the common core [-1,1]^2 in the xy projection.
"""

import csv
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNTS = (2, 4, 8, 16, 32, 64, 128, 256, 512, 1000)
WALKS = (1, 100, 10000, 1000000)


def smt_positive(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"(/ {value.numerator} {value.denominator})"


def smt_negative(value: Fraction) -> str:
    return f"(- {smt_positive(value)})"


def make_formula(count: int, walk: int, denominator: int | None = None,
                 spread: int = 100, short_radius: int = 1) -> str:
    denominator = denominator if denominator is not None else 20 * walk
    lines = ["(set-logic QF_LRA)",
             "(declare-fun x0 () Real)",
             "(declare-fun x1 () Real)",
             "(declare-fun x2 () Real)",
             "(assert (or"]
    for i in range(count):
        x = Fraction(short_radius * (count - 1) + (spread - short_radius) * i,
                     count - 1)
        y = Fraction(spread * short_radius, 1) / x
        z = Fraction(1, denominator)
        lines.extend(["  (and",
                      f"    (<= {smt_negative(x)} x0)",
                      f"    (<= x0 {smt_positive(x)})",
                      f"    (<= {smt_negative(y)} x1)",
                      f"    (<= x1 {smt_positive(y)})",
                      f"    (<= {smt_negative(z)} x2)",
                      f"    (<= x2 {smt_positive(z)})",
                      "  )"])
    return "\n".join(lines + ["))", "(check-sat)"]) + "\n"


def exact_union_quadrant(count: int, spread: int = 100,
                         short_radius: int = 1) -> Fraction:
    # Sorted x radii increase, while y radii decrease. The union is a
    # staircase: sum (x_i-x_(i-1))*y_i, with x_(-1)=0.
    area = Fraction(0)
    previous_x = Fraction(0)
    for i in range(count):
        x = Fraction(short_radius * (count - 1) + (spread - short_radius) * i,
                     count - 1)
        y = Fraction(spread * short_radius, 1) / x
        area += (x - previous_x) * y
        previous_x = x
    return area


def main() -> None:
    benchmarks = HERE / "benchmarks"
    benchmarks.mkdir(parents=True, exist_ok=True)
    rows = []
    for walk in WALKS:
        denominator = 20 * walk
        for count in COUNTS:
            name = f"SAMPONLY_fixed3_m{count:04d}_target{walk:07d}"
            filename = f"{name}.smt2"
            (benchmarks / filename).write_text(make_formula(count, walk))
            body_volume = Fraction(800, denominator)
            union_volume = Fraction(8, denominator) * exact_union_quadrant(count)
            rows.append({
                "name": name, "bodies": count, "dimension": 3,
                "target_walk_length": walk, "thin_denominator": denominator,
                "exact_body_volume": f"{float(body_volume):.17g}",
                "exact_union_volume": f"{float(union_volume):.17g}",
                "union_to_body_ratio": f"{float(union_volume / body_volume):.17g}",
                "predicted_trap_error": f"{float(1 - body_volume / union_volume):.17g}",
                "max_walk_displacement": f"{4 / (20 * 3**0.5):.17g}",
                "sampler": "ball", "seed": 42, "epsilon": 0.2, "delta": 0.2,
                "volume_oracle": "exact-box-volumes",
                "path": f"benchmarks/{filename}",
            })
    with (HERE / "manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Generated {len(rows)} cases across {len(COUNTS)} body counts and {len(WALKS)} walk lengths")


if __name__ == "__main__":
    main()
