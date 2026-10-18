#!/usr/bin/env python3
"""Generate analytic-volume Kakeya boxes for exact-oracle sampler tests."""

import csv
from decimal import Decimal
from pathlib import Path

HERE = Path(__file__).resolve().parent
CASES = [
    (2, 32, 2, "0"),
    (3, 20, 2, "0"),
    (3, 30, 2, "0"),
    (3, 30, 3, "0"),
    (3, 34, 2, "0"),
    (4, 10, 3, "0"),
    (4, 10, 3, "0.01"),
    (4, 10, 4, "0"),
    (4, 10, 4, "0.01"),
    (4, 16, 3, "0"),
    (12, 1000, 12, "0", "ball"),
    (20, 1000, 20, "0", "ball"),
    (3, 1000, 2, "0", "billiard", 1000),
]


def number(value: Decimal) -> str:
    return f"(- {abs(value)})" if value < 0 else str(value)


def formula(n: int, width: int, m: int, shift: Decimal,
            thin_denominator: int = 1) -> str:
    if not (n >= m >= 2 and width >= 4 and width % 2 == 0):
        raise ValueError("Expected 2 <= bodies <= dimension and positive even width")
    if thin_denominator < 1 or (thin_denominator > 1 and m >= n):
        raise ValueError("A thin shared coordinate requires bodies < dimension")
    lines = ["(set-logic QF_LRA)"]
    lines.extend(f"(declare-fun x{i} () Real)" for i in range(n))
    lines.append("(assert (or")
    for axis in range(m):
        lines.append("  (and")
        for i in range(n):
            if i == n-1 and thin_denominator > 1:
                lines.extend((f"    (<= (- (/ 1 {thin_denominator})) x{i})",
                              f"    (<= x{i} (/ 1 {thin_denominator}))"))
                continue
            radius = Decimal(width // 2 if i == axis else 1)
            offset = shift if i == axis == m - 1 else Decimal(0)
            lower = (str(-radius+offset) if thin_denominator > 1
                     else number(-radius+offset))
            lines.extend((f"    (<= {lower} x{i})",
                          f"    (<= x{i} {number(radius+offset)})"))
        lines.append("  )")
    return "\n".join(lines + ["))", "(check-sat)"]) + "\n"


def main() -> None:
    out = HERE / "benchmarks"
    out.mkdir(exist_ok=True)
    rows = []
    for item in CASES:
        n, width, m, shift_text = item[:4]
        sampler = item[4] if len(item) > 4 else "billiard"
        thin_denominator = item[5] if len(item) > 5 else 1
        shift = Decimal(shift_text)
        name = f"SAMPONLY_m{m}_n{n}_w{width}"
        if shift:
            name += "_shiftp001"
        if sampler != "billiard":
            name += f"_{sampler}"
        if thin_denominator > 1:
            name += f"_thin{thin_denominator}"
        path = out / f"{name}.smt2"
        path.write_text(formula(n, width, m, shift, thin_denominator))
        core = Decimal(2**n) / Decimal(thin_denominator)
        rows.append(dict(name=name, dimension=n, bodies=m, width=width,
                         last_spoke_shift=shift_text,
                         thin_denominator=thin_denominator,
                         exact_body_volume=core*width/2,
                         exact_union_volume=core*(1+m*(width//2-1)),
                         seed=42, sampler=sampler,
                         path=f"benchmarks/{path.name}"))
    with (HERE / "manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Generated {len(rows)} exact-volume box unions")


if __name__ == "__main__":
    main()
