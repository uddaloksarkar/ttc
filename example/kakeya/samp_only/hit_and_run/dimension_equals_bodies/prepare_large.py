#!/usr/bin/env python3
"""Prepare larger dimension=body-count box unions without changing the matrix."""

import csv
from pathlib import Path

import probe

HERE = Path(__file__).resolve().parent
BENCHMARKS = HERE / "large_counts" / "benchmarks"
COUNTS = (64, 100)


def main() -> None:
    BENCHMARKS.mkdir(parents=True, exist_ok=True)
    rows = []
    for count in COUNTS:
        path = BENCHMARKS / f"RDHR_dimeqbody_m{count:04d}_d{probe.DENOMINATOR:08d}.smt2"
        formula = probe.make_formula(count)
        if path.exists() and path.read_text() != formula:
            raise ValueError(f"existing benchmark differs: {path}")
        if not path.exists():
            path.write_text(formula)
        body, union = probe.verify_written(path, count)
        rows.append(dict(benchmark=path.name, bodies=count, dimension=count,
                         thin_denominator=probe.DENOMINATOR,
                         exact_body_volume=str(body), exact_union_volume=str(union),
                         path=str(path.relative_to(HERE))))
        print(f"independently verified {count} bodies in dimension {count}", flush=True)
    manifest = HERE / "large_counts" / "geometry.csv"
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0])
        writer.writeheader()
        writer.writerows(rows)
    print(manifest, flush=True)


if __name__ == "__main__":
    main()
