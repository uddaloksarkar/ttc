#!/usr/bin/env python3
"""Rebuild measurements for boxes run with TTC's estimated volume oracle."""

import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def parse(case: dict, log: Path) -> dict:
    raw = log.read_text()
    if "run_status=0" not in raw:
        raise ValueError(f"Incomplete run: {log}")
    def one(pattern: str) -> str:
        match = re.search(pattern, raw, re.M)
        if not match:
            raise ValueError(f"Missing {pattern!r}: {log}")
        return match.group(1)
    m = int(case["bodies"])
    n = int(case["dimension"])
    body = float(case["exact_body_volume"])
    union = float(case["exact_union_volume"])
    walk = int(one(r"^c sample walk length: (\d+)$"))
    if int(one(r"^c polytopes: (\d+)$")) != m:
        raise ValueError(f"Body count mismatch: {log}")
    if int(one(r"^c volume walk length: (\d+)$")) != int(case["volume_walk_length"]):
        raise ValueError(f"Volume walk mismatch: {log}")
    vols = [float(x) for x in re.findall(
        r"^c\s+\d+\s+([\d.]+)\s+\d+\s+\d+$", raw, re.M)]
    if len(vols) != m:
        raise ValueError(f"Expected {m} body estimates: {log}")
    returned = float(one(r"^s volume ([\d.eE+-]+)$"))
    err_body = max(abs(v/body-1) for v in vols)
    err_union = abs(returned/union-1)
    return dict(name=case["name"], bodies=m, dimension=n, width=20,
                volume_walk_length=int(case["volume_walk_length"]),
                sample_walk_length=walk, seed=42,
                exact_body_volume=body, exact_union_volume=union,
                estimated_body_sum=sum(vols), estimated_union_volume=returned,
                max_body_relative_error=err_body, union_relative_error=err_union,
                apparent_overlap=1-returned/sum(vols),
                true_overlap=1-union/(m*body),
                body_oracle_fails=err_body > .2/12,
                ttc_fails=err_union > .2,
                log_path=str(log.relative_to(HERE)))


def main() -> None:
    cases = list(csv.DictReader((HERE / "manifest.csv").open(newline="")))
    rows = []
    for case in cases:
        for log in (HERE / "benchmarks").glob(f"{case['name']}_samp*.log"):
            rows.append(parse(case, log))
    rows.sort(key=lambda r: (r["volume_walk_length"], r["bodies"], r["sample_walk_length"]))
    if not rows:
        raise SystemExit("No completed runs")
    with (HERE / "results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['name']} walk={row['sample_walk_length']}: "
              f"union error={row['union_relative_error']:.2%}, "
              f"body error={row['max_body_relative_error']:.2%}")


if __name__ == "__main__":
    main()
