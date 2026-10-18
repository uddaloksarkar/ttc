#!/usr/bin/env python3
"""Validate compressed cluster logs, write one CSV, and regenerate figures.

Run from anywhere: python3 example/kakeya/samp_only/accumulate.py
"""

from __future__ import annotations

import argparse
import csv
import lzma
import re
import sys
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True
from plot_all import draw
from sunflower import exact_box_volumes

HERE = Path(__file__).resolve().parent
BODY_ROW = re.compile(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$", re.M)
ESTIMATE = re.compile(r"^s volume ([0-9.eE+-]+)\s*$", re.M)


def read(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(newline="")))


def validate(folder: Path, case: dict[str, str], meta: dict[str, str],
             body_volumes: list[Fraction], sampler: str) -> tuple[str, str, str]:
    log = folder / case["log_file"].replace(".log", ".out.xz")
    if not log.exists():
        return "pending", "", ""
    try:
        with lzma.open(log, "rt", errors="replace") as stream:
            raw = stream.read()
    except (OSError, lzma.LZMAError):
        return "invalid_xz", "", ""
    expected = (f"c sample walk length: {case['walk_length']}",
                f"c sampler walk: {sampler}",
                "c body volume oracle: exact axis-aligned box",
                f"c polytopes: {meta['bodies']}",
                f"c real variables: {meta['dimension']}",
                "c exit 0")
    rows = BODY_ROW.findall(raw)
    estimates = ESTIMATE.findall(raw)
    marker = re.findall(r"^run_status=(\S+)\s*$", raw, re.M)
    if (not marker or marker[-1] != "0" or any(item not in raw for item in expected)
            or len(rows) != len(body_volumes) or len(estimates) != 1
            or [int(row[0]) for row in rows] != list(range(1, len(body_volumes) + 1))):
        return "invalid_output", "", ""
    for row, exact in zip(rows, body_volumes):
        measured = float(row[1])
        if abs(measured - float(exact)) > max(5.1e-7, 1e-6 * float(exact)):
            return "wrong_body_volume", "", ""
    estimate = float(estimates[0])
    if estimate <= 0:
        return "invalid_estimate", "", ""
    error = abs(estimate / float(Fraction(case["exact_union_volume"])) - 1)
    return "0", f"{estimate:.17g}", f"{error:.17g}"


def accumulate(sampler_folder: str) -> None:
    folder = HERE / sampler_folder
    sampler = "ball" if sampler_folder == "ball_walk" else "rdhr"
    catalog = {row["benchmark"]: row for row in read(folder / "catalog.csv")}
    prior = {(row["benchmark"], row["walk_length"]): row
             for row in read(folder / "measurements.csv")}
    volumes = {}
    for name, meta in catalog.items():
        path = folder / "benchmarks" / name
        if not path.is_file():
            raise FileNotFoundError(path)
        if len(exact_box_volumes(path)) != int(meta["bodies"]):
            raise ValueError(f"body count mismatch: {path}")
        volumes[name] = exact_box_volumes(path)
    results = []
    for case in read(folder / "cases.csv"):
        name = case["benchmark"]
        meta = catalog[name]
        status, estimate, error = validate(folder, case, meta, volumes[name], sampler)
        baseline = prior.get((name, case["walk_length"]))
        if status != "0" and baseline is not None:
            # A prior validated completion remains usable if a later cluster
            # rerun is absent, interrupted, or malformed.
            error = baseline["relative_error"]
            status_out = "0"
            evidence = "prior_validated_csv"
        else:
            status_out = status
            evidence = "compressed_ttc_log" if status == "0" else ""
        results.append(dict(benchmark=name, bodies=meta["bodies"],
                            dimension=meta["dimension"],
                            walk_length=case["walk_length"],
                            exact_union_volume=case["exact_union_volume"],
                            estimated_union_volume=estimate,
                            relative_error=error,
                            fails_20pct=(int(float(error) > 0.2) if error else ""),
                            run_status=status_out, evidence=evidence,
                            cluster_log_status=status,
                            compressed_log=(case["log_file"].replace(".log", ".out.xz")
                                            if status != "pending" else "")))
    target = folder / "combined_results.csv"
    temporary = folder / "combined_results.csv.tmp"
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)
    temporary.replace(target)
    draw(sampler_folder)
    print(f"{sampler_folder}: {sum(row['run_status']=='0' for row in results)}/{len(results)} "
          "validated cases")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler", choices=("ball_walk", "hit_and_run", "both"),
                        default="both")
    arguments = parser.parse_args()
    for name in ("ball_walk", "hit_and_run"):
        if arguments.sampler in ("both", name):
            accumulate(name)
