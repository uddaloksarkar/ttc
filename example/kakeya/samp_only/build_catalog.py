#!/usr/bin/env python3
"""Collect validated sampler-only failures into flat sampler directories.

The old experiment trees are input only.  This script never deletes or edits
them, which matters while TTC processes still write to their original logs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import re
import shutil
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
SOURCES = {
    "ball_walk": (
        "results.csv", "target_results.csv",
        "body_count_walk_extension/results.csv",
        "body_scaling/results.csv",
        "body_scaling/dimension_scaling/results.csv",
        "body_scaling/thickness_results.csv",
        "body_scaling/same_geometry_results.csv",
        "body_scaling/side_refinement.csv",
    ),
    "hit_and_run": (
        "hit_and_run/results.csv",
        "hit_and_run/dimension_results.csv",
        "hit_and_run/dimension_pilot_results.csv",
    ),
}


def read_csv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(newline=""))) if path.exists() else []


def input_file(row: dict[str, str], table: Path, root: Path) -> Path | None:
    name = row.get("benchmark") or row.get("path")
    if not name:
        name = f"benchmarks/{row.get('name', '')}.smt2"
    candidates = [root / "hit_and_run" / "benchmarks" / Path(name).name,
                  root / "hit_and_run" / "unselected_benchmarks" / Path(name).name,
                  table.parent / name, root / name]
    for candidate in candidates:
        if candidate.is_file() and candidate.suffix == ".smt2":
            return candidate
    matches = list(root.rglob(Path(name).name))
    return matches[0] if len(matches) == 1 else None


def measurement(row: dict[str, str], sampler: str) -> tuple[int, int, int, float, bool] | None:
    if sampler == "ball_walk" and row.get("sampler", "ball") not in ("ball", ""):
        return None
    if row.get("run_status") not in (None, "", "0"):
        return None
    if row.get("ttc_exit_code") not in (None, "", "0"):
        return None
    if row.get("exact_body_oracle") not in (None, "", "1", "True", "true"):
        return None
    error = (row.get("relative_error") or row.get("union_relative_error"))
    walk = (row.get("walk_length") or row.get("sample_walk_length")
            or row.get("target_walk_length"))
    if not error or not walk:
        return None
    try:
        bodies = int(row["bodies"])
        dimension = int(row.get("dimension") or (3 if sampler == "hit_and_run" else 0))
        value = float(error)
        return bodies, dimension, int(walk), value, value > 0.2
    except (ValueError, KeyError):
        return None


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def build(sampler: str, root: Path) -> None:
    destination = HERE / sampler
    inputs = destination / "benchmarks"
    inputs.mkdir(parents=True, exist_ok=True)
    (destination / "figures").mkdir(exist_ok=True)
    observed: dict[Path, list[dict[str, object]]] = defaultdict(list)
    planned: dict[Path, dict[int, str]] = defaultdict(dict)
    tables = list(SOURCES[sampler])
    for relative in dict.fromkeys(tables):
        table = root / relative
        for row in read_csv(table):
            item = measurement(row, sampler)
            path = input_file(row, table, root)
            if path is not None:
                walk_text = (row.get("walk_length") or row.get("sample_walk_length")
                             or row.get("target_walk_length"))
                if walk_text and walk_text.isdigit():
                    planned[path.resolve()][int(walk_text)] = row.get("exact_union_volume", "")
            if item is None or path is None:
                continue
            bodies, dimension, walk, error, failed = item
            if dimension == 0:
                dimension = len(re.findall(r"^\(declare-fun\s+\S+\s+\(\)\s+Real\)",
                                           path.read_text(), re.M))
                if dimension == 0:
                    raise ValueError(f"cannot determine dimension: {path}")
            observed[path.resolve()].append(dict(bodies=bodies, dimension=dimension,
                                                walk_length=walk, relative_error=error,
                                                fails_20pct=int(failed),
                                                source_table=relative))
    names: dict[str, str] = {}
    catalog = []
    results = []
    cases_to_run = []
    for source in sorted(observed):
        cases = observed[source]
        if not any(case["fails_20pct"] for case in cases):
            continue
        digest = hashlib.sha256(source.read_bytes()).hexdigest()[:12]
        name = f"{source.stem}_{digest}.smt2"
        target = inputs / name
        if sampler == "ball_walk":
            if target.exists() and target.read_bytes() != source.read_bytes():
                raise ValueError(f"collision: {target}")
            shutil.copy2(source, target)
        names[str(source)] = name
        catalog.append(dict(benchmark=name, bodies=cases[0]["bodies"],
                            dimension=cases[0]["dimension"],
                            sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                            original_path=os.path.relpath(source, root.resolve())))
        for walk, exact_union in sorted(planned[source].items()):
            cases_to_run.append(dict(benchmark=name, walk_length=walk,
                                     exact_union_volume=exact_union,
                                     log_file=f"{Path(name).stem}_walk{walk:07d}.log"))
        unique = {}
        for case in cases:
            key = case["walk_length"]
            if key in unique and abs(unique[key]["relative_error"] - case["relative_error"]) > 1e-10:
                raise ValueError(f"conflicting measured walk: {source}, {key}")
            unique[key] = case
        for case in sorted(unique.values(), key=lambda item: item["walk_length"]):
            results.append(dict(benchmark=name, **case))
    write_csv(destination / "catalog.csv", catalog,
              ["benchmark", "bodies", "dimension", "sha256", "original_path"])
    write_csv(destination / "measurements.csv", results,
              ["benchmark", "bodies", "dimension", "walk_length", "relative_error",
               "fails_20pct", "source_table"])
    write_csv(destination / "cases.csv", cases_to_run,
              ["benchmark", "walk_length", "exact_union_volume", "log_file"])
    print(f"{sampler}: {len(catalog)} failing geometries; {len(results)} validated measurements; "
          f"{len(cases_to_run)} TTC cases")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path,
                        help="optional override for both legacy source trees")
    args = parser.parse_args()
    for sampler in SOURCES:
        if args.source_root:
            root = args.source_root.resolve()
        elif sampler == "ball_walk":
            legacy = HERE.parent.parent / "kakeya_previous" / "samp_only_legacy"
            root = legacy if legacy.exists() else HERE
        else:
            root = HERE
        build(sampler, root)


if __name__ == "__main__":
    main()
