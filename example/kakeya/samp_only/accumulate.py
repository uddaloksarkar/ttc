#!/usr/bin/env python3
"""Validate compressed cluster logs, write one CSV, and regenerate figures.

Figures (in <campaign>/figures):
  walk_length_vs_relative_error_n3     observed epsilon against walk length w
  walk_length_vs_relative_error_n_eq_m for n = 3 and n = m; hue is the family,
                                       lightness the number of bodies m
  boxes3d_early_<benchmark>            3D views of the 3D benchmarks with the
  boxes3d_late_<benchmark>             most bodies: the earliest to reach
                                       epsilon <= TARGET, and one that never
                                       does; both share one stretched x2 scale

Run from anywhere: python3 example/kakeya/samp_only/accumulate.py [directory]
"""

from __future__ import annotations

import argparse
import csv
import lzma
import re
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

sys.dont_write_bytecode = True
from sunflower import PALETTE, STYLE, box_bounds, draw_boxes_3d, exact_box_volumes

HERE = Path(__file__).resolve().parent
TARGET = 0.2
INK, GRID = PALETTE["ink"], PALETTE["grid"]
# Light -> dark ramp per family; a body count m gets the same step in each.
FAMILIES = {name: PALETTE[name] for name in ("n3", "n_eq_m")}
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


def accumulate(folder: Path, sampler_name: str, logs_folder: Path | None = None) -> None:
    logs_folder = logs_folder or folder
    sampler = "ball" if sampler_name == "ball_walk" else "rdhr"
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
    cases = read(folder / "cases.csv")
    if logs_folder != folder:
        volume_by_benchmark = {case["benchmark"]: case["exact_union_volume"]
                               for case in cases}
        known_logs = {case["log_file"].replace(".log", ".out.xz")
                      for case in cases}
        for log in sorted(logs_folder.glob("*.out.xz")):
            if log.name in known_logs:
                continue
            match = re.fullmatch(r"(.+)_walk(\d+)\.out\.xz", log.name)
            benchmark = f"{match.group(1)}.smt2" if match else ""
            if benchmark not in catalog or benchmark not in volume_by_benchmark:
                raise ValueError(f"log has no catalog/case metadata: {log}")
            cases.append(dict(benchmark=benchmark, walk_length=str(int(match.group(2))),
                              exact_union_volume=volume_by_benchmark[benchmark],
                              log_file=log.name.removesuffix(".out.xz") + ".log"))
    results = []
    for case in cases:
        name = case["benchmark"]
        meta = catalog[name]
        status, estimate, error = validate(logs_folder, case, meta, volumes[name], sampler)
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
    draw(folder)
    print(f"{folder}: {sum(row['run_status']=='0' for row in results)}/{len(results)} "
          "validated cases")


def family(bodies: int, dimension: int) -> str | None:
    if dimension == 3 and bodies != 3:
        return "n3"
    return "n_eq_m" if dimension == bodies else None


def walk_to_target(curve: list[tuple[int, float]]) -> int | None:
    """Smallest tested w with epsilon <= TARGET there and at every longer w."""
    reached = None
    for walk, error in reversed(curve):
        if error > TARGET:
            break
        reached = walk
    return reached


def save(fig, output: Path) -> None:
    output.parent.mkdir(exist_ok=True)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), bbox_inches="tight", dpi=220)
    plt.close(fig)


def draw(folder: Path) -> None:
    plt.rcParams.update(STYLE)
    curves, shape = defaultdict(list), {}
    for row in read(folder / "combined_results.csv"):
        if row["run_status"] == "0" and row["relative_error"]:
            curves[row["benchmark"]].append((int(row["walk_length"]),
                                             float(row["relative_error"])))
            shape[row["benchmark"]] = (int(row["bodies"]), int(row["dimension"]))
    for curve in curves.values():
        curve.sort()
    levels = sorted({bodies for bodies, _ in shape.values()})
    walks = [walk for curve in curves.values() for walk, _ in curve]
    top = 1.05 * max(error for curve in curves.values() for _, error in curve)

    # One figure per family, with identical axes so they can be set side by side.
    for name, ramp in FAMILIES.items():
        members = sorted((b for b in curves if family(*shape[b]) == name),
                         key=lambda b: shape[b])
        if not members:
            continue
        colormap = LinearSegmentedColormap.from_list(name, ramp)
        color = {m: colormap(levels.index(m) / max(len(levels) - 1, 1)) for m in levels}
        fig, ax = plt.subplots(figsize=(5.2, 3.6))
        for benchmark in members:
            xs, ys = zip(*curves[benchmark])
            ax.plot(xs, ys, color=color[shape[benchmark][0]], linewidth=1.4, alpha=0.9)
        ax.axhspan(0, TARGET, color=ramp[2], alpha=0.07, lw=0)
        ax.axhline(TARGET, color=INK, linestyle=(0, (4, 3)), linewidth=0.9)
        ax.text(1.0, TARGET, rf"$\varepsilon = {TARGET:g}$ ", color=INK, ha="right",
                va="bottom", fontsize=8.5, transform=ax.get_yaxis_transform())
        ax.set_xscale("log")
        ax.set_xlim(min(walks) / 1.3, max(walks) * 1.3)
        ax.set_ylim(0, top)
        ax.grid(True, which="major", color=GRID, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.set_xlabel(r"walk length $w$")
        ax.set_ylabel(r"observed $\varepsilon$")
        present = sorted({shape[b][0] for b in members})
        ax.legend(handles=[Line2D([0], [0], color=color[m], linewidth=2,
                                  label=rf"$m = {m}$") for m in present],
                  loc="upper right", fontsize=8, ncol=1 if len(present) <= 5 else 2,
                  handlelength=1.4, labelspacing=0.3, borderaxespad=0.2)
        save(fig, folder / "figures" / f"walk_length_vs_relative_error_{name}")

    # 3D views: among 3D benchmarks with the most bodies where one reaches
    # the target and another never does, the earliest to reach it and the one
    # that never does (tested to the longest walk; ties: highest final error).
    # Without a never-reaching one, the latest to reach the target stands in.
    reach = {b: walk_to_target(curves[b]) for b in curves if shape[b][1] == 3}
    counts = sorted({shape[b][0] for b in reach}, reverse=True)
    counts = ([m for m in counts if any(reach[b] for b in reach if shape[b][0] == m)
               and any(not reach[b] for b in reach if shape[b][0] == m)]
              or [m for m in counts if sum(bool(reach[b]) for b in reach
                                           if shape[b][0] == m) > 1])
    if not counts:
        return
    group = [b for b in reach if shape[b][0] == counts[0]]
    met = sorted((reach[b], b) for b in group if reach[b])
    never = [b for b in group if not reach[b]]
    early = met[0][1]
    late = (max(never, key=lambda b: curves[b][-1]) if never else met[-1][1])
    for pattern in ("sunflower_early_*", "sunflower_late_*", "boxes3d_*"):
        for stale in (folder / "figures").glob(pattern):
            stale.unlink()
    z_reach = max(abs(value) for benchmark in (early, late)
                  for box in box_bounds(folder / "benchmarks" / benchmark)
                  for value in box[2])
    for label, benchmark in (("early", early), ("late", late)):
        draw_boxes_3d(folder / "benchmarks" / benchmark, folder / "figures"
                      / f"boxes3d_{label}_{Path(benchmark).stem}", z_reach)
        walk = reach[benchmark]
        print(f"{folder}: boxes3d_{label}: {benchmark} "
              + (f"reaches epsilon <= {TARGET:g} at w = {walk}" if walk else
                 f"ends at epsilon = {curves[benchmark][-1][1]:.3f}"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path,
                        help="campaign or out_<sampler> directory (relative to the current working directory)")
    parser.add_argument("--sampler", choices=("ball_walk", "hit_and_run", "both"),
                        default="both")
    arguments = parser.parse_args()
    if arguments.directory is not None:
        input_folder = arguments.directory.resolve()
        if not input_folder.is_dir():
            parser.error(f"not a directory: {arguments.directory}")
        is_output_folder = input_folder.name.startswith("out_")
        inferred_name = (input_folder.name[4:] if is_output_folder
                         else input_folder.name)
        sampler_name = inferred_name if arguments.sampler == "both" else arguments.sampler
        if sampler_name not in ("ball_walk", "hit_and_run"):
            parser.error("cannot infer sampler from directory name; pass --sampler")
        folder = input_folder.parent / sampler_name if is_output_folder else input_folder
        if not (folder / "catalog.csv").is_file():
            parser.error(f"campaign metadata not found in: {folder}")
        accumulate(folder, sampler_name, input_folder)
    else:
        for name in ("ball_walk", "hit_and_run"):
            if arguments.sampler in ("both", name):
                accumulate(HERE / name, name)
