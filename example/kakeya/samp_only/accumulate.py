#!/usr/bin/env python3
"""Validate compressed cluster logs and plot walk length versus relative error."""

from __future__ import annotations

import argparse
import csv
import hashlib
import lzma
import re
import sys
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, LogNorm, Normalize
from matplotlib.lines import Line2D

sys.dont_write_bytecode = True
from sunflower import PALETTE, STYLE
from hit_and_run.utility import exact_box_volumes

HERE = Path(__file__).resolve().parent
TARGET = 0.2
EXCLUDE_BOX_DENOMINATOR_FROM = 4_000_000
INK = PALETTE["ink"]
FAMILIES = {"n3": ("3D boxes", PALETTE["n3"][2]),
            "n_eq_m": ("n=m boxes", PALETTE["n_eq_m"][2]),
            "radial": ("radial urchins", PALETTE["radial"][2])}
BODY_ROW = re.compile(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*$", re.M)
ESTIMATE = re.compile(r"^s volume ([0-9.eE+-]+)\s*$", re.M)


def read(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(newline="")))


def output_suffix(folder: Path, logs_folder: Path, sampler_name: str) -> str:
    """Name artifacts after the cluster output directory's job ID."""
    if logs_folder == folder:
        return ""
    prefix = f"out_{sampler_name}_"
    if not logs_folder.name.startswith(prefix):
        return ""
    suffix = logs_folder.name[len(prefix):]
    if not re.fullmatch(r"[A-Za-z0-9_-]+", suffix):
        raise ValueError(f"invalid output job ID: {logs_folder.name}")
    return f"_{suffix}"


def benchmark_path(folder: Path, name: str, meta: dict[str, str]) -> Path:
    staged = folder / "benchmarks" / name
    return staged if staged.is_file() else HERE / meta["original_path"]


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


def validate_radial(folder: Path, meta: dict[str, str], walk: int) -> tuple[str, str, str]:
    log = folder / f"{Path(meta['benchmark']).stem}_walk{walk:07d}.out.xz"
    if not log.exists():
        return "pending", "", ""
    try:
        with lzma.open(log, "rt", errors="replace") as stream:
            raw = stream.read()
    except (OSError, lzma.LZMAError):
        return "invalid_xz", "", ""
    marker = re.findall(r"^run_status=(\S+)\s*$", raw, re.M)
    if not marker or marker[-1] != "0":
        return marker[-1] if marker else "invalid_output", "", ""
    n = int(meta["bodies"])
    expected = (f"c sample walk length: {walk}", "c sampler walk: rdhr",
                "c body volume oracle: exact parallelepiped",
                f"c polytopes: {n}", f"c real variables: {n}", "c exit 0")
    rows, estimates = BODY_ROW.findall(raw), ESTIMATE.findall(raw)
    if (any(item not in raw for item in expected) or len(rows) != n
            or [int(row[0]) for row in rows] != list(range(1, n + 1))
            or len(estimates) != 1):
        return "invalid_output", "", ""
    exact_body = float(Fraction(meta["exact_body_volume"]))
    if any(abs(float(row[1]) - exact_body) > max(1e-6, abs(exact_body) * 1e-6)
           for row in rows):
        return "wrong_body_volume", "", ""
    estimate = float(estimates[0])
    if estimate <= 0:
        return "invalid_estimate", "", ""
    error = abs(estimate / float(Fraction(meta["exact_union_volume"])) - 1)
    return "0", f"{estimate:.17g}", f"{error:.17g}"


def accumulate(folder: Path, sampler_name: str, logs_folder: Path | None = None) -> None:
    logs_folder = logs_folder or folder
    sampler = "ball" if sampler_name == "ball_walk" else "rdhr"
    catalog = {row["benchmark"]: row for row in read(folder / "catalog.csv")}
    prior = {(row["benchmark"], row["walk_length"]): row
             for row in read(folder / "measurements.csv")}
    volumes = {}
    for name, meta in catalog.items():
        path = benchmark_path(folder, name, meta)
        if not path.is_file():
            raise FileNotFoundError(path)
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta["sha256"]:
            raise ValueError(f"benchmark content differs from catalog: {path}")
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
            if sampler_name == "hit_and_run" and log.name.startswith("RDHR_radial_urchin_"):
                # The radial family has its own exact-parallelepiped validator.
                continue
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
    if sampler_name == "hit_and_run" and logs_folder != folder:
        radial_geometry = folder / "radial_geometry.csv"
        radial = read(radial_geometry)
        walks = sorted(int(row["target_walk_length"]) for row in radial)
        for meta in radial:
            benchmark = meta["benchmark"]
            if not (folder / "benchmarks" / benchmark).is_file():
                raise FileNotFoundError(folder / "benchmarks" / benchmark)
            for walk in walks:
                status, estimate, error = validate_radial(logs_folder, meta, walk)
                compressed = f"{Path(benchmark).stem}_walk{walk:07d}.out.xz"
                results.append(dict(benchmark=benchmark, bodies=meta["bodies"],
                                    dimension=meta["dimension"], walk_length=walk,
                                    exact_union_volume=meta["exact_union_volume"],
                                    estimated_union_volume=estimate,
                                    relative_error=error,
                                    fails_20pct=(int(float(error) > TARGET) if error else ""),
                                    run_status=status,
                                    evidence="compressed_ttc_log" if status == "0" else "",
                                    cluster_log_status=status,
                                    compressed_log=compressed if status != "pending" else ""))
    suffix = output_suffix(folder, logs_folder, sampler_name)
    target = folder / f"combined_results{suffix}.csv"
    temporary = folder / f"combined_results{suffix}.csv.tmp"
    with temporary.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=results[0])
        writer.writeheader()
        writer.writerows(results)
    temporary.replace(target)
    draw(folder, target, suffix, current_output_only=logs_folder != folder)
    print(f"{target}: {sum(row['run_status']=='0' for row in results)}/{len(results)} "
          "validated cases")


def family(benchmark: str, bodies: int, dimension: int) -> str | None:
    if benchmark.startswith("RDHR_radial_urchin_"):
        return "radial"
    if dimension == 3 and bodies != 3:
        return "n3"
    return "n_eq_m" if dimension == bodies else None


def save(fig, output: Path) -> None:
    output.parent.mkdir(exist_ok=True)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), bbox_inches="tight", dpi=220)
    plt.close(fig)


def draw_hit_and_run(folder: Path, curves: dict[str, list[tuple[int, float]]],
                     shape: dict[str, tuple[int, int]], suffix: str) -> None:
    plt.rcParams.update({"font.size": 24, "axes.labelsize": 24})
    radial = {name: points for name, points in curves.items()
              if name.startswith("RDHR_radial_urchin_")}
    boxes = {}
    widths = {}
    for name, points in curves.items():
        if not name.startswith(("RDHR_box_", "RDHR_dimeqbody_", "RDHR_factored_")):
            continue
        if shape[name][1] < 3:
            continue  # Two-dimensional boxes have no thin coordinate.
        match = re.search(r"_d(\d+)_", name)
        if match is None:
            raise ValueError(f"box denominator missing from {name}")
        denominator = int(match.group(1))
        if denominator >= EXCLUDE_BOX_DENOMINATOR_FROM:
            continue
        boxes[name] = points
        widths[name] = 2 / denominator
    plotted = list(radial.values()) + list(boxes.values())
    if not plotted:
        return
    walks = [walk for points in plotted for walk, _ in points]
    xlimits = (min(walks), max(walks) * 1.1)
    ymax = max(0.86, 1.05 * max(error for points in plotted for _, error in points))
    output = folder / f"figures{suffix}"

    def setup_axis(ax) -> None:
        ax.axhspan(0, TARGET, facecolor="#edf7f2", zorder=0)
        ax.axhline(TARGET, color="#303030", linewidth=1.8,
                   linestyle=(0, (4, 3)), zorder=2)
        ax.text(0.98, 0.95, rf"$\varepsilon={TARGET:g}$",
                ha="right", va="top", color="#303030",
                fontsize=24, fontweight="bold",
                transform=ax.transAxes)
        ax.set_xscale("log")
        ax.set_xlim(*xlimits)
        ax.set_ylim(0, ymax)
        ax.set_xlabel("walk length")
        ax.set_ylabel("relative error")
        ax.tick_params(labelsize=24)
        ax.grid(False)

    if radial:
        counts = {name: shape[name][0] for name in radial}
        ramp = plt.get_cmap("YlOrBr")
        cmap = LinearSegmentedColormap.from_list(
            "radial_bodies", [ramp(0.35), ramp(0.94)])
        norm = Normalize(vmin=min(counts.values()), vmax=max(counts.values()))
        fig, ax = plt.subplots(figsize=(7.2, 5.1), layout="constrained")
        for name in sorted(radial, key=lambda item: (counts[item], item)):
            xs, ys = zip(*radial[name])
            ax.plot(xs, ys, color=cmap(norm(counts[name])),
                    linewidth=1.4, alpha=0.9, solid_capstyle="round", zorder=3)
        setup_axis(ax)
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                                fraction=0.046, pad=0.035)
        colorbar.set_label("bodies = dimensions", fontsize=24)
        colorbar.ax.tick_params(labelsize=24)
        colorbar.set_ticks([n for n in sorted(set(counts.values())) if n in
                            (4, 8, 16, 24, 32, 40, 48, 56, 64, 72)])
        save(fig, output / f"radial_walk_vs_relative_error{suffix}")

    if boxes:
        ramp = plt.get_cmap("Blues")
        cmap = LinearSegmentedColormap.from_list(
            "box_width", [ramp(0.96), ramp(0.31)])
        norm = LogNorm(vmin=min(widths.values()), vmax=max(widths.values()))
        fig, ax = plt.subplots(figsize=(7.2, 5.1), layout="constrained")
        for name in sorted(boxes, key=lambda item: (widths[item], item)):
            xs, ys = zip(*boxes[name])
            ax.plot(xs, ys, color=cmap(norm(widths[name])),
                    linewidth=1.4, alpha=0.9, solid_capstyle="round", zorder=3)
        setup_axis(ax)
        colorbar = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), ax=ax,
                                fraction=0.046, pad=0.035)
        colorbar.set_label("box width", fontsize=24)
        colorbar.ax.tick_params(labelsize=24)
        ticks = [tick for tick in (0.002, 0.02, 0.2, 2)
                 if norm.vmin <= tick <= norm.vmax]
        colorbar.set_ticks(ticks)
        colorbar.set_ticklabels([f"{tick:g}" for tick in ticks])
        save(fig, output / f"boxes_walk_vs_relative_error{suffix}")


def draw(folder: Path, results_path: Path, suffix: str = "",
         current_output_only: bool = False) -> None:
    plt.rcParams.update(STYLE)
    curves, shape = defaultdict(list), {}
    for row in read(results_path):
        if (row["run_status"] == "0" and row["relative_error"]
                and (not current_output_only or row["evidence"] == "compressed_ttc_log")):
            curves[row["benchmark"]].append((int(row["walk_length"]),
                                             float(row["relative_error"])))
            shape[row["benchmark"]] = (int(row["bodies"]), int(row["dimension"]))
    for curve in curves.values():
        curve.sort()
    if not curves:
        return
    if folder.name == "hit_and_run":
        draw_hit_and_run(folder, curves, shape, suffix)
        return
    walks = [walk for curve in curves.values() for walk, _ in curve]
    top = 1.05 * max(error for curve in curves.values() for _, error in curve)
    fig, ax = plt.subplots(figsize=(7.0, 4.2))
    present = set()
    for benchmark in sorted(curves):
        group = family(benchmark, *shape[benchmark])
        if group not in FAMILIES:
            continue
        present.add(group)
        xs, ys = zip(*curves[benchmark])
        ax.plot(xs, ys, color=FAMILIES[group][1], linewidth=1.05,
                alpha=0.55 if group != "radial" else 0.65)
    ax.axhline(TARGET, color=INK, linestyle=(0, (4, 3)), linewidth=0.9)
    ax.text(1.0, TARGET, rf"$\varepsilon = {TARGET:g}$ ", color=INK,
            ha="right", va="bottom", fontsize=8.5,
            transform=ax.get_yaxis_transform())
    ax.set_xscale("log")
    ax.set_xlim(min(walks) / 1.3, max(walks) * 1.3)
    ax.set_ylim(0, top)
    ax.grid(False)
    ax.set_xlabel(r"walk length $w$")
    ax.set_ylabel("relative error")
    ax.legend(handles=[Line2D([0], [0], color=FAMILIES[group][1], linewidth=2,
                              label=FAMILIES[group][0]) for group in FAMILIES if group in present],
              loc="upper center", bbox_to_anchor=(0.5, 1.11), ncol=len(present), fontsize=8)
    save(fig, folder / f"figures{suffix}" / f"walk_length_vs_relative_error{suffix}")


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
        output_match = re.fullmatch(r"out_(hit_and_run|ball_walk)_(?:[A-Za-z0-9_-]+)",
                                    input_folder.name)
        inferred_name = (output_match.group(1) if output_match else
                         input_folder.name[4:] if is_output_folder else
                         input_folder.name)
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
