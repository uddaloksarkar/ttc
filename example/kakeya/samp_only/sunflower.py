#!/usr/bin/env python3
"""Draw a sunflower-like x0/x1 projection of one box-union SMT2 benchmark.

Example: python3 sunflower.py --sampler hit_and_run --benchmark NAME.smt2
"""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from hit_and_run.utility import boxes, number, parse_sexprs

HERE = Path(__file__).resolve().parent
# The one palette for every figure in this folder (accumulate.py imports it).
# Each family is a light -> dark ramp; its middle step is the family's base.
PALETTE = {
    "ink": "#0b0b0b", "muted": "#52514e", "grid": "#e4e3df",
    "n3": ("#c4bcf4", "#9384e8", "#6250d6", "#3a2a9a", "#170d4a"),
    "n_eq_m": ("#8bd9bb", "#4fc496", "#1baf7a", "#11815a", "#074a33"),
    "radial": ("#f6c98b", "#eda34d", "#d97824", "#a95312", "#713209"),
}
STYLE = {
    "font.size": 9, "axes.labelsize": 10, "axes.edgecolor": PALETTE["muted"],
    "axes.labelcolor": PALETTE["ink"], "axes.linewidth": 0.8,
    "xtick.color": PALETTE["muted"], "ytick.color": PALETTE["muted"],
    "text.color": PALETTE["ink"], "axes.spines.top": False,
    "axes.spines.right": False, "legend.frameon": False,
    "pgf.texsystem": "pdflatex",
}


def project(path: Path) -> list[tuple[float, float, float, float]]:
    exprs = parse_sexprs(path.read_text())
    assertions = [expr[1] for expr in exprs if expr[0] == "assert"]
    if not assertions:
        raise ValueError("no SMT2 assertions")
    all_boxes = boxes(["and", *assertions])
    result = []
    for box in all_boxes:
        coords = [box.get((kind, f"x{axis}")) for axis in (0, 1) for kind in ("l", "u")]
        if any(item is None for item in coords):
            raise ValueError(f"missing x0/x1 bound in {path}")
        x0, x1, y0, y1 = map(float, coords)
        if x0 >= x1 or y0 >= y1:
            raise ValueError(f"empty projected body in {path}")
        result.append((x0, x1, y0, y1))
    return result


def box_bounds(path: Path, axes: int = 3) -> list[list[tuple[float, float]]]:
    """(lower, upper) along x0..x{axes-1} for each box of the benchmark."""
    exprs = parse_sexprs(path.read_text())
    assertions = [expr[1] for expr in exprs if expr[0] == "assert"]
    result = []
    for box in boxes(["and", *assertions]):
        bounds = [(box.get(("l", f"x{axis}")), box.get(("u", f"x{axis}")))
                  for axis in range(axes)]
        if any(None in pair for pair in bounds):
            raise ValueError(f"missing bound in {path}")
        result.append([(float(low), float(high)) for low, high in bounds])
    return result


def draw_projection(path: Path, output: Path, color: str | None = None) -> None:
    """Save the x0/x1 projection of a box-union benchmark as PDF, PGF and PNG."""
    rectangles = project(path)
    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(5.4, 5.4), facecolor="white")
    # Every body in one color: the n = 3 base (the projections are 3D inputs).
    color = color or PALETTE["n3"][2]
    for x0, x1, y0, y1 in rectangles:
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                               facecolor=to_rgba(color, 0.18),
                               edgecolor=to_rgba(color, 0.75), linewidth=1.05))
    reach = 1.05 * max(abs(value) for rectangle in rectangles for value in rectangle)
    ax.set_xlim(-reach, reach)
    ax.set_ylim(-reach, reach)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$x_0$", fontsize=11)
    ax.set_ylabel(r"$x_1$", fontsize=11)
    ax.tick_params(labelsize=9, direction="out", length=3)
    fig.tight_layout(pad=0.5)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def draw_boxes_3d(path: Path, output: Path, z_reach: float | None = None) -> None:
    """Save a 3D view of the boxes as PDF, PGF and PNG.

    The x2 axis spans +-z_reach (default: this benchmark's own extent) and is
    drawn stretched, so two benchmarks that differ only in thickness can be
    compared by passing both the same z_reach.
    """
    bounds = box_bounds(path)
    extent = max(abs(v) for box in bounds for pair in box[:2] for v in pair)
    reach = 1.05 * extent
    z_reach = z_reach or max(abs(v) for box in bounds for v in box[2])
    color = PALETTE["n3"][2]
    plt.rcParams.update({**STYLE, "font.size": 24, "axes.labelsize": 24})
    fig = plt.figure(figsize=(7.2, 6.0), facecolor="white")
    ax = fig.add_subplot(projection="3d", computed_zorder=False)
    faces = []
    for (a0, a1), (b0, b1), (c0, c1) in bounds:
        corner = lambda i, j, k: ((a0, a1)[i], (b0, b1)[j], (c0, c1)[k])
        for axis in range(3):
            for side in (0, 1):
                square = [(0, 0), (1, 0), (1, 1), (0, 1)]
                face = []
                for u, v in square:
                    index = [u, v]
                    index.insert(axis, side)
                    face.append(corner(*index))
                faces.append(face)
    ax.add_collection3d(Poly3DCollection(faces, facecolor=to_rgba(color, 0.07),
                                         edgecolor=to_rgba(color, 0.55),
                                         linewidth=0.5))
    ax.set_xlim(-reach, reach)
    ax.set_ylim(-reach, reach)
    ax.set_zlim(-1.1 * z_reach, 1.1 * z_reach)
    ax.set_box_aspect((1, 1, 0.45), zoom=1.5)
    ax.view_init(elev=24, azim=-58)
    ax.set_axis_off()
    fig.subplots_adjust(left=0.01, right=0.86, bottom=0.02, top=0.98)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler", required=True, choices=("ball_walk", "hit_and_run"))
    parser.add_argument("--benchmark", required=True,
                        help="Filename from the selected sampler's catalog.csv")
    args = parser.parse_args()
    folder = HERE / args.sampler
    catalog = {row["benchmark"]: row
               for row in csv.DictReader((folder / "catalog.csv").open())}
    if args.benchmark not in catalog:
        parser.error("benchmark is absent from catalog.csv")
    path = folder / "benchmarks" / args.benchmark
    if not path.is_file():
        path = HERE / catalog[args.benchmark]["original_path"]
    if hashlib.sha256(path.read_bytes()).hexdigest() != catalog[args.benchmark]["sha256"]:
        parser.error(f"benchmark content differs from catalog: {path}")
    output = folder / "figures" / f"sunflower_{path.stem}"
    draw_projection(path, output)
    print(output)


if __name__ == "__main__":
    main()
