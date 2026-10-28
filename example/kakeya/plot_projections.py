#!/usr/bin/env python3
"""Plot selected Kakeya box and radial benchmark geometry."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import sys
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon, Rectangle
from matplotlib.colors import to_rgba
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

sys.dont_write_bytecode = True

from utility import boxes, number, parse_sexprs

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


def figure_directory(folder: Path, suffix: str = "") -> Path:
    """Keep figures inside their sampler campaign directory."""
    return folder / f"figures{suffix}"


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


def projection_polygon(vectors: list[tuple[Fraction, Fraction]],
                       extended: int, length: int) -> list[tuple[float, float]]:
    """Vertices of the 2D zonotope, without enumerating 2^n box corners."""
    center = [0.0, 0.0]
    generators = []
    for index, (a, b) in enumerate(vectors):
        lo, hi = -1, length if index == extended else 1
        center[0] += float(a) * (lo + hi) / 2
        center[1] += float(b) * (lo + hi) / 2
        dx, dy = float(a) * (hi - lo) / 2, float(b) * (hi - lo) / 2
        if dy < 0 or (dy == 0 and dx < 0):
            dx, dy = -dx, -dy
        generators.append((dx, dy))
    generators.sort(key=lambda g: math.atan2(g[1], g[0]))
    x = center[0] - sum(g[0] for g in generators)
    y = center[1] - sum(g[1] for g in generators)
    vertices = [(x, y)]
    for sign in (1, -1):
        for dx, dy in generators:
            x += sign * 2 * dx
            y += sign * 2 * dy
            vertices.append((x, y))
    return vertices[:-1]


def plot(path: Path, vectors: list[tuple[Fraction, Fraction]], length: int) -> None:
    plt.rcParams.update({**STYLE, "font.size": 24, "axes.labelsize": 24})
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    palette = PALETTE["radial"]
    extent = 0.0
    for index in range(len(vectors)):
        vertices = projection_polygon(vectors, index, length)
        extent = max(extent, *(abs(c) for vertex in vertices for c in vertex))
        color = palette[round(4 * index / max(len(vectors) - 1, 1))]
        ax.add_patch(Polygon(vertices, closed=True,
                             facecolor=to_rgba(color, 0.13),
                             edgecolor=to_rgba(color, 0.75), linewidth=0.8))
    limit = extent * 1.04
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()
    fig.tight_layout(pad=0.5)
    for suffix in ("pdf", "pgf", "png"):
        fig.savefig(path.with_suffix(f".{suffix}"), dpi=300)
    plt.close(fig)


def plot_box_width(path: Path, output: Path, x2_limit: float,
                   horizontal_axis: int) -> None:
    """Project 3D boxes onto x0/x2 or x1/x2 with a shared thin-axis scale."""
    bounds = box_bounds(path)
    plt.rcParams.update(STYLE)
    fig, ax = plt.subplots(figsize=(5.4, 3.7), facecolor="white")
    color = PALETTE["n3"][2]
    for box in bounds:
        lower, upper = box[horizontal_axis]
        z0, z1 = box[2]
        ax.add_patch(Rectangle((lower, z0), upper - lower, z1 - z0,
                               facecolor=to_rgba(color, 0.16),
                               edgecolor=to_rgba(color, 0.8), linewidth=1.0))
    reach = 1.05 * max(abs(x) for box in bounds for x in box[horizontal_axis])
    ax.set_xlim(-reach, reach)
    ax.set_ylim(-x2_limit, x2_limit)
    ax.set_xlabel(rf"$x_{horizontal_axis}$")
    ax.set_ylabel(r"$x_2$")
    ax.grid(False)
    ax.tick_params(labelsize=9, direction="out", length=3)
    fig.tight_layout(pad=0.6)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), dpi=300)
    plt.close(fig)


def linear(expr: object) -> dict[int, Fraction]:
    if isinstance(expr, str) and expr.startswith("x") and expr[1:].isdigit():
        return {int(expr[1:]): Fraction(1)}
    if not isinstance(expr, list) or not expr:
        raise ValueError(f"unsupported linear form: {expr}")
    if expr[0] == "+":
        result: dict[int, Fraction] = {}
        for term in expr[1:]:
            for index, value in linear(term).items():
                result[index] = result.get(index, Fraction()) + value
        return result
    if expr[0] == "-" and len(expr) == 2:
        return {index: -value for index, value in linear(expr[1]).items()}
    if expr[0] == "*" and len(expr) == 3:
        scale = number(expr[1])
        return {index: scale * value for index, value in linear(expr[2]).items()}
    raise ValueError(f"unsupported linear form: {expr}")


def projection_geometry(path: Path) -> tuple[list[tuple[Fraction, Fraction]], int]:
    exprs = parse_sexprs(path.read_text())
    declarations = [expr[1] for expr in exprs if expr[0] == "declare-fun"]
    n = len(declarations)
    if n < 3 or declarations != [f"x{index}" for index in range(n)]:
        raise ValueError(f"invalid variable declarations: {path}")
    assertions = [expr[1] for expr in exprs if expr[0] == "assert"]
    if len(assertions) != 1 or assertions[0][0] != "or":
        raise ValueError(f"expected one union assertion: {path}")
    bodies = assertions[0][1:]
    if len(bodies) != n:
        raise ValueError(f"expected {n} bodies: {path}")
    forms = None
    extensions = []
    length = None
    for body in bodies:
        if body[0] != "and" or len(body[1:]) != 2 * n:
            raise ValueError(f"invalid body bounds: {path}")
        current_forms = []
        upper = []
        for index in range(n):
            lower_bound, upper_bound = body[1 + 2 * index:3 + 2 * index]
            if (lower_bound[0] != "<=" or upper_bound[0] != "<="
                    or number(lower_bound[1]) != -1
                    or lower_bound[2] != upper_bound[1]):
                raise ValueError(f"inconsistent facet bounds: {path}")
            current_forms.append(linear(upper_bound[1]))
            upper.append(number(upper_bound[2]))
        if forms is None:
            forms = current_forms
        elif current_forms != forms:
            raise ValueError(f"body facets differ: {path}")
        extended = [index for index, value in enumerate(upper) if value > 1]
        if len(extended) != 1 or any(value < 1 for value in upper):
            raise ValueError(f"expected one extended coordinate per body: {path}")
        extension = extended[0]
        if any(value != 1 for index, value in enumerate(upper) if index != extension):
            raise ValueError(f"unexpected coordinate bound: {path}")
        if length is None:
            length = upper[extension]
        elif upper[extension] != length:
            raise ValueError(f"unequal arm lengths: {path}")
        extensions.append(extension)
    if sorted(extensions) != list(range(n)) or length.denominator != 1:
        raise ValueError(f"incomplete radial arms or noninteger arm length: {path}")
    for index in range(2, n):
        if forms[index] != {index: Fraction(1)}:
            raise ValueError(f"unexpected higher-dimensional facet: {path}")
    a, b = forms[:2]
    determinant = a.get(0, 0) * b.get(1, 0) - a.get(1, 0) * b.get(0, 0)
    if determinant == 0:
        raise ValueError(f"singular planar transform: {path}")
    vectors = []
    for index in range(n):
        first = Fraction(index == 0) - a.get(index, 0) if index >= 2 else Fraction(index == 0)
        second = Fraction(index == 1) - b.get(index, 0) if index >= 2 else Fraction(index == 1)
        vectors.append(((b.get(1, 0) * first - a.get(1, 0) * second) / determinant,
                        (-b.get(0, 0) * first + a.get(0, 0) * second) / determinant))
    return vectors, int(length)



def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", nargs="?", type=Path,
                        help="optional out_<sampler>_<jobid> directory")
    parser.add_argument("--sampler", choices=("hit_and_run", "ball_walk"),
                        help="required when no output directory is given")
    parser.add_argument("--benchmark", action="append", required=True,
                        help="selected SMT2 filename; repeat to plot more")
    parser.add_argument("--plane", choices=("x0x1", "x0x2", "x1x2", "3d"),
                        default="x0x1")
    parser.add_argument("--x2-limit", type=float,
                        help="shared positive x2 limit for box x0x2, x1x2, or 3d views")
    args = parser.parse_args()
    if args.x2_limit is not None and args.x2_limit <= 0:
        parser.error("--x2-limit must be positive")
    if args.plane == "x0x1" and args.x2_limit is not None:
        parser.error("--x2-limit requires --plane x0x2, x1x2, or 3d")
    output_dir = args.output_dir.resolve() if args.output_dir else None
    if output_dir:
        if not output_dir.is_dir():
            parser.error(f"not a directory: {output_dir}")
        match = re.fullmatch(r"out_(hit_and_run|ball_walk)_([A-Za-z0-9_-]+)",
                             output_dir.name)
        if not match:
            parser.error("expected out_hit_and_run_<jobid> or out_ball_walk_<jobid>")
        sampler, job_id = match.groups()
        if args.sampler and args.sampler != sampler:
            parser.error("--sampler differs from output directory")
        suffix = f"_{job_id}"
    else:
        if not args.sampler:
            parser.error("--sampler is required without an output directory")
        sampler, suffix = args.sampler, ""
    folder = HERE / sampler
    with (folder / "catalog.csv").open(newline="") as stream:
        catalog = {row["benchmark"]: row for row in csv.DictReader(stream)}
    radial = {}
    if sampler == "hit_and_run":
        with (folder / "radial_geometry.csv").open(newline="") as stream:
            radial = {row["benchmark"]: row for row in csv.DictReader(stream)}
    figure_dir = figure_directory(folder, suffix)
    figure_dir.mkdir(exist_ok=True)
    for raw_name in args.benchmark:
        name = Path(raw_name).name
        if name not in catalog and name not in radial:
            parser.error(f"benchmark absent from selected catalog: {name}")
        if output_dir and not any(output_dir.glob(f"{Path(name).stem}_walk*.out.xz")):
            parser.error(f"no compressed output for benchmark: {name}")
        path = folder / "benchmarks" / name
        if name in radial:
            if args.plane != "x0x1":
                parser.error("radial projections support only x0x1")
            output = figure_dir / f"{path.stem}_projection{suffix}"
            vectors, length = projection_geometry(path)
            plot(output, vectors, length)
        else:
            meta = catalog[name]
            if not path.is_file():
                path = HERE / meta["original_path"]
            if hashlib.sha256(path.read_bytes()).hexdigest() != meta["sha256"]:
                raise ValueError(f"benchmark differs from catalog: {path}")
            if args.plane != "x0x1" and int(meta["dimension"]) != 3:
                parser.error(f"{args.plane} requires a 3D box")
            plane_suffix = f"_{args.plane}" if args.plane != "x0x1" else ""
            output = figure_dir / f"{path.stem}{plane_suffix}_projection{suffix}"
            if args.plane in ("x0x2", "x1x2", "3d"):
                bounds = box_bounds(path)
                thin_radius = max(abs(z) for box in bounds for z in box[2])
                if args.x2_limit is not None and args.x2_limit < thin_radius:
                    parser.error("--x2-limit would crop a box")
                limit = args.x2_limit or 1.05 * thin_radius
                if args.plane == "3d":
                    draw_boxes_3d(path, output, z_reach=limit)
                else:
                    plot_box_width(path, output, limit,
                                   horizontal_axis=0 if args.plane == "x0x2" else 1)
            else:
                color = PALETTE["n3" if int(meta["dimension"]) == 3 else "n_eq_m"][2]
                draw_projection(path, output, color=color)
        print(output)


if __name__ == "__main__":
    main()
