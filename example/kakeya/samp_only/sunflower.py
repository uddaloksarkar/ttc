#!/usr/bin/env python3
"""Draw a sunflower-like x0/x1 projection of one box-union SMT2 benchmark.

Example: python3 sunflower.py --sampler hit_and_run --benchmark NAME.smt2
"""

from __future__ import annotations

import argparse
import csv
import re
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent


def parse_sexprs(text: str) -> list:
    tokens = re.findall(r"\(|\)|[^\s()]+", re.sub(r";[^\n]*", "", text))
    stack: list[list] = []
    result = []
    for token in tokens:
        if token == "(":
            item: list = []
            (stack[-1] if stack else result).append(item)
            stack.append(item)
        elif token == ")":
            if not stack:
                raise ValueError("unbalanced SMT2 closing parenthesis")
            stack.pop()
        else:
            if not stack:
                raise ValueError("token outside SMT2 expression")
            stack[-1].append(token)
    if stack:
        raise ValueError("unbalanced SMT2 opening parenthesis")
    return result


def number(expr) -> Fraction:
    if isinstance(expr, str):
        return Fraction(expr)
    if expr[0] == "-" and len(expr) == 2:
        return -number(expr[1])
    if expr[0] == "/" and len(expr) == 3:
        return number(expr[1]) / number(expr[2])
    raise ValueError(f"unsupported numeric SMT2 expression: {expr}")


def merge(left: dict, right: dict) -> dict:
    output = dict(left)
    for key, value in right.items():
        if key in output:
            output[key] = max(output[key], value) if key[0] == "l" else min(output[key], value)
        else:
            output[key] = value
    return output


def boxes(expr) -> list[dict]:
    operator = expr[0]
    if operator == "or":
        return [box for child in expr[1:] for box in boxes(child)]
    if operator == "and":
        result = [{}]
        for child in expr[1:]:
            result = [merge(existing, addition) for existing in result for addition in boxes(child)]
        return result
    if operator == "<=":
        left, right = expr[1:]
        if isinstance(left, str) and left.startswith("x"):
            return [{("u", left): number(right)}]
        if isinstance(right, str) and right.startswith("x"):
            return [{("l", right): number(left)}]
    raise ValueError(f"benchmark is not a supported union of boxes: {expr}")


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


def exact_box_volumes(path: Path) -> list[Fraction]:
    """Calculate each box volume directly from the SMT2 inequalities."""
    exprs = parse_sexprs(path.read_text())
    variables = [expr[1] for expr in exprs if expr[0] == "declare-fun"
                 and expr[-1] == "Real"]
    assertions = [expr[1] for expr in exprs if expr[0] == "assert"]
    output = []
    for box in boxes(["and", *assertions]):
        volume = Fraction(1)
        for variable in variables:
            lower = box.get(("l", variable))
            upper = box.get(("u", variable))
            if lower is None or upper is None or lower >= upper:
                raise ValueError(f"non-box or missing {variable} bound: {path}")
            volume *= upper - lower
        output.append(volume)
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sampler", required=True, choices=("ball_walk", "hit_and_run"))
    parser.add_argument("--benchmark", required=True,
                        help="Filename from the selected sampler's catalog.csv")
    args = parser.parse_args()
    folder = HERE / args.sampler
    catalog = {row["benchmark"] for row in csv.DictReader((folder / "catalog.csv").open())}
    if args.benchmark not in catalog:
        parser.error("benchmark is absent from catalog.csv")
    path = folder / "benchmarks" / args.benchmark
    rectangles = project(path)
    fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    colors = plt.get_cmap("hsv", len(rectangles) + 1)
    for index, (x0, x1, y0, y1) in enumerate(rectangles):
        color = colors(index)
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                               facecolor=color, edgecolor=color,
                               linewidth=1.2, alpha=0.12))
        ax.plot([0, x1], [0, y1], color=color, alpha=0.5, linewidth=0.8)
        ax.plot([0, x1], [0, y0], color=color, alpha=0.5, linewidth=0.8)
        ax.plot([0, x0], [0, y1], color=color, alpha=0.5, linewidth=0.8)
        ax.plot([0, x0], [0, y0], color=color, alpha=0.5, linewidth=0.8)
    ax.axhline(0, color="black", linewidth=0.5)
    ax.axvline(0, color="black", linewidth=0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.autoscale_view()
    ax.set_xlabel("x0")
    ax.set_ylabel("x1")
    ax.set_title(f"Projected box union: {len(rectangles)} bodies")
    output = folder / "figures" / f"sunflower_{path.stem}"
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), bbox_inches="tight", dpi=220)
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()
