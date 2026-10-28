"""Exact box and radial geometry for the Kakeya sampler benchmarks."""

from __future__ import annotations

import re
from fractions import Fraction
from functools import lru_cache
from pathlib import Path


def parse_sexprs(text: str) -> list:
    """Parse the small SMT2 subset used by the box benchmarks."""
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
            result = [merge(existing, addition)
                      for existing in result for addition in boxes(child)]
        return result
    if operator == "<=":
        left, right = expr[1:]
        if isinstance(left, str) and left.startswith("x"):
            return [{("u", left): number(right)}]
        if isinstance(right, str) and right.startswith("x"):
            return [{("l", right): number(left)}]
    raise ValueError(f"benchmark is not a supported union of boxes: {expr}")


def box_body_volume(half_widths: list[Fraction]) -> Fraction:
    """Volume of a centered axis-aligned box with rational half-widths."""
    if not half_widths or any(width <= 0 for width in half_widths):
        raise ValueError("box half-widths must be positive")
    volume = Fraction(1)
    for width in half_widths:
        volume *= 2 * width
    return volume


def box_radius_ratio_squared(half_widths: list[Fraction]) -> Fraction:
    """Exact (circumradius / inradius)^2 for one centered box."""
    if not half_widths or any(width <= 0 for width in half_widths):
        raise ValueError("box half-widths must be positive")
    return sum(width * width for width in half_widths) / min(half_widths) ** 2


def max_box_radius_ratio_squared(bodies: list[list[Fraction]]) -> Fraction:
    """Maximum squared radius ratio across individual boxes."""
    if not bodies:
        raise ValueError("at least one box is required")
    return max(box_radius_ratio_squared(body) for body in bodies)


def radial_body_volume(determinant: Fraction, length: int,
                       dimension: int) -> Fraction:
    """Exact volume of one sea urchin body after its linear transform."""
    if determinant <= 0 or length <= 0 or dimension < 2:
        raise ValueError("invalid radial body geometry")
    return determinant * (length + 1) * 2 ** (dimension - 1)


def exact_box_volumes(path: Path) -> list[Fraction]:
    """Parse SMT2 inequalities and return the exact volume of each box."""
    expressions = parse_sexprs(path.read_text())
    variables = [expr[1] for expr in expressions if expr[0] == "declare-fun"
                 and expr[-1] == "Real"]
    assertions = [expr[1] for expr in expressions if expr[0] == "assert"]
    if not variables or not assertions:
        raise ValueError(f"missing variables or assertion in {path}")
    output = []
    for body in boxes(["and", *assertions]):
        half_widths = []
        for variable in variables:
            lower = body.get(("l", variable))
            upper = body.get(("u", variable))
            if lower is None or upper is None or lower >= upper:
                raise ValueError(f"non-box or missing {variable} bound: {path}")
            half_widths.append((upper - lower) / 2)
        output.append(box_body_volume(half_widths))
    return output


def exact_box_geometry(path: Path) -> tuple[list[Fraction], Fraction, int]:
    """Return exact body volumes, union volume, and dimension from box SMT2."""
    expressions = parse_sexprs(path.read_text())
    variables = [expr[1] for expr in expressions if expr[0] == "declare-fun"
                 and expr[-1] == "Real"]
    assertions = [expr[1] for expr in expressions if expr[0] == "assert"]
    if not variables or not assertions:
        raise ValueError(f"missing variables or assertion in {path}")
    bodies_expr = boxes(["and", *assertions])
    bounds = tuple(tuple((body.get(("l", variable)), body.get(("u", variable)))
                         for variable in variables) for body in bodies_expr)
    if any(low is None or high is None or low >= high
           for body in bounds for low, high in body):
        raise ValueError(f"invalid box bounds in {path}")
    body_volumes = []
    for body in bounds:
        volume = Fraction(1)
        for low, high in body:
            volume *= high - low
        body_volumes.append(volume)

    @lru_cache(None)
    def union(active: tuple[int, ...], axis: int) -> Fraction:
        if axis == len(variables):
            return Fraction(1)
        endpoints = sorted({value for index in active for value in bounds[index][axis]})
        volume = Fraction(0)
        for low, high in zip(endpoints, endpoints[1:]):
            covering = tuple(index for index in active
                             if bounds[index][axis][0] <= low
                             and bounds[index][axis][1] >= high)
            if covering:
                volume += (high - low) * union(covering, axis + 1)
        return volume

    return body_volumes, union(tuple(range(len(bounds))), 0), len(variables)
