#!/usr/bin/env python3
"""Prepare and verify the selected Kakeya sampler benchmark suites."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import re
import sys
import tempfile
from collections import defaultdict
from fractions import Fraction
from pathlib import Path

sys.dont_write_bytecode = True

from utility import (box_body_volume, exact_box_geometry, max_box_radius_ratio_squared,
                     radial_body_volume)

ROOT = Path(__file__).resolve().parent
HERE = ROOT / "hit_and_run"
BENCHMARKS = HERE / "benchmarks"
SPREAD = 100
SHORT = 1
RADIAL_SCALE = 10000
RADIAL_COUNTS = (4, 6, 8, 10, 12, 14, 16, 20, 24, 28,
                 32, 36, 40, 44, 48, 52, 56, 60, 64, 72)
WALK_LENGTHS = (5, 10, 18, 34, 65, 124, 236, 449, 853, 1622,
                3083, 5861, 11143, 21183, 40271, 76559, 145544,
                276692, 526015, 1000000)
BOX_INVENTORY = {
    2: (10, 100, 1000),
    4: (1, 10, 100, 1000),
    8: (1, 10, 55, 56, 57, 59, 75, 100, 1000),
    16: (1, 10, 100, 1000),
}
STAIRCASE_INVENTORY = {
    8: (8, 9, 10, 25, 50, 500),
}
GEOMETRY_FIELDS = ("name", "bodies", "dimension", "thin_denominator",
                   "exact_body_volume", "exact_union_volume",
                   "max_body_R_over_r", "max_body_R_over_r_approx", "benchmark")


def radius_ratio_fields(squared: Fraction) -> dict[str, str]:
    radicand = (str(squared.numerator) if squared.denominator == 1
                else str(squared))
    return {"max_body_R_over_r": f"sqrt({radicand})",
            "max_body_R_over_r_approx": f"{math.sqrt(float(squared)):.10g}"}


def write_box_geometry(rows: list[dict]) -> None:
    manifest = HERE / "geometry.csv"
    existing = ({row["name"]: row for row in csv.DictReader(manifest.open())
                 if row["benchmark"].startswith("benchmarks/")
                 and int(row["thin_denominator"]) < 4_000_000}
                if manifest.exists() else {})
    existing.update({row["name"]: row for row in rows})
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=GEOMETRY_FIELDS)
        writer.writeheader()
        writer.writerows(existing[name] for name in sorted(existing))


# An independent parse of the written fixed-3D SMT2 bounds.
BOX_NAME = re.compile(r"SAMPONLY_thickness_m(\d+)_walk(\d+)_d(\d+)"
                      r"(?:_r(\d+))?(?:_s(\d+))?\.smt2$")
BOX_SIDE_NAME = re.compile(r"SAMPONLY_side_m(\d+)_walk(\d+)_d(\d+)_r(\d+)_"
                           r"s(\d+)over(\d+)\.smt2$")


def box_radius(token):
    token = token.strip()
    fraction = re.fullmatch(r"\(/\s+(\d+)\s+(\d+)\)", token)
    if fraction:
        return Fraction(int(fraction.group(1)), int(fraction.group(2)))
    return Fraction(int(token))


def box_bounds(raw, coordinate):
    upper = re.findall(rf"^\s*\(<= x{coordinate} (.+)\)\s*$", raw, re.M)
    lower = re.findall(rf"^\s*\(<= \(- (.+)\) x{coordinate}\)\s*$", raw, re.M)
    upper, lower = [box_radius(x) for x in upper], [box_radius(x) for x in lower]
    if lower != upper:
        raise ValueError(f"asymmetric x{coordinate} bounds")
    return upper


def verify_box_file(path):
    match = BOX_NAME.fullmatch(path.name)
    side_match = BOX_SIDE_NAME.fullmatch(path.name) if not match else None
    if not match and not side_match:
        raise ValueError(f"unexpected benchmark name: {path.name}")
    if match:
        count, walk, denominator, spread, short = match.groups()
        count, walk, denominator = map(int, (count, walk, denominator))
        spread, short = int(spread or 100), Fraction(int(short or 1))
    else:
        count, walk, denominator, spread, numerator, side_denominator = map(
            int, side_match.groups())
        spread, short = int(spread), Fraction(numerator, side_denominator)
    raw = path.read_text()
    xs, ys, zs = (box_bounds(raw, coordinate) for coordinate in range(3))
    if not (len(xs) == len(ys) == len(zs) == count):
        raise ValueError(f"wrong body count: {path.name}")
    expected_x = [(short * (count - 1) + (spread - short) * i) / (count - 1)
                  for i in range(count)]
    if xs != expected_x:
        raise ValueError(f"wrong x staircase: {path.name}")
    if ys != [spread * short / x for x in xs]:
        raise ValueError(f"wrong y staircase: {path.name}")
    if zs != [Fraction(1, denominator)] * count:
        raise ValueError(f"wrong z thickness: {path.name}")
    bodies = [[x, y, z] for x, y, z in zip(xs, ys, zs)]
    volumes = [box_body_volume(body) for body in bodies]
    expected_body = Fraction(8 * spread, denominator) * short
    if volumes != [expected_body] * count:
        raise ValueError(f"wrong body volume: {path.name}")
    max_ratio_squared = max_box_radius_ratio_squared(bodies)
    # Integrate the parsed boxes over x slabs, independent of the generator's
    # exact_union_quadrant implementation.
    endpoints = [Fraction(0)] + sorted(set(xs))
    quadrant = sum((right - left) * max(y for x, y in zip(xs, ys) if x >= right)
                   for left, right in zip(endpoints, endpoints[1:]))
    union = 8 * zs[0] * quadrant
    return {"name": path.stem, "bodies": count, "walk_length": walk,
            "denominator": denominator, "spread": spread,
            "short_radius": str(short),
            "exact_body_volume": str(expected_body),
            "exact_union_volume": str(union),
            "max_body_radius_ratio_squared": str(max_ratio_squared),
            "verified": 1}


def benchmark_path(name: str) -> Path:
    return BENCHMARKS / name


def smt_positive(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"(/ {value.numerator} {value.denominator})"

def smt_negative(value: Fraction) -> str:
    return f"(- {smt_positive(value)})"

def box_formula(count: int, walk: int, denominator: int | None = None,
                 spread: int = 100, short_radius: int = 1) -> str:
    denominator = denominator if denominator is not None else 20 * walk
    lines = ["(set-logic QF_LRA)",
             "(declare-fun x0 () Real)",
             "(declare-fun x1 () Real)",
             "(declare-fun x2 () Real)",
             "(assert (or"]
    for i in range(count):
        x = Fraction(short_radius * (count - 1) + (spread - short_radius) * i,
                     count - 1)
        y = Fraction(spread * short_radius, 1) / x
        z = Fraction(1, denominator)
        lines.extend(["  (and",
                      f"    (<= {smt_negative(x)} x0)",
                      f"    (<= x0 {smt_positive(x)})",
                      f"    (<= {smt_negative(y)} x1)",
                      f"    (<= x1 {smt_positive(y)})",
                      f"    (<= {smt_negative(z)} x2)",
                      f"    (<= x2 {smt_positive(z)})",
                      "  )"])
    return "\n".join(lines + ["))", "(check-sat)"]) + "\n"

def exact_union_quadrant(count: int, spread: int = 100,
                         short_radius: int = 1) -> Fraction:
    # Sorted x radii increase, while y radii decrease. The union is a
    # staircase: sum (x_i-x_(i-1))*y_i, with x_(-1)=0.
    area = Fraction(0)
    previous_x = Fraction(0)
    for i in range(count):
        x = Fraction(short_radius * (count - 1) + (spread - short_radius) * i,
                     count - 1)
        y = Fraction(spread * short_radius, 1) / x
        area += (x - previous_x) * y
        previous_x = x
    return area

def stem(count: int, denominator: int) -> str:
    return f"RDHR_box_m{count:04d}_d{denominator:08d}_s{SHORT:04d}"

def box_geometry(count: int, denominator: int) -> tuple[Fraction, Fraction]:
    return (box_body_volume([Fraction(SHORT), Fraction(SPREAD),
                             Fraction(1, denominator)]),
            Fraction(8, denominator) * exact_union_quadrant(count, SPREAD, SHORT))

def generate_boxes(counts: list[int], denominators: list[int]) -> None:
    if any(denominator >= 4_000_000 for denominator in denominators):
        raise ValueError("box denominator must be below 4,000,000")
    if any(denominator not in BOX_INVENTORY.get(count, ())
           for count in counts for denominator in denominators):
        raise ValueError("box count/denominator is outside the selected suite")
    BENCHMARKS.mkdir(parents=True, exist_ok=True)
    rows = []
    for denominator in sorted(set(denominators)):
        for count in sorted(set(counts)):
            if count < 2:
                raise ValueError("at least two bodies are needed")
            path = benchmark_path(f"{stem(count, denominator)}.smt2")
            path.parent.mkdir(parents=True, exist_ok=True)
            formula = box_formula(count, 5, denominator, SPREAD, SHORT)
            if path.exists():
                if path.read_text() != formula:
                    raise ValueError(f"existing benchmark differs: {path}")
            else:
                path.write_text(formula)
            # The independent verifier checks a historical filename grammar.
            # Feed it the exact bytes written for this RDHR case under that
            # grammar, without changing or duplicating the benchmark itself.
            with tempfile.TemporaryDirectory() as temporary:
                verification_path = (Path(temporary) /
                                     f"SAMPONLY_thickness_m{count:04d}_walk0000005_"
                                     f"d{denominator:08d}.smt2")
                verification_path.write_bytes(path.read_bytes())
                parsed = verify_box_file(verification_path)
            body, union = box_geometry(count, denominator)
            max_ratio_squared = Fraction(parsed["max_body_radius_ratio_squared"])
            if (int(parsed["bodies"]) != count
                    or Fraction(parsed["exact_body_volume"]) != body
                    or Fraction(parsed["exact_union_volume"]) != union):
                raise ValueError(f"independent geometry check failed: {path}")
            rows.append(dict(name=stem(count, denominator), bodies=count,
                             dimension=3, thin_denominator=denominator,
                             exact_body_volume=str(body), exact_union_volume=str(union),
                             **radius_ratio_fields(max_ratio_squared),
                             benchmark=str(path.relative_to(HERE))))
    write_box_geometry(rows)
    print(f"Prepared and independently verified {len(rows)} box unions", flush=True)

def common_factor(count: int, denominator: int) -> Fraction:
    return (Fraction(1) if count == 2 else
            Fraction(2, denominator) * 2 ** (count - 3))

def exact_volumes(count: int, denominator: int) -> tuple[Fraction, Fraction]:
    factor = common_factor(count, denominator)
    radii = [Fraction(1), Fraction(100)]
    if count >= 3:
        radii += [Fraction(1, denominator)] + [Fraction(1)] * (count - 3)
    return (box_body_volume(radii), 4 * exact_union_quadrant(count) * factor)

def staircase_formula(count: int, denominator: int) -> str:
    lines = ["(set-logic QF_LRA)"]
    lines += [f"(declare-fun x{j} () Real)" for j in range(count)]
    lines.append("(assert (or")
    for i in range(count):
        x = Fraction(count - 1 + 99 * i, count - 1)
        y = Fraction(100) / x
        radii = [x, y] + ([Fraction(1, denominator)] +
                           [Fraction(1)] * (count - 3) if count >= 3 else [])
        lines.append("  (and")
        for j, radius in enumerate(radii):
            lines.append(f"    (<= {smt_negative(radius)} x{j})")
            lines.append(f"    (<= x{j} {smt_positive(radius)})")
        lines.append("  )")
    return "\n".join(lines + ["))", "(check-sat)"]) + "\n"

def parse_fraction(token: str) -> Fraction:
    match = re.fullmatch(r"\(/\s+(\d+)\s+(\d+)\)", token.strip())
    return Fraction(int(match[1]), int(match[2])) if match else Fraction(token)

def verify_staircase(path: Path, count: int, denominator: int) -> tuple[Fraction, Fraction, Fraction]:
    raw = path.read_text()
    declared = re.findall(r"^\(declare-fun x(\d+) \(\) Real\)$", raw, re.M)
    if declared != [str(j) for j in range(count)]:
        raise ValueError(f"variable count differs: {path}")
    bounds = []
    for j in range(count):
        upper = [parse_fraction(x) for x in re.findall(
            rf"^\s*\(<= x{j} (.+)\)$", raw, re.M)]
        lower = [parse_fraction(x) for x in re.findall(
            rf"^\s*\(<= \(- (.+)\) x{j}\)$", raw, re.M)]
        if len(upper) != count or upper != lower:
            raise ValueError(f"body bounds differ: {path}, coordinate {j}")
        bounds.append(upper)
    xs, ys = bounds[:2]
    if xs != [Fraction(count - 1 + 99 * i, count - 1)
              for i in range(count)] or ys != [Fraction(100) / x for x in xs]:
        raise ValueError(f"staircase geometry differs: {path}")
    for j in range(2, count):
        expected = Fraction(1, denominator) if j == 2 else Fraction(1)
        if bounds[j] != [expected] * count:
            raise ValueError(f"common coordinates differ: {path}")
    factor = Fraction(1)
    for j in range(2, count):
        factor *= 2 * bounds[j][0]
    bodies = [[bounds[j][i] for j in range(count)] for i in range(count)]
    body = box_body_volume(bodies[0])
    if any(box_body_volume(item) != body for item in bodies):
        raise ValueError(f"individual body volumes differ: {path}")
    endpoints = [Fraction(0)] + sorted(set(xs))
    quadrant = sum((right - left) * max(y for x, y in zip(xs, ys) if x >= right)
                   for left, right in zip(endpoints, endpoints[1:]))
    union = 4 * quadrant * factor
    if (body, union) != exact_volumes(count, denominator):
        raise ValueError(f"exact integral differs: {path}")
    return body, union, max_box_radius_ratio_squared(bodies)

def walk_lengths() -> tuple[int, ...]:
    if len(WALK_LENGTHS) != len(RADIAL_COUNTS):
        raise ValueError("radial counts and target walks differ")
    return WALK_LENGTHS

def smt_number(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator) if value >= 0 else f"(- {-value.numerator})"
    quotient = f"(/ {abs(value.numerator)} {value.denominator})"
    return quotient if value >= 0 else f"(- {quotient})"

def linear_form(coefficients: list[Fraction]) -> str:
    terms = []
    for index, coefficient in enumerate(coefficients):
        if coefficient == 0:
            continue
        variable = f"x{index}"
        if coefficient == 1:
            terms.append(variable)
        elif coefficient == -1:
            terms.append(f"(- {variable})")
        else:
            terms.append(f"(* {smt_number(coefficient)} {variable})")
    return terms[0] if len(terms) == 1 else f"(+ {' '.join(terms)})"

def directions(n: int) -> list[tuple[Fraction, Fraction]]:
    """Rational approximations to evenly spaced directions in the plane."""
    result = []
    for index in range(n):
        angle = 2 * math.pi * (index + 0.5) / n
        result.append((Fraction(round(RADIAL_SCALE * math.cos(angle)), RADIAL_SCALE),
                       Fraction(round(RADIAL_SCALE * math.sin(angle)), RADIAL_SCALE)))
    return result

def inverse_forms(vectors: list[tuple[Fraction, Fraction]]) -> tuple[list[str], Fraction]:
    """x=P y, with x_j=y_j for j>=2; derive y as exact linear forms in x."""
    n = len(vectors)
    (a0, b0), (a1, b1) = vectors[:2]
    determinant = a0 * b1 - a1 * b0
    if determinant <= 0:
        raise ValueError("the planar change of coordinates is singular")
    y0 = [b1 / determinant, -a1 / determinant]
    y1 = [-b0 / determinant, a0 / determinant]
    for a, b in vectors[2:]:
        y0.append((-b1 * a + a1 * b) / determinant)
        y1.append((b0 * a - a0 * b) / determinant)
    forms = [linear_form(y0), linear_form(y1)]
    for index in range(2, n):
        forms.append(f"x{index}")
    return forms, determinant

def benchmark_text(n: int, length: int, forms: list[str]) -> str:
    lines = ["(set-logic QF_LRA)",
             f"; {n} dimensions, {n} slanted bodies; exact y-to-x determinant in geometry.csv"]
    lines += [f"(declare-fun x{index} () Real)" for index in range(n)]
    lines.append("(assert (or")
    for extended in range(n):
        lines.append("  (and")
        for index, form in enumerate(forms):
            lines.append(f"    (<= (- 1) {form})")
            lines.append(f"    (<= {form} {length if index == extended else 1})")
        lines.append("  )")
    lines += ["))", "(check-sat)"]
    return "\n".join(lines) + "\n"

def generate_radial(counts: tuple[int, ...]) -> None:
    benchmark_dir = HERE / "benchmarks"
    benchmark_dir.mkdir(exist_ok=True)
    targets = dict(zip(RADIAL_COUNTS, walk_lengths()))
    records = []
    for n in counts:
        if n < 3:
            raise ValueError("n must be at least 3")
        length = 2 * n + 1
        vectors = directions(n)
        forms, determinant = inverse_forms(vectors)
        name = f"RDHR_radial_urchin_m{n:04d}_n{n:04d}.smt2"
        path = benchmark_dir / name
        expected = benchmark_text(n, length, forms)
        if path.exists() and path.read_text() != expected:
            raise ValueError(f"existing benchmark differs: {path}")
        if not path.exists():
            path.write_text(expected)
        core = determinant * 2**n
        body = radial_body_volume(determinant, length, n)
        union = core + n * determinant * (length - 1) * 2**(n - 1)
        if determinant <= 0 or body <= 0 or union <= body:
            raise ValueError(f"invalid exact geometry for n={n}")
        records.append(dict(benchmark=name, bodies=n, dimension=n, length=length,
                            target_walk_length=targets.get(n, ""),
                            determinant=determinant, exact_core_volume=core,
                            exact_body_volume=body, exact_union_volume=union))
        print(path)
    manifest = HERE / "radial_geometry.csv"
    existing = {}
    if manifest.exists():
        with manifest.open(newline="") as stream:
            existing = {row["benchmark"]: row for row in csv.DictReader(stream)}
    existing.update({row["benchmark"]: row for row in records})
    with manifest.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(sorted(existing.values(), key=lambda row: int(row["bodies"])))

def parse_counts(raw: str) -> tuple[int, ...]:
    counts = tuple(int(item.strip()) for item in raw.split(","))
    if not counts or any(count < 2 for count in counts) or len(counts) != len(set(counts)):
        raise ValueError("counts must be distinct integers of at least 2")
    return counts


def write_verified(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_text() != content:
            raise ValueError(f"existing benchmark differs: {path}")
    else:
        path.write_text(content)


def generate_staircase(counts: tuple[int, ...], denominator: int) -> None:
    if denominator >= 4_000_000:
        raise ValueError("box denominator must be below 4,000,000")
    if any(denominator not in STAIRCASE_INVENTORY.get(count, ()) for count in counts):
        raise ValueError("staircase count/denominator is outside the selected suite")
    rows = []
    for count in counts:
        path = benchmark_path(f"RDHR_dimeqbody_m{count:04d}_d{denominator:08d}.smt2")
        write_verified(path, staircase_formula(count, denominator))
        body, union, ratio_squared = verify_staircase(path, count, denominator)
        rows.append(dict(name=path.stem, bodies=count, dimension=count,
                         thin_denominator=denominator,
                         exact_body_volume=str(body), exact_union_volume=str(union),
                         **radius_ratio_fields(ratio_squared),
                         benchmark=str(path.relative_to(HERE))))
        print(f"{path}: n=m={count}, exact body={body}, exact union={union}")
    write_box_geometry(rows)


def regenerate_inventory() -> None:
    by_denominator: dict[int, list[int]] = defaultdict(list)
    for count, denominators in BOX_INVENTORY.items():
        for denominator in denominators:
            by_denominator[denominator].append(count)
    for denominator, counts in sorted(by_denominator.items()):
        generate_boxes(counts, [denominator])
    for count, denominators in STAIRCASE_INVENTORY.items():
        for denominator in denominators:
            generate_staircase((count,), denominator)
    generate_radial(RADIAL_COUNTS)


def prepare_ball_walk() -> None:
    """Verify curated SMT2 inputs and rebuild their exact geometry metadata."""
    folder = ROOT / "ball_walk"
    with (folder / "catalog.csv").open(newline="") as stream:
        catalog = list(csv.DictReader(stream))
    with (folder / "cases.csv").open(newline="") as stream:
        cases = list(csv.DictReader(stream))
    expected = defaultdict(set)
    for case in cases:
        expected[case["benchmark"]].add(Fraction(case["exact_union_volume"]))
    records = []
    for meta in catalog:
        path = folder / "benchmarks" / meta["benchmark"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != meta["sha256"]:
            raise ValueError(f"benchmark differs from catalog: {path}")
        body_volumes, union, dimension = exact_box_geometry(path)
        if (len(body_volumes) != int(meta["bodies"])
                or dimension != int(meta["dimension"])
                or expected[meta["benchmark"]] != {union}):
            raise ValueError(f"exact geometry differs from catalog or cases: {path}")
        records.append(dict(benchmark=meta["benchmark"], bodies=len(body_volumes),
                            dimension=dimension,
                            exact_body_volumes=";".join(map(str, body_volumes)),
                            exact_union_volume=str(union)))
    if set(expected) != {meta["benchmark"] for meta in catalog}:
        raise ValueError("ball-walk cases and catalog differ")
    with (folder / "geometry.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=records[0])
        writer.writeheader()
        writer.writerows(records)
    print(f"Verified {len(records)} ball-walk benchmarks and exact volumes")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", nargs="?", default="all",
                        choices=("all", "box", "staircase", "radial"),
                        help="benchmark family (default: all)")
    parser.add_argument("--sampler", choices=("both", "hit_and_run", "ball_walk"),
                        default="both", help="sampler suite (default: both)")
    parser.add_argument("--counts", help="comma-separated body counts for one family")
    parser.add_argument("--denominator", type=int, default=500,
                        help="thin-coordinate denominator for staircase")
    args = parser.parse_args()
    if not 1 <= args.denominator < 4_000_000:
        parser.error("denominator must be between 1 and 3,999,999")
    if args.mode == "all" and args.counts:
        parser.error("--counts requires one specific mode")
    if args.mode in ("box", "radial", "all") and args.denominator != 500:
        parser.error("--denominator applies only to staircase")
    if args.sampler == "ball_walk" and (args.mode != "all" or args.counts):
        parser.error("ball-walk preparation uses the complete curated suite")
    try:
        counts = parse_counts(args.counts) if args.counts else None
    except ValueError as error:
        parser.error(str(error))
    if args.sampler == "ball_walk":
        prepare_ball_walk()
    elif args.mode == "all":
        regenerate_inventory()
    elif args.mode == "box":
        selected = set(counts) if counts else set(BOX_INVENTORY)
        if selected - BOX_INVENTORY.keys():
            parser.error("box counts must be among 2,4,8,16")
        by_denominator: dict[int, list[int]] = defaultdict(list)
        for count in sorted(selected):
            for denominator in BOX_INVENTORY[count]:
                by_denominator[denominator].append(count)
        for denominator, members in sorted(by_denominator.items()):
            generate_boxes(members, [denominator])
    elif args.mode == "staircase":
        selected = counts or tuple(STAIRCASE_INVENTORY)
        if any(args.denominator not in STAIRCASE_INVENTORY.get(count, ())
               for count in selected):
            parser.error("staircase count/denominator is outside the selected suite")
        generate_staircase(selected, args.denominator)
    else:
        selected = counts or RADIAL_COUNTS
        if set(selected) - set(RADIAL_COUNTS):
            parser.error("radial counts must come from the selected suite")
        generate_radial(selected)
    if args.sampler == "both":
        prepare_ball_walk()


if __name__ == "__main__":
    main()
