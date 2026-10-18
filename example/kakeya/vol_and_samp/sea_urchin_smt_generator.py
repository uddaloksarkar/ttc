#!/usr/bin/env python3
"""
Generate a convex "sea-urchin" SMT-LIB benchmark for QF_LRA.

Construction in y-coordinates:
    C = [-1,1]^n.
For each coordinate i and sign +/- there is one convex body K_{i,sign}.
K_{i,+} extends C from y_i=1 out to y_i=L.
K_{i,-} extends C from y_i=-1 out to y_i=-L.
All other coordinates remain in [-1,1].

Thus there are exactly 2n convex polytopes, every pair intersects exactly
in the common core C, and the petals are pairwise disjoint.

beta is the stationary volume fraction of the petal in each body:
    beta = (L-1)/(L+1),  so  L = (1+beta)/(1-beta).

The exact union volume is
    2^n * (1 + 2n*beta/(1-beta)).

With --rotate, a fixed rational orthogonal transformation is applied,
so volumes, Euclidean R/r, and hit-and-run geometry are unchanged while
the inequalities are no longer axis-aligned.

--subset writes an extra union file over a chosen set of bodies.  Since the
petals are pairwise disjoint and every body contains the core, a union of m
distinct bodies has exact volume 2^n * (1 + m*beta/(1-beta)), which is
recorded in metadata.txt next to the file it belongs to.
"""

import argparse
from fractions import Fraction
from pathlib import Path


def frac(s):
    return Fraction(str(s))


def smt_num(q: Fraction) -> str:
    if q.denominator == 1:
        return str(q.numerator)
    if q.numerator < 0:
        return f"(- (/ {-q.numerator} {q.denominator}))"
    return f"(/ {q.numerator} {q.denominator})"


def identity(n):
    return [[Fraction(int(i == j), 1) for j in range(n)] for i in range(n)]


def apply_givens_rows(M, i, j):
    """Left-multiply rows i,j by [[3/5,4/5],[-4/5,3/5]]."""
    a, b = Fraction(3, 5), Fraction(4, 5)
    ri, rj = M[i][:], M[j][:]
    M[i] = [a*x + b*y for x, y in zip(ri, rj)]
    M[j] = [-b*x + a*y for x, y in zip(ri, rj)]


def rational_rotation(n):
    # Two exact 3-4-5 Givens layers, enough to destroy simple axis alignment.
    M = identity(n)
    for i in range(0, n - 1, 2):
        apply_givens_rows(M, i, i + 1)
    for i in range(1, n - 1, 2):
        apply_givens_rows(M, i, i + 1)
    return M


def linexpr(coeffs, vars_):
    terms = []
    for c, x in zip(coeffs, vars_):
        if c == 0:
            continue
        if c == 1:
            terms.append(x)
        elif c == -1:
            terms.append(f"(- {x})")
        else:
            terms.append(f"(* {smt_num(c)} {x})")
    if not terms:
        return "0"
    if len(terms) == 1:
        return terms[0]
    return "(+ " + " ".join(terms) + ")"


def bounds_for_body(n, axis, sign, L, forms):
    cons = []
    one = Fraction(1, 1)
    for j in range(n):
        y = forms[j]
        if j != axis:
            cons.append(f"(<= (- 1) {y})")
            cons.append(f"(<= {y} 1)")
        elif sign == "+":
            cons.append(f"(<= (- 1) {y})")
            cons.append(f"(<= {y} {smt_num(L)})")
        else:
            cons.append(f"(<= {smt_num(-L)} {y})")
            cons.append(f"(<= {y} 1)")
    return "(and\n      " + "\n      ".join(cons) + "\n    )"


def preamble(n):
    lines = ["(set-logic QF_LRA)"]
    for i in range(n):
        lines.append(f"(declare-fun x{i} () Real)")
    return "\n".join(lines) + "\n"


def write_single(path, n, body):
    txt = preamble(n)
    txt += "(assert\n  " + body + "\n)\n(check-sat)\n"
    path.write_text(txt)


def write_union(path, n, bodies, header=None):
    """Assert the disjunction of `bodies`.  A single body is written as a bare
    assert -- SMT-LIB gives `or` arity >= 2, and not every solver is lenient."""
    txt = preamble(n)
    if header:
        txt += "".join(f"; {line}\n" for line in header.splitlines())
    if len(bodies) == 1:
        txt += "(assert\n  " + bodies[0] + "\n)\n(check-sat)\n"
    else:
        txt += "(assert\n  (or\n"
        txt += "\n".join("    " + b.replace("\n", "\n    ") for b in bodies)
        txt += "\n  )\n)\n(check-sat)\n"
    path.write_text(txt)


def parse_subset(spec, n):
    """Parse a body-selection spec into an ordered, de-duplicated list of
    (axis, sign) pairs.

    Comma-separated selectors.  Each is an axis part -- N, N..M (inclusive) or
    * -- optionally suffixed with + or - to take one sign; without a suffix
    both signs of each named axis are taken.  Ranges spell their separator
    ".." so that the "-" suffix stays unambiguous.
    """
    sel, seen = [], set()
    for tok in spec.split(","):
        tok = tok.strip()
        if not tok:
            continue
        if tok[-1] in "+-":
            signs, axpart = (tok[-1],), tok[:-1].strip()
        else:
            signs, axpart = ("+", "-"), tok
        if axpart in ("*", "all"):
            axes = range(n)
        elif ".." in axpart:
            lo, _, hi = axpart.partition("..")
            try:
                lo, hi = int(lo), int(hi)
            except ValueError:
                raise SystemExit(f"Bad subset selector {tok!r}: expected N..M.")
            if lo > hi:
                raise SystemExit(f"Bad subset selector {tok!r}: empty range.")
            axes = range(lo, hi + 1)
        else:
            try:
                axes = [int(axpart)]
            except ValueError:
                raise SystemExit(f"Bad subset selector {tok!r}.")
        for a in axes:
            if not 0 <= a < n:
                raise SystemExit(
                    f"Bad subset selector {tok!r}: axis {a} outside 0..{n-1}.")
            for s in signs:
                if (a, s) not in seen:
                    seen.add((a, s))
                    sel.append((a, s))
    if not sel:
        raise SystemExit(f"Subset spec {spec!r} selects no bodies.")
    return sel


def subset_stem(spec):
    """Filename stem for a spec given without an explicit =NAME."""
    out = (spec.replace("..", "to").replace("+", "p").replace("-", "m")
              .replace("*", "all").replace(",", "_"))
    out = "".join(ch for ch in out if ch.isalnum() or ch == "_")
    return f"subset_{out}" if out else "subset"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--beta", required=True,
                    help="Petal volume fraction, e.g. 0.9 or 9/10; must lie in (0,1).")
    ap.add_argument("--outdir", default="sea_urchin")
    ap.add_argument("--rotate", action="store_true",
                    help="Apply an exact rational orthogonal rotation.")
    ap.add_argument("--subset", action="append", metavar="SPEC[=NAME]",
                    help="Also write a union over just these bodies, e.g. "
                         "'0..4' (5 axes, both signs), '0+,1+,2+', '*+' (all "
                         "plus bodies), '3-=one_petal'. Repeatable.")
    args = ap.parse_args()

    n = args.n
    beta = frac(args.beta)
    if n < 2:
        raise SystemExit("Use n >= 2.")
    if not (0 < beta < 1):
        raise SystemExit("beta must be in (0,1).")

    L = (1 + beta) / (1 - beta)
    ext = L - 1
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    vars_ = [f"x{i}" for i in range(n)]
    Q = rational_rotation(n) if args.rotate else identity(n)
    forms = [linexpr(row, vars_) for row in Q]

    bodies = []
    names = []
    for i in range(n):
        for sign in ("+", "-"):
            bodies.append(bounds_for_body(n, i, sign, L, forms))
            names.append(f"K_{i:03d}_{'plus' if sign == '+' else 'minus'}")

    write_union(out / "union.smt2", n, bodies)

    for name, body in zip(names, bodies):
        write_single(out / f"{name}.smt2", n, body)

    core_vol = Fraction(2**n, 1)
    petal_vol = core_vol * beta / (1 - beta)
    body_vol = core_vol / (1 - beta)
    union_vol = core_vol * (1 + Fraction(2*n, 1) * beta / (1 - beta))
    half_long = Fraction(1, 1) / (1 - beta)
    r2 = Fraction(n - 1, 1) + half_long * half_long

    meta = f"""n = {n}
number_of_bodies = {2*n}
beta = {beta}
L = {L}
extension_length = L-1 = {ext}

exact_core_volume = {core_vol}
exact_petal_volume_each = {petal_vol}
exact_body_volume_each = {body_vol}
exact_union_volume = {union_vol}

exact_inradius_each_body = 1
circumradius_squared_each_body = {r2}
R_over_r = sqrt({r2})
longitudinal_side_length = L+1 = {L+1}
transverse_side_length = 2
body_axis_aspect_ratio = (L+1)/2 = {half_long}
petal_length_to_transverse_width = (L-1)/2 = {ext/2}

rotated = {args.rotate}
"""
    subsets = []
    for raw in (args.subset or []):
        spec, sep, name = raw.partition("=")
        spec = spec.strip()
        name = name.strip() if sep else subset_stem(spec)
        if not name:
            raise SystemExit(f"Empty output name in --subset {raw!r}.")
        sel = parse_subset(spec, n)
        picked = [names[2*a + (0 if s == "+" else 1)] for a, s in sel]
        idx = [2*a + (0 if s == "+" else 1) for a, s in sel]
        # Every body holds the core and one petal, and the petals are pairwise
        # disjoint, so m distinct bodies union to core + m petals exactly.
        vol = core_vol + len(sel) * petal_vol
        write_union(out / f"{name}.smt2", n, [bodies[i] for i in idx],
                    header=f"subset spec: {spec}\n"
                           f"{len(sel)} of {2*n} bodies: {', '.join(picked)}\n"
                           f"exact union volume: {vol}")
        subsets.append((name, spec, len(sel), vol))

    if subsets:
        meta += "\nsubsets\n"
        for name, spec, count, vol in subsets:
            meta += (f"  {name}.smt2: spec = {spec}, bodies = {count}, "
                     f"exact_union_volume = {vol}\n")

    (out / "metadata.txt").write_text(meta)
    print(meta)
    print(f"Wrote {out / 'union.smt2'} and {2*n} individual polytope files.")
    for name, spec, count, vol in subsets:
        print(f"Wrote {out / (name + '.smt2')}: {count} bodies ({spec}), "
              f"exact volume {vol}.")


if __name__ == "__main__":
    main()
