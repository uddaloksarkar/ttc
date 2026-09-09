#!/usr/bin/env python3
"""
har_mixing.py
=============

Empirical mixing-time analysis of a hit-and-run random walk sampler on
convex polytopes given in H-representation (Ax <= b).

Everything that touches the random walk itself -- direction sampling,
chord-length computation along a line, the Metropolis-free hit-and-run
update -- is implemented from scratch on top of numpy. No library ever
computes a polytope's volume, and no library ever draws directly from
the target distribution; the walk is the only sampling mechanism.

scipy is used for two auxiliary, non-volumetric primitives that make
the geometry controllable and checkable:
  * `linprog` to find a Chebyshev center / inscribed-ball radius (a
    single linear program -- this is just "find the deepest feasible
    point", not a volume computation), and
  * `minimize` (SLSQP) to find the minimal enclosing ball of an
    explicit, analytically-known vertex set (again, a distance/radius
    computation, not a volume computation).

CORE PARAMETERS
----------------
  n            : ambient dimension of the polytope
  num_steps    : number of hit-and-run steps per chain
  R_r_ratio    : target "fatness" index rho = R/r (circumradius over
                 inradius) of the synthetic test polytope. rho = 1 is
                 a ball; rho grows as the body becomes thinner /
                 more ill-conditioned.

Two synthetic polytope families are provided, both with rho set
*exactly* (up to numerical tolerance) via closed-form or calibrated
construction:

  box     : an axis-aligned hyperrectangle stretched along one axis.
            A box is itself a zonotope (Minkowski sum of orthogonal
            segments), so this doubles as the "zonotope" example.
            rho is achievable exactly in closed form; minimum
            achievable rho for a box in R^n is sqrt(n) (the cube).

  simplex : a right-angled simplex with vertices {0, kappa*e_1, e_2,
            ..., e_n}, calibrated by 1-D bisection on kappa so that
            R/r hits the target exactly. Minimum achievable rho is a
            slightly larger, n-dependent constant (a simplex is never
            as "round" as a box of the same dimension); the script
            reports and clips to it if the target is infeasible.

Usage
-----
  Single run:
    python har_mixing.py --n 10 --num_steps 500 --R_r_ratio 8 --kind box

  Step-count sweep (finds an empirical mixing step for one polytope):
    python har_mixing.py --mode step_sweep --n 10 --R_r_ratio 8 --kind box --plot

  Bottleneck isolation (varies n and R/r independently, fits log-log
  scaling exponents, and compares them to the commonly cited
  O*(n^2 (R/r)^2) hit-and-run reference bound):
    python har_mixing.py --mode bottleneck --kind box --plot

REAL BENCHMARK POLYTOPES
-------------------------
The same analysis runs on the polytopes ttc extracts from an SMT2 LRA
benchmark. ttc canonicalizes the formula into a union of polytopes and
`--dump-ine PREFIX` writes each one as PREFIX_cubeN.ine in the cdd /
Avis-Fukuda H-representation; this script shells out to it, parses the
files, and treats every cube as a test body:

    python har_mixing.py --mode benchmark --smt2 bench.smt2.xz --plot
    python har_mixing.py --mode benchmark --smt2 dim9.smt2.xz dim14.smt2.xz dim20.smt2.xz
    python har_mixing.py --mode benchmark --ine dumped_ines/      # skip re-running ttc

Two things necessarily change for a benchmark body:
  * there is no closed-form marginal, so convergence is measured against a
    long-run reference sample (two-sample KS) rather than exact truth --
    see section 4b for why that is weaker and how it is guarded, and
  * R has no cheap exact value, so it is bracketed by LPs and the reported
    R/r is an UPPER bound.
All polytopes of one benchmark share a dimension, so pass SEVERAL files of
different dim if you want the n exponent to be identifiable.
"""

import argparse
import glob
import itertools
import json
import lzma
import os
import shutil
import subprocess
import sys
import tempfile
import time

import numpy as np
from scipy.optimize import linprog, minimize, brentq
from scipy import stats


# ===========================================================================
# 1. Geometry utilities (LP / small optimization -- NOT volume computation)
# ===========================================================================

def chebyshev_center(A, b):
    """Largest inscribed ball of {x : Ax <= b} via a single linear program.

    maximize r  s.t.  A_i x + r*||A_i||_2 <= b_i for all i,  r >= 0.

    Returns (center, radius). This is the standard LP formulation of the
    Chebyshev center; it never touches volume.
    """
    n = A.shape[1]
    row_norms = np.linalg.norm(A, axis=1)
    c = np.zeros(n + 1)
    c[-1] = -1.0  # maximize r  <=>  minimize -r
    A_ub = np.hstack([A, row_norms.reshape(-1, 1)])
    res = linprog(c, A_ub=A_ub, b_ub=b,
                   bounds=[(None, None)] * n + [(0, None)], method='highs')
    if not res.success:
        raise RuntimeError(f"Chebyshev-center LP failed: {res.message}")
    return res.x[:-1], float(res.x[-1])


def smallest_enclosing_ball(vertices):
    """Minimal enclosing ball of a finite, explicitly known point set.

    Formulated as the convex feasibility program
        minimize   R2
        subject to R2 >= ||v_i - c||^2   for every vertex v_i
    and solved with SLSQP from a centroid-based feasible start. Convex
    program with a unique global optimum, so any local solver that
    converges finds THE minimum -- this is a radius computation over a
    known vertex list, not a volume estimate.
    """
    V = np.asarray(vertices, dtype=float)
    c0 = V.mean(axis=0)
    R0_sq = np.max(np.sum((V - c0) ** 2, axis=1))
    x0 = np.concatenate([c0, [R0_sq]])

    def objective(x):
        return x[-1]

    def cons_fun(x):
        c, R2 = x[:-1], x[-1]
        return R2 - np.sum((V - c) ** 2, axis=1)

    res = minimize(objective, x0, constraints=[{'type': 'ineq', 'fun': cons_fun}],
                    method='SLSQP', options={'maxiter': 2000, 'ftol': 1e-14})
    c = res.x[:-1]
    R = float(np.sqrt(max(res.x[-1], 0.0)))
    return c, R


# ===========================================================================
# 2. Synthetic polytopes with an exactly-controlled R/r ("fatness") ratio
# ===========================================================================

class Polytope:
    def __init__(self, A, b, center, r, R, kind, params):
        self.A = np.asarray(A, dtype=float)
        self.b = np.asarray(b, dtype=float)
        self.center = np.asarray(center, dtype=float)
        self.r = r
        self.R = R
        self.kind = kind
        self.params = params

    @property
    def n(self):
        return self.A.shape[1]

    @property
    def ratio(self):
        return self.R / self.r

    def __repr__(self):
        return (f"Polytope(kind={self.kind}, n={self.n}, r={self.r:.5g}, "
                f"R={self.R:.5g}, R/r={self.ratio:.5g}, params={self.params})")


def make_box(n, R_r_ratio):
    """Axis-aligned box (zonotope): side kappa along axis 0, side 1 along
    every other axis, centered at the origin.

        r = 1/2
        R = 0.5 * sqrt(kappa^2 + (n-1))
        R/r = sqrt(kappa^2 + (n-1))

    Solving for kappa given a target ratio rho is closed form:
        kappa = sqrt(rho^2 - (n-1))
    which requires rho >= sqrt(n) (the cube, kappa=1, is the roundest
    achievable box in R^n).
    """
    rho_min = np.sqrt(n)
    rho = R_r_ratio
    if rho < rho_min - 1e-9:
        print(f"[warn] R_r_ratio={rho:.4f} is below the minimum achievable "
              f"for a box in R^{n} (sqrt(n)={rho_min:.4f}); clipping to sqrt(n).")
        rho = rho_min
    kappa = float(np.sqrt(max(rho ** 2 - (n - 1), 1.0)))
    L = np.ones(n)
    L[0] = kappa
    A = np.vstack([np.eye(n), -np.eye(n)])
    b = np.concatenate([L / 2, L / 2])
    center = np.zeros(n)
    r = 0.5
    R = 0.5 * np.sqrt(kappa ** 2 + (n - 1))
    return Polytope(A, b, center, r, R, 'box', {'kappa': kappa, 'L': L.tolist()})


def _build_simplex_Ab(n, kappa):
    A = -np.eye(n)  # x_i >= 0
    coef = np.ones(n)
    coef[0] = 1.0 / kappa
    A = np.vstack([A, coef.reshape(1, -1)])  # sum x_i/kappa_i <= 1
    b = np.concatenate([np.zeros(n), [1.0]])
    return A, b


def _simplex_vertices(n, kappa):
    V = np.zeros((n + 1, n))
    V[1, 0] = kappa
    for i in range(1, n):
        V[i + 1, i] = 1.0
    return V


def _simplex_ratio(n, kappa):
    A, b = _build_simplex_Ab(n, kappa)
    _, r = chebyshev_center(A, b)
    V = _simplex_vertices(n, kappa)
    _, R = smallest_enclosing_ball(V)
    return R / r, A, b, r, R


def make_simplex(n, R_r_ratio, xtol=1e-7, max_kappa=1e8):
    """Right-angled simplex with vertices {0, kappa*e_1, e_2, ..., e_n}.

    Unlike the box, R and r have no simple closed form once kappa != 1,
    so kappa is calibrated by 1-D bisection (scipy.optimize.brentq) on
    the ratio R(kappa)/r(kappa), which is monotonically increasing in
    kappa for kappa >= 1. r comes from the Chebyshev LP on the exact
    H-representation; R from the minimal enclosing ball of the exactly
    known vertex list.
    """
    rho_min, A0, b0, r0, R0 = _simplex_ratio(n, 1.0)
    if R_r_ratio < rho_min - 1e-6:
        print(f"[warn] R_r_ratio={R_r_ratio:.4f} is below the minimum achievable "
              f"for a right simplex in R^{n} ({rho_min:.4f}); clipping.")
        center, _ = chebyshev_center(A0, b0)
        return Polytope(A0, b0, center, r0, R0, 'simplex', {'kappa': 1.0})

    def f(kappa):
        return _simplex_ratio(n, kappa)[0] - R_r_ratio

    kappa_hi = 2.0
    while f(kappa_hi) < 0:
        kappa_hi *= 2
        if kappa_hi > max_kappa:
            raise RuntimeError("Could not bracket target R_r_ratio for simplex; "
                                "ratio requested is implausibly large.")
    kappa_star = brentq(f, 1.0, kappa_hi, xtol=xtol)
    _, A, b, r, R = _simplex_ratio(n, kappa_star)
    center, _ = chebyshev_center(A, b)
    return Polytope(A, b, center, r, R, 'simplex', {'kappa': kappa_star})


def make_polytope(n, R_r_ratio, kind):
    if kind == 'box':
        return make_box(n, R_r_ratio)
    elif kind == 'simplex':
        return make_simplex(n, R_r_ratio)
    else:
        raise ValueError(f"unknown polytope kind: {kind}")


# ===========================================================================
# 2b. Polytopes from the SMT benchmarks, via ttc --dump-ine
# ===========================================================================
#
# ttc canonicalizes an LRA formula into a union of polytopes and can write each
# one out with
#
#     ./ttc -v 2 --dump-ine <PREFIX> file.smt2
#
# producing <PREFIX>_cube1.ine, _cube2.ine, ... in the cdd / Avis-Fukuda
# H-representation (see writePolytopeIne in src/volume/lra_volume.cpp).  The
# body of that file is
#
#     begin
#      <m> <n+1> real
#      b_i  -a_i0  -a_i1  ...        (m rows, each encoding a_i . x <= b_i)
#     end
#
# so column 0 is b and the remaining columns are -A.
#
# Unlike the synthetic families, a benchmark polytope has:
#   * no closed-form marginal, so convergence is measured against a long-run
#     reference sample instead of an analytic CDF (see section 4b), and
#   * no cheap exact circumradius, so R is bracketed by LPs (below).

def read_ine(path):
    """Parse a cdd H-representation file into (A, b) with A x <= b."""
    with open(path) as fh:
        tokens = fh.read().split("\n")

    rows, i = None, 0
    while i < len(tokens):
        line = tokens[i].strip()
        if line.lower().startswith("linearity"):
            raise ValueError(f"{path}: 'linearity' (equality) rows are not "
                             "supported -- the body is lower-dimensional and a "
                             "hit-and-run chord would be degenerate")
        if line == "begin":
            rows = i
            break
        i += 1
    if rows is None:
        raise ValueError(f"{path}: no 'begin' section")

    header = tokens[rows + 1].split()
    m, ncols = int(header[0]), int(header[1])
    data = []
    for line in tokens[rows + 2: rows + 2 + m]:
        vals = [float(v) for v in line.split()]
        if len(vals) != ncols:
            raise ValueError(f"{path}: expected {ncols} columns, got {len(vals)}")
        data.append(vals)
    if len(data) != m:
        raise ValueError(f"{path}: expected {m} rows, got {len(data)}")

    M = np.asarray(data, dtype=float)
    b = M[:, 0].copy()
    A = -M[:, 1:].copy()
    return A, b


def centered_chebyshev(A, b, lo, hi, slack=1e-9):
    """Chebyshev center with a canonical tie-break, plus its radius.

    The largest inscribed ball is often not unique -- on a box it slides freely
    along the long axis, and HiGHS then returns an arbitrary point of that
    optimal face.  R measured about such a point is far too large: a box with
    true R/r = 6 came out at 10.6 before this tie-break existed.

    So: solve for the radius r* first, then among all centers that still admit
    a ball of radius r* pick the one closest (in max-coordinate distance) to
    the bounding-box center, by a second LP over (c, t):

        minimize t
        s.t.  a_i.c + r* ||a_i|| <= b_i          (keeps the inscribed radius)
              |c_j - m_j| <= t   for all j,       m = (lo + hi)/2

    On a simplex, where the optimum is already unique, this changes nothing.
    """
    n = A.shape[1]
    c0, r_star = chebyshev_center(A, b)
    if r_star <= 0:
        return c0, r_star
    m = 0.5 * (lo + hi)
    row_norms = np.linalg.norm(A, axis=1)

    # variables: [c (n), t (1)]
    A_ub = [np.hstack([A, np.zeros((A.shape[0], 1))])]
    b_ub = [b - r_star * (1.0 - slack) * row_norms]
    eye = np.eye(n)
    ones = np.ones((n, 1))
    A_ub += [np.hstack([eye, -ones]), np.hstack([-eye, -ones])]
    b_ub += [m, -m]

    obj = np.zeros(n + 1)
    obj[-1] = 1.0
    res = linprog(obj, A_ub=np.vstack(A_ub), b_ub=np.concatenate(b_ub),
                  bounds=[(None, None)] * n + [(0, None)], method='highs')
    if not res.success:
        return c0, r_star
    return res.x[:n], r_star


def inscribed_radius_at(A, b, c):
    """Radius of the largest ball centered exactly at c inside {Ax <= b}."""
    return float(np.min((b - A @ c) / np.linalg.norm(A, axis=1)))


def bounding_box(A, b):
    """Per-coordinate min/max over {Ax <= b} by 2n linear programs.

    Returns (lo, hi, vertices) where `vertices` are the LP optima themselves --
    genuine points of the polytope, used below for a lower bound on R and as
    over-dispersed chain starts.  Returns None if the body is unbounded or
    infeasible in any coordinate.
    """
    n = A.shape[1]
    lo, hi, verts = np.empty(n), np.empty(n), []
    for j in range(n):
        c = np.zeros(n)
        c[j] = 1.0
        for sense, sink in ((1.0, lo), (-1.0, hi)):
            res = linprog(sense * c, A_ub=A, b_ub=b,
                          bounds=[(None, None)] * n, method='highs')
            if not res.success:
                return None
            sink[j] = res.x[j]
            verts.append(res.x.copy())
    return lo, hi, np.asarray(verts)


def polytope_from_ine(path, min_inradius=1e-9):
    """Build a Polytope from a .ine file, or return (None, reason).

    r is exact (Chebyshev LP).  R follows the SAME convention as the synthetic
    families -- the radius of the minimal enclosing ball, about its own center,
    not concentric with r -- so the exponents this mode fits are comparable
    with the ones `bottleneck` mode fits.  Exact vertex enumeration is
    exponential, so R is bracketed instead:

        R_lo = MEB radius of the 2n LP optima -- a subset of K, so its
               enclosing ball is no larger than K's
        R_hi = 0.5 * ||hi - lo||, the ball around the bounding box, which
               certainly contains K

    We report R_hi, so `poly.ratio` is an UPPER bound on the true fatness; the
    gap to ratio_lower says how loose the bracket is on that body.
    """
    A, b = read_ine(path)
    box = bounding_box(A, b)
    if box is None:
        return None, "unbounded or infeasible"
    lo, hi, verts = box

    # the center is only used to start chains and to disperse them, so take the
    # canonically-placed Chebyshev optimum (see centered_chebyshev)
    center, r = centered_chebyshev(A, b, lo, hi)
    if r < min_inradius:
        return None, f"degenerate (inradius {r:.3g} < {min_inradius:g})"

    R_hi = 0.5 * float(np.linalg.norm(hi - lo))
    _, R_lo = smallest_enclosing_ball(verts)
    R_lo = float(min(R_lo, R_hi))
    poly = Polytope(A, b, center, r, R_hi, 'ine',
                    {'path': path, 'facets': int(A.shape[0]),
                     'R_lower': R_lo, 'R_upper': R_hi,
                     'ratio_lower': R_lo / r,
                     'lo': lo.tolist(), 'hi': hi.tolist()})
    poly.lp_vertices = verts
    return poly, None


def decompress_if_xz(path, workdir):
    """Return a plain .smt2 path, decompressing <name>.smt2.xz if needed."""
    if not path.endswith(".xz"):
        return path
    out = os.path.join(workdir, os.path.basename(path)[:-3])
    with lzma.open(path, "rb") as src, open(out, "wb") as dst:
        shutil.copyfileobj(src, dst)
    return out


def dump_ine_from_smt2(smt2_path, ttc_bin, workdir, extra_args=(), verbose=True,
                       return_output=False):
    """Run `ttc --dump-ine` on one benchmark and return the .ine paths, sorted
    by the cube index ttc assigned (cube1, cube2, ... -- ttc numbers from 1).

    With return_output=True the (paths, stdout) pair comes back instead; ttc
    continues into the normal volume computation after dumping, so its stdout
    carries the per-polytope volume table and the sample-count totals.
    """
    ttc = shutil.which(ttc_bin) or os.path.abspath(ttc_bin)
    if not os.path.exists(ttc):
        sys.exit(f"ttc binary not found: {ttc_bin} (pass --ttc_bin)")
    smt2 = decompress_if_xz(smt2_path, workdir)
    prefix = os.path.join(workdir, "p")
    cmd = [ttc, "-v", "2", "--dump-ine", prefix, *extra_args, smt2]
    if verbose:
        print(f"$ {' '.join(cmd)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    files = sorted(glob.glob(prefix + "_cube*.ine"),
                   key=lambda p: int(p.rsplit("cube", 1)[1].split(".")[0]))
    if not files:
        sys.exit("ttc produced no .ine files; last output:\n"
                 + (proc.stdout or proc.stderr)[-2000:])
    if verbose:
        print(f"  ttc exit {proc.returncode}, {len(files)} polytopes dumped")
    if return_output:
        return files, (proc.stdout or "")
    return files


def load_benchmark_polytopes(args, workdir):
    """Resolve --smt2 / --ine into a list of Polytope objects, reporting every
    polytope that had to be skipped rather than dropping it silently."""
    sources = {}
    if args.smt2:
        paths = []
        for k, smt2 in enumerate(args.smt2):
            sub = os.path.join(workdir, f"b{k}")
            os.makedirs(sub, exist_ok=True)
            got = dump_ine_from_smt2(smt2, args.ttc_bin, sub,
                                     extra_args=args.ttc_args.split() if args.ttc_args else ())
            for p in got:
                sources[p] = os.path.basename(smt2)
            paths.extend(got)
    elif args.ine:
        paths = []
        for pat in args.ine:
            paths.extend(sorted(glob.glob(os.path.join(pat, "*.ine"))
                                if os.path.isdir(pat) else glob.glob(pat)))
        if not paths:
            sys.exit(f"no .ine files matched {args.ine}")
    else:
        return []

    polys, skipped = [], []
    for p in paths:
        poly, reason = polytope_from_ine(p)
        if poly is None:
            skipped.append((os.path.basename(p), reason))
            continue
        poly.params['source'] = sources.get(
            p, os.path.basename(os.path.dirname(os.path.abspath(p))))
        polys.append(poly)
    if skipped:
        print(f"[warn] skipped {len(skipped)}/{len(paths)} polytopes:")
        for name, reason in skipped[:10]:
            print(f"         {name}: {reason}")
        if len(skipped) > 10:
            print(f"         ... and {len(skipped) - 10} more")
    if args.max_polytopes and len(polys) > args.max_polytopes:
        print(f"[note] using the first {args.max_polytopes} of {len(polys)} polytopes")
        polys = polys[:args.max_polytopes]
    return polys


# ===========================================================================
# 3. Hit-and-run random walk -- implemented from scratch
# ===========================================================================

def sample_direction(n, rng):
    """Uniform direction on the unit sphere S^{n-1} via a normalized
    isotropic Gaussian vector."""
    d = rng.normal(size=n)
    norm = np.linalg.norm(d)
    while norm < 1e-12:
        d = rng.normal(size=n)
        norm = np.linalg.norm(d)
    return d / norm


def chord_interval(x, d, A, b, eps=1e-12):
    """t_min, t_max such that x + t*d stays inside {Ax <= b}.

    For each row i: t*(A_i . d) <= b_i - A_i.x =: slack_i (>= 0 since x
    is feasible). Rows with A_i.d > 0 give upper bounds on t, rows with
    A_i.d < 0 give lower bounds (sign flip on division), rows with
    A_i.d ~= 0 are parallel to the chord and impose no bound.
    """
    Ad = A @ d
    slack = b - A @ x
    pos = Ad > eps
    neg = Ad < -eps
    t_max = np.min(slack[pos] / Ad[pos]) if np.any(pos) else np.inf
    t_min = np.max(slack[neg] / Ad[neg]) if np.any(neg) else -np.inf
    return t_min, t_max


def hit_and_run_step(x, A, b, rng, max_retries=25):
    """One hit-and-run update: pick a uniformly random direction through
    x, find the chord it cuts through the polytope, move to a point
    drawn uniformly along that chord."""
    n = len(x)
    for _ in range(max_retries):
        d = sample_direction(n, rng)
        t_min, t_max = chord_interval(x, d, A, b)
        if np.isfinite(t_min) and np.isfinite(t_max) and t_max > t_min:
            t = rng.uniform(t_min, t_max)
            return x + t * d
    raise RuntimeError("hit_and_run_step: no valid chord found after "
                        f"{max_retries} retries (numerical degeneracy?).")


def run_chain(x0, A, b, num_steps, rng):
    """Run one hit-and-run chain for num_steps steps, return the final point."""
    x = np.array(x0, dtype=float)
    for _ in range(num_steps):
        x = hit_and_run_step(x, A, b, rng)
    return x


def run_parallel_chains(poly, num_steps, num_chains, rng, x0=None, starts=None):
    """Run num_chains INDEPENDENT hit-and-run chains, each of length
    num_steps. Returns an (num_chains, n) array of final states -- this is
    the batch of "samples after num_steps" whose distribution we test
    for convergence to uniform.

    By default every chain starts at the same point x0 (the polytope center).
    Pass `starts` -- an (num_chains, n) array, e.g. from dispersed_starts -- to
    launch them from over-dispersed points instead, which is what makes the
    reference sample and any between-chain diagnostic meaningful.
    """
    if starts is None:
        if x0 is None:
            x0 = poly.center
        starts = np.tile(np.asarray(x0, dtype=float), (num_chains, 1))
    samples = np.empty((num_chains, poly.n))
    for m in range(num_chains):
        chain_rng = np.random.default_rng(rng.integers(0, 2 ** 63 - 1))
        samples[m] = run_chain(starts[m % len(starts)], poly.A, poly.b,
                               num_steps, chain_rng)
    return samples


def dispersed_starts(poly, num_chains, rng, shrink=1e-6):
    """Over-dispersed feasible starting points: optima of random linear
    objectives (i.e. random vertices of the polytope), pulled a hair toward the
    Chebyshev center so they stay strictly interior.

    This is an LP per point -- a "find the extreme point in this direction"
    computation, not a volume or a draw from the target.
    """
    n = poly.n
    pool = list(getattr(poly, 'lp_vertices', []))
    need = max(0, min(num_chains, 64) - len(pool))
    for _ in range(need):
        c = sample_direction(n, rng)
        res = linprog(-c, A_ub=poly.A, b_ub=poly.b,
                      bounds=[(None, None)] * n, method='highs')
        if res.success:
            pool.append(res.x.copy())
    if not pool:
        return np.tile(poly.center, (num_chains, 1))
    P = np.asarray(pool)
    idx = rng.integers(0, len(P), size=num_chains)
    return poly.center + (1.0 - shrink) * (P[idx] - poly.center)


# ===========================================================================
# 4b. Convergence without a closed form: long-run reference sample
# ===========================================================================
#
# The synthetic families have analytic marginals, so convergence is measured
# against exact truth.  A benchmark polytope has none, so the target is instead
# approximated by a long-run reference sample and the comparison becomes a
# TWO-sample test.  That is strictly weaker: if the walk mixes so slowly that
# even the reference is not converged, the test compares two equally-wrong
# distributions and reports convergence too early.  Two things guard against
# it -- the reference runs from over-dispersed near-vertex starts (so it does
# not inherit the tested chains' central start), and reference_selfcheck splits
# it in half and KS-tests the halves against each other, which fails loudly
# when ref_steps is too small.

def projection_direction(poly, mode, coord, rng):
    """Unit vector whose projection of the samples is the tested statistic.

    'coord'   : the requested coordinate axis (default for the synthetic bodies)
    'longest' : the longest axis of the LP bounding box -- the direction a
                stretched benchmark polytope mixes slowest along
    'random'  : one fixed random direction
    """
    n = poly.n
    if mode == 'random':
        return sample_direction(n, rng)
    if mode == 'longest' and 'lo' in poly.params:
        extent = np.asarray(poly.params['hi']) - np.asarray(poly.params['lo'])
        w = np.zeros(n)
        w[int(np.argmax(extent))] = 1.0
        return w
    w = np.zeros(n)
    w[coord] = 1.0
    return w


def reference_sample(poly, rng, ref_chains, ref_steps, verbose=True):
    """Long-run reference batch from over-dispersed starts."""
    if verbose:
        print(f"  building reference sample: {ref_chains} chains x {ref_steps} steps")
    starts = dispersed_starts(poly, ref_chains, rng)
    return run_parallel_chains(poly, ref_steps, ref_chains, rng, starts=starts)


def reference_selfcheck(ref, w, alpha=0.01):
    """Split the reference in half and two-sample KS the halves.

    A small p-value means the reference itself has not settled, so any
    convergence verdict measured against it is unreliable.
    """
    proj = ref @ w
    half = len(proj) // 2
    if half < 30:
        return {'ok': None, 'pvalue': None}
    p = float(stats.ks_2samp(proj[:half], proj[half:]).pvalue)
    return {'ok': bool(p > alpha), 'pvalue': p}


def evaluate_convergence_empirical(samples, ref, w):
    """Two-sample comparison of a projected marginal against the reference."""
    x, y = samples @ w, ref @ w
    ks = stats.ks_2samp(x, y)
    ref_var = float(np.var(y, ddof=1))
    emp_var = float(np.var(x, ddof=1))
    return {
        'n_samples': int(len(x)),
        'empirical_variance': emp_var,
        'theoretical_variance': ref_var,   # reference, not exact truth
        'variance_ratio': emp_var / ref_var if ref_var > 0 else float('nan'),
        'ks_statistic': float(ks.statistic),
        'ks_pvalue': float(ks.pvalue),
        'reference': True,
    }


# ===========================================================================
# 4. Convergence diagnostics
# ===========================================================================

def theoretical_marginal(poly, coord):
    """Exact closed-form CDF and variance of coordinate `coord` under the
    TRUE uniform distribution on the polytope. Used as ground truth for
    the convergence test -- this is analytic, not simulated.
    """
    n = poly.n
    if poly.kind == 'box':
        L = poly.params['L'][coord]
        lo, hi = -L / 2, L / 2
        cdf = lambda t, lo=lo, hi=hi: np.clip((t - lo) / (hi - lo), 0, 1)
        var = L ** 2 / 12
        return cdf, var
    elif poly.kind == 'simplex':
        # y_coord = kappa_coord * Beta(1, n),  kappa_0 = kappa, kappa_{i>0} = 1
        kappa = poly.params['kappa'] if coord == 0 else 1.0
        cdf = lambda t, kappa=kappa, n=n: 1 - (1 - np.clip(t / kappa, 0, 1)) ** n
        var = kappa ** 2 * n / ((n + 1) ** 2 * (n + 2))
        return cdf, var
    else:
        raise ValueError(poly.kind)


def evaluate_convergence(poly, samples, coord=0, ref=None, w=None):
    """Empirical variance + Kolmogorov-Smirnov test of one marginal.

    Against the exact analytic CDF for the synthetic families; against the
    long-run reference sample (two-sample KS on the projection `w`) when one is
    supplied, which is the only option for a benchmark polytope.
    """
    if ref is not None:
        if w is None:
            w = np.zeros(poly.n)
            w[coord] = 1.0
        return evaluate_convergence_empirical(samples, ref, w)
    cdf, true_var = theoretical_marginal(poly, coord)
    x = samples[:, coord]
    emp_var = float(np.var(x, ddof=1))
    ks_stat, ks_p = stats.kstest(x, cdf)
    return {
        'n_samples': int(len(x)),
        'empirical_variance': emp_var,
        'theoretical_variance': true_var,
        'variance_ratio': emp_var / true_var,
        'ks_statistic': float(ks_stat),
        'ks_pvalue': float(ks_p),
    }


# ===========================================================================
# 5. Mixing-time sweeps and bottleneck isolation
# ===========================================================================

def theoretical_bound(n, R_r_ratio, const=1.0):
    """A commonly cited O*(n^2 (R/r)^2) reference for hit-and-run mixing
    time (Lovasz-type bound, up to log factors and an unspecified
    constant). Provided ONLY as a qualitative reference curve to compare
    empirical scaling against -- not a rigorous prediction, and `const`
    is not fit to anything by default.
    """
    return const * (n ** 2) * (R_r_ratio ** 2)


def step_sweep(poly, step_grid, num_chains, rng, coord=0, var_tol=0.05, ks_alpha=0.05,
                verbose=True, ref=None, w=None):
    """Run the sampler at each step count in step_grid, evaluate
    convergence, and report the first step count at which BOTH the
    variance is within var_tol of the true value AND the KS test does
    not reject uniformity at level ks_alpha. That step count is the
    script's empirical mixing-time estimate.
    """
    rows = []
    mixing_step = None
    for steps in step_grid:
        samples = run_parallel_chains(poly, steps, num_chains, rng)
        row = evaluate_convergence(poly, samples, coord, ref=ref, w=w)
        row['num_steps'] = steps
        converged = (abs(row['variance_ratio'] - 1) < var_tol) and (row['ks_pvalue'] > ks_alpha)
        row['converged'] = converged
        if converged and mixing_step is None:
            mixing_step = steps
        rows.append(row)
        if verbose:
            print(f"  steps={steps:6d}  var_ratio={row['variance_ratio']:.3f}  "
                  f"KS_p={row['ks_pvalue']:.3g}  converged={converged}")
    return rows, mixing_step


def find_mixing_step_adaptive(poly, num_chains, rng, coord=0, var_tol=0.05, ks_alpha=0.05,
                               start_steps=10, max_steps=20000, verbose=True,
                               ref=None, w=None):
    """Geometrically double the step count from start_steps until BOTH the
    variance and KS convergence criteria are met, or max_steps is
    exceeded. Used by bottleneck_analysis instead of a fixed step_grid so
    that slow-mixing configurations are never silently dropped from the
    log-log fit (a fixed grid would right-censor them and bias the fitted
    exponent -- e.g. it can even flip its sign).

    Returns (mixing_step or None, last_row). None means the walk still
    hadn't converged by max_steps.
    """
    steps = start_steps
    last_row = None
    while steps <= max_steps:
        samples = run_parallel_chains(poly, steps, num_chains, rng)
        row = evaluate_convergence(poly, samples, coord, ref=ref, w=w)
        row['num_steps'] = steps
        converged = (abs(row['variance_ratio'] - 1) < var_tol) and (row['ks_pvalue'] > ks_alpha)
        row['converged'] = converged
        last_row = row
        if verbose:
            print(f"    steps={steps:6d}  var_ratio={row['variance_ratio']:.3f}  "
                  f"KS_p={row['ks_pvalue']:.3g}  converged={converged}")
        if converged:
            return steps, last_row
        steps *= 2
    print(f"  [warn] did not converge by max_steps={max_steps}; "
          f"excluded from the scaling fit (not silently ignored -- flagged here).")
    return None, last_row


def _fit_loglog_exponent(xs, ys):
    """Fit log(y) = a*log(x) + c by least squares; return (a, c, r_squared)."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    lx, ly = np.log(xs), np.log(ys)
    a, c = np.polyfit(lx, ly, 1)
    pred = a * lx + c
    ss_res = np.sum((ly - pred) ** 2)
    ss_tot = np.sum((ly - ly.mean()) ** 2)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')
    return float(a), float(c), float(r2)


def bottleneck_analysis(kind, n_grid, ratio_grid, base_n, base_ratio,
                         num_chains, seed, var_tol=0.05, ks_alpha=0.05,
                         start_steps=10, max_steps=200000):
    """Vary n (at fixed R/r = base_ratio) and vary R/r (at fixed n =
    base_n) independently, estimate the empirical mixing step for each
    configuration via an ADAPTIVE geometric search (see
    find_mixing_step_adaptive), and fit a log-log scaling exponent for
    each axis. A fixed step grid would right-censor slow-mixing
    configurations and bias the fitted exponent -- adaptivity avoids
    that. Comparing the fitted exponents to the O*(n^2 (R/r)^2)
    reference's exponents (2 and 2) is the "which factor is the
    bottleneck" answer: whichever fitted exponent is larger is scaling
    mixing time faster in this empirical range.
    """
    rng = np.random.default_rng(seed)
    results = {'vs_n': [], 'vs_ratio': []}

    print(f"\n=== Sweep 1/2: mixing step vs n  (R/r fixed at {base_ratio:.3f}) ===")
    for n in n_grid:
        poly = make_polytope(n, base_ratio, kind)
        print(f"n={n}  (built {poly.kind}, actual R/r={poly.ratio:.3f}, kappa={poly.params.get('kappa'):.3f})")
        mstep, _ = find_mixing_step_adaptive(poly, num_chains, rng, var_tol=var_tol,
                                              ks_alpha=ks_alpha, start_steps=start_steps,
                                              max_steps=max_steps)
        print(f"  -> empirical mixing step: {mstep}")
        results['vs_n'].append({'n': n, 'R_r_ratio': poly.ratio, 'mixing_step': mstep})

    print(f"\n=== Sweep 2/2: mixing step vs R/r  (n fixed at {base_n}) ===")
    for rho in ratio_grid:
        poly = make_polytope(base_n, rho, kind)
        print(f"R/r target={rho:.3f}  (built {poly.kind}, actual R/r={poly.ratio:.3f})")
        mstep, _ = find_mixing_step_adaptive(poly, num_chains, rng, var_tol=var_tol,
                                              ks_alpha=ks_alpha, start_steps=start_steps,
                                              max_steps=max_steps)
        print(f"  -> empirical mixing step: {mstep}")
        results['vs_ratio'].append({'n': base_n, 'R_r_ratio': poly.ratio, 'mixing_step': mstep})

    summary = {}
    valid_n = [(r['n'], r['mixing_step']) for r in results['vs_n'] if r['mixing_step'] is not None]
    valid_r = [(r['R_r_ratio'], r['mixing_step']) for r in results['vs_ratio'] if r['mixing_step'] is not None]

    if len(valid_n) >= 2:
        ns, ms = zip(*valid_n)
        a_n, c_n, r2_n = _fit_loglog_exponent(ns, ms)
        summary['n_exponent'] = a_n
        summary['n_fit_r2'] = r2_n
    else:
        summary['n_exponent'] = None
        summary['n_fit_r2'] = None

    if len(valid_r) >= 2:
        rs, ms = zip(*valid_r)
        a_r, c_r, r2_r = _fit_loglog_exponent(rs, ms)
        summary['ratio_exponent'] = a_r
        summary['ratio_fit_r2'] = r2_r
    else:
        summary['ratio_exponent'] = None
        summary['ratio_fit_r2'] = None

    print("\n=== Bottleneck summary ===")
    print("Reference O*(n^2 (R/r)^2) bound implies exponent 2 on EACH axis.")
    if summary['n_exponent'] is not None:
        print(f"  Empirical fit:  mixing_step ~ n^{summary['n_exponent']:.2f}   "
              f"(R^2={summary['n_fit_r2']:.3f})")
    else:
        print("  Empirical fit vs n: not enough converged points to fit.")
    if summary['ratio_exponent'] is not None:
        print(f"  Empirical fit:  mixing_step ~ (R/r)^{summary['ratio_exponent']:.2f}   "
              f"(R^2={summary['ratio_fit_r2']:.3f})")
    else:
        print("  Empirical fit vs R/r: not enough converged points to fit.")

    if summary['n_exponent'] is not None and summary['ratio_exponent'] is not None:
        if summary['ratio_exponent'] > summary['n_exponent']:
            print("  -> In this range, R/r (fatness) scales mixing time faster than n: "
                  "the FATNESS INDEX looks like the primary bottleneck.")
        elif summary['n_exponent'] > summary['ratio_exponent']:
            print("  -> In this range, n scales mixing time faster than R/r: "
                  "DIMENSION looks like the primary bottleneck.")
        else:
            print("  -> The two exponents are comparable in this range.")

    return results, summary


def benchmark_analysis(polys, num_chains, seed, ref_chains, ref_steps,
                        projection='longest', coord=0, var_tol=0.05, ks_alpha=0.05,
                        start_steps=10, max_steps=20000, verbose=True):
    """Empirical mixing step for every polytope of an SMT benchmark, then a fit
    of how it scales with n and with R/r.

    The synthetic `bottleneck` mode varies one factor at a time; here the
    geometry is whatever the benchmark produced, so n and R/r co-vary and the
    marginal log-log slopes confound each other.  The primary result is
    therefore a JOINT regression

        log(mixing_step) = a*log(n) + b*log(R/r) + c

    whose partial exponents are directly comparable to the 2 and 2 of the
    O*(n^2 (R/r)^2) reference.  The marginal fits are printed alongside, plus
    the correlation between log n and log R/r -- when that is close to +-1 the
    two exponents are not separately identifiable and the printout says so.

    Note that every R/r here is an UPPER bound (see polytope_from_ine), so the
    fitted R/r exponent inherits that looseness.
    """
    rng = np.random.default_rng(seed)
    rows = []
    for i, poly in enumerate(polys):
        src = poly.params.get('source', '?')
        print(f"\n[{i + 1}/{len(polys)}] {src} :: {os.path.basename(poly.params['path'])}  "
              f"n={poly.n}, facets={poly.params['facets']}, r={poly.r:.4g}, "
              f"R/r in [{poly.params['ratio_lower']:.3g}, {poly.ratio:.3g}]")
        w = projection_direction(poly, projection, coord, rng)
        ref = reference_sample(poly, rng, ref_chains, ref_steps, verbose=verbose)
        check = reference_selfcheck(ref, w)
        if check['ok'] is False:
            print(f"  [warn] reference halves disagree (KS p={check['pvalue']:.3g}): "
                  f"ref_steps={ref_steps} looks too small for this body; "
                  f"the mixing step below is a LOWER bound")
        mstep, last = find_mixing_step_adaptive(
            poly, num_chains, rng, coord=coord, var_tol=var_tol, ks_alpha=ks_alpha,
            start_steps=start_steps, max_steps=max_steps, verbose=verbose,
            ref=ref, w=w)
        print(f"  -> empirical mixing step: {mstep}")
        rows.append({
            'index': i,
            'source': src,
            'ine': os.path.basename(poly.params['path']),
            'n': poly.n,
            'facets': poly.params['facets'],
            'r': poly.r,
            'R_upper': poly.params['R_upper'],
            'R_lower': poly.params['R_lower'],
            'ratio_upper': poly.ratio,
            'ratio_lower': poly.params['ratio_lower'],
            'mixing_step': mstep,
            'reference_ok': check['ok'],
            'reference_pvalue': check['pvalue'],
            'last_variance_ratio': None if last is None else last['variance_ratio'],
            'last_ks_pvalue': None if last is None else last['ks_pvalue'],
        })

    valid = [r for r in rows if r['mixing_step'] is not None]
    print(f"\n=== Benchmark summary: {len(valid)}/{len(rows)} polytopes converged "
          f"by max_steps={max_steps} ===")
    head = f"{'#':>3} {'n':>4} {'facets':>7} {'R/r(up)':>9} {'mix':>7}  source"
    print(head)
    print("-" * (len(head) + 12))
    for r in rows:
        print(f"{r['index']:>3} {r['n']:>4} {r['facets']:>7} {r['ratio_upper']:>9.3g} "
              f"{str(r['mixing_step']):>7}  {r['source']}")

    summary = {'n_polytopes': len(rows), 'n_converged': len(valid)}
    if len(valid) < 3:
        print("\nnot enough converged polytopes to fit a scaling law.")
        return rows, summary

    ns = np.array([r['n'] for r in valid], dtype=float)
    rs = np.array([r['ratio_upper'] for r in valid], dtype=float)
    ms = np.array([r['mixing_step'] for r in valid], dtype=float)
    ln, lr, lm = np.log(ns), np.log(rs), np.log(ms)

    print("\n=== Scaling fit ===")
    print("Reference O*(n^2 (R/r)^2) implies exponent 2 on each axis.")
    if ln.std() < 1e-12:
        print(f"  n is constant at {int(ns[0])} across these polytopes "
              f"(one benchmark = one dimension) -- the n exponent is NOT "
              f"identifiable; pass several --smt2 files of different dim to fit it.")
        a_r, _, r2_r = _fit_loglog_exponent(rs, ms)
        summary.update({'ratio_exponent': a_r, 'ratio_fit_r2': r2_r,
                        'n_exponent': None})
        print(f"  mixing_step ~ (R/r)^{a_r:.2f}  (R^2={r2_r:.3f}, n fixed)")
        return rows, summary

    corr = float(np.corrcoef(ln, lr)[0, 1])
    M = np.column_stack([ln, lr, np.ones_like(ln)])
    coef, *_ = np.linalg.lstsq(M, lm, rcond=None)
    pred = M @ coef
    ss_res = float(np.sum((lm - pred) ** 2))
    ss_tot = float(np.sum((lm - lm.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

    a_n_m, _, r2_n_m = _fit_loglog_exponent(ns, ms)
    a_r_m, _, r2_r_m = _fit_loglog_exponent(rs, ms)
    summary.update({'n_exponent': float(coef[0]), 'ratio_exponent': float(coef[1]),
                    'joint_r2': r2, 'loglog_corr': corr,
                    'n_exponent_marginal': a_n_m, 'ratio_exponent_marginal': a_r_m})

    print(f"  joint:    mixing_step ~ n^{coef[0]:.2f} (R/r)^{coef[1]:.2f}  (R^2={r2:.3f})")
    print(f"  marginal: mixing_step ~ n^{a_n_m:.2f}          (R^2={r2_n_m:.3f})")
    print(f"  marginal: mixing_step ~ (R/r)^{a_r_m:.2f}      (R^2={r2_r_m:.3f})")
    print(f"  corr(log n, log R/r) = {corr:+.3f}")
    if abs(corr) > 0.9:
        print("  [warn] n and R/r are nearly collinear in this set; the two joint "
              "exponents are not separately identifiable. Widen the benchmark mix.")
    elif coef[1] > coef[0]:
        print("  -> R/r (fatness) carries the larger exponent: FATNESS looks like "
              "the bottleneck on these benchmarks.")
    elif coef[0] > coef[1]:
        print("  -> n carries the larger exponent: DIMENSION looks like the "
              "bottleneck on these benchmarks.")
    return rows, summary


# ===========================================================================
# 6. Plotting (optional)
# ===========================================================================

def plot_step_sweep(rows, poly, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    steps = [r['num_steps'] for r in rows]
    var_ratio = [r['variance_ratio'] for r in rows]
    ks_p = [r['ks_pvalue'] for r in rows]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    axes[0].axhline(1.0, color='gray', linestyle='--', linewidth=1)
    axes[0].plot(steps, var_ratio, marker='o')
    axes[0].set_xscale('log')
    axes[0].set_xlabel('num_steps')
    axes[0].set_ylabel('empirical / theoretical variance')
    axes[0].set_title('Variance convergence')

    axes[1].axhline(0.05, color='gray', linestyle='--', linewidth=1)
    axes[1].plot(steps, ks_p, marker='o', color='tab:orange')
    axes[1].set_xscale('log')
    axes[1].set_yscale('log')
    axes[1].set_xlabel('num_steps')
    axes[1].set_ylabel('KS p-value')
    axes[1].set_title('KS test vs uniform')

    fig.suptitle(f"{poly.kind}, n={poly.n}, R/r={poly.ratio:.2f}")
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_bottleneck(results, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))

    ns = [r['n'] for r in results['vs_n'] if r['mixing_step'] is not None]
    ms_n = [r['mixing_step'] for r in results['vs_n'] if r['mixing_step'] is not None]
    axes[0].loglog(ns, ms_n, marker='o')
    axes[0].set_xlabel('n')
    axes[0].set_ylabel('empirical mixing step')
    axes[0].set_title('Mixing step vs dimension n')

    rs = [r['R_r_ratio'] for r in results['vs_ratio'] if r['mixing_step'] is not None]
    ms_r = [r['mixing_step'] for r in results['vs_ratio'] if r['mixing_step'] is not None]
    axes[1].loglog(rs, ms_r, marker='o', color='tab:orange')
    axes[1].set_xlabel('R / r')
    axes[1].set_ylabel('empirical mixing step')
    axes[1].set_title('Mixing step vs fatness R/r')

    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_benchmark(rows, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    valid = [r for r in rows if r['mixing_step'] is not None]
    if len(valid) < 2:
        print("[warn] too few converged polytopes to plot")
        return
    ns = [r['n'] for r in valid]
    rs = [r['ratio_upper'] for r in valid]
    ms = [r['mixing_step'] for r in valid]

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, xs, xlabel, color in ((axes[0], ns, 'n', '#2a78d6'),
                                   (axes[1], rs, 'R / r (upper bound)', '#eb6834')):
        ax.scatter(xs, ms, s=34, color=color, zorder=3,
                   edgecolor='#fcfcfb', linewidth=0.8)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(xlabel)
        ax.set_ylabel('empirical mixing step')
        ax.grid(True, which='major', color='#e6e5e2', lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        if len(set(xs)) > 1:
            a, c, r2 = _fit_loglog_exponent(xs, ms)
            xf = np.linspace(min(xs), max(xs), 50)
            ax.plot(xf, np.exp(c) * xf ** a, color='#52514e', lw=1.2,
                    ls=(0, (4, 3)), zorder=2)
            ax.set_title(f"exponent {a:.2f} (R^2={r2:.2f})", fontsize=10, loc='left')
        else:
            ax.set_xlim(xs[0] * 0.5, xs[0] * 2.0)
            ax.set_title(f"{xlabel} constant at {xs[0]:g} -- no exponent to fit",
                         fontsize=10, loc='left')

    fig.suptitle(f"Mixing step across {len(valid)} benchmark polytopes", x=0.01,
                 ha='left')
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


# ===========================================================================
# 7. Self-test: cross-check the construction against independent methods
# ===========================================================================

def selftest():
    ok = True
    print("Self-test 1: box r/R closed form vs numeric vertex-based check (n=3)")
    poly = make_box(3, 6.0)
    verts = np.array(list(itertools.product(*[(-l / 2, l / 2) for l in poly.params['L']])))
    _, R_numeric = smallest_enclosing_ball(verts)
    _, r_numeric = chebyshev_center(poly.A, poly.b)
    pass1 = abs(R_numeric - poly.R) < 1e-4 and abs(r_numeric - poly.r) < 1e-4
    print(f"  closed-form (r={poly.r:.5f}, R={poly.R:.5f}) vs "
          f"numeric (r={r_numeric:.5f}, R={R_numeric:.5f}): {'PASS' if pass1 else 'FAIL'}")
    ok &= pass1

    print("Self-test 2: simplex marginal formula vs an INDEPENDENT uniform sampler "
          "(Dirichlet / normalized-exponential trick, not hit-and-run)")
    n = 6
    rng = np.random.default_rng(0)
    poly2 = make_simplex(n, 12.0)
    kappa = poly2.params['kappa']
    N = 200_000
    E = rng.exponential(size=(N, n + 1))
    D = E / E.sum(axis=1, keepdims=True)
    X = D[:, :n]
    Y0 = kappa * X[:, 0]
    _, true_var = theoretical_marginal(poly2, 0)
    emp_var = Y0.var()
    cdf, _ = theoretical_marginal(poly2, 0)
    ks = stats.kstest(Y0, cdf)
    pass2 = abs(emp_var - true_var) / true_var < 0.02 and ks.pvalue > 0.01
    print(f"  theoretical var={true_var:.5f}, Dirichlet-sampled var={emp_var:.5f}, "
          f"KS p={ks.pvalue:.3f}: {'PASS' if pass2 else 'FAIL'}")
    ok &= pass2

    print("Self-test 3: hit-and-run itself recovers known variance on a stretched box")
    poly3 = make_box(3, 10.0)
    rng2 = np.random.default_rng(3)
    samples = run_parallel_chains(poly3, num_steps=400, num_chains=2000, rng=rng2)
    row = evaluate_convergence(poly3, samples, coord=0)
    pass3 = abs(row['variance_ratio'] - 1) < 0.1 and row['ks_pvalue'] > 0.01
    print(f"  variance_ratio={row['variance_ratio']:.3f}, KS p={row['ks_pvalue']:.3f}: "
          f"{'PASS' if pass3 else 'FAIL'}")
    ok &= pass3

    print("Self-test 4: .ine round-trip -- write a synthetic body in cdd format, read "
          "it back through the benchmark path, check r exactly and R's bracket")
    pass4 = True
    tmpd = tempfile.mkdtemp(prefix="har_selftest_")
    try:
        for label, p in (("box n=4 rho=6", make_box(4, 6.0)),
                          ("simplex n=5 rho=9", make_simplex(5, 9.0))):
            path = os.path.join(tmpd, "t_cube1.ine")
            with open(path, "w") as fh:
                fh.write(f"t\nH-representation\nbegin\n {p.A.shape[0]} {p.n + 1} real\n")
                for i in range(p.A.shape[0]):
                    fh.write(" " + " ".join(["%.17g" % p.b[i]]
                                            + ["%.17g" % -a for a in p.A[i]]) + "\n")
                fh.write("end\n")
            q, reason = polytope_from_ine(path)
            if q is None:
                print(f"  {label}: FAIL ({reason})")
                pass4 = False
                continue
            r_ok = abs(q.r - p.r) < 1e-6 * max(1.0, p.r)
            brack = (q.params['R_lower'] <= p.R * (1 + 1e-6)
                     and p.R <= q.params['R_upper'] * (1 + 1e-6))
            print(f"  {label}: r {q.r:.6f} vs {p.r:.6f}, R in "
                  f"[{q.params['R_lower']:.4f}, {q.params['R_upper']:.4f}] vs exact "
                  f"{p.R:.4f}: {'PASS' if (r_ok and brack) else 'FAIL'}")
            pass4 &= r_ok and brack
    finally:
        shutil.rmtree(tmpd, ignore_errors=True)
    ok &= pass4

    print(f"\nSelf-test {'PASSED' if ok else 'FAILED'}")
    return ok


# ===========================================================================
# 8. CLI
# ===========================================================================

def build_argparser():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mode',
                    choices=['single', 'step_sweep', 'bottleneck', 'benchmark', 'selftest'],
                    default='single')
    p.add_argument('--smt2', nargs='+', default=None,
                    help='SMT2 benchmark(s), plain or .xz. Each is passed through '
                         '`ttc --dump-ine` and every polytope it emits is analysed. '
                         'Give files of DIFFERENT dim to make the n exponent '
                         'identifiable -- one file is one dimension.')
    p.add_argument('--ine', nargs='+', default=None,
                    help='use existing .ine files (paths, globs, or a directory) '
                         'instead of re-running ttc')
    p.add_argument('--ttc_bin', default='./ttc', help='ttc binary (default: ./ttc)')
    p.add_argument('--ttc_args', default=None,
                    help='extra args forwarded to ttc, quoted, e.g. "--no-cdd-simp"')
    p.add_argument('--keep_ine', default=None,
                    help='directory to copy the dumped .ine files into')
    p.add_argument('--max_polytopes', type=int, default=0,
                    help='cap the number of polytopes analysed (0 = all)')
    p.add_argument('--poly_index', type=int, default=0,
                    help='which benchmark polytope to use in single/step_sweep mode')
    p.add_argument('--ref_chains', type=int, default=2000,
                    help='chains in the long-run reference sample (benchmark modes)')
    p.add_argument('--ref_steps', type=int, default=5000,
                    help='steps per reference chain; must comfortably exceed the '
                         'mixing step being measured or the test is optimistic')
    p.add_argument('--projection', choices=['coord', 'longest', 'random'],
                    default='longest',
                    help='marginal tested for benchmark polytopes (default: the '
                         'longest bounding-box axis)')
    p.add_argument('--n', type=int, default=10, help='ambient dimension')
    p.add_argument('--num_steps', type=int, default=500, help='hit-and-run steps per chain')
    p.add_argument('--R_r_ratio', type=float, default=8.0, help='target fatness index R/r')
    p.add_argument('--kind', choices=['box', 'simplex'], default='box')
    p.add_argument('--num_chains', type=int, default=500,
                    help='independent chains per step count (batch size for the test)')
    p.add_argument('--coord', type=int, default=0, help='coordinate to evaluate (0 = stretched axis)')
    p.add_argument('--var_tol', type=float, default=0.05, help='variance-ratio convergence tolerance')
    p.add_argument('--ks_alpha', type=float, default=0.05, help='KS test significance level')
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--plot', action='store_true', help='save a PNG plot of the sweep')
    p.add_argument('--out_prefix', type=str, default='har_mixing', help='output file prefix')
    p.add_argument('--step_grid', type=str, default=None,
                    help='comma-separated step counts, e.g. "10,25,50,100,250,500,1000"')
    p.add_argument('--n_grid', type=str, default='4,8,16,32',
                    help='(bottleneck mode) comma-separated n values')
    p.add_argument('--ratio_grid', type=str, default=None,
                    help='(bottleneck mode) comma-separated R/r values')
    p.add_argument('--start_steps', type=int, default=10,
                    help='(bottleneck mode) starting step count for the adaptive search')
    p.add_argument('--max_steps', type=int, default=20000,
                    help='(bottleneck mode) cap on the adaptive search; configs that have not '
                         'converged by this point are reported and excluded from the fit')
    return p


def default_step_grid():
    return [10, 25, 50, 100, 250, 500, 1000, 2000]


def main(argv=None):
    args = build_argparser().parse_args(argv)
    t0 = time.time()

    if args.mode == 'selftest':
        ok = selftest()
        sys.exit(0 if ok else 1)

    step_grid = ([int(s) for s in args.step_grid.split(',')]
                 if args.step_grid else default_step_grid())

    # A benchmark source replaces the synthetic construction everywhere.
    workdir = tempfile.mkdtemp(prefix="har_ine_")
    bench_polys = []
    try:
        if args.smt2 or args.ine:
            bench_polys = load_benchmark_polytopes(args, workdir)
            if not bench_polys:
                sys.exit("no usable polytopes from the benchmark(s)")
            if args.keep_ine:
                os.makedirs(args.keep_ine, exist_ok=True)
                for p in bench_polys:
                    shutil.copy(p.params['path'], args.keep_ine)
                print(f"copied {len(bench_polys)} .ine files to {args.keep_ine}")

        def pick_polytope():
            """Benchmark polytope if one was loaded, else the synthetic body."""
            if bench_polys:
                idx = min(args.poly_index, len(bench_polys) - 1)
                return bench_polys[idx]
            return make_polytope(args.n, args.R_r_ratio, args.kind)

        def reference_for(poly, rng):
            """Reference sample + projection, or (None, None) when the body has
            an analytic marginal to test against."""
            if poly.kind not in ('box', 'simplex'):
                w = projection_direction(poly, args.projection, args.coord, rng)
                ref = reference_sample(poly, rng, args.ref_chains, args.ref_steps)
                check = reference_selfcheck(ref, w)
                if check['ok'] is False:
                    print(f"[warn] reference halves disagree (KS p={check['pvalue']:.3g}); "
                          f"raise --ref_steps above {args.ref_steps}")
                return ref, w
            return None, None

        if args.mode == 'single':
            poly = pick_polytope()
            print(poly)
            rng = np.random.default_rng(args.seed)
            ref, w = reference_for(poly, rng)
            samples = run_parallel_chains(poly, args.num_steps, args.num_chains, rng)
            row = evaluate_convergence(poly, samples, coord=args.coord, ref=ref, w=w)
            row['num_steps'] = args.num_steps
            print(json.dumps(row, indent=2))

        elif args.mode == 'step_sweep':
            poly = pick_polytope()
            print(poly)
            rng = np.random.default_rng(args.seed)
            ref, w = reference_for(poly, rng)
            rows, mixing_step = step_sweep(poly, step_grid, args.num_chains, rng,
                                            coord=args.coord, var_tol=args.var_tol,
                                            ks_alpha=args.ks_alpha, ref=ref, w=w)
            print(f"\nEstimated empirical mixing step: {mixing_step}")
            bound = theoretical_bound(poly.n, poly.ratio)
            print(f"Reference O*(n^2 (R/r)^2) value (unnormalized): {bound:.3g}")
            if args.plot:
                out_path = f"{args.out_prefix}_step_sweep.png"
                plot_step_sweep(rows, poly, out_path)
                print(f"Saved plot to {out_path}")

        elif args.mode == 'benchmark':
            if not bench_polys:
                sys.exit("--mode benchmark needs --smt2 or --ine")
            rows, summary = benchmark_analysis(
                bench_polys, num_chains=args.num_chains, seed=args.seed,
                ref_chains=args.ref_chains, ref_steps=args.ref_steps,
                projection=args.projection, coord=args.coord,
                var_tol=args.var_tol, ks_alpha=args.ks_alpha,
                start_steps=args.start_steps, max_steps=args.max_steps)
            out_json = f"{args.out_prefix}_benchmark.json"
            with open(out_json, "w") as fh:
                json.dump({'rows': rows, 'summary': summary}, fh, indent=2)
            print(f"\nSaved results to {out_json}")
            if args.plot:
                out_path = f"{args.out_prefix}_benchmark.png"
                plot_benchmark(rows, out_path)
                print(f"Saved plot to {out_path}")

        elif args.mode == 'bottleneck':
            n_grid = [int(x) for x in args.n_grid.split(',')]
            ratio_grid = ([float(x) for x in args.ratio_grid.split(',')] if args.ratio_grid
                          else [np.sqrt(args.n) * m for m in (1.5, 3, 6, 12)])
            results, summary = bottleneck_analysis(
                args.kind, n_grid, ratio_grid, base_n=args.n, base_ratio=args.R_r_ratio,
                num_chains=args.num_chains, seed=args.seed,
                var_tol=args.var_tol, ks_alpha=args.ks_alpha,
                start_steps=args.start_steps, max_steps=args.max_steps)
            if args.plot:
                out_path = f"{args.out_prefix}_bottleneck.png"
                plot_bottleneck(results, out_path)
                print(f"Saved plot to {out_path}")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    print(f"\n[done in {time.time() - t0:.1f}s]")


if __name__ == '__main__':
    main()
