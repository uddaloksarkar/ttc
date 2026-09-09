#!/usr/bin/env python3
"""Report the fatness index R/r of every polytope in an SMT2 LRA benchmark.

Runs `ttc --dump-ine` on each input, parses the cdd H-representation files it
writes, and prints one row per polytope.  The geometry code is imported from
har_mixing.py so both tools agree by construction.

r is exact (Chebyshev LP).  R -- the minimal enclosing ball radius -- is
bracketed rather than computed exactly, because exact vertex enumeration is
exponential in the dimension:

    R_lo : enclosing ball of the 2n LP optima, which are points of the body
    R_hi : ball around the LP bounding box, which certainly contains the body

so the true R/r lies in [ratio_lo, ratio_hi].  Rows print both; the CSV keeps
every column.

Usage
-----
    ./polytope_ratios.py bench.smt2.xz
    ./polytope_ratios.py *.smt2.xz --ttc_bin ./build/ttc --csv ratios.csv
    ./polytope_ratios.py --ine dumped/          # reuse .ine files already dumped
"""

import argparse
import csv
import math
import os
import re
import shutil
import statistics
import sys
import tempfile

from har_mixing import (dump_ine_from_smt2, polytope_from_ine)

FIELDS = ["source", "ine", "index", "n", "facets", "r",
          "R_lower", "R_upper", "ratio_lower", "ratio_upper",
          "volume", "expected_samples", "samples_deleted", "store_after"]

# Constants of the union-volume estimator, from lra_volume.cpp:818-838.
K_EPSILON = 0.8
K_DELTA = 0.2

# "c        1       12.345678                3             100"
RE_VOLROW = re.compile(r"^c\s+(\d+)\s+([\d.eE+-]+)\s+(\d+)\s+(\d+)\s*$")
RE_TOTAL_GEN = re.compile(r"^c total samples generated:\s*(\d+)")
RE_TOTAL_DEL = re.compile(r"^c total samples deleted:\s*(\d+)")


def sample_threshold(num_polytopes):
    """thresh from lra_volume.cpp:836 -- the (eps, delta) sample-size target
    that bounds how many points the estimator keeps at once."""
    eps = K_EPSILON / 2.0
    return max(12.0 * math.log(24.0 / K_DELTA) / (eps * eps),
               6.0 * (math.log(6.0 / K_DELTA) + math.log(max(1.0, num_polytopes))))


def parse_ttc_volumes(output):
    """(rows, totals) from ttc's stdout: the per-polytope volume table and the
    run's sample counters.  rows are in the order the union algorithm processed
    them, which is the dump order minus any cube dropped after being written."""
    rows, totals = [], {}
    for line in output.splitlines():
        m = RE_VOLROW.match(line)
        if m:
            rows.append({'index': int(m.group(1)), 'volume': float(m.group(2)),
                         'samples_deleted': int(m.group(3)),
                         'store_after': int(m.group(4))})
            continue
        m = RE_TOTAL_GEN.match(line)
        if m:
            totals['generated'] = int(m.group(1))
            continue
        m = RE_TOTAL_DEL.match(line)
        if m:
            totals['deleted'] = int(m.group(1))
    return rows, totals


def expected_draws(vol_rows, thresh):
    """Replay the estimator's sampling probability to predict E[n_i].

    The algorithm draws n_i ~ Poisson(p * vol_i) and, before each polytope,
    halves p (down-sampling the stored set with it) while p * vol_i > thresh
    (lra_volume.cpp:704-708).  p never increases, so that schedule follows from
    the volumes alone.  A second rule then re-draws while n + |store| > thresh
    (lines 713-719), which caps the draw at the room left in the store; the
    stored sizes come from ttc's own table, so we apply that cap too.

    Even so this over-predicts -- roughly 1.4x on a dim-9 benchmark -- because
    each re-draw also halves p for every later polytope.  Read it as an upper
    estimate on the sampling cost, not a forecast.
    """
    p, out, prev_store = 1.0, [], 0
    for row in vol_rows:
        v = row['volume']
        while p * v > thresh:
            p /= 2.0
        room = max(0.0, thresh - max(0, prev_store - row['samples_deleted']))
        out.append(min(p * v, room))
        prev_store = row['store_after']
    return out


def collect(paths, source_of):
    """(rows, skipped) for a list of .ine files."""
    rows, skipped = [], []
    for i, path in enumerate(paths):
        poly, reason = polytope_from_ine(path)
        if poly is None:
            skipped.append((os.path.basename(path), reason))
            continue
        rows.append({
            "source": source_of(path),
            "ine": os.path.basename(path),
            "index": i,
            "n": poly.n,
            "facets": poly.params["facets"],
            "r": poly.r,
            "R_lower": poly.params["R_lower"],
            "R_upper": poly.params["R_upper"],
            "ratio_lower": poly.params["ratio_lower"],
            "ratio_upper": poly.ratio,
        })
    return rows, skipped


RUN_TOTALS = []


def attach_sample_counts(rows, ttc_output, source):
    """Join ttc's volume table onto our per-polytope rows.

    ttc writes a cube's .ine *before* the inner-ball check that can still drop
    it (lra_volume.cpp:880-892), so the dumped files can outnumber the volume
    rows.  Both sequences are in processing order, so when the counts agree the
    join is positional; when they do not, per-polytope counts are left blank
    rather than silently mis-aligned -- the run totals are still exact.
    """
    vol_rows, totals = parse_ttc_volumes(ttc_output)
    thresh = sample_threshold(len(vol_rows))
    totals['thresh'] = thresh
    totals['source'] = source
    totals['polytopes'] = len(vol_rows)

    if vol_rows and len(vol_rows) == len(rows):
        exp = expected_draws(vol_rows, thresh)
        for row, vr, e in zip(rows, vol_rows, exp):
            row['volume'] = vr['volume']
            row['samples_deleted'] = vr['samples_deleted']
            row['store_after'] = vr['store_after']
            row['expected_samples'] = e
        totals['expected_total'] = sum(exp)
    elif vol_rows:
        totals['misaligned'] = (len(rows), len(vol_rows))
        totals['expected_total'] = sum(expected_draws(vol_rows, thresh))
    RUN_TOTALS.append(totals)


def print_sample_summary():
    if not RUN_TOTALS:
        return
    print("\nSampling cost (from ttc's own run)")
    for t in RUN_TOTALS:
        print(f"  {t['source']}")
        print(f"    (eps={K_EPSILON}, delta={K_DELTA}) sample threshold: "
              f"{t['thresh']:.0f} points held at once")
        if 'expected_total' in t:
            print(f"    predicted total draws   : {t['expected_total']:,.0f}  "
                  f"(upper estimate over {t['polytopes']} polytopes)")
        if 'generated' in t:
            print(f"    actual samples generated: {t['generated']:,}")
        if 'deleted' in t:
            print(f"    actual samples deleted  : {t['deleted']:,}")
        if 'misaligned' in t:
            got, want = t['misaligned']
            print(f"    [note] {got} .ine files vs {want} volume rows -- some cubes "
                  f"were dumped then dropped, so per-polytope counts are omitted")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("smt2", nargs="*", help="SMT2 benchmark(s), plain or .xz")
    ap.add_argument("--ine", nargs="+", default=None,
                    help="use existing .ine files (paths, globs or a directory) "
                         "instead of running ttc")
    ap.add_argument("--ttc_bin", default="./ttc", help="ttc binary (default: ./ttc)")
    ap.add_argument("--ttc_args", default=None,
                    help='extra args forwarded to ttc, quoted, e.g. "--no-cdd-simp"')
    ap.add_argument("--csv", default=None, help="also write the rows to this CSV")
    ap.add_argument("--keep_ine", default=None,
                    help="directory to keep the dumped .ine files in")
    ap.add_argument("--quiet", action="store_true", help="suppress the ttc command line")
    args = ap.parse_args()

    if not args.smt2 and not args.ine:
        ap.error("give one or more SMT2 files, or --ine")

    workdir = tempfile.mkdtemp(prefix="ratios_")
    try:
        rows, skipped = [], []
        if args.ine:
            import glob
            paths = []
            for pat in args.ine:
                paths.extend(sorted(glob.glob(os.path.join(pat, "*.ine"))
                                    if os.path.isdir(pat) else glob.glob(pat)))
            if not paths:
                sys.exit(f"no .ine files matched {args.ine}")
            r, s = collect(paths, lambda p: os.path.basename(os.path.dirname(
                os.path.abspath(p))))
            rows += r
            skipped += s
        else:
            for k, smt2 in enumerate(args.smt2):
                sub = os.path.join(workdir, f"b{k}")
                os.makedirs(sub, exist_ok=True)
                paths, out = dump_ine_from_smt2(
                    smt2, args.ttc_bin, sub,
                    extra_args=args.ttc_args.split() if args.ttc_args else (),
                    verbose=not args.quiet, return_output=True)
                name = os.path.basename(smt2)
                r, s = collect(paths, lambda p, name=name: name)
                attach_sample_counts(r, out, name)
                rows += r
                skipped += s
                if args.keep_ine:
                    os.makedirs(args.keep_ine, exist_ok=True)
                    for p in paths:
                        shutil.copy(p, args.keep_ine)

        if skipped:
            print(f"\n[warn] skipped {len(skipped)} polytope(s):")
            for name, reason in skipped:
                print(f"         {name}: {reason}")
        if not rows:
            sys.exit("no usable polytopes")

        has_samples = any('expected_samples' in r for r in rows)
        extra = f" {'volume':>13} {'E[samples]':>11} {'deleted':>8}" if has_samples else ""
        head = (f"\n{'#':>4} {'n':>4} {'facets':>7} {'r':>12} "
                f"{'R_lo':>12} {'R_hi':>12} {'R/r lo':>10} {'R/r hi':>10}"
                f"{extra}  polytope")
        print(head)
        print("-" * (len(head) + 10))
        for row in rows:
            line = (f"{row['index']:>4} {row['n']:>4} {row['facets']:>7} "
                    f"{row['r']:>12.5g} {row['R_lower']:>12.5g} {row['R_upper']:>12.5g} "
                    f"{row['ratio_lower']:>10.4g} {row['ratio_upper']:>10.4g}")
            if has_samples:
                if 'expected_samples' in row:
                    line += (f" {row['volume']:>13.6g} {row['expected_samples']:>11.1f}"
                             f" {row['samples_deleted']:>8}")
                else:
                    line += f" {'-':>13} {'-':>11} {'-':>8}"
            print(line + f"  {row['ine']}")

        hi = [r["ratio_upper"] for r in rows]
        lo = [r["ratio_lower"] for r in rows]
        print(f"\n{len(rows)} polytopes, dim {min(r['n'] for r in rows)}"
              f"-{max(r['n'] for r in rows)}")
        print(f"  R/r upper bound: min {min(hi):.4g}, median {statistics.median(hi):.4g}, "
              f"max {max(hi):.4g}")
        print(f"  R/r lower bound: min {min(lo):.4g}, median {statistics.median(lo):.4g}, "
              f"max {max(lo):.4g}")
        widest = max(rows, key=lambda r: r["ratio_upper"] / r["ratio_lower"])
        print(f"  widest bracket:  {widest['ine']} "
              f"[{widest['ratio_lower']:.4g}, {widest['ratio_upper']:.4g}] "
              f"({widest['ratio_upper'] / widest['ratio_lower']:.2f}x)")

        print_sample_summary()

        if args.csv:
            with open(args.csv, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=FIELDS)
                w.writeheader()
                w.writerows(rows)
            print(f"\nwrote {len(rows)} rows -> {args.csv}")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)


if __name__ == "__main__":
    main()
