#!/usr/bin/env python3
"""Compare volume estimates for one benchmark across samplers and walk lengths.

Reads the CSV produced by scripts/collect_ttc.py and prints, for a single
benchmark:

  1. sampler x walklen table at the two walk lengths of interest (default
     10 and 100), each cell averaged over seeds;
  2. per-sampler walk-length sweep over every walk length present.

Volumes are averaged over seeds and reported both absolutely and as a
percentage deviation, since the absolute numbers (~1e22) are unreadable.
Because there is no ground truth, deviations are relative to the grand mean of
all solved runs of that benchmark -- a spread, not an error.

Usage:
    ./scripts/compare_bench.py results-2001477.csv dim10_numpol10_seed638
    ./scripts/compare_bench.py results-2001477.csv --list
    ./scripts/compare_bench.py results-2001477.csv <bench> --walklens 10 1000 --seeds
"""

import argparse
import csv
import statistics
import sys
from collections import defaultdict


def load(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def fnum(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return None


def pick_benchmark(rows, needle):
    names = sorted({r["benchmark"] for r in rows})
    if needle in names:
        return needle
    matches = [n for n in names if needle in n]
    if not matches:
        sys.exit(f"no benchmark matching {needle!r} "
                 f"({len(names)} available; use --list)")
    if len(matches) > 1:
        head = "\n  ".join(matches[:10])
        sys.exit(f"{needle!r} is ambiguous, {len(matches)} matches:\n  {head}"
                 + ("\n  …" if len(matches) > 10 else ""))
    return matches[0]


class Cell:
    """The seed replicates of one (sampler, walklen) pair."""

    def __init__(self):
        self.vols = []
        self.samp = []
        self.wall = []
        self.seeds = []
        self.unsolved = 0

    def add(self, row):
        if row["solved"] != "1":
            self.unsolved += 1
            return
        v = fnum(row["volume"])
        if v is None:
            self.unsolved += 1
            return
        self.vols.append(v)
        self.seeds.append((row["seed"], v))
        for key, sink in (("t_sampling", self.samp), ("wall_s", self.wall)):
            x = fnum(row[key])
            if x is not None:
                sink.append(x)

    @property
    def n(self):
        return len(self.vols)

    @property
    def mean(self):
        return statistics.fmean(self.vols) if self.vols else None

    @property
    def cv(self):
        """Seed-to-seed coefficient of variation, in percent."""
        if len(self.vols) < 2:
            return None
        m = statistics.fmean(self.vols)
        if m == 0:
            return None
        return 100.0 * statistics.stdev(self.vols) / m

    def mean_of(self, which):
        xs = self.samp if which == "samp" else self.wall
        return statistics.fmean(xs) if xs else None


def fmt(x, spec="{:.4g}", width=12, dash="-"):
    return f"{dash:>{width}}" if x is None else f"{x:>{width}{spec[2:-1]}}"


def pct(x, width=9):
    return f"{'-':>{width}}" if x is None else f"{x:>+{width}.2f}"


def rel(value, ref):
    if value is None or ref in (None, 0):
        return None
    return 100.0 * (value / ref - 1.0)


# ---------------------------------------------------------------------------
# sweep mode: every benchmark, every sampler, as a ratio to the denominator
# ---------------------------------------------------------------------------

# Categorical slots 1-3 of the validated reference palette.  Three series is the
# documented cap for an all-pairs chart form; worst all-pairs CVD dE 9.2.
SERIES_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#8a8a85"
SURFACE = "#fcfcfb"


def build_cells(rows):
    """(benchmark, sampler, walklen) -> Cell."""
    cells = defaultdict(Cell)
    for r in rows:
        cells[(r["benchmark"], r["sampler"], r["walklen"])].add(r)
    return cells


def sweep_ratios(rows, walklens, denom):
    """{walklen: {sampler: [(benchmark, ratio), ...]}}, seed-averaged.

    A benchmark contributes to a (walklen, sampler) series only when both that
    sampler and the denominator produced a seed-averaged volume at that walk
    length, so every ratio compares like with like.
    """
    cells = build_cells(rows)
    benches = sorted({b for b, _, _ in cells})
    samplers = sorted({s for _, s, _ in cells})
    if denom not in samplers:
        sys.exit(f"denominator sampler {denom!r} not in data: {samplers}")

    out, skipped = {}, defaultdict(int)
    for wl in walklens:
        per_sampler = {}
        for s in samplers:
            if s == denom:
                continue
            pairs = []
            for b in benches:
                num = cells[(b, s, wl)].mean if (b, s, wl) in cells else None
                den = cells[(b, denom, wl)].mean if (b, denom, wl) in cells else None
                if num is None or den in (None, 0):
                    skipped[(wl, s)] += 1
                    continue
                pairs.append((b, num / den))
            per_sampler[s] = pairs
        out[wl] = per_sampler
    return out, skipped


def geomean(xs):
    return statistics.geometric_mean(xs) if xs else None


def print_sweep_table(ratios, skipped, denom):
    """The table view -- also the accessibility relief for the low-contrast hue."""
    print(f"\nratio of seed-averaged volume to {denom}, over all benchmarks\n")
    head = (f"{'walklen':>8}  {'sampler':<10}{'n':>5}{'geomean':>10}{'median':>10}"
            f"{'p10':>9}{'p90':>9}{'max|dev|%':>11}{'frac>1':>9}{'skipped':>9}")
    print(head)
    print("-" * len(head))
    for wl, per_sampler in ratios.items():
        for s, pairs in per_sampler.items():
            vals = [r for _, r in pairs]
            if not vals:
                print(f"{wl:>8}  {s:<10}{0:>5}")
                continue
            q = statistics.quantiles(vals, n=10) if len(vals) >= 10 else None
            worst = max(abs(100.0 * (v - 1.0)) for v in vals)
            frac = sum(1 for v in vals if v > 1.0) / len(vals)
            print(f"{wl:>8}  {s:<10}{len(vals):>5}"
                  f"{geomean(vals):>10.4f}{statistics.median(vals):>10.4f}"
                  f"{(q[0] if q else float('nan')):>9.4f}"
                  f"{(q[-1] if q else float('nan')):>9.4f}"
                  f"{worst:>11.2f}{frac:>9.3f}{skipped[(wl, s)]:>9}")
        print()


def write_sweep_csv(ratios, path, denom):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["walklen", "sampler", "denominator", "benchmark", "ratio",
                    "dev_pct"])
        for wl, per_sampler in ratios.items():
            for s, pairs in per_sampler.items():
                for b, r in sorted(pairs):
                    w.writerow([wl, s, denom, b, f"{r:.6g}",
                                f"{100.0 * (r - 1.0):.4f}"])
    print(f"per-benchmark ratios -> {path}")


def plot_sweep(ratios, out_path, denom, nbench):
    """One panel per walk length; per sampler, the sorted ratio curve.

    Each sampler's ratios are sorted ascending and plotted against rank, so the
    x axis runs 1..nbench and every benchmark is one point on every curve; the
    midpoint is the median ratio.  Flat-and-on-1.0 means "agrees with the
    denominator everywhere"; a steep curve means benchmark-to-benchmark
    disagreement.  A curve that stops short of nbench had benchmarks skipped.
    """
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.ticker import FuncFormatter
    except ImportError:
        print("matplotlib not installed; skipping the plot "
              "(the CSV and the table above are unaffected)", file=sys.stderr)
        return

    wls = list(ratios)
    allvals = [r for per in ratios.values() for pairs in per.values()
               for _, r in pairs]
    if not allvals:
        print("nothing to plot", file=sys.stderr)
        return
    lo, hi = min(allvals), max(allvals)
    pad = 0.04 * (hi - lo) or 0.02
    ylim = (lo - pad, hi + pad)
    # a log axis over a ~0.8-1.2 span gets one decade tick from the default
    # locator, so place ratio ticks by hand at a readable step
    span = ylim[1] - ylim[0]
    step = next(s for s in (0.02, 0.05, 0.1, 0.2, 0.5) if span / s <= 8)
    ticks = []
    t = round(1.0 - step * int((1.0 - ylim[0]) / step + 1), 10)
    while t <= ylim[1]:
        if t >= ylim[0]:
            ticks.append(round(t, 10))
        t = round(t + step, 10)

    fig, axes = plt.subplots(1, len(wls), figsize=(5.6 * len(wls), 4.6),
                             sharey=True, squeeze=False)
    fig.patch.set_facecolor(SURFACE)

    for ax, wl in zip(axes[0], wls):
        ax.set_facecolor(SURFACE)
        ax.axhline(1.0, color=INK_MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
        ax.annotate(f"{denom} = 1.00", xy=(0.015, 1.0),
                    xycoords=("axes fraction", "data"), xytext=(0, 3),
                    textcoords="offset points",
                    va="bottom", ha="left", fontsize=8, color=INK_SECONDARY)

        ends = []
        for i, (s, pairs) in enumerate(sorted(ratios[wl].items())):
            vals = sorted(r for _, r in pairs)
            if not vals:
                continue
            xs = list(range(1, len(vals) + 1))
            color = SERIES_COLORS[i % len(SERIES_COLORS)]
            ax.plot(xs, vals, color=color, lw=2.0, solid_capstyle="round",
                    label=s, zorder=3)
            ends.append((vals[-1], s))

        # direct labels in ink (the curve end beside each carries identity),
        # nudged apart so the three ends do not overprint
        ax.set_yscale("log")
        ax.set_ylim(*ylim)
        import math
        f_lo, f_hi = math.log(ylim[0]), math.log(ylim[1])
        placed = []
        for v, s in sorted(ends):
            frac = (math.log(v) - f_lo) / (f_hi - f_lo)
            if placed and frac - placed[-1] < 0.055:
                frac = placed[-1] + 0.055
            placed.append(frac)
            ax.annotate(s, xy=(1.01, frac), xycoords="axes fraction",
                        va="center", ha="left", fontsize=9,
                        color=INK_SECONDARY, zorder=4, annotation_clip=False)

        ax.set_yticks(ticks)
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:.2f}"))
        ax.yaxis.set_minor_locator(plt.NullLocator())
        ax.set_title(f"walklen = {wl}", fontsize=11, color=INK_PRIMARY,
                     pad=8, loc="left")
        ax.set_xlabel(f"benchmarks, sorted by ratio (1-{nbench})", fontsize=9,
                      color=INK_SECONDARY)
        ax.set_xlim(0, nbench)
        ax.grid(axis="y", color="#e6e5e2", lw=0.8, zorder=0)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        for side in ("left", "bottom"):
            ax.spines[side].set_color("#d9d8d4")
        ax.tick_params(colors=INK_SECONDARY, labelsize=9, length=0)

    axes[0][0].set_ylabel(f"volume ratio  (sampler / {denom})", fontsize=9,
                          color=INK_SECONDARY)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", frameon=False, fontsize=9,
               labelcolor=INK_SECONDARY, ncol=len(labels),
               bbox_to_anchor=(0.99, 0.99))
    fig.suptitle(f"Volume agreement with {denom} across {nbench} benchmarks",
                 fontsize=13, color=INK_PRIMARY, x=0.008, ha="left", y=0.985)
    fig.tight_layout(rect=(0, 0, 0.97, 0.93))
    fig.savefig(out_path, dpi=160, facecolor=SURFACE)
    print(f"plot -> {out_path}")


def run_sweep(rows, args):
    ratios, skipped = sweep_ratios(rows, args.walklens, args.denom)
    nbench = len({r["benchmark"] for r in rows})
    print_sweep_table(ratios, skipped, args.denom)
    write_sweep_csv(ratios, args.csv_out, args.denom)
    plot_sweep(ratios, args.out, args.denom, nbench)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv", help="results CSV from collect_ttc.py")
    ap.add_argument("benchmark", nargs="?",
                    help="benchmark name or unique substring")
    ap.add_argument("--walklens", nargs="+", default=["10", "100"],
                    metavar="WL",
                    help="walk lengths to compare (default: 10 100); "
                         "single-benchmark mode uses the first two")
    ap.add_argument("--sweep", action="store_true",
                    help="sweep every benchmark and plot each sampler's volume "
                         "as a ratio to --denom, one panel per walk length")
    ap.add_argument("--denom", default="ball",
                    help="denominator sampler for --sweep (default: ball)")
    ap.add_argument("--out", default="ratio-sweep.png",
                    help="plot file for --sweep (default: ratio-sweep.png)")
    ap.add_argument("--csv-out", default="ratio-sweep.csv",
                    help="per-benchmark ratio CSV for --sweep "
                         "(default: ratio-sweep.csv)")
    ap.add_argument("--seeds", action="store_true",
                    help="also print the individual per-seed volumes")
    ap.add_argument("--list", action="store_true",
                    help="list available benchmarks and exit")
    args = ap.parse_args()

    rows = load(args.csv)
    if args.list:
        for n in sorted({r["benchmark"] for r in rows}):
            print(n)
        return
    if args.sweep:
        run_sweep(rows, args)
        return
    if not args.benchmark:
        ap.error("a benchmark is required (or pass --sweep / --list)")
    if len(args.walklens) < 2:
        ap.error("single-benchmark mode needs two --walklens")

    bench = pick_benchmark(rows, args.benchmark)
    sub = [r for r in rows if r["benchmark"] == bench]

    cells = defaultdict(Cell)
    for r in sub:
        cells[(r["sampler"], r["walklen"])].add(r)

    samplers = sorted({s for s, _ in cells})
    walklens = sorted({w for _, w in cells}, key=lambda w: int(w))

    solved = [v for c in cells.values() for v in c.vols]
    if not solved:
        sys.exit(f"{bench}: no solved runs to compare")
    grand = statistics.fmean(solved)

    print(f"benchmark : {bench}")
    print(f"runs      : {len(solved)} solved / {len(sub)} total")
    print(f"grand mean volume (all solved runs, the reference) : {grand:.6g}")

    # ---- 1. samplers side by side at the two walk lengths ----------------
    wl_a, wl_b = args.walklens[:2]
    print(f"\n[1] sampler comparison at walklen {wl_a} vs {wl_b} "
          f"(volume averaged over seeds)\n")
    head = (f"{'sampler':<10}"
            f"{'wl=' + wl_a:>14}{'dev%':>9}{'cv%':>8}{'n':>4}"
            f"{'wl=' + wl_b:>14}{'dev%':>9}{'cv%':>8}{'n':>4}"
            f"{'B/A-1 %':>10}")
    print(head)
    print("-" * len(head))
    for s in samplers:
        ca, cb = cells.get((s, wl_a), Cell()), cells.get((s, wl_b), Cell())
        line = f"{s:<10}"
        for c in (ca, cb):
            line += (fmt(c.mean, "{:.6g}", 14) + pct(rel(c.mean, grand))
                     + fmt(c.cv, "{:.2f}", 8) + f"{c.n:>4}")
        line += pct(rel(cb.mean, ca.mean), 10)
        print(line)

    for wl in (wl_a, wl_b):
        means = [cells[(s, wl)].mean for s in samplers if cells.get((s, wl))]
        means = [m for m in means if m is not None]
        if len(means) > 1 and min(means) > 0:
            print(f"  spread across samplers at wl={wl}: "
                  f"max/min = {max(means) / min(means):.4f} "
                  f"({100 * (max(means) / min(means) - 1):+.2f}%)")

    # ---- 2. walk-length sweep, one block per sampler ---------------------
    print(f"\n[2] walk-length comparison per sampler "
          f"(dev% vs that sampler's own wl={walklens[-1]})\n")
    for s in samplers:
        base = cells[(s, walklens[-1])].mean if (s, walklens[-1]) in cells else None
        print(f"{s}:")
        head = (f"  {'walklen':>8}{'mean volume':>16}{'dev%':>9}{'cv%':>8}"
                f"{'n':>4}{'t_samp':>9}{'wall':>9}{'unsolved':>10}")
        print(head)
        print("  " + "-" * (len(head) - 2))
        for wl in walklens:
            c = cells.get((s, wl))
            if c is None:
                continue
            print(f"  {wl:>8}" + fmt(c.mean, "{:.6g}", 16)
                  + pct(rel(c.mean, base)) + fmt(c.cv, "{:.2f}", 8)
                  + f"{c.n:>4}" + fmt(c.mean_of('samp'), "{:.2f}", 9)
                  + fmt(c.mean_of('wall'), "{:.2f}", 9)
                  + f"{c.unsolved:>10}")
            if args.seeds:
                for seed, v in sorted(c.seeds):
                    print(f"{'seed ' + seed:>18}{v:>16.6g}")
        print()

    print("note: dev% is spread, not error -- there is no ground-truth volume.")
    print("      cv% is the seed-to-seed spread within a cell; a walk-length")
    print("      effect is only real if it is large next to cv%.")


if __name__ == "__main__":
    main()
