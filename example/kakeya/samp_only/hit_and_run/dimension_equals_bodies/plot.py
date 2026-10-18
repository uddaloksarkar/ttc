#!/usr/bin/env python3
"""Compare completed fixed-3D and dimension-matched TTC pilot results."""

import csv
from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"


def save(fig, name):
    FIGURES.mkdir(exist_ok=True)
    fig.tight_layout()
    for suffix in ("pdf", "pgf", "png"):
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    plt.rcParams.update({"font.size": 10, "pgf.texsystem": "pdflatex"})
    old = {(int(r["bodies"]), int(r["walk_length"])): r
           for r in csv.DictReader((HERE.parent / "results.csv").open())
           if r["run_status"] == "0" and int(r["thin_denominator"]) == 5000000}
    new = {(int(r["bodies"]), int(r["walk_length"])): r
           for r in csv.DictReader((HERE / "results.csv").open())
           if r["run_status"] == "0"}
    all_rows = list(csv.DictReader((HERE / "results.csv").open()))
    fig, ax = plt.subplots(figsize=(7.0, 4.3))
    for count in (2, 4, 8, 16, 32):
        series = sorted((r for r in all_rows if int(r["bodies"]) == count
                         and r["run_status"] == "0"),
                        key=lambda r: int(r["walk_length"]))
        if series:
            ax.plot([int(r["walk_length"]) for r in series],
                    [100 * float(r["relative_error"]) for r in series],
                    marker="o", linewidth=1.7, markersize=4,
                    label=f"{count} bodies")
    ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
               label="20% error limit")
    ax.set_xscale("log")
    ax.set_xlim(4, 1300000)
    ax.set_xlabel("Sampler walk length")
    ax.set_ylabel("Union volume relative error (%)")
    ax.set_title("Dimension = body count: exact-volume RDHR cases")
    ax.grid(alpha=0.2)
    ax.legend(ncol=3, frameon=False)
    save(fig, "matrix_error_vs_walk")

    fig, ax = plt.subplots(figsize=(7.2, 3.7))
    counts_all = (2, 4, 8, 16, 32)
    for status, label, color, marker in (("0", "Validated TTC run", "#328455", "o"),
                                         ("other", "Not yet validated", "0.65", "x")):
        points = [(int(r["walk_length"]), counts_all.index(int(r["bodies"])) + 1)
                  for r in all_rows if (r["run_status"] == "0") == (status == "0")]
        ax.scatter([p[0] for p in points], [p[1] for p in points],
                   s=34, color=color, marker=marker, label=label)
    ax.set_xscale("log")
    ax.set_xlim(4, 1300000)
    ax.set_yticks(range(1, 6), [f"{count} bodies = dimension {count}"
                                  for count in counts_all])
    ax.set_xlabel("Sampler walk length")
    ax.set_title("Coverage of 100 planned dimension-matched TTC cases")
    ax.grid(axis="x", alpha=0.2)
    ax.legend(ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.26),
              frameon=False)
    save(fig, "matrix_coverage")

    xs = [float(Fraction(7 + 99 * i, 7)) for i in range(8)]
    ys = [100 / x for x in xs]
    fig, ax = plt.subplots(figsize=(5.3, 4.8))
    for x, y in zip(xs, ys):
        ax.plot([0.5, x, x], [y, y, 0.5], color="0.7", linewidth=0.8)
    boundary_x, boundary_y = [0.5, xs[0]], [ys[0], ys[0]]
    for i in range(1, 8):
        boundary_x.extend([xs[i - 1], xs[i]])
        boundary_y.extend([ys[i], ys[i]])
    ax.plot(boundary_x, boundary_y, color="#1764a1", linewidth=2.3,
            label="Union boundary (8 bodies)")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 110)
    ax.set_ylim(0.5, 110)
    ax.set_xlabel("$x_0$ radius (positive quadrant)")
    ax.set_ylabel("$x_1$ radius (positive quadrant)")
    ax.set_title("Projection of the eight-body box union")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    save(fig, "projection_m8_xy")
    counts = [count for count in (2, 4, 8, 16, 32)
              if (count, 5) in old and (count, 5) in new]
    if counts:
        fig, ax = plt.subplots(figsize=(6.5, 4.0))
        xs = list(range(len(counts)))
        ax.bar([x - 0.18 for x in xs],
               [100 * float(old[count, 5]["relative_error"]) for count in counts],
               width=0.36, color="#1764a1", label="Dimension 3")
        ax.bar([x + 0.18 for x in xs],
               [100 * float(new[count, 5]["relative_error"]) for count in counts],
               width=0.36, color="#d76027", label="Dimension = bodies")
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        ax.set_xticks(xs, [str(count) for count in counts])
        ax.set_xlabel("Number of bodies")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Same x-y boxes, exact body volumes, walk 5")
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)
        ax.legend(frameon=False)
        save(fig, "dimension_comparison_walk5")

    candidate = HERE / "fatness_d00000050" / "result.csv"
    baseline = [r for r in csv.DictReader((HERE.parent / "results.csv").open())
                if r["bodies"] == "8" and r["thin_denominator"] == "50"
                and r["walk_length"] == "1000" and r["run_status"] == "0"]
    if candidate.exists() and baseline:
        comparison = list(csv.DictReader(candidate.open()))[0]
        if comparison["run_status"] == "0":
            fig, ax = plt.subplots(figsize=(5.5, 3.9))
            ax.bar(["Dimension 3", "Dimension 8"],
                   [100 * float(baseline[0]["relative_error"]),
                    100 * float(comparison["relative_error"])],
                   color=["#1764a1", "#d76027"])
            ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                       label="20% error limit")
            ax.set_ylabel("Union volume relative error (%)")
            ax.set_title("Eight bodies, z half-width 1/50, walk 1000")
            ax.grid(axis="y", alpha=0.2)
            ax.set_axisbelow(True)
            ax.legend(frameon=False)
            save(fig, "fatness_comparison_m8_walk1000_d50")

    paired = []
    fatness_complete = []
    old_by_d = {int(r["thin_denominator"]): r
                for r in csv.DictReader((HERE.parent / "results.csv").open())
                if r["bodies"] == "8" and r["walk_length"] == "1000"
                and r["run_status"] == "0"}
    for path in sorted(HERE.glob("fatness_d*/result.csv")):
        row = list(csv.DictReader(path.open()))[0]
        denominator = int(row["thin_denominator"])
        if row["run_status"] == "0":
            fatness_complete.append((denominator, row))
            if denominator in old_by_d:
                paired.append((denominator, old_by_d[denominator], row))
    if len(fatness_complete) >= 2:
        fatness_complete.sort(key=lambda item: 1 / item[0])
        fig, ax = plt.subplots(figsize=(6.8, 4.2))
        ax.plot([1 / d for d, row in fatness_complete],
                [100 * float(row["relative_error"]) for d, row in fatness_complete],
                marker="o", linewidth=1.8, label="8 bodies, dimension 8")
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        failures = [(d, row) for d, row in fatness_complete
                    if float(row["relative_error"]) > 0.2]
        if failures:
            fail_d, fail_row = min(failures, key=lambda item: item[0])
            ax.annotate(f"fattest measured fail: 1/{fail_d}",
                        (1 / fail_d, 100 * float(fail_row["relative_error"])),
                        xytext=(12, 12), textcoords="offset points", fontsize=8,
                        arrowprops=dict(arrowstyle="->", color="0.3"))
            nearer_passes = [(d, row) for d, row in fatness_complete
                             if d < fail_d and float(row["relative_error"]) <= 0.2]
            if nearer_passes:
                pass_d, pass_row = max(nearer_passes, key=lambda item: item[0])
                ax.annotate(f"thicker pass: 1/{pass_d}",
                            (1 / pass_d, 100 * float(pass_row["relative_error"])),
                            xytext=(-105, -20), textcoords="offset points", fontsize=8,
                            arrowprops=dict(arrowstyle="->", color="0.3"))
        ax.set_xscale("log")
        ax.set_xlabel("z half-width (fatter to the right)")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Dimension-eight RDHR error versus thickness at walk 1000")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
        save(fig, "fatness_boundary_dimension8_m8_walk1000")
    if len(paired) >= 2:
        paired.sort(key=lambda item: 1 / item[0])
        fig, ax = plt.subplots(figsize=(6.6, 4.1))
        xs = [1 / item[0] for item in paired]
        ax.plot(xs, [100 * float(item[1]["relative_error"]) for item in paired],
                marker="o", linewidth=1.8, label="Dimension 3")
        ax.plot(xs, [100 * float(item[2]["relative_error"]) for item in paired],
                marker="o", linewidth=1.8, label="Dimension 8")
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        ax.set_xlabel("z half-width (fatter to the right)")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Same eight x-y boxes, exact volumes, walk 1000")
        ax.grid(alpha=0.2)
        ax.legend(frameon=False)
        save(fig, "fatness_comparison_m8_walk1000")

    matched = []
    for count, directory in ((4, "walk0010000_m0004_d00000500"),
                             (8, "walk0010000_d00000500")):
        result = HERE / directory / "result.csv"
        if not result.exists():
            break
        row = list(csv.DictReader(result.open()))[0]
        if (row["run_status"] != "0" or int(row["bodies"]) != count
                or int(row["dimension"]) != count
                or int(row["walk_length"]) != 10000
                or int(row["thin_denominator"]) != 500):
            break
        matched.append(row)
    if len(matched) == 2:
        fig, ax = plt.subplots(figsize=(5.4, 3.9))
        ax.bar(["4 bodies", "8 bodies"],
               [100 * float(row["relative_error"]) for row in matched],
               color=["#d76027", "#328455"])
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Same walk 10000 and half-width 1/500")
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)
        ax.legend(frameon=False)
        save(fig, "body_count_effect_walk10000_d500")


if __name__ == "__main__":
    main()
