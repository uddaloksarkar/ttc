#!/usr/bin/env python3
"""Plot validated factored 8/64/100-body pilots and their x-y projections."""

from fractions import Fraction
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from prepare_factored import COUNTS, HERE
from run_factored import parse

FIGURES = HERE / "figures"


def save(fig, stem: str) -> None:
    FIGURES.mkdir(exist_ok=True)
    fig.tight_layout()
    for suffix in ("pdf", "pgf", "png"):
        fig.savefig(FIGURES / f"{stem}.{suffix}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    plt.rcParams.update({"font.size": 10, "pgf.texsystem": "pdflatex"})
    rows = {count: parse(count, 5) for count in COUNTS}
    validated = [count for count in COUNTS if rows[count]["run_status"] == "0"]

    fig, ax = plt.subplots(figsize=(6.6, 4.1))
    if validated:
        ax.plot(validated,
                [100 * float(rows[count]["relative_error"]) for count in validated],
                marker="o", linewidth=1.8, color="#1764a1", label="Validated TTC")
    ax.axhline(20, linestyle="--", color="0.35", label="20% error limit")
    ax.set_xscale("log")
    ax.set_xlim(7, 115)
    ax.set_xticks(COUNTS, [str(count) for count in COUNTS])
    ax.set_xlabel("Number of bodies (= dimension)")
    ax.set_ylabel("Union volume relative error (%)")
    ax.set_title("Factored exact-volume RDHR pilots at walk 5")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    save(fig, "error_vs_body_count_walk5")

    fig, ax = plt.subplots(figsize=(6.6, 2.8))
    for count in COUNTS:
        complete = count in validated
        ax.scatter(count, 1, s=130, marker="o" if complete else "x",
                   color="#328455" if complete else "0.55")
        ax.annotate("validated" if complete else "not validated", (count, 1),
                    xytext=(0, 12), textcoords="offset points", ha="center")
    ax.set_xscale("log")
    ax.set_xlim(7, 115)
    ax.set_ylim(0.75, 1.45)
    ax.set_xticks(COUNTS, [str(count) for count in COUNTS])
    ax.set_yticks([])
    ax.set_xlabel("Number of bodies (= dimension)")
    ax.set_title("TTC coverage of the same walk-5 pilots")
    save(fig, "coverage_walk5")

    fig, ax = plt.subplots(figsize=(6.2, 4.8))
    for count in COUNTS:
        xs = [float(Fraction(count - 1 + 99 * i, count - 1))
              for i in range(count)]
        ys = [100 / x for x in xs]
        boundary_x = [0.5, xs[0]]
        boundary_y = [ys[0], ys[0]]
        for i in range(1, count):
            boundary_x.extend((xs[i - 1], xs[i]))
            boundary_y.extend((ys[i], ys[i]))
        ax.plot(boundary_x, boundary_y, linewidth=1.7,
                label=f"{count} bodies")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.5, 110)
    ax.set_ylim(0.5, 110)
    ax.set_xlabel("x0 radius (positive quadrant)")
    ax.set_ylabel("x1 radius (positive quadrant)")
    ax.set_title("Projection of the factored box unions")
    ax.grid(alpha=0.2)
    ax.legend(frameon=False)
    save(fig, "xy_projection")


if __name__ == "__main__":
    main()
