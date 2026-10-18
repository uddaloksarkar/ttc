#!/usr/bin/env python3
"""Plot every validated run for every failing benchmark in each sampler folder."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent


def draw(sampler: str) -> None:
    folder = HERE / sampler
    table = folder / "combined_results.csv"
    if not table.exists():
        table = folder / "measurements.csv"
    rows = [row for row in csv.DictReader(table.open())
            if row.get("run_status", "0") == "0" and row.get("relative_error")]
    by_input = defaultdict(list)
    for row in rows:
        by_input[row["benchmark"]].append(row)
    groups = sorted({(int(row["bodies"]), int(row["dimension"])) for row in rows})
    palette = plt.get_cmap("tab20", max(len(groups), 1))
    colors = {group: palette(i) for i, group in enumerate(groups)}
    fig, ax = plt.subplots(figsize=(11, 7.2))
    for benchmark, cases in sorted(by_input.items()):
        cases.sort(key=lambda row: int(row["walk_length"]))
        group = (int(cases[0]["bodies"]), int(cases[0]["dimension"]))
        color = colors[group]
        walks = [int(row["walk_length"]) for row in cases]
        errors = [100 * float(row["relative_error"]) for row in cases]
        ax.plot(walks, errors, color=color, alpha=0.35, linewidth=1.1)
        ax.scatter(walks, errors, color=[color], alpha=0.65, s=14)
    ax.axhline(20, color="black", linestyle="--", linewidth=1,
               label="20% failure threshold")
    ax.set_xscale("log")
    ax.set_xlabel("Sampler walk length")
    ax.set_ylabel("Absolute relative union-volume error (%)")
    ax.set_title(f"{sampler.replace('_', ' ').title()}: all validated failing-benchmark measurements")
    ax.grid(True, which="both", alpha=0.2)
    handles = [Line2D([0], [0], color=colors[group], marker="o", markersize=4,
                      label=f"{group[0]} bodies, {group[1]} dimensions")
               for group in groups]
    handles.append(Line2D([0], [0], color="black", linestyle="--",
                          label="20% failure threshold"))
    fig.legend(handles=handles, fontsize=8, ncol=3, loc="lower center",
               bbox_to_anchor=(0.5, 0.015), frameon=True)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.92, bottom=0.24)
    output = folder / "figures" / "walk_length_vs_relative_error"
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(output.with_suffix(f".{extension}"), bbox_inches="tight", dpi=220)
    plt.close(fig)
    print(f"{sampler}: plotted {len(rows)} measurements from {len(by_input)} benchmarks")


if __name__ == "__main__":
    for name in ("ball_walk", "hit_and_run"):
        draw(name)
