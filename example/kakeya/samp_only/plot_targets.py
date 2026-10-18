#!/usr/bin/env python3
"""Plot measured TTC target failures, clearly marking unrun targets."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
plt.rcParams.update({"font.size": 10, "pgf.texsystem": "pdflatex"})


def save(fig, root: Path, stem: str) -> None:
    out = root / "figures"
    out.mkdir(exist_ok=True)
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(out / f"{stem}.{extension}", dpi=220, bbox_inches="tight")
    plt.close(fig)


def plot_category(root: Path) -> None:
    rows = list(csv.DictReader((root / "target_results.csv").open(newline="")))
    complete = [r for r in rows if r["run_status"] == "0"]
    targets = sorted({int(r["target_walk_length"]) for r in rows})
    by_target = defaultdict(list)
    for row in complete:
        by_target[int(row["target_walk_length"])].append(row)

    fig, ax = plt.subplots(figsize=(8.2, 4.5))
    for i, (target, group) in enumerate(sorted(by_target.items())):
        ax.scatter([target] * len(group),
                   [100 * float(r["union_relative_error"]) for r in group],
                   s=32, color="#2166ac", alpha=.75,
                   label="2 bodies" if i == 0 else None)
    ax.axhline(20, color="#b2182b", linestyle="--", label="TTC tolerance (20%)")
    ax.set(xlabel="Target sampler walk length", ylabel="Measured union relative error (%)",
           xlim=(.8, 1.3e5), ylim=(0, 105))
    ax.set_xscale("log")
    ax.grid(alpha=.25)
    ax.legend()
    save(fig, root, "target_union_error_measured")

    fig, ax = plt.subplots(figsize=(8.2, 3.9))
    counts = [len(by_target[t]) for t in targets]
    ax.scatter(targets, [5] * len(targets), marker="_", s=160,
               color="#bdbdbd", label="Generated: 5 two-body cases")
    ax.plot(targets, counts, "o-", color="#2166ac", label="Completed TTC runs")
    ax.set(xlabel="Target sampler walk length", ylabel="Cases per target",
           xlim=(.8, 1.3e5), ylim=(-.4, 5.8), yticks=range(6))
    ax.set_xscale("log")
    ax.grid(alpha=.25)
    ax.legend()
    save(fig, root, "target_run_coverage")

    fig, (ax_top, ax_bottom) = plt.subplots(2, 1, sharex=True,
                                            figsize=(8.2, 5.2),
                                            gridspec_kw={"height_ratios": [2, 1]})
    for target, group in sorted(by_target.items()):
        ax_top.scatter([target] * len(group),
                       [100 * float(r["sampler_deleted_fraction"]) for r in group],
                       s=26, color="#2166ac", alpha=.7)
        ax_bottom.scatter([target] * len(group),
                          [200 / float(r["width"]) for r in group],
                          s=26, color="#b2182b", alpha=.7)
    ax_top.set(ylabel="TTC deleted samples (%)", ylim=(-4, 105))
    ax_bottom.set(xlabel="Target sampler walk length",
                  ylabel="True body\nintersection (%)", ylim=(0, 2.4),
                  xlim=(.8, 1.3e5))
    ax_bottom.set_xscale("log")
    for ax in (ax_top, ax_bottom):
        ax.grid(alpha=.25)
    save(fig, root, "target_overlap_mismatch_measured")

    fig, ax = plt.subplots(figsize=(6.0, 5.0))
    widths = sorted({int(r["width"]) for r in rows}, reverse=True)
    for width in widths:
        half = width / 2
        ax.add_patch(Rectangle((-half, -1), width, 2, fill=False,
                               edgecolor="#2166ac", linewidth=1.3, alpha=.65))
        ax.add_patch(Rectangle((-1, -half), 2, width, fill=False,
                               edgecolor="#ef8a62", linewidth=1.3, alpha=.65))
    ax.add_patch(Rectangle((-1, -1), 2, 2, fill=False,
                           edgecolor="#b2182b", linewidth=2,
                           label="Shared overlap"))
    ax.plot([], [], color="#2166ac", label="Body 1 projections")
    ax.plot([], [], color="#ef8a62", label="Body 2 projections")
    limit = max(widths) / 2 + 3
    ax.set(xlabel="$x_0$", ylabel="$x_1$", xlim=(-limit, limit),
           ylim=(-limit, limit), title="Projection of all target boxes")
    ax.set_aspect("equal")
    ax.grid(alpha=.2)
    ax.legend()
    save(fig, root, "target_projection_x0_x1")


if __name__ == "__main__":
    plot_category(HERE)
