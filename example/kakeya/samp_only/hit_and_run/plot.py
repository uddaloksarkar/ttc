#!/usr/bin/env python3
"""Plot only validated random-direction hit-and-run TTC measurements."""

import csv
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
FIGURES = HERE / "figures"


def save(fig, name: str) -> None:
    FIGURES.mkdir(exist_ok=True)
    fig.tight_layout()
    for extension in ("pdf", "pgf", "png"):
        fig.savefig(FIGURES / f"{name}.{extension}", dpi=220,
                    bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    rows = [row for row in csv.DictReader((HERE / "results.csv").open())
            if row["run_status"] == "0"]
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 10, "pgf.texsystem": "pdflatex"})
    by_denominator = defaultdict(list)
    for row in rows:
        assert int(row["body_rows"]) == int(row["bodies"])
        by_denominator[int(row["thin_denominator"])].append(row)

    for denominator, group in sorted(by_denominator.items()):
        if len({int(row["walk_length"]) for row in group}) < 2:
            continue
        by_count = defaultdict(list)
        for row in group:
            by_count[int(row["bodies"])].append(row)
        fig, ax = plt.subplots(figsize=(6.7, 4.3))
        for count, series in sorted(by_count.items()):
            series.sort(key=lambda row: int(row["walk_length"]))
            ax.plot([int(row["walk_length"]) for row in series],
                    [100 * float(row["relative_error"]) for row in series],
                    marker="o", linewidth=1.7, markersize=4.5,
                    label=f"{count} bodies")
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        ax.set_xscale("log")
        ax.set_xlabel("Sampler walk length")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Random-direction hit-and-run on exact-volume boxes")
        ax.grid(alpha=0.22)
        ax.legend(ncol=2, frameon=False)
        save(fig, f"union_error_vs_walk_d{denominator:08d}")

        overlap = sorted((row for row in group if int(row["bodies"]) == 2
                          and row.get("deletion_ratio")),
                         key=lambda row: int(row["walk_length"]))
        if len(overlap) >= 2:
            fig, ax = plt.subplots(figsize=(6.7, 4.3))
            walks = [int(row["walk_length"]) for row in overlap]
            observed = [100 * float(row["deletion_ratio"]) for row in overlap]
            true = 100 * float(overlap[0]["true_overlap_fraction"])
            ax.plot(walks, observed, marker="o", linewidth=1.8,
                    label="TTC deletion ratio")
            ax.axhline(true, color="0.35", linestyle="--", linewidth=1.1,
                       label="Exact intersection fraction")
            ax.set_xscale("log")
            ax.set_ylim(0, 105)
            ax.set_xlabel("Sampler walk length")
            ax.set_ylabel("Fraction of first-body samples (%)")
            ax.set_title("RDHR sample deletion versus exact body overlap")
            ax.grid(alpha=0.22)
            ax.legend(frameon=False)
            save(fig, f"overlap_mismatch_d{denominator:08d}")

    for denominator, group in sorted(by_denominator.items()):
        if denominator == 5000000:
            continue
        for walk in sorted({int(row["walk_length"]) for row in group}):
            body_rows = [row for row in group
                         if int(row["walk_length"]) == walk]
            if len(body_rows) < 2:
                continue
            body_rows.sort(key=lambda row: int(row["bodies"]))
            fig, ax = plt.subplots(figsize=(6.2, 4.0))
            xs = [int(row["bodies"]) for row in body_rows]
            ys = [100 * float(row["relative_error"]) for row in body_rows]
            ax.plot(xs, ys, marker="o", linewidth=1.8, color="#1764a1")
            ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                       label="20% error limit")
            ax.set_xticks(xs, [str(x) for x in xs])
            ax.set_xlabel("Number of bodies")
            ax.set_ylabel("Union volume relative error (%)")
            ax.set_title(f"RDHR sampling error by body count at walk {walk}")
            ax.grid(alpha=0.22)
            ax.legend(frameon=False)
            save(fig, f"body_count_effect_d{denominator:08d}_walk{walk}")

    thickness = sorted((row for row in rows
                        if int(row["bodies"]) == 8
                        and int(row["walk_length"]) == 1000
                        and int(row["thin_denominator"]) <= 100),
                       key=lambda row: 1 / int(row["thin_denominator"]))
    if len(thickness) >= 2:
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.1))
        xs = [1 / int(row["thin_denominator"]) for row in thickness]
        ys = [100 * float(row["relative_error"]) for row in thickness]
        for ax in axes:
            ax.plot(xs, ys, marker="o", linewidth=1.8, label="8 bodies")
            ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                       label="20% error limit")
            ax.set_xlabel("z half-width (fatter to the right)")
            ax.grid(alpha=0.22)
        axes[0].set_xscale("log")
        axes[0].set_ylabel("Union volume relative error (%)")
        axes[0].set_title("Full thickness sweep")
        axes[0].legend(frameon=False)
        axes[1].set_xlim(0.017, 0.0202)
        axes[1].set_ylim(18.7, 20.6)
        axes[1].set_title("Passing/failing boundary")
        by_d = {int(row["thin_denominator"]): row for row in thickness}
        if 54 in by_d and 55 in by_d:
            axes[1].annotate("1/54: pass", (1 / 54, 100 * float(by_d[54]["relative_error"])),
                             xytext=(0.0187, 19.30), fontsize=8,
                             arrowprops=dict(arrowstyle="->", color="0.3"))
            axes[1].annotate("1/55: fail", (1 / 55, 100 * float(by_d[55]["relative_error"])),
                             xytext=(0.0171, 20.35), fontsize=8,
                             arrowprops=dict(arrowstyle="->", color="0.3"))
        fig.suptitle("Random-direction hit-and-run: 8 bodies, walk 1000")
        save(fig, "fatness_boundary_m8_walk1000")

    fig, axes = plt.subplots(1, 2, figsize=(8.3, 3.9))
    for ax, limit, title in zip(axes, (105, 2.2),
                                ("Full xy projection", "Shared core enlarged")):
        ax.add_patch(Rectangle((-1, -100), 2, 200, fill=False,
                               edgecolor="#1764a1", linewidth=2,
                               label="Body 1"))
        ax.add_patch(Rectangle((-100, -1), 200, 2, fill=False,
                               edgecolor="#d76027", linewidth=2,
                               label="Body 2"))
        ax.add_patch(Rectangle((-1, -1), 2, 2, facecolor="#399d59",
                               edgecolor="none", alpha=0.45,
                               label="Intersection"))
        ax.set(xlim=(-limit, limit), ylim=(-limit, limit),
               xlabel="x", ylabel="y", title=title)
        ax.set_aspect("equal")
        ax.grid(alpha=0.18)
    axes[1].legend(loc="upper left", frameon=False, fontsize=8)
    save(fig, "two_box_xy_projection")

    target_path = HERE / "target_status.csv"
    if target_path.exists():
        targets = list(csv.DictReader(target_path.open()))
        fig, ax = plt.subplots(figsize=(7.0, 2.6))
        styles = (("0", "Completed valid TTC run", "#328455", "o"),
                  ("running", "TTC running", "#d9831f", "o"),
                  ("pending", "Pending", "0.55", "x"))
        for status, label, color, marker in styles:
            walks = [int(row["walk_length"]) for row in targets
                     if row["run_status"] == status]
            if walks:
                ax.scatter(walks, [1] * len(walks), s=58, marker=marker,
                           color=color, label=label, zorder=3)
        other = [int(row["walk_length"]) for row in targets
                 if row["run_status"] not in {"0", "running", "pending"}]
        if other:
            ax.scatter(other, [1] * len(other), s=58, marker="x",
                       color="#a02c2c", label="Interrupted or unreported",
                       zorder=3)
        ax.set_xscale("log")
        ax.set_xlim(4, 1300000)
        ax.set_ylim(0.75, 1.25)
        ax.set_yticks([])
        ax.set_xlabel("Target sampler walk length")
        ax.set_title("RDHR target coverage (one fixed two-box benchmark)")
        ax.grid(axis="x", alpha=0.22)
        ax.legend(ncol=2, loc="lower center", frameon=False, fontsize=8)
        save(fig, "target_run_coverage")

    matrix_path = HERE / "benchmark_matrix.csv"
    if matrix_path.exists():
        matrix = list(csv.DictReader(matrix_path.open()))
        walks = sorted({int(row["walk_length"]) for row in matrix})
        variants = sorted({int(row["bodies"]) for row in matrix})
        fig, ax = plt.subplots(figsize=(7.2, 4.3))
        colors = plt.cm.viridis([0.05, 0.27, 0.50, 0.73, 0.95])
        for count, color in zip(variants, colors):
            completed = sorted((row for row in matrix
                                if int(row["bodies"]) == count
                                and row["run_status"] == "0"),
                               key=lambda row: int(row["walk_length"]))
            if completed:
                ax.plot([int(row["walk_length"]) for row in completed],
                        [100 * float(row["relative_error"]) for row in completed],
                        marker="o", markersize=4, color=color, linewidth=1.7,
                        label=f"{count} bodies")
        ax.axhline(20, color="0.35", linestyle="--", linewidth=1.1,
                   label="20% error limit")
        ax.set_xscale("log")
        ax.set_xlim(4, 1300000)
        ax.set_xlabel("Sampler walk length")
        ax.set_ylabel("Union volume relative error (%)")
        ax.set_title("Exact-volume box benchmarks by body count")
        ax.grid(alpha=0.22)
        ax.legend(ncol=3, frameon=False)
        save(fig, "benchmark_matrix_union_error")

        fig, ax = plt.subplots(figsize=(7.3, 3.8))
        symbols = (("0", "Completed TTC run", "#328455", "o"),
                   ("running", "TTC running", "#d9831f", "o"),
                   ("pending", "Pending", "0.65", "x"))
        for status, label, color, marker in symbols:
            points = [(int(row["walk_length"]), variants.index(int(row["bodies"])) + 1)
                      for row in matrix if row["run_status"] == status]
            if points:
                ax.scatter([p[0] for p in points], [p[1] for p in points],
                           s=36, marker=marker, color=color, label=label)
        ax.set_xscale("log")
        ax.set_xlim(4, 1300000)
        ax.set_yticks(range(1, len(variants) + 1),
                      [f"{count} bodies" for count in variants])
        ax.set_xlabel("Sampler walk length")
        ax.set_title("Status of 100 planned hit-and-run TTC cases")
        ax.grid(axis="x", alpha=0.22)
        ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.27),
                  frameon=False, fontsize=8)
        save(fig, "benchmark_matrix_coverage")

    control = HERE / "controls" / "CDHR_control_m0002_d05000000_walk0000005.log"
    random_log = HERE / "logs" / "RDHR_box_m0002_d05000000_s0001_walk0000005.log"
    if control.exists() and random_log.exists():
        values = []
        for sampler, path in (("rdhr", random_log), ("cdhr", control)):
            raw = path.read_text()
            required = (f"c sampler walk: {sampler}", "c sample walk length: 5",
                        "c body volume oracle: exact axis-aligned box",
                        "c exit 0", "run_status=0")
            if any(item not in raw for item in required):
                raise ValueError(f"invalid hit-and-run control log: {path}")
            bodies = re.findall(r"^c\s+(\d+)\s+([0-9.eE+-]+)\s+"
                                r"(\d+)\s+(\d+)\s*$", raw, re.M)
            estimate = re.findall(r"^s volume ([0-9.eE+-]+)\s*$", raw, re.M)
            if ([int(row[0]) for row in bodies] != [1, 2]
                    or any(abs(float(row[1]) - 0.00016) > 5.1e-7 for row in bodies)
                    or len(estimate) != 1):
                raise ValueError(f"invalid exact-volume control result: {path}")
            error = 100 * abs(float(estimate[0]) / 0.0003184 - 1)
            deletion = 100 * int(bodies[1][2]) / int(bodies[0][3])
            values.append((error, deletion))
        fig, axes = plt.subplots(1, 2, figsize=(8.1, 3.7))
        labels = ["Random directions", "Coordinate directions"]
        colors = ["#d76027", "#1764a1"]
        axes[0].bar(labels, [v[0] for v in values], color=colors)
        axes[0].axhline(20, color="0.35", linestyle="--", linewidth=1.1)
        axes[0].set_ylabel("Union volume relative error (%)")
        axes[0].set_title("TTC error")
        axes[1].bar(labels, [v[1] for v in values], color=colors)
        axes[1].axhline(1, color="0.35", linestyle="--", linewidth=1.1)
        axes[1].set_ylabel("First-body samples deleted (%)")
        axes[1].set_title("Logged deletion ratio")
        for ax in axes:
            ax.tick_params(axis="x", labelrotation=15)
            ax.grid(axis="y", alpha=0.18)
            ax.set_axisbelow(True)
        fig.suptitle("Same two boxes, exact body volumes, and walk length 5")
        save(fig, "rdhr_vs_cdhr_control")


if __name__ == "__main__":
    main()
