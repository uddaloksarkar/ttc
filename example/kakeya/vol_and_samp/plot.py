#!/usr/bin/env python3
"""Draw body-count-only plots, separated by the fixed volume-walk setting."""

import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
OUT = HERE / "figures"
plt.rcParams.update({"font.size": 11, "pgf.texsystem": "pdflatex"})


def save(fig, stem: str) -> None:
    for suffix in ("pdf", "pgf", "png"):
        fig.savefig(OUT / f"{stem}.{suffix}", bbox_inches="tight", dpi=220)
    plt.close(fig)


def main() -> None:
    rows = list(csv.DictReader((HERE / "results.csv").open(newline="")))
    by_vol = defaultdict(list)
    for row in rows:
        by_vol[int(row["volume_walk_length"])].append(row)
    for vol, subset in sorted(by_vol.items()):
        by_count = defaultdict(list)
        for row in subset:
            by_count[int(row["bodies"])].append(row)
        walks = sorted({int(row["sample_walk_length"]) for row in subset})
        plots = [
            ("union_error", "Union relative error (%)",
             lambda r: 100*float(r["union_relative_error"]), 20, "TTC tolerance"),
            ("overlap_mismatch", "Apparent minus exact overlap (percentage points)",
             lambda r: 100*(float(r["apparent_overlap"])-float(r["true_overlap"])),
             0, "Exact overlap"),
            ("body_error", "Largest body relative error (%)",
             lambda r: 100*float(r["max_body_relative_error"]),
             100*.2/12, "Body budget"),
        ]
        for stem, ylabel, value, reference, ref_label in plots:
            fig, ax = plt.subplots(figsize=(6.6, 4.4))
            for m, values in sorted(by_count.items()):
                values.sort(key=lambda r: int(r["sample_walk_length"]))
                ax.plot([int(r["sample_walk_length"]) for r in values],
                        [value(r) for r in values], "o-", lw=1.8, label=f"{m} bodies")
            ax.axhline(reference, color="firebrick", ls="--", label=ref_label)
            ax.set(xlabel="Sampling walk length", ylabel=ylabel,
                   title=f"Fixed volume walk length {vol}")
            ax.set_xscale("log")
            ax.set_xticks(walks, labels=[str(w) for w in walks])
            ax.grid(alpha=.25)
            ax.legend(ncol=2, fontsize=9)
            save(fig, f"VOLANDSAMP_v{vol}_{stem}_vs_walklen")


if __name__ == "__main__":
    main()
