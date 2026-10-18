#!/usr/bin/env python3
"""Create the 20-by-5 bounded-step ball-walk target suite.

Both categories use the same five x-y widths at every target walk. At walks
through 127 the inputs are identical except for the volume oracle. At higher
walks, the ordinary-volume suite uses a different thin denominator so volesti
can return positive body volumes; comparisons there involve different boxes.
"""

import csv
import math
from decimal import Decimal
from pathlib import Path

from generate import formula

WALKS = [1, 2, 3, 6, 11, 21, 38, 70, 127, 234,
         428, 785, 1438, 2637, 4833, 8859, 16238, 29764, 54556, 100000]
WIDTHS = [100, 120, 140, 160, 180]


def generate_category(root: Path, prefix: str) -> None:
    benchmark_dir = root / "benchmarks"
    benchmark_dir.mkdir(exist_ok=True)
    rows = []
    for target_index, walk in enumerate(WALKS):
        # The original 20*walk thinning made volesti stall or return an
        # invalid body volume before the higher ordinary-volume runs reached
        # sampling. Keep a moderate, distinct thickness at those targets.
        denominator = (1000 + target_index
                       if prefix == "VOLANDSAMP" and walk >= 234
                       else 20 * walk)
        for replica, width in enumerate(WIDTHS, 1):
            name = f"{prefix}_target{walk:06d}_r{replica:02d}_m2_n3_w{width}"
            if prefix == "VOLANDSAMP" and walk >= 234:
                name += f"_thin{denominator}"
            path = benchmark_dir / f"{name}.smt2"
            path.write_text(formula(3, width, 2, Decimal(0), denominator))
            rows.append({
                "name": name, "target_walk_length": walk, "replica": replica,
                "dimension": 3, "bodies": 2, "width": width,
                "thin_denominator": denominator,
                "exact_body_volume": format(Decimal(4 * width) / Decimal(denominator), "f"),
                "exact_union_volume": format(Decimal(8 * (width - 1)) / Decimal(denominator), "f"),
                "true_overlap_fraction": format(Decimal(1) / Decimal(width), "f"),
                "max_walk_displacement": f"{4 * walk / (denominator * math.sqrt(3)):.12f}",
                "seed": 42, "sampler": "ball", "epsilon": "0.2", "delta": "0.2",
                "volume_walk_length": 1 if prefix == "VOLANDSAMP" else "exact",
                "path": f"benchmarks/{path.name}",
            })
    with (root / "target_manifest.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print(f"Created {len(rows)} {prefix} targets at {len(WALKS)} walk lengths")


if __name__ == "__main__":
    generate_category(Path(__file__).resolve().parent, "SAMPONLY")
