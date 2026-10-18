#!/usr/bin/env python3
"""Join all target geometries with completed TTC measurements and run status."""

import csv
import re
from pathlib import Path

HERE = Path(__file__).resolve().parent


def analyze_category(root: Path) -> None:
    cases = list(csv.DictReader((root / "target_manifest.csv").open(newline="")))
    rows = []
    for case in cases:
        row = dict(case)
        log = root / "benchmarks" / f"{case['name']}.log"
        row.update(run_status="unrun", estimated_body_sum="",
                   estimated_union_volume="", union_relative_error="",
                   apparent_overlap_fraction="", sampler_deleted_fraction="",
                   max_body_relative_error="", wall_seconds="")
        if log.exists():
            raw = log.read_text()
            status = re.search(r"^run_status=(\S+)$", raw, re.M)
            row["run_status"] = status.group(1) if status else "incomplete"
            elapsed = re.search(r"^wall_seconds=([\d.]+)$", raw, re.M)
            if elapsed:
                row["wall_seconds"] = elapsed.group(1)
            if row["run_status"] == "0":
                result = re.search(r"^s volume ([\d.eE+-]+)$", raw, re.M)
                volumes = [float(x) for x in re.findall(
                    r"^c\s+\d+\s+([\d.eE+-]+)\s+\d+\s+\d+$", raw, re.M)]
                deletion = re.findall(
                    r"^c\s+\d+\s+[\d.eE+-]+\s+(\d+)\s+(\d+)$", raw, re.M)
                if not result:
                    row["run_status"] = "invalid_no_result"
                    rows.append(row)
                    continue
                if len(volumes) != int(case["bodies"]) or len(deletion) != int(case["bodies"]):
                    row["run_status"] = "invalid_body_rows"
                    rows.append(row)
                    continue
                expected_oracle = ("exact axis-aligned box" if case["volume_walk_length"] == "exact"
                                   else "volesti estimate")
                if f"c body volume oracle: {expected_oracle}" not in raw:
                    raise ValueError(f"Wrong body volume oracle: {log}")
                if f"c sample walk length: {case['target_walk_length']}" not in raw:
                    raise ValueError(f"Wrong sample walk length: {log}")
                if any(int(total) == 0 for _, total in deletion):
                    row["run_status"] = "invalid_zero_sample"
                    rows.append(row)
                    continue
                returned = float(result.group(1))
                if returned <= 0:
                    row["run_status"] = "invalid_nonpositive_result"
                    rows.append(row)
                    continue
                exact_union = float(case["exact_union_volume"])
                body_sum = sum(volumes)
                exact_body = float(case["exact_body_volume"])
                row["estimated_body_sum"] = body_sum
                row["estimated_union_volume"] = returned
                row["union_relative_error"] = abs(returned / exact_union - 1)
                # TTC currently prints body rows at six fixed decimals. Very
                # small positive estimates can therefore appear as 0.000000.
                if body_sum > 0:
                    row["apparent_overlap_fraction"] = 1 - returned / body_sum
                row["sampler_deleted_fraction"] = int(deletion[1][0]) / int(deletion[0][1])
                if body_sum > 0:
                    row["max_body_relative_error"] = max(abs(v / exact_body - 1) for v in volumes)
        rows.append(row)
    with (root / "target_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    complete = [r for r in rows if r["run_status"] == "0"]
    failed = sum(float(r["union_relative_error"]) > .2 for r in complete)
    print(f"{root.name}: {len(rows)} targets, {len(complete)} complete, "
          f"{failed} measured TTC union failures")


if __name__ == "__main__":
    analyze_category(HERE)
