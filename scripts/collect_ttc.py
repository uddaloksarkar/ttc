#!/usr/bin/env python3
"""Accumulate the per-run outputs of a script_ttc.sh sweep into one CSV.

Reads $SCRATCH/outfiles_ttc/manifest-<jobid>.txt to learn which out-<jobid>-<N>
directory holds which config, then parses every <bench>.out.xz / .timeout.xz
pair in those directories.

Usage:
    ./scripts/collect_ttc.py 2001463                     # dir defaults to $SCRATCH/outfiles_ttc
    ./scripts/collect_ttc.py 2001463 -d /path/to/outfiles_ttc -o sweep.csv
"""

import argparse
import csv
import lzma
import os
import re
import sys
from pathlib import Path

# "s volume 1.2345e+06"
RE_VOLUME = re.compile(r"^s volume (\S+)")
# "c      12.34    45.67% sampling"
RE_PROFILE = re.compile(r"^c\s+([\d.]+)\s+([\d.]+)%\s+(.+?)\s*$")
RE_SAMPLES_GEN = re.compile(r"^c total samples generated: (\d+)")
RE_SAMPLES_DEL = re.compile(r"^c total samples deleted: (\d+)")
RE_POLYTOPES = re.compile(r"^c polytopes: (\d+)")
RE_REALVARS = re.compile(r"^c Real variables: (\d+)")

# /usr/bin/time --verbose
# the label contains "(h:mm:ss or m:ss)" and the value contains colons too, so
# anchor on the closing paren of the label rather than on any colon
RE_WALL = re.compile(r"Elapsed \(wall clock\) time[^)]*\):\s*(\S+)")
RE_MAXRSS = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
RE_EXIT = re.compile(r"Exit status:\s*(\d+)")

# config flags out of the manifest's command line
RE_SAMPLER = re.compile(r"--sampler\s+(\S+)")
RE_WALKLEN = re.compile(r"--walklen-samp\s+(\S+)")
RE_SEED = re.compile(r"--seed\s+(\S+)")

PROFILE_KEYS = {
    "parse": "t_parse",
    "DNFization": "t_dnfization",
    "polytope volume": "t_polytope_volume",
    "sampling": "t_sampling",
    "solve": "t_solve",
}


def read_maybe_xz(path):
    """Read a file that may or may not have been compressed by the job."""
    if path.suffix == ".xz":
        with lzma.open(path, "rt", errors="replace") as f:
            return f.read()
    with open(path, "r", errors="replace") as f:
        return f.read()


def parse_out(text):
    row = {"volume": "", "exit0": 0}
    for line in text.splitlines():
        if line.startswith("s volume "):
            m = RE_VOLUME.match(line)
            if m:
                row["volume"] = m.group(1)
        elif line == "c exit 0":
            row["exit0"] = 1
        elif line.startswith("c "):
            m = RE_PROFILE.match(line)
            if m and m.group(3) in PROFILE_KEYS:
                row[PROFILE_KEYS[m.group(3)]] = m.group(1)
                continue
            for regex, key in (
                (RE_SAMPLES_GEN, "samples_generated"),
                (RE_SAMPLES_DEL, "samples_deleted"),
                (RE_POLYTOPES, "polytopes"),
                (RE_REALVARS, "real_vars"),
            ):
                m = regex.match(line)
                if m:
                    row[key] = m.group(1)
                    break
    return row


def hms_to_seconds(s):
    """'1:02:03.45' or '2:03.45' -> seconds."""
    parts = s.split(":")
    try:
        secs = float(parts[-1])
        if len(parts) > 1:
            secs += 60.0 * int(parts[-2])
        if len(parts) > 2:
            secs += 3600.0 * int(parts[-3])
    except ValueError:
        return ""
    return f"{secs:.2f}"


def parse_timeout(text):
    row = {}
    m = RE_WALL.search(text)
    if m:
        row["wall_s"] = hms_to_seconds(m.group(1))
    m = RE_MAXRSS.search(text)
    if m:
        row["maxrss_kb"] = m.group(1)
    m = RE_EXIT.search(text)
    if m:
        row["exit_status"] = m.group(1)
    return row


def load_manifest(path):
    """dir name -> (raw config string, sampler, walklen, seed)."""
    configs = {}
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            dirname, _, opts = line.partition("\t")
            entry = {"config": opts.strip()}
            for key, regex in (("sampler", RE_SAMPLER),
                               ("walklen", RE_WALKLEN),
                               ("seed", RE_SEED)):
                m = regex.search(opts)
                entry[key] = m.group(1) if m else ""
            configs[dirname] = entry
    return configs


FIELDS = [
    "config_dir",
    "sampler",
    "walklen",
    "seed",
    "benchmark",
    "solved",
    "volume",
    "t_sampling",
    "t_polytope_volume",
    "t_parse",
    "t_dnfization",
    "t_solve",
    "wall_s",
    "maxrss_kb",
    "exit_status",
    "samples_generated",
    "samples_deleted",
    "polytopes",
    "real_vars",
    "config",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("jobid", help="SLURM job id, e.g. 2001463")
    ap.add_argument("-d", "--dir",
                    default=os.path.join(os.environ.get("SCRATCH", "."), "outfiles_ttc"),
                    help="directory holding manifest-<jobid>.txt and out-<jobid>-* "
                         "(default: $SCRATCH/outfiles_ttc)")
    ap.add_argument("-o", "--out", default=None,
                    help="output CSV (default: results-<jobid>.csv)")
    args = ap.parse_args()

    root = Path(args.dir)
    manifest = root / f"manifest-{args.jobid}.txt"
    if not manifest.is_file():
        sys.exit(f"no manifest at {manifest}")
    configs = load_manifest(manifest)

    out_path = Path(args.out or f"results-{args.jobid}.csv")
    nrows = 0
    missing_dirs = []

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()

        for dirname in sorted(configs, key=lambda d: int(d.rsplit("-", 1)[1])):
            cfg = configs[dirname]
            cdir = root / dirname
            if not cdir.is_dir():
                missing_dirs.append(dirname)
                continue

            # .out may be .out.xz (normal) or bare .out (if xz was skipped)
            outs = sorted(list(cdir.glob("*.out.xz")) + list(cdir.glob("*.out")))
            for outfile in outs:
                bench = outfile.name
                for suffix in (".out.xz", ".out"):
                    if bench.endswith(suffix):
                        bench = bench[: -len(suffix)]
                        break

                row = {
                    "config_dir": dirname,
                    "sampler": cfg["sampler"],
                    "walklen": cfg["walklen"],
                    "seed": cfg["seed"],
                    "benchmark": bench,
                    "config": cfg["config"],
                }
                try:
                    row.update(parse_out(read_maybe_xz(outfile)))
                except (lzma.LZMAError, OSError) as e:
                    print(f"warn: unreadable {outfile}: {e}", file=sys.stderr)

                for cand in (cdir / f"{bench}.timeout.xz", cdir / f"{bench}.timeout"):
                    if cand.exists():
                        try:
                            row.update(parse_timeout(read_maybe_xz(cand)))
                        except (lzma.LZMAError, OSError) as e:
                            print(f"warn: unreadable {cand}: {e}", file=sys.stderr)
                        break

                row["solved"] = row.get("exit0", 0)
                writer.writerow(row)
                nrows += 1

    print(f"wrote {nrows} rows from {len(configs) - len(missing_dirs)} configs -> {out_path}")
    if missing_dirs:
        print(f"note: {len(missing_dirs)} config dirs from the manifest are absent: "
              f"{', '.join(missing_dirs[:5])}"
              f"{' …' if len(missing_dirs) > 5 else ''}", file=sys.stderr)


if __name__ == "__main__":
    main()
