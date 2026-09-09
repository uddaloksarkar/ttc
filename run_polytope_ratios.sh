#!/bin/bash

# Standalone follow-up to script_ttc.sh: compute the R/r fatness ratios
# (polytope_ratios.py) for every benchmark in a directory, independent of the
# sampler/walklen/seed sweep. Kept out of script_ttc.sh's hot path because
# `ttc --dump-ine` still runs the full volume computation rather than a cheap
# LP-only dump, so this is roughly as expensive as one more full sweep config
# -- run it separately, not on the critical path of the MPI job.
#
# Each benchmark's invocation is capped with doalarm (the same tool
# script_ttc.sh uses for its sampler runs), and written to its own CSV, so one
# hung or zero-volume benchmark (e.g. an infeasible polytope, which
# polytope_ratios.py hard-exits on) never takes the rest of the batch down
# with it. Re-running the script skips benchmarks already done, so it is safe
# to resume after an interruption.
#
# Usage:
#   ./run_polytope_ratios.sh [filespos] [output_dir]
#
# Env overrides:
#   TTC_BIN   ttc binary (default ./ttc)
#   TLIMIT    per-benchmark wall-clock cap in seconds, enforced via doalarm
#             (default 300)
#   DOALARM   path to the doalarm binary (default: first existing one of
#             ./doalarm, bins/doalarm, $SLURM_SUBMIT_DIR/bins/doalarm)
#   NFILES    if >0, cap the number of benchmarks processed (default 0 = all)

set -u

filespos="${1:-syntheticLRA}"
outdir="${2:-polytope_ratios_out}"
ttc_bin="${TTC_BIN:-./ttc}"
tlimit="${TLIMIT:-300}"
nfiles="${NFILES:-0}"

doalarm="${DOALARM:-}"
if [[ -z "$doalarm" ]]; then
    for cand in ./doalarm bins/doalarm "${SLURM_SUBMIT_DIR:-}/bins/doalarm"; do
        if [[ -n "$cand" && -x "$cand" ]]; then
            doalarm="$cand"
            break
        fi
    done
fi
if [[ -z "$doalarm" ]]; then
    echo "error: doalarm binary not found (looked in ./doalarm, bins/doalarm," \
         "\$SLURM_SUBMIT_DIR/bins/doalarm); set DOALARM=<path> to override" >&2
    exit 1
fi
if [[ ! -x "$ttc_bin" ]]; then
    echo "error: ttc binary not found or not executable: $ttc_bin" >&2
    exit 1
fi
if [[ ! -d "$filespos" ]]; then
    echo "error: benchmark directory not found: $filespos" >&2
    exit 1
fi

mkdir -p "${outdir}/per_file" "${outdir}/logs"

shopt -s nullglob
files=("${filespos}"/*.smt2.xz)
if [[ ${#files[@]} -eq 0 ]]; then
    echo "error: no *.smt2.xz files in ${filespos}" >&2
    exit 1
fi
if [[ ${nfiles} -gt 0 && ${nfiles} -lt ${#files[@]} ]]; then
    files=("${files[@]:0:${nfiles}}")
fi

echo "processing ${#files[@]} benchmark(s) from ${filespos}, tlimit=${tlimit}s per file"

skipped=()
for f in "${files[@]}"; do
    name=$(basename "$f")
    csv="${outdir}/per_file/${name%.smt2.xz}.csv"
    log="${outdir}/logs/${name%.smt2.xz}.log"
    if [[ -s "$csv" ]]; then
        continue   # already done -- lets the script be safely re-run to resume
    fi
    "$doalarm" -t real "$tlimit" \
        ./polytope_ratios.py "$f" --ttc_bin "$ttc_bin" --csv "$csv" --quiet \
        > "$log" 2>&1
    rc=$?
    if [[ $rc -ne 0 ]]; then
        echo "  [skip] ${name} (exit ${rc}, see ${log})"
        skipped+=("$name")
        rm -f "$csv"   # a failed run must not look like a completed one
    fi
done

# merge every per-file CSV into one, keeping a single header
combined="${outdir}/polytope_ratios_all.csv"
: > "$combined"
header_written=0
for csv in "${outdir}"/per_file/*.csv; do
    [[ -e "$csv" ]] || continue
    if [[ $header_written -eq 0 ]]; then
        cat "$csv" >> "$combined"
        header_written=1
    else
        tail -n +2 "$csv" >> "$combined"
    fi
done

echo "done: $((${#files[@]} - ${#skipped[@]}))/${#files[@]} benchmarks written -> ${combined}"
if [[ ${#skipped[@]} -gt 0 ]]; then
    printf '  skipped: %s\n' "${skipped[@]}"
fi
