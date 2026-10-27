#!/bin/bash

# Run a Kakeya campaign from the repository root (hit_and_run by default).
# This block is self-contained
# so the per-campaign creator.sh files are not needed on the cluster.
case "${1:-hit_and_run}" in
    ball_walk|hit_and_run|hit_and_run_radial)
        requested_campaign="${1:-hit_and_run}"
        campaign="$requested_campaign"
        radial_only=false
        if [[ "$requested_campaign" == hit_and_run_radial ]]; then
            campaign=hit_and_run
            radial_only=true
        fi
        ulimit -t unlimited
        project_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
        campaign_dir="${project_dir}/example/kakeya/samp_only/${campaign}"
        # /project is read-only on compute nodes, so results go to a writable
        # directory under the ttc root (project_dir keeps the /scratch path).
        job_id=${SLURM_JOB_ID:-local}
        out_dir="${OUT_DIR:-${project_dir}/out_${campaign}}_${job_id}"
        sampler=ball
        [[ "$campaign" == hit_and_run ]] && sampler=rdhr
        walklens=(5 10 18 34 65 124 236 449 853 1622 3083 5861 11143 21183 40271 76559 145544 276692 526015 1000000)
        ttc_bin=${TTC_BIN:-${project_dir}/ttc}
        tlimit=${TLIMIT_SECONDS:-1209600}
        rank=${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-0}}
        numthreads=${OMPI_COMM_WORLD_SIZE:-${SLURM_NTASKS:-1}}
        workdir="${SCRATCH:-/tmp}/scratch/${job_id}_${rank}"

        [[ -x "$ttc_bin" ]] || { echo "TTC executable not found: $ttc_bin" >&2; exit 1; }
        mkdir -p "$workdir" "$out_dir"
        cd "$workdir" || exit 1
        cp -a "${campaign_dir}/benchmarks" .
        if [[ "$campaign" == hit_and_run && "$radial_only" == false ]]; then
            # The catalog names the outputs; stage selected originals with those
            # names only in scratch, without keeping duplicate repository files.
            while IFS=, read -r benchmark bodies dimension sha original_path; do
                [[ "$benchmark" == benchmark ]] && continue
                original_path=${original_path%$'\r'}
                source="${project_dir}/example/kakeya/samp_only/${original_path}"
                [[ -f "$source" ]] || { echo "Benchmark not found: $source" >&2; exit 1; }
                actual=$(sha256sum "$source")
                [[ "${actual%% *}" == "$sha" ]] || {
                    echo "Benchmark differs from catalog: $source" >&2; exit 1;
                }
                cp "$source" "benchmarks/$benchmark"
            done < "${campaign_dir}/catalog.csv"
        fi
        shopt -s nullglob
        rm -f todo
        touch todo

        numlines=0
        if [[ "$campaign" == hit_and_run ]]; then
            help_text=$("$ttc_bin" --help 2>&1)
            [[ "$help_text" == *"exact-parallelotope-volumes"* ]] || {
                echo "TTC binary lacks --exact-parallelotope-volumes; rebuild it" >&2
                exit 1
            }
        fi
        for file in benchmarks/*.smt2; do
            benchmark=${file#benchmarks/}
            if [[ "$benchmark" == RDHR_radial_urchin_*.smt2 ]]; then
                oracle=--exact-parallelotope-volumes
            else
                [[ "$radial_only" == false ]] || continue
                # Only selected catalog inputs are included alongside radial inputs.
                [[ "$benchmark" =~ _[0-9a-f]{12}\.smt2$ ]] || continue
                oracle=--exact-box-volumes
            fi
            for walk in "${walklens[@]}"; do
                base="${benchmark%.smt2}_walk$(printf '%07d' "$walk")"
                printf 'if [[ ! -f %q ]]; then /usr/bin/time --verbose -o %q timeout -k 2 %q %q %q --seed 42 --epsilon 0.2 --delta 0.2 --sampler %q --walklen-samp %q %q > %q 2>&1; status=$?; printf "\\nrun_status=%%s\\n" "$status" >> %q; xz -f %q %q; cp %q %q %q; fi\n' \
                    "${out_dir}/${base}.out.xz" "${base}.timeout" "$tlimit" \
                    "$ttc_bin" "$file" "$sampler" "$walk" "$oracle" \
                    "${base}.out" "${base}.out" "${base}.out" "${base}.timeout" \
                    "${base}.out.xz" "${base}.timeout.xz" "$out_dir" >> todo
                numlines=$((numlines+1))
            done
        done

        if [[ "${DRY_RUN:-0}" == 1 ]]; then
            echo "Prepared ${numlines} cases for ${requested_campaign}; first command:"
            head -n 1 todo
            exit 0
        fi

        numper=$(((numlines + numthreads - 1) / numthreads))
        for ((i=0; i<numthreads; i++)); do
            start=$((i * numper + 1))
            sed -n "${start},$((start + numper - 1))p" todo > "todo_${i}.sh"
            sed -i '1i#!/bin/bash' "todo_${i}.sh"
            chmod +x "todo_${i}.sh"
        done
        echo "This is MPI exec number ${rank}; ${numlines} ${requested_campaign} cases across ${numthreads} ranks"
        ./"todo_${rank}.sh" > "rank_${rank}.log" 2>&1
        echo "Finished waiting rank ${rank}"
        exit 0
        ;;
    legacy) shift ;;
    *) echo "Usage: $0 [hit_and_run|hit_and_run_radial|ball_walk|legacy]" >&2; exit 2 ;;
esac

# Legacy syntheticLRA sweep; reached only with the explicit "legacy" argument.
# source $SCRATCH/ttc-venv/bin/activate

# filespos="lra_benchmarks"
# filespos="bounded-set-copy"
# filespos="DAC25_benchmarks"
filespos="syntheticLRA"

ulimit -t unlimited
shopt -s nullglob
rm -f todo
touch todo

# Sampling sweep for the LRA volume engine: sampler x walk length x seed.
# Grid design, from measured behaviour on a dim-15 instance:
#  * volume_cooling_gaussians is ~90% of each run and is identical across every
#    config here, so walk length is nearly free (billiard wl=1000 costs 3.2s
#    total vs 1.9s at wl=10).  No reason to trim the top of the range.
#  * Single-seed cells are noise-dominated: the same instance moved +-11%
#    between adjacent walk lengths, non-monotonically.  Replicating over seeds
#    is what makes the walk-length axis readable, so the budget goes there
#    rather than into more walk lengths.  The volume estimator's RNG is
#    seed-independent (fixed compile-time seed), so --seed moves only the
#    sampling/union half -- exactly the half under study.
#  * Walk lengths are log-uniform (ratio 10) so each point carries equal
#    information on a log-x plot; wl=1 anchors "does mixing matter at all".
#  * Equal --walklen-samp is NOT equal work: per step, billiard ~4.7x rdhr
#    ~34x cdhr.  Plot quality against measured sampling time, not against wl.
#  * --walklen-samp only sizes the sampling walk; --walklen-vol is left at its
#    default so the two walks stay separate axes.
# ttc is the node-local copy made below -- do NOT run it from $SLURM_SUBMIT_DIR.
samplers=(billiard ball rdhr cdhr)
walklens=(1 10 100 1000)
seeds=(41 42 43)

opts_arr=()
for smp in "${samplers[@]}"
do
    for wl in "${walklens[@]}"
    do
        for sd in "${seeds[@]}"
        do
            opts_arr+=("./ttc -v 2 --seed ${sd} --sampler ${smp} --walklen-samp ${wl} ")
        done
    done
done
# Output dir index is
#   at_opt = ((sampler idx)*${#walklens[@]} + (walklen idx))*${#seeds[@]} + seed idx
# so seeds vary fastest and each block of ${#seeds[@]} consecutive dirs is one
# (sampler, walklen) cell.  Read the exact mapping from
# ${outputdir}/manifest-${SLURM_JOB_ID}.txt rather than recomputing it.

# Reference points -- append (never insert) so existing indices do not shift.
# opts_arr+=("./ttc -v 2 --seed 42 --nogmp ")     # all-double baseline
# opts_arr+=("./ttc -v 2 --seed 42 --fullgmp ")   # walk itself in GMP
# opts_arr+=("./ttc -v 2 --seed 42 --no-cdd-simp ")
# opts_arr+=("./cvc5 -S -m --incremental --bv-sat-solver=cryptominisat --hashsm=xor  ")
# opts_arr+=("./sharpSMT-CG -p ")
# opts_arr+=("./sharp_smt ")
# opts_arr+=("./sharpSMT -l --ge 1 --bunch 1 --verb 2 --scale 10 --canon 1 ")

output="out"
tlimit="3600"
#tlimit="3600"
#6GB mem limitq
memlimit="6000000"
numthreads=$((OMPI_COMM_WORLD_SIZE))



SERVER=$SLURM_SUBMIT_HOST
WORKDIR="$SCRATCH/scratch/${SLURM_JOB_ID}_${OMPI_COMM_WORLD_RANK}"
sleep 0
output="${output}-${SLURM_JOB_ID}"

# May comment out the below echo commands later
echo ------------------------------------------------------
echo "Job is running on node ${SLURM_JOB_NODELIST}"
echo ------------------------------------------------------
echo "Rank is: ${OMPI_COMM_WORLD_RANK}"
echo "SLURM: sbatch is running on $SLURM_SUBMIT_HOST"
echo "SLURM: working directory is $SLURM_SUBMIT_DIR"
echo "SLURM: job identifier is ${SLURM_JOB_ID}"
echo "SLURM: job name is $SLURM_JOB_NAME"
echo "SLURM: node file is $SLURM_JOB_NODELIST"
echo "SLURM: current home directory is $HOME"
echo "SLURM: PATH = $SLURM_O_PATH"
echo "server      is ${SERVER}"
echo "workdir     is ${WORKDIR}"
echo "servpermdir is ${SERVPERMDIR}"
echo "Output dir  is ${output}"


# echo "Transferring files from server to compute node"
mkdir -p "${WORKDIR}"
cd "${WORKDIR}" || exit

files=$(find ${SLURM_SUBMIT_DIR}/${filespos}/ -maxdepth 1 -name "*.smt2.xz" | shuf  --random-source=${SLURM_SUBMIT_DIR}/bins/myrnd)
# Benchmark subset for the sweep.  The shuf above is seeded from bins/myrnd, so
# the first N files are a fixed random sample -- the same N across every config
# and across reruns.  Set to 0 to use the whole set.  With 80 configs, 1131
# files is ~90k runs; 300 files keeps the sweep near the budget of a 20-config
# full-set run.  Confirm the winning config on the full set afterwards.
nfiles=300
if [[ ${nfiles} -gt 0 ]]; then
    files=$(echo "${files}" | head -n ${nfiles})
fi
# files=$(ls ${SLURM_SUBMIT_DIR}/${filespos}/*.smt2.xz | shuf --random-source=${SLURM_SUBMIT_DIR}/bins/myrnd)
outputdir="${SCRATCH}/outfiles_ttc"
mkdir ./bin
cp ${SLURM_SUBMIT_DIR}/bins/doalarm .
cp -L ${SLURM_SUBMIT_DIR}/ttc .

# create todo
rm -f todo
mkdir -p ${output}
at_opt=0
numlines=0
# Record which config each output dir holds; rank 0 only, so the file is written
# once per job rather than once per MPI rank.
if [[ ${OMPI_COMM_WORLD_RANK} -eq 0 ]]; then
    mkdir -p "${outputdir}"
    : > "${outputdir}/manifest-${SLURM_JOB_ID}.txt"
fi
for opts in "${opts_arr[@]}"
do
    fin_out_dir="${output}-${at_opt}"
    mkdir -p "${fin_out_dir}" || exit
    if [[ ${OMPI_COMM_WORLD_RANK} -eq 0 ]]; then
        printf '%s\t%s\n' "${fin_out_dir}" "${opts}" \
            >> "${outputdir}/manifest-${SLURM_JOB_ID}.txt"
    fi
    for file in $files
    do
        filename=$(basename "$file")
        filenameunzipped=${filename%.xz}

        # create dir
        echo "mkdir -p ${outputdir}/${fin_out_dir}" >> todo
        echo "cp ${SLURM_SUBMIT_DIR}/${filespos}/${filename} ." >> todo
        echo "unxz ${filename}" >> todo
        baseout="${fin_out_dir}/${filename}"

        # run
        echo "/usr/bin/time --verbose -o ${baseout}.timeout ./doalarm -t real ${tlimit} ${opts} ${filenameunzipped} > ${baseout}.out 2>&1" >> todo

        #copy back result
        echo "xz ${baseout}.out*" >> todo
        echo "xz ${baseout}.timeout*" >> todo
        echo "rm -f core.*" >> todo
        echo "rm -f ${baseout}_*" >> todo

        echo "mv ${baseout}.out* ${outputdir}/${fin_out_dir}/" >> todo
        echo "mv ${baseout}.timeout* ${outputdir}/${fin_out_dir}/"  >> todo

        #lines:
        # 3+1+4+2 = 10

        numlines=$((numlines+1))
    done
    let at_opt=at_opt+1
done
mylinesper=10

# create per-core todos
numper=$((numlines/numthreads))
remain=$((numlines-numper*numthreads))
if [[ $remain -ge 1 ]]; then
    numper=$((numper+1))
fi
remain=$((numlines-numper*(numthreads-1)))

moretime=$((tlimit+20))
mystart=0
for ((myi=0; myi < numthreads ; myi++))
do
    rm -f todo_$myi.sh
    if [[ $myi -eq $OMPI_COMM_WORLD_RANK ]]; then
        touch todo_$myi.sh
        echo "#!/bin/bash" > todo_$myi.sh
        echo "ulimit -t $moretime" >> todo_$myi.sh
        echo "ulimit -v $memlimit" >> todo_$myi.sh
        echo "ulimit -c 0" >> todo_$myi.sh
        echo "set -x" >> todo_$myi.sh
    fi
    typeset -i myi
    typeset -i numper
    typeset -i mystart
    mystart=$((mystart + numper))
    if [[ $myi -lt $((numthreads-1)) ]]; then
        if [[ $mystart -gt $((numlines+numper)) ]]; then
            sleep 0
        else
            if [[ $mystart -lt $numlines ]]; then
                myp=$((numper*mylinesper))
                mys=$((mystart*mylinesper))
                if [[ $myi -eq $OMPI_COMM_WORLD_RANK ]]; then
                    head -n $mys todo | tail -n $myp >> todo_$myi.sh
                fi
            else
                #we are at boundary, e.g. numlines is 100, numper is 3, mystart is 102
                #we must only print the last numper-(mystart-numlines) = 3-2 = 1
                mys=$((mystart*mylinesper))
                p=$(( numper-mystart+numlines ))
                if [[ $p -gt 0 ]]; then
                    myp=$((p*mylinesper))
                    if [[ $myi -eq $OMPI_COMM_WORLD_RANK ]]; then
                        head -n $mys todo | tail -n $myp >> todo_$myi.sh
                    fi
                fi
            fi
        fi
    else
        if [[ $remain -gt 0 ]]; then
            mys=$((mystart*mylinesper))
            mr=$((remain*mylinesper))
            if [[ $myi -eq $OMPI_COMM_WORLD_RANK ]]; then
                head -n $mys todo | tail -n $mr >> todo_$myi.sh
            fi
        fi
    fi
    if [[ $myi -eq $OMPI_COMM_WORLD_RANK ]]; then
        echo "exit 0" >> todo_$myi.sh
        chmod +x todo_$myi.sh
    fi
done
# echo "Done."

# Execute todos
echo "This is MPI exec number $OMPI_COMM_WORLD_RANK"
rm -f ${output}/out_${OMPI_COMM_WORLD_RANK}
./todo_${OMPI_COMM_WORLD_RANK}.sh > ${output}/out_${OMPI_COMM_WORLD_RANK}
echo "Finished waiting rank $OMPI_COMM_WORLD_RANK"

rm -rf *
rm -f doalarm
rm -f cvc*
exit 0
