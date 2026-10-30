#!/bin/bash -el
#SBATCH --job-name=post_process
#SBATCH --account=e3sm
#SBATCH --nodes=1
#SBATCH --output=post_process.o%j
#SBATCH --exclusive
#SBATCH --time=24:00:00
#SBATCH --qos=regular
#SBATCH --constraint=cpu

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

run1_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
run2_dir="/global/cfs/cdirs/e3sm/S2S2D/simulation"
work_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVarOcn"

max_jobs=1

#leads=(1980050100 1980110100)
#leads=(1981050100 1981110100)
#leads=(1982050100 1982110100)
#leads=(1983050100 1983110100)
#leads=(1984050100 1984110100)
#leads=(1985050100 1985110100)
#leads=(1986050100 1986110100)
#leads=(1987050100 1987110100)
#leads=(1988050100 1988110100)
#leads=(1989050100 1989110100) 
#leads=(1990050100 1990110100)
#leads=(1991050100 1991110100)
#leads=(1992050100 1992110100)
#leads=(1993050100 1993110100)
#leads=(1994050100 1994110100)
#leads=(1995050100 1995110100)
#leads=(1996050100 1996110100)
#leads=(1997050100 1997110100)
#leads=(1998050100 1998110100)
#leads=(1999050100 1999110100)
#leads=(2000050100 2000110100)
#leads=(2001050100 2001110100)
#leads=(2002050100 2002110100)
#leads=(2004050100 2004110100)
#leads=(2010050100 2010110100)
#leads=(2011050100 2011110100)
leads=(2012050100 2012110100)
#leads=(2013050100 2013110100)
#leads=(2014050100 2014110100)
#leads=(2015050100 2015110100)
#leads=(2016050100 2016110100) 
#leads=(2017050100 2017110100)
#leads=(2018050100 2018110100)

log_dir="${work_dir}/post_process_logs"
mkdir -p "${log_dir}"

main_log="${log_dir}/post_process_${SLURM_JOB_ID:-manual}.log"
fail_log="${log_dir}/post_process_${SLURM_JOB_ID:-manual}.failed"

: > "${main_log}"
: > "${fail_log}"

timestamp() {
    date "+%Y-%m-%d %H:%M:%S"
}

log_msg() {
    echo "[$(timestamp)] $*" | tee -a "${main_log}"
}

trap 'log_msg "Received termination signal. Killing child jobs..."; jobs -rp | xargs -r kill || true; wait || true' SIGINT SIGTERM

wait_for_slot() {
    while (( $(jobs -rp | wc -l) >= max_jobs )); do
        wait -n || true
    done
}

wait_for_all() {
    while (( $(jobs -rp | wc -l) > 0 )); do
        wait -n || true
    done
}

# ------------------------------------------------------------------
# Copy zstash script to stable work directory for safer execution
# ------------------------------------------------------------------
zstash_scp="${work_dir}/zstash_S2D.sh"
if [[ ! -f "${zstash_scp}" ]]; then
    echo "ERROR: Source zstash script not found: ${zstash_scp}"
    exit 1
fi
chmod +x "${zstash_scp}"
log_msg "Using zstash script: ${zstash_scp}"

log_msg "Starting post-process workflow"
log_msg "run1_dir=${run1_dir}"
log_msg "run2_dir=${run2_dir}"
log_msg "work_dir=${work_dir}"
log_msg "max_jobs=${max_jobs}"

# =========================================================
# Stage 1:
# Create zstash archive from run1_dir
# =========================================================
log_msg "Starting Stage 1: zstash create from run1_dir"

cd "${run1_dir}"

for lead in "${leads[@]}"; do
    case_name="${exp}_${lead}"

    if [[ -d "${case_name}" ]]; then
        wait_for_slot
        (
            log_msg "Stage 1 start: ${case_name}"

            if [[ -f "${case_name}/zstash/index.db" ]]; then
                log_msg "Existing zstash metadata found for ${case_name} → running update"
                "${zstash_scp}" update "${case_name}"
            else
                log_msg "No existing zstash metadata found for ${case_name} → running create"
                "${zstash_scp}" create "${case_name}"
            fi

            log_msg "Stage 1 complete: ${case_name}"
        ) || echo "Stage 1 failed: ${case_name}" >> "${fail_log}" &
    else
        log_msg "Stage 1 missing directory: ${run1_dir}/${case_name}"
        echo "Stage 1 missing directory: ${run1_dir}/${case_name}" >> "${fail_log}"
    fi

done

wait_for_all

if [[ -s "${fail_log}" ]]; then
    log_msg "Stage 1 completed with failures. See ${fail_log}"
    exit 1
else
    log_msg "Stage 1 completed successfully."
fi

# Reset failure log
: > "${fail_log}"

# =========================================================
# Stage 2:
# Copy zstash metadata from Stage 1 and continue archive
# from run2_dir to the same HPSS target
# =========================================================
log_msg "Starting Stage 2: continue zstash archive from run2_dir"

cd "${run2_dir}"

for lead in "${leads[@]}"; do
    case_name="${exp}_${lead}"

    if [[ -d "${case_name}" ]]; then
        wait_for_slot
        (
            log_msg "Stage 2 start: ${case_name}"

            mkdir -p "${case_name}/zstash"

            if [[ -f "${case_name}/zstash/index.db" ]]; then
                log_msg "Existing local zstash metadata found for ${case_name} -> running update"
                "${zstash_scp}" update "${case_name}"
            else
                if [[ -d "${run1_dir}/${case_name}/zstash" ]]; then
                    log_msg "Copying zstash metadata from Stage 1 for ${case_name}"
                    cp -rp "${run1_dir}/${case_name}/zstash/"* \
                        "${case_name}/zstash/" 2>/dev/null || true
                else
                    log_msg "Stage 2 warning: missing zstash metadata in ${run1_dir}/${case_name}"
                fi

                if [[ -f "${case_name}/zstash/index.db" ]]; then
                    log_msg "Copied zstash metadata found for ${case_name} -> running update"
                    "${zstash_scp}" update "${case_name}"
                else
                    log_msg "No zstash index.db found for ${case_name} -> running create"
                    "${zstash_scp}" create "${case_name}"
                fi
            fi

            log_msg "Stage 2 complete: ${case_name}"
        ) || echo "Stage 2 failed: ${case_name}" >> "${fail_log}" &
    else
        log_msg "Stage 2 missing directory: ${run2_dir}/${case_name}"
        echo "Stage 2 missing directory: ${run2_dir}/${case_name}" >> "${fail_log}"
    fi
done

wait_for_all

if [[ -s "${fail_log}" ]]; then
    log_msg "Stage 2 completed with failures. See ${fail_log}"
    exit 1
else
    log_msg "Stage 2 completed successfully."
fi

# Reset failure log
: > "${fail_log}"

# =========================================================
# Stage 3:
# Final verification of the combined archive using zstash check
# =========================================================
log_msg "Starting Stage 3: final zstash check"

cd "${run2_dir}"

for lead in "${leads[@]}"; do
    case_name="${exp}_${lead}"

    if [[ -d "${case_name}" ]]; then
        wait_for_slot
        (
            log_msg "Stage 3 check start: ${case_name}"
            "${zstash_scp}" check "${case_name}"
            log_msg "Stage 3 check complete: ${case_name}"
        ) || echo "Stage 3 check failed: ${case_name}" >> "${fail_log}" &
    else
        log_msg "Stage 3 missing directory: ${run2_dir}/${case_name}"
        echo "Stage 3 missing directory: ${run2_dir}/${case_name}" >> "${fail_log}"
    fi
done

wait_for_all

if [[ -s "${fail_log}" ]]; then
    log_msg "Stage 3 completed with check failures. See ${fail_log}"
    exit 1
else
    log_msg "Stage 3 completed successfully."
    rm -f "${fail_log}"
fi

log_msg "All post-processing completed successfully."

exit 0
