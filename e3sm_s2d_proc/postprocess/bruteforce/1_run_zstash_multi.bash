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

run_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
#work_dir="/global/cfs/cdirs/e3sm/S2S2D/s2d_diag"
work_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_BruteForce"
max_jobs=1

leads=(2000020100 2000080100)
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
#leads=(2003050100 2003110100)
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
log_msg "run_dir=${run_dir}"
log_msg "work_dir=${work_dir}"
log_msg "max_jobs=${max_jobs}"

# =========================================================
# Stage 1:
# Create zstash archive from run_dir
# =========================================================
log_msg "Starting Stage 1: zstash create from run_dir"

cd "${run_dir}"

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
        log_msg "Stage 1 missing directory: ${run_dir}/${case_name}"
        echo "Stage 1 missing directory: ${run_dir}/${case_name}" >> "${fail_log}"
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
# Final verification of the combined archive using zstash check
# =========================================================
log_msg "Starting Stage 2: final zstash check"

cd "${run_dir}"

for lead in "${leads[@]}"; do
    case_name="${exp}_${lead}"

    if [[ -d "${case_name}" ]]; then
        wait_for_slot
        (
            log_msg "Stage 2 check start: ${case_name}"
            "${zstash_scp}" check "${case_name}"
            log_msg "Stage 2 check complete: ${case_name}"
        ) || echo "Stage 2 check failed: ${case_name}" >> "${fail_log}" &
    else
        log_msg "Stage 2 missing directory: ${run_dir}/${case_name}"
        echo "Stage 2 missing directory: ${run_dir}/${case_name}" >> "${fail_log}"
    fi
done

wait_for_all

if [[ -s "${fail_log}" ]]; then
    log_msg "Stage 2 completed with check failures. See ${fail_log}"
    exit 1
else
    log_msg "Stage 2 completed successfully."
    rm -f "${fail_log}"
fi

log_msg "All post-processing completed successfully."

exit 0
