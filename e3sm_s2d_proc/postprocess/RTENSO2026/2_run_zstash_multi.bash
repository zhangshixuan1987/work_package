#!/bin/bash -el
#SBATCH --job-name=rtenso_zstash
#SBATCH --qos=xfer
#SBATCH --time=2-00:00:00
#SBATCH --output=/global/cfs/cdirs/e3sm/zhan391/s2d_diag/post_process_logs/zstash_%x_%j.out
#
# Archive RTENSO2026 cases to HPSS (HPSS_PTH/<case>) in three stages:
#   1. zstash create/update from run1_dir (pscratch: restarts, logs, run/, case_scripts/)
#   2. copy the index to run2_dir (CFS) and zstash update with the moved history
#   3. zstash check of the combined archive
# Resumable: a case with an index.db is updated, so a killed job can be resubmitted as is.
# Runs on the xfer queue (login nodes, 2-day limit) so a login-node restart does not kill it.
#
#   IC_METHOD=obs LEADS="2015020100 2015050100" sbatch 2_run_zstash_multi.bash
#
# Run stage 1 only after the case's history has been moved off pscratch.

source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

run1_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
run2_dir="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/simulation"   # 1_run_mov_archive.bash out_dir
#work_dir="/global/cfs/cdirs/e3sm/S2S2D/s2d_diag"
work_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag"
export HPSS_PTH="/home/projects/e3sm/www/RTENSO2026"   # read by zstash_S2D.sh

ic_method="${IC_METHOD:?set IC_METHOD (clim | obs | obs_ice_nsidc)}"
exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_RTENSO2026_${ic_method}"
max_jobs=1

leads=(${LEADS:?set LEADS, space-separated start dates, e.g. "2026020100 2026050100"})

log_dir="${work_dir}/post_process_logs"
mkdir -p "${log_dir}"

tag="zstash_${ic_method}_${SLURM_JOB_ID:-manual$$}"     # one log per job, so parallel jobs do not share it
main_log="${log_dir}/${tag}.log"
fail_log="${log_dir}/${tag}.failed"

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
log_msg "HPSS_PTH=${HPSS_PTH}"
log_msg "max_jobs=${max_jobs}"
log_msg "ic_method=${ic_method} leads=${leads[*]}"

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
