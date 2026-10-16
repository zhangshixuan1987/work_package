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
set -u 

###############################################################################
# USER SETTINGS
###############################################################################

run_dir="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
work_dir="/global/cfs/cdirs/e3sm/zhan391/s2d_diag"

exp="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_BruteForce"

#leads=(1987050100)

leads=(1987110100)
ensembles=(EN00 EN01 EN02 EN03 EN04 EN05 EN06 EN07 EN08 EN09)

hpss_root="/home/projects/e3sm/www/S2S2D"

###############################################################################
# MAIN LOOP
###############################################################################

mkdir -p "${run_dir}"
mkdir -p "${work_dir}"

for lead in "${leads[@]}"; do

    case_name="${exp}_${lead}"
    case_dir="${run_dir}/${case_name}"
    zstash_cache="${case_dir}/zstash"

    echo "======================================================================"
    echo "Case:       ${case_name}"
    echo "Case dir:   ${case_dir}"
    echo "HPSS path:  ${hpss_root}/${case_name}"
    echo "Cache dir:  ${zstash_cache}"
    echo "======================================================================"

    mkdir -p "${case_dir}"

    # Start with a clean local zstash cache.
    rm -rf "${zstash_cache}"
    mkdir -p "${zstash_cache}"

    cd "${case_dir}" || {
        echo "ERROR: Cannot enter ${case_dir}" >&2
        continue
    }

    declare -A expected_hist_counts=()
    for ensemble in "${ensembles[@]}"; do
        expected_hist_counts["${ensemble}"]=0
    done

    if ! remote_hist_listing=$(zstash ls \
        --hpss="${hpss_root}/${case_name}" \
        "EN0[0-9]/archive/*/hist/*"); then
        echo "ERROR: Failed to list remote hist files for ${case_name}" >&2
        continue
    fi

    while IFS= read -r remote_path; do
        if [[ "${remote_path}" =~ ^(EN[0-9][0-9])/archive/.*/hist/ ]]; then
            ensemble="${BASH_REMATCH[1]}"
            ((expected_hist_counts["${ensemble}"]+=1))
        fi
    done <<< "${remote_hist_listing}"

    pids=()
    started_ensembles=()
    failed=0

    for ensemble in "${ensembles[@]}"; do
        local_hist_count=$(find "${case_dir}/${ensemble}/archive" -type f -path '*/hist/*' 2>/dev/null | wc -l)
        expected_hist_count=${expected_hist_counts["${ensemble}"]}

        if [[ "${expected_hist_count}" -eq 0 ]]; then
            echo "Skipping ${case_name} ${ensemble}; no remote hist files found"
            continue
        fi

        if [[ "${local_hist_count}" -ge "${expected_hist_count}" ]]; then
            echo "Skipping ${case_name} ${ensemble}; found ${local_hist_count}/${expected_hist_count} hist files locally"
            continue
        fi

        if [[ "${local_hist_count}" -gt 0 ]]; then
            echo "Re-extracting ${case_name} ${ensemble}; found ${local_hist_count}/${expected_hist_count} hist files locally"
        else
            echo "Starting extract for ${case_name} ${ensemble}"
        fi

        zstash extract \
            -v \
            --hpss="${hpss_root}/${case_name}" \
            "${ensemble}/archive/*/hist/*" &
        pids+=("$!")
        started_ensembles+=("${ensemble}")
    done

    for index in "${!pids[@]}"; do
        if ! wait "${pids[$index]}"; then
            echo "ERROR: zstash extraction failed for ${case_name} ${started_ensembles[$index]}" >&2
            failed=1
        fi
    done

    if [[ "${#started_ensembles[@]}" -eq 0 ]]; then
        echo "No extraction needed for ${case_name}; all ensemble files already exist"
    elif [[ "${failed}" -eq 0 ]]; then
        echo "Successfully extracted ${case_name} ensembles ${started_ensembles[*]}"
    else
        echo "ERROR: One or more ensemble extracts failed for ${case_name}" >&2
    fi

done
