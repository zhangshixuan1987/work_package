#!/usr/bin/env bash

set -uo pipefail

src_root="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
dst_root="/global/cfs/cdirs/e3sm/S2S2D/post_process"
case="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_BruteForce"
ens=10
dry_run="${DRY_RUN:-false}"

if [[ "${1:-}" == "--dry-run" ]]; then
    dry_run=true
fi

# Maximum number of simultaneous rsync jobs.
max_jobs=4

if [[ "${dry_run}" == true ]]; then
    log_root="${LOG_ROOT:-${PWD}/rsync_dry_run_logs}"
else
    log_root="${LOG_ROOT:-${dst_root}}"
    mkdir -p "${dst_root}"
fi

mkdir -p "${log_root}"

shopt -s nullglob
lead_dirs=("${src_root}/${case}"_*)
shopt -u nullglob

if (( ${#lead_dirs[@]} == 0 )); then
    echo "No lead-date directories matching:"
    echo "  ${src_root}/${case}_*"
    exit 0
fi

tasks=()
for lead_dir in "${lead_dirs[@]}"; do
    case_date=$(basename "${lead_dir}")
    lead_date="${case_date#${case}_}"

    if [[ "${lead_date}" == "${case_date}" ]]; then
        echo "Skipping unexpected directory name: ${lead_dir}"
        continue
    fi

    for (( en=0; en<ens; en++ )); do
        ENS=$(printf "EN%02d" "${en}")
        src_post="${src_root}/${case}_${lead_date}/${ENS}/post"

        if [[ -d "${src_post}" ]]; then
            tasks+=("${case}_${lead_date}/${ENS}")
        else
            echo "Skipping missing source post dir: ${src_post}"
        fi
    done
done

if (( ${#tasks[@]} == 0 )); then
    echo "No ensemble post directories found under:"
    echo "  ${src_root}/${case}_*/EN00-EN$(printf "%02d" "$((ens - 1))")/post"
    exit 0
fi

echo "Found ${#tasks[@]} ensemble post directories."
echo "Source root: ${src_root}"
echo "Destination root: ${dst_root}"
echo "Ensembles: EN00-EN$(printf "%02d" "$((ens - 1))")"
echo "Parallel jobs: ${max_jobs}"
echo "Dry run: ${dry_run}"
echo "Log root: ${log_root}"
echo

move_one_post_directory() {
    local rel="$1"
    local src="${src_root}/${rel}/post"
    local dst="${dst_root}/${rel}/post"
    local safe_name="${rel//\//_}"
    local log="${log_root}/rsync_${safe_name}_post.log"
    local verify_log="${log_root}/verify_${safe_name}_post.log"
    local rsync_opts=(-aH --partial --info=progress2)

    echo "[$(date '+%F %T')] Starting ${rel}/post"

    if [[ "${dry_run}" == true ]]; then
        rsync_opts+=(-n)
    else
        mkdir -p "${dst}"
    fi

    if ! rsync "${rsync_opts[@]}" \
        "${src}/" \
        "${dst}/" >"${log}" 2>&1; then

        echo "[$(date '+%F %T')] COPY FAILED: ${rel}/post"
        echo "  Log: ${log}"
        return 1
    fi

    if [[ "${dry_run}" == true ]]; then
        echo "[$(date '+%F %T')] DRY RUN passed: ${rel}/post"
        echo "[$(date '+%F %T')] DRY RUN would remove source after verified copy: ${src}"
        return 0
    fi

    # A dry-run comparison should produce no file-list differences.
    if ! rsync -aHnci \
        --delete \
        "${src}/" \
        "${dst}/" >"${verify_log}" 2>&1; then

        echo "[$(date '+%F %T')] VERIFICATION COMMAND FAILED: ${rel}/post"
        echo "  Log: ${verify_log}"
        return 1
    fi

    if [[ -s "${verify_log}" ]]; then
        echo "[$(date '+%F %T')] VERIFICATION FAILED: ${rel}/post"
        echo "  Differences are listed in: ${verify_log}"
        return 1
    fi

    echo "[$(date '+%F %T')] Verification passed: ${rel}/post"
    echo "[$(date '+%F %T')] Removing source: ${src}"

    if ! rm -rf -- "${src}"; then
        echo "[$(date '+%F %T')] SOURCE REMOVAL FAILED: ${rel}/post"
        return 1
    fi

    rm -f "${verify_log}"

    echo "[$(date '+%F %T')] Successfully moved: ${rel}/post"
}

export -f move_one_post_directory
export src_root
export dst_root
export dry_run
export log_root

printf '%s\0' "${tasks[@]}" |
    xargs -0 -n 1 -P "${max_jobs}" bash -c '
        move_one_post_directory "$1"
    ' _

status=$?

echo
if (( status == 0 )); then
    echo "All parallel transfer jobs completed successfully."
else
    echo "One or more transfer jobs failed."
    echo "Check the logs under:"
    echo "  ${dst_root}"
    echo
    echo "Sources are removed only when their individual verification passes."
fi

exit "${status}"
