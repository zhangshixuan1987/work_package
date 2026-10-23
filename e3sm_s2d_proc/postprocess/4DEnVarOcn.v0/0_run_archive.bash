#!/bin/bash

set -euo pipefail

# ==============================================================================
# User settings
# ==============================================================================

ROOT_DIR="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"

CASE_PATTERN="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVarOcn"

# Empty        => all dates
# "1980*"      => all 1980 starts
# "1980050100" => one specific start date
CASE_DATE="200*"

MAX_JOBS=10
WAIT_TIME=10

LOGDIR="archive_logs"

# ==============================================================================
# Setup
# ==============================================================================

mkdir -p "${LOGDIR}"
LOGDIR="$(realpath "${LOGDIR}")"

if [ -z "${CASE_DATE}" ]; then
    SEARCH_PATTERN="${ROOT_DIR}/${CASE_PATTERN}_*/*/case_scripts"
else
    SEARCH_PATTERN="${ROOT_DIR}/${CASE_PATTERN}_${CASE_DATE}/*/case_scripts"
fi

echo "========================================================="
echo "Root directory  : ${ROOT_DIR}"
echo "Case pattern    : ${CASE_PATTERN}"
echo "Case date       : ${CASE_DATE:-ALL}"
echo "Search pattern  : ${SEARCH_PATTERN}"
echo "Max jobs        : ${MAX_JOBS}"
echo "Log directory   : ${LOGDIR}"
echo "Start time      : $(date)"
echo "========================================================="

# ==============================================================================
# Main loop
# ==============================================================================

for dir in ${SEARCH_PATTERN}
do
    [ -d "${dir}" ] || continue

    (
        case_root=$(basename "$(dirname "$(dirname "${dir}")")")
        ens_name=$(basename "$(dirname "${dir}")")
        case_name="${case_root}_${ens_name}"

        echo "[$(date)] Starting ${case_name}"

        cd "${dir}" || exit 1

        ./case.st_archive \
            > "${LOGDIR}/${case_name}.out" \
            2> "${LOGDIR}/${case_name}.err"

        rc=$?

        echo "[$(date)] Finished ${case_name} with rc=${rc}"

        exit "${rc}"
    ) &

    while [ "$(jobs -rp | wc -l)" -ge "${MAX_JOBS}" ]
    do
        sleep "${WAIT_TIME}"
    done
done

wait

echo "========================================================="
echo "All archive jobs completed."
echo "Finish time: $(date)"
echo "========================================================="
