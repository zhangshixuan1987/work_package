#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=post
#SBATCH  --nodes=1
#SBATCH  --output=runtmp/logs/post.%j
#SBATCH  --exclusive
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --partition=short
#------------------------------------------------------------------------------
# Average the regridded 6-hourly series (eam_6hourly from Step 4) to daily.
# Files already in ts/daily (eam_daily from eam.h1) are kept. Requires the
# eam_6hourly completion record of each member; writes
# postprocess_complete.eam_daily_from_6hourly.ENxx.<range>.
#
#   sbatch 5_post_6hr_to_daily.sh [members]
#------------------------------------------------------------------------------

# User settings (may also be given in the environment)
MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"
POST_START="${POST_START:-}"                        # same range as Step 4
POST_END="${POST_END:-}"
OVERWRITE_EXISTING="${OVERWRITE_EXISTING:-FALSE}"
MAX_CONCURRENT_WORKERS="${MAX_CONCURRENT_WORKERS:-8}"

set -Eeo pipefail
if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}")
else
  WORK_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")
fi
[[ -r "${WORK_DIR}/workflow_lib/common.sh" ]] || { echo "ERROR: ${WORK_DIR} is not a hindcast workflow directory; submit from it" >&2; exit 1; }
WF_SCRIPT_NAME=$(basename "${BASH_SOURCE[0]}")
source "${WORK_DIR}/workflow_lib/common.sh"
trap wf_exit_summary EXIT
wf_init "${WORK_DIR}"
source "${my_workflow_lib}/post/post_lib.sh"

printf -v OVERWRITE_EXISTING '%s' "$(normalize_bool OVERWRITE_EXISTING "$OVERWRITE_EXISTING")"
[[ $MAX_CONCURRENT_WORKERS =~ ^[1-9][0-9]*$ ]] || die "MAX_CONCURRENT_WORKERS must be a positive integer"
POST_START=${POST_START:-$my_post_start}
POST_END=${POST_END:-$my_post_end}
export POST_START POST_END
if [[ $OVERWRITE_EXISTING == TRUE ]]; then export FORCE=1; fi
select_members "$MEMBER_SELECTION"
summary_header
stage_lock step4_convert_daily
require_own_case
source_env "$my_e3sm_unified"
command -v ncra >/dev/null || die "ncra not on PATH after ${my_e3sm_unified}"

missing=()
for m in "${MEMBERS[@]}"; do
  record_valid "$(post_record eam_6hourly "$m")" || missing+=("$m")
done
(( ${#missing[@]} == 0 )) || die "eam_6hourly not complete for $(tags "${missing[@]}") (range ${POST_START}_${POST_END}); run ./4_post_process_all.sh first"

convert_one() { post_job eam_daily_from_6hourly "$1"; }
JOBS=$MAX_CONCURRENT_WORKERS
run_pool convert_daily convert_one "${MEMBERS[@]}"
(( ${#POOL_FAIL[@]} == 0 )) || die "daily conversion failed for $(tags "${POOL_FAIL[@]}")"
log "daily conversion finished"
