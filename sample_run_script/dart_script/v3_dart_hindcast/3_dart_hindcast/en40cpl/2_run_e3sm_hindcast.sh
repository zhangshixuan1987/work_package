#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --account=esmd
#SBATCH  --time=08:00:00
#SBATCH  --partition=slurm
#SBATCH  --job-name=e3sm_hindcast_run
#SBATCH  --nodes=1
#SBATCH  --output=runtmp/logs/e3sm_hindcast_run.%j
#SBATCH  --exclusive
#SBATCH  --no-kill
#SBATCH  --requeue
#------------------------------------------------------------------------------
# Step 2 -- submit the hindcast of each member.
# Requires setup_complete.ENxx from Step 1. For each selected member: points
# it at its DA initial conditions, verifies the generated namelists / MPAS
# streams / restart pointers, then case.submit. Writes hindcast_submitted.ENxx
# and invalidates that member's downstream records.
#
#   sbatch 2_run_e3sm_hindcast.sh [members]     (or ./2_run_e3sm_hindcast.sh)
#   sbatch 2_run_e3sm_hindcast.sh failed        resubmit only failed members
#
# Members that are queued, running, finished or complete are skipped, so the
# script is safe to resubmit. Check progress with ./3_check_hindcast.sh.
#------------------------------------------------------------------------------

# User settings (may also be given in the environment)
MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"   # all | 1,3,7-9 | EN05 | failed | unfinished
RESUBMIT_COMPLETED="${RESUBMIT_COMPLETED:-FALSE}"   # also rerun FINISHED/COMPLETE members
DRY_RUN="${DRY_RUN:-FALSE}"                         # only show what would be submitted

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

for v in RESUBMIT_COMPLETED DRY_RUN; do printf -v "$v" '%s' "$(normalize_bool "$v" "${!v}")"; done
select_members "$MEMBER_SELECTION"
summary_header
stage_lock step2_submit
check_inputs

todo=() skipped=() missing=()
active_jobs >/dev/null
for m in "${MEMBERS[@]}"; do
  st=$(member_status "$m")
  case $st in
    NO_CASE)          missing+=("$m") ;;
    PENDING|RUNNING)  skipped+=("$(mtag "$m"):${st}") ;;
    FINISHED|COMPLETE)
      if [[ $RESUBMIT_COMPLETED == TRUE ]]; then todo+=("$m"); else skipped+=("$(mtag "$m"):${st}"); fi ;;
    *)                todo+=("$m") ;;    # READY FAILED FAILED_ARCHIVE INVALID
  esac
done
(( ${#missing[@]} == 0 )) || die "members not set up: $(tags "${missing[@]}"); run: sbatch 1_run_e3sm_ensemble_setup.sh $(csv "${missing[@]}")"
(( ${#skipped[@]} == 0 )) || log "skipping ${skipped[*]}"
(( ${#todo[@]} )) || { log "nothing to submit"; exit 0; }

if [[ $DRY_RUN == TRUE ]]; then
  log "DRY_RUN: would configure and submit $(tags "${todo[@]}")"
  for m in "${todo[@]}"; do log "  $(mtag "$m"): $(mcase "$m") <- $(ref_restdir "$m")"; done
  exit 0
fi

configure_and_submit() {
  configure_member "$1"
  submit_member "$1"
}
JOBS=${my_max_parallel_setup:-4}
run_pool run configure_and_submit "${todo[@]}"
log "submitted=[$(tags "${POOL_OK[@]}")] failed=[$(tags "${POOL_FAIL[@]}")]"
(( ${#POOL_FAIL[@]} == 0 )) || die "fix the failures above, then: sbatch 2_run_e3sm_hindcast.sh $(csv "${POOL_FAIL[@]}")"
log "next: ./3_check_hindcast.sh (until all members are COMPLETE)"
