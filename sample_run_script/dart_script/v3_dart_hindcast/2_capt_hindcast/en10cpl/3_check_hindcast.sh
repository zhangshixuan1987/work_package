#!/bin/bash
#------------------------------------------------------------------------------
# Step 3 -- check and validate the hindcast (run on a login node).
# Prints the state of every member. A member whose run has ended successfully
# is validated (final restart + expected history files, see
# my_expected_*_streams); it is then recorded as COMPLETE
# (hindcast_complete.ENxx), which Steps 4 and 5 require, or INVALID.
#
#   ./3_check_hindcast.sh [members]
#
#   NO_CASE   not set up              READY     set up, never submitted
#   PENDING   queued                  RUNNING   running
#   FINISHED  run ended, not yet validated
#   COMPLETE  outputs validated       INVALID   outputs missing/unreadable
#   FAILED    run failed/cancelled    FAILED_ARCHIVE  run OK, st_archive failed
#------------------------------------------------------------------------------

MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"

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

select_members "$MEMBER_SELECTION"
stage_lock step3_check
command -v ncdump >/dev/null || die "ncdump not on PATH"
active_jobs >/dev/null

declare -A count=()
failed=() problems=()
printf '%-6s %-15s %-22s %-10s %s\n' MEMBER STATUS LAST_ARCHIVED_RESTART JOB NOTE
for m in "${MEMBERS[@]}"; do
  t=$(mtag "$m"); note=""
  st=$(member_status "$m")
  if [[ $st == FINISHED || $st == INVALID ]]; then
    if out=$(validate_member_outputs "$m"); then
      record_rm "hindcast_invalid.${t}"
      record_write "hindcast_complete.${t}" "member=${t}" "end=${END_STAMP}" "archive=$(marchive "$m")"
      st=COMPLETE; note="validated"
    else
      record_write "hindcast_invalid.${t}" "member=${t}" "reason=$(echo "$out" | head -1)"
      st=INVALID; note=$(echo "$out" | head -1)
    fi
  fi
  count[$st]=$(( ${count[$st]:-0} + 1 ))
  case $st in FAILED*|INVALID) failed+=("$m") ;; esac
  rest=$(ls -1 "$(marchive "$m")/rest" 2>/dev/null | tail -1 || true)
  job=$(record_get "hindcast_submitted.${t}" job_id)
  printf '%-6s %-15s %-22s %-10s %s\n' "$t" "$st" "${rest:--}" "${job:--}" "$note"
done

summary=""
for st in "${!count[@]}"; do summary+="${st}=${count[$st]} "; done
log "summary: ${summary}"
if (( ${#failed[@]} )); then
  log "resubmit failed/invalid members: sbatch 2_run_e3sm_hindcast.sh $(csv "${failed[@]}")"
  log "FAILED_ARCHIVE only: run ./case.st_archive in that member's case_scripts, then rerun this check"
fi
if [[ ${count[COMPLETE]:-0} == "${#MEMBERS[@]}" ]]; then
  log "all selected members COMPLETE; next: ./4_post_process_all.sh"
fi
