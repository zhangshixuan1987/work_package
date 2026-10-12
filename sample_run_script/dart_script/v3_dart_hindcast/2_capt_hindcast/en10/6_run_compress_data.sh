#!/bin/bash -l
#SBATCH --account=esmd
#SBATCH --time=02:00:00
#SBATCH --partition=short
#SBATCH --job-name=e3sm_hindcast_compress
#SBATCH --nodes=1
#SBATCH --output=runtmp/logs/e3sm_hindcast_compress.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue
#------------------------------------------------------------------------------
# Step 5 -- compress the member archives (optional).
# Requires hindcast_complete.ENxx (Step 3). Each member is compressed by
# workflow_lib/compress/compress_member.sh: netCDF-4 classic, deflate level 5,
# temp copy + ncdump check + atomic rename; already-compressed files are
# skipped. Writes compress_complete.ENxx; completed members are skipped.
#
#   sbatch 6_run_compress_data.sh [members]
#------------------------------------------------------------------------------

# User settings (may also be given in the environment)
MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"
COMPRESS_HISTORY="${COMPRESS_HISTORY:-TRUE}"
COMPRESS_RESTARTS="${COMPRESS_RESTARTS:-TRUE}"
COMPRESS_MAX_PARALLEL="${COMPRESS_MAX_PARALLEL:-40}"     # nccopy processes
REQUIRE_UPSTREAM_RECORDS="${REQUIRE_UPSTREAM_RECORDS:-TRUE}"
SKIP_COMPLETED="${SKIP_COMPLETED:-TRUE}"

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

for v in COMPRESS_HISTORY COMPRESS_RESTARTS REQUIRE_UPSTREAM_RECORDS SKIP_COMPLETED; do
  printf -v "$v" '%s' "$(normalize_bool "$v" "${!v}")"
done
[[ $COMPRESS_MAX_PARALLEL =~ ^[1-9][0-9]*$ ]] || die "COMPRESS_MAX_PARALLEL must be a positive integer"
[[ $COMPRESS_HISTORY == TRUE || $COMPRESS_RESTARTS == TRUE ]] || die "both compression modes are disabled"
select_members "$MEMBER_SELECTION"
summary_header
stage_lock step5_compress
require_own_case
[[ -r $my_machine_env_file ]] || die "machine environment not readable: ${my_machine_env_file}"
source_env "$my_machine_env_file"
for c in nccopy ncdump xargs; do command -v "$c" >/dev/null || die "$c not on PATH"; done

problems=()
active_jobs >/dev/null
for m in "${MEMBERS[@]}"; do
  t=$(mtag "$m")
  case $(member_status "$m") in PENDING|RUNNING) problems+=("${t}: model job still active") ;; esac
  if [[ $REQUIRE_UPSTREAM_RECORDS == TRUE ]] && ! record_valid "hindcast_complete.${t}"; then
    problems+=("${t}: no hindcast_complete record (run ./3_check_hindcast.sh)")
  fi
done
if (( ${#problems[@]} )); then
  for p in "${problems[@]}"; do printf '  PROBLEM: %s\n' "$p" >&2; done
  die "compression preflight failed with ${#problems[@]} problem(s)"
fi

compress_one() {
  local m=$1 t; t=$(mtag "$m")
  if [[ $SKIP_COMPLETED == TRUE ]] && record_valid "compress_complete.${t}"; then
    log "${t}: compressed earlier, skipping"; return 0
  fi
  record_rm "compress_complete.${t}"
  COMPRESS_DRIVER_ACTIVE=TRUE WORKFLOW_ROOT=$my_workflow_root COMPRESS_MEMBER=$m \
    COMPRESS_HISTORY=$COMPRESS_HISTORY COMPRESS_RESTARTS=$COMPRESS_RESTARTS \
    COMPRESS_MAX_PARALLEL=$COMPRESS_MAX_PARALLEL bash "${my_workflow_lib}/compress/compress_member.sh"
  record_write "compress_complete.${t}" "member=${t}" "history=${COMPRESS_HISTORY}" "restarts=${COMPRESS_RESTARTS}"
}

# members one after another; files within a member in parallel
JOBS=1
run_pool compress compress_one "${MEMBERS[@]}"
(( ${#POOL_FAIL[@]} == 0 )) || die "compression failed for $(tags "${POOL_FAIL[@]}"); resubmit to retry (finished files are skipped)"
log "compression finished for $(tags "${MEMBERS[@]}")"
