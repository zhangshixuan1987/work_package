#!/bin/bash -l
#------------------------------------------------------------------------------
# Batch system directives
#------------------------------------------------------------------------------
#SBATCH  --job-name=e3sm_hindcast_setup
#SBATCH  --nodes=1
#SBATCH  --output=runtmp/logs/e3sm_hindcast_setup.%j
#SBATCH  --exclusive
#SBATCH  --account=esmd
#SBATCH  --time=02:00:00
#SBATCH  --partition=short
#------------------------------------------------------------------------------
# Step 1 -- ensemble setup.
# Builds e3sm.exe once (base case) with the template
# workflow_lib/run_template/compile_and_setup_e3sm.${my_compset}.sh, then
# creates every member case with that executable and configures it for its DA
# initial conditions (verified with preview_namelists).
#
#   sbatch 1_run_e3sm_ensemble_setup.sh [members]
#
# Writes runtmp/status/build_complete, setup_complete.ENxx for each member and
# setup_complete once every member is ready. Members that are already set up
# are skipped, so a failed or partial setup is simply resubmitted.
#------------------------------------------------------------------------------

# User settings (may also be given in the environment)
MEMBER_SELECTION="${1:-${MEMBER_SELECTION:-all}}"   # all | 1,3,7-9 | EN05 | unfinished
REBUILD="${REBUILD:-FALSE}"                         # recompile e3sm.exe even if built
RECREATE_EXISTING="${RECREATE_EXISTING:-FALSE}"     # recreate members already set up
ADOPT_EXISTING_CASE="${ADOPT_EXISTING_CASE:-FALSE}" # take over a case made by the old scripts
PREFLIGHT_ONLY="${PREFLIGHT_ONLY:-FALSE}"           # only check inputs, change nothing

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

for v in REBUILD RECREATE_EXISTING ADOPT_EXISTING_CASE PREFLIGHT_ONLY; do
  printf -v "$v" '%s' "$(normalize_bool "$v" "${!v}")"
done
select_members "$MEMBER_SELECTION"
summary_header
stage_lock step1_setup
check_inputs
if [[ $PREFLIGHT_ONLY == TRUE ]]; then log "PREFLIGHT_ONLY=TRUE: nothing changed"; exit 0; fi

active=()
active_jobs >/dev/null
for m in "${MEMBERS[@]}"; do
  case $(member_status "$m") in PENDING|RUNNING) active+=("$m") ;; esac
done
(( ${#active[@]} == 0 )) || die "members with queued/running jobs, refusing to set them up: $(tags "${active[@]}")"
if [[ -n $(ownership_problem) && $ADOPT_EXISTING_CASE != TRUE ]]; then die "$(ownership_problem)"; fi
claim_case_dir

record_rm setup_complete
record_write setup_in_progress "members=$(csv "${MEMBERS[@]}")"
materialize_run_template

#--- base case: build e3sm.exe once ---------------------------------------------
if [[ $REBUILD == TRUE ]] || ! record_valid build_complete || [[ ! -x $(base_exe) ]]; then
  log "building e3sm.exe (log: ${LOG_DIR}/build.log)"
  set +e
  ( set -Eeo pipefail; build_base ) > "${LOG_DIR}/build.log" 2>&1
  rc=$?
  set -e
  if (( rc != 0 )); then tail -20 "${LOG_DIR}/build.log" >&2; die "build failed, see ${LOG_DIR}/build.log"; fi
else
  log "reusing $(base_exe) (REBUILD=TRUE to recompile)"
fi

#--- member cases ---------------------------------------------------------------
todo=()
for m in "${MEMBERS[@]}"; do
  t=$(mtag "$m")
  if record_valid "setup_complete.${t}" && [[ $RECREATE_EXISTING != TRUE ]]; then
    if ! cmp -s "$(base_exe)" "$(mdir "$m")/build/e3sm.exe"; then
      cp -p "$(base_exe)" "$(mdir "$m")/build/e3sm.exe"
      log "${t}: already set up; refreshed e3sm.exe from the new build"
    else
      log "${t}: already set up, skipping (RECREATE_EXISTING=TRUE to recreate)"
    fi
  else
    todo+=("$m")
  fi
done
JOBS=${my_max_parallel_setup:-4}
if (( ${#todo[@]} )); then
  run_pool setup setup_member "${todo[@]}"
  log "setup: ok=[$(tags "${POOL_OK[@]}")] failed=[$(tags "${POOL_FAIL[@]}")]"
  (( ${#POOL_FAIL[@]} == 0 )) || die "setup failed for $(tags "${POOL_FAIL[@]}"); fix and resubmit: sbatch 1_run_e3sm_ensemble_setup.sh $(csv "${POOL_FAIL[@]}")"
fi

#--- overall completion record ------------------------------------------------------
missing=()
for m in $(seq "$ENS_START" "$ENS_SIZE"); do
  record_valid "setup_complete.$(mtag "$m")" || missing+=("$m")
done
record_rm setup_in_progress
if (( ${#missing[@]} == 0 )); then
  record_write setup_complete "exe=$(base_exe)"
  log "all members ready; next: sbatch 2_run_e3sm_hindcast.sh"
else
  log "selected members ready; not yet set up: $(tags "${missing[@]}")"
fi
