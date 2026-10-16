#!/bin/bash
#------------------------------------------------------------------------------
# Step 4: continuous free run (replaces DA cycling for this case).
# Run on a login node after Steps 1-3 (do not sbatch this script):
#
#   ./4_run_e3sm_continuous.sh            # configure and submit pending members
#   ./4_run_e3sm_continuous.sh status     # show each member's progress and jobs
#   ./4_run_e3sm_continuous.sh finalize   # after all members finish
#
# Every member starts once (hybrid) from its Step 3 perturbed eam.i file at
# my_e3sm_start_date and runs to my_e3sm_end_date in one model run; the only
# restart is written at the end time. my_free_layout selects how members run:
#   member_jobs  each member is its own CIME batch job (case.run, st_archive)
#                on my_free_member_nodes nodes in my_free_job_queue
#   single_job   one my_job_nnodes-node allocation in my_free_single_job_queue
#                runs all members side by side, my_nodes_per_member nodes each,
#                and finalizes itself when every member is done
# Rerunning submit resumes any unfinished member from its latest archived
# restart. finalize verifies history, deletes intermediate restart sets, and
# writes the completion record. No eam.i (inithist) files are written.
#------------------------------------------------------------------------------
set -Eeuo pipefail

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

MODE="${1:-submit}"
[[ "${MODE}" =~ ^(submit|status|finalize|run)$ ]] || fail "usage: $0 [submit|status|finalize]"
# "run" is the single_job allocation's entry point, submitted by "submit".
if [[ "${MODE}" == "run" ]]; then
  [[ -n "${SLURM_JOB_ID:-}" ]] || fail "run mode is internal to the single_job allocation; use: $0 submit"
  my_wkdir=$(readlink -f "${2:-}") || fail "run mode needs the case script directory"
else
  [[ -z "${SLURM_JOB_ID:-}" ]] || fail "run this script on a login node, not inside a Slurm job"
  my_wkdir=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")") || fail "could not resolve script path"
fi
SCRIPT_PATH="${my_wkdir}/4_run_e3sm_continuous.sh"
CONFIG_FILE="${my_wkdir}/create_and_setup_case.sh"
[[ -r "${SCRIPT_PATH}" ]] || fail "missing driver: ${SCRIPT_PATH}"
[[ -r "${CONFIG_FILE}" ]] || fail "missing configuration: ${CONFIG_FILE}"
source "${CONFIG_FILE}"
[[ -r "${my_dart_env_file:-}" ]] || fail "configured machine environment is not readable: ${my_dart_env_file:-unset}"
source "${my_dart_env_file}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"
LOG_DIR="${my_log_dir}"

for cmd in awk cksum cp date ex find flock ncdump sbatch squeue srun; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
# status only reads, so it can run while a submit or allocation is active.
if [[ "${MODE}" != "status" ]]; then
  exec 9>"${my_lock_dir}/e3sm_free_run.lock"
  flock -n 9 || fail "another free-run driver is active for ${my_wkdir}"
fi

validate_positive_int() {
  [[ "$2" =~ ^[1-9][0-9]*$ ]] || fail "$1 must be a positive integer, got: $2"
}

LAYOUT="${my_free_layout:-member_jobs}"
KEEP_INTERMEDIATE="${my_free_keep_intermediate_restarts:-FALSE}"
KEEP_INTERMEDIATE="${KEEP_INTERMEDIATE^^}"
case "${LAYOUT}" in
  member_jobs)
    MEMBER_NODES="${my_free_member_nodes:-20}"
    JOB_QUEUE="${my_free_job_queue:-short}"
    JOB_WALLTIME="${my_free_job_walltime:-02:00:00}"
    ;;
  single_job)
    MEMBER_NODES="${my_nodes_per_member:-2}"
    JOB_QUEUE="${my_free_single_job_queue:-slurm}"
    JOB_WALLTIME="${my_free_single_job_walltime:-24:00:00}"
    MAX_PARALLEL_SETUP="${my_max_parallel_setup:-8}"
    validate_positive_int "my_nodes_per_member" "${MEMBER_NODES}"
    validate_positive_int "my_job_nnodes" "${my_job_nnodes}"
    validate_positive_int "my_max_parallel_setup" "${MAX_PARALLEL_SETUP}"
    (( my_job_nnodes % MEMBER_NODES == 0 )) || fail "my_job_nnodes (${my_job_nnodes}) must be divisible by my_nodes_per_member (${MEMBER_NODES})"
    BATCH_SIZE=$((my_job_nnodes / MEMBER_NODES))
    ;;
  *) fail "my_free_layout must be member_jobs or single_job, got: ${LAYOUT}" ;;
esac
validate_positive_int "member nodes" "${MEMBER_NODES}"
validate_positive_int "my_task_per_node" "${my_task_per_node}"
validate_positive_int "my_ensnum" "${my_ensnum}"
[[ "${JOB_WALLTIME}" =~ ^[0-9]{1,3}:[0-5][0-9]:[0-5][0-9]$ ]] || fail "job walltime must be HH:MM:SS, got: ${JOB_WALLTIME}"
[[ "${KEEP_INTERMEDIATE}" == "TRUE" || "${KEEP_INTERMEDIATE}" == "FALSE" ]] || fail "my_free_keep_intermediate_restarts must be TRUE or FALSE"
[[ "${my_runtype}" == "AMIP" || "${my_runtype}" == "Full-CPL" ]] || fail "my_runtype must be AMIP or Full-CPL"
MEMBER_TASKS=$((MEMBER_NODES * my_task_per_node))
DRIVER_JOB_NAME="e3sm-free-$(printf '%s' "${my_wkdir}" | cksum | awk '{print $1}')"

#------------------------------------------------------------------------------
# Time helpers. Stamps use the E3SM form YYYY-MM-DD-SSSSS.
#------------------------------------------------------------------------------
stamp_to_epoch() {
  local stamp="$1" ymd tod
  [[ "${stamp}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]] || return 1
  ymd="${stamp%-*}"; tod=$((10#${stamp##*-}))
  (( tod < 86400 )) || return 1
  echo $(( $(date -u -d "${ymd}" +%s) + tod ))
}
epoch_to_stamp() {
  local ymd h m s
  read -r ymd h m s < <(date -u -d "@$1" '+%Y-%m-%d %H %M %S')
  printf '%s-%05d\n' "${ymd}" "$((10#${h} * 3600 + 10#${m} * 60 + 10#${s}))"
}

START_STAMP="${my_e3sm_start_date}-${my_e3sm_start_tod}"
END_STAMP="${my_e3sm_end_date}-${my_e3sm_end_tod}"
START_EPOCH=$(stamp_to_epoch "${START_STAMP}") || fail "invalid start time: ${START_STAMP}"
END_EPOCH=$(stamp_to_epoch "${END_STAMP}") || fail "invalid end time: ${END_STAMP}"
(( START_EPOCH % 86400 == 0 && END_EPOCH % 86400 == 0 && END_EPOCH > START_EPOCH )) || fail "start and end must be 00Z with end after start"
[[ "${START_STAMP}" == "${my_refdate}-${my_reftod}" ]] || fail "start time ${START_STAMP} differs from the Step 3 reference time ${my_refdate}-${my_reftod}"
COMPLETE_MARKER="${my_status_dir}/free_run_complete.${START_STAMP}-${END_STAMP}"

marker_field() { awk -F= -v key="$2" '$1 == key {sub(/^[^=]*=/, ""); print; found=1; exit} END {if (!found) exit 1}' "$1"; }
perturb_marker="${my_status_dir}/perturb_complete.${my_refdate}-${my_reftod}"
[[ -s "${perturb_marker}" ]] || fail "missing Step 3 completion record: ${perturb_marker}"
[[ "$(marker_field "${perturb_marker}" case)" == "${my_casename}" ]] || fail "Step 3 record is for a different case: ${perturb_marker}"
[[ "$(marker_field "${perturb_marker}" ensemble_size)" == "${my_ensnum}" ]] || fail "Step 3 record has a different ensemble size: ${perturb_marker}"

#------------------------------------------------------------------------------
# Member helpers
#------------------------------------------------------------------------------
member_name() { printf 'EN%02d' "$1"; }
member_case_dir() { echo "${my_modeldir}/$1/case_scripts"; }
member_run_dir() { echo "${my_modeldir}/$1/run"; }
member_archive() { echo "${my_modeldir}/$1/archive"; }

# Latest archived time with a real (non-empty) EAM restart, else the start time.
# Step 2 leaves a zero-byte eam.r placeholder at the start time.
member_state_stamp() {
  local case_name="${my_casename}.$1" archive dir stamp epoch best_epoch="${START_EPOCH}"
  archive=$(member_archive "$1")
  for dir in "${archive}"/rest/????-??-??-?????; do
    [[ -d "${dir}" ]] || continue
    stamp=$(basename "${dir}")
    epoch=$(stamp_to_epoch "${stamp}") || continue
    (( epoch > best_epoch && epoch <= END_EPOCH )) || continue
    [[ -s "${dir}/${case_name}.eam.r.${stamp}.nc" ]] || continue
    best_epoch=${epoch}
  done
  epoch_to_stamp "${best_epoch}"
}

# Slurm jobs (member_jobs case.run or st_archive) belonging to this member's case.
member_jobs() {
  squeue -h -u "${USER}" -o '%i %j %T' | awk -v c="${my_casename}.$1" 'index($2, c) {print $1 ":" $3}' | paste -sd, -
}
# The single_job allocation for this case, if queued or running.
driver_jobs() {
  squeue -h -u "${USER}" -n "${DRIVER_JOB_NAME}" -o '%i:%T' | paste -sd, -
}

user_eam_nl() {
  ex user_nl_eam <<EOF
g/^ *ncdata *=/s@=.*@= "$1"@
g/^ *inithist *=/s@=.*@= 'NONE'@
wq
EOF
}

user_elm_nl() {
  ex user_nl_elm <<EOF
g/^ *finidat *=/s@=.*@= "$1"@
g/^ *check_finidat_year_consistency *=/s@=.*@= .false.@
g/^ *check_dynpft_consistency *=/s@=.*@= .false.@
g/^ *check_finidat_fsurdat_consistency *=/s@=.*@= .false.@
g/^ *check_finidat_pct_consistency *=/s@=.*@= .false.@
wq
EOF
}

user_mosart_nl() {
  ex user_nl_mosart <<EOF
g/^ *finidat_rtm *=/s@=.*@= "$1"@
wq
EOF
}

user_mpassi_nl() {
  ex user_nl_mpassi <<EOF
g/^ *config_start_time *=/s@=.*@= '$1_$2'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_initial_condition_type *=/s@=.*@= 'restart'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ice'@
wq
EOF
}

user_mpaso_nl() {
  ex user_nl_mpaso <<EOF
g/^ *config_start_time *=/s@=.*@= '$1_$2'@
g/^ *config_calendar_type *=/s@=.*@= 'gregorian'@
g/^ *config_do_restart *=/s@=.*@= .true.@
g/^ *config_output_reference_time *=/s@=.*@= '$3'@
g/^ *config_restart_timestamp_name *=/s@=.*@= 'rpointer.ocn'@
wq
EOF
}

# Point the MPAS initial-condition input stream at the staged restart file
# (same edit as the DA cycle worker). It is read only at the hybrid start.
patch_mpas_streams() {
  local run_dir="$1" streams="$2" src_dir="$3" ic_file="$4" case_dir="$5" line
  for line in "${@:6}"; do
    sed -i "${line}s#filename_template=.*#filename_template=\"${ic_file}\"#" "${run_dir}/${streams}"
  done
  cp "${run_dir}/${streams}" "${case_dir}/SourceMods/${src_dir}/"
}

stage_file_atomic() {
  local source_file="$1" destination_file="$2" tmp_file="$2.tmp.$$"
  [[ -s "${source_file}" ]] || fail "restart source is missing or empty: ${source_file}"
  if [[ "${source_file}" == *.nc ]]; then
    ncdump -h "${source_file}" >/dev/null 2>&1 || fail "restart source is invalid NetCDF: ${source_file}"
  fi
  cp -p "${source_file}" "${tmp_file}"
  mv -f "${tmp_file}" "${destination_file}"
}

# PE layout, run length, and queue for this member. The run goes to the end
# time in one job; the only restart is written there.
configure_member_case() {
  local archive="$1" run_hours="$2"
  ./xmlchange NTASKS="${MEMBER_TASKS}",NTHRDS=1,ROOTPE=0
  ./xmlchange MAX_MPITASKS_PER_NODE="${my_task_per_node}",MAX_TASKS_PER_NODE="${my_task_per_node}"
  ./case.setup --reset
  # case.setup --reset always clears BUILD_COMPLETE. The executable has no
  # compile-time task count and is reused, as Step 1 does for every member.
  [[ -x "$(./xmlquery EXEROOT --value)/e3sm.exe" ]] || fail "missing member executable; rerun Step 1"
  ./xmlchange BUILD_COMPLETE=TRUE
  ./xmlchange STOP_OPTION=nhours,STOP_N="${run_hours}"
  ./xmlchange REST_OPTION=nhours,REST_N="${run_hours}"
  ./xmlchange RESUBMIT=0
  ./xmlchange DOUT_S=TRUE,DOUT_S_ROOT="${archive}"
  ./xmlchange JOB_QUEUE="${JOB_QUEUE}",JOB_WALLCLOCK_TIME="${JOB_WALLTIME}" --subgroup case.run
  ./xmlchange JOB_QUEUE="${JOB_QUEUE}" --subgroup case.st_archive
  ./xmlchange PROJECT="${my_project}",CHARGE_ACCOUNT="${my_project}"
}

# Hybrid start from the Step 3 perturbed eam.i file and Step 2 restarts.
prepare_member_start() (
  local enstr="$1" run_hours="$2"
  local case_name="${my_casename}.${enstr}" case_dir run_dir archive ref_dir ymd tod hour ocn_time file
  case_dir=$(member_case_dir "${enstr}"); run_dir=$(member_run_dir "${enstr}"); archive=$(member_archive "${enstr}")
  ymd="${my_e3sm_start_date}"; tod="${my_e3sm_start_tod}"
  hour=$(printf '%02d' $((10#${tod} / 3600)))
  ocn_time="${ymd}_${hour}:00:00"
  ref_dir="${archive}/rest/${ymd}-${tod}"
  set -e

  cd "${case_dir}"
  # Reference-case settings first, so case.setup validates the member's own
  # restart directory rather than a stale RUN_REFDIR.
  ./xmlchange RUN_TYPE=hybrid,CONTINUE_RUN=FALSE,GET_REFCASE=FALSE
  ./xmlchange RUN_STARTDATE="${ymd}",START_TOD="${tod}"
  ./xmlchange RUN_REFCASE="${case_name}",RUN_REFDATE="${ymd}",RUN_REFTOD="${tod}",RUN_REFDIR="${ref_dir}"

  # Stage restarts before case.setup: MPAS-SI buildnml requires its hybrid
  # restart file in the run directory when namelists are generated.
  cd "${run_dir}"
  rm -f rpointer*
  echo "${case_name}.eam.r.${ymd}-${tod}.nc" > rpointer.atm
  echo "./${case_name}.elm.r.${ymd}-${tod}.nc" > rpointer.lnd
  echo "./${case_name}.mosart.r.${ymd}-${tod}.nc" > rpointer.rof
  echo "${case_name}.cpl.r.${ymd}-${tod}.nc" > rpointer.drv
  echo "${ocn_time}" > rpointer.ice
  if [[ "${my_runtype}" == "AMIP" ]]; then
    printf '%s\n%s\n' "${case_name}.docn.r.${ymd}_${tod}.nc" "${case_name}.docn.rs1.${ymd}_${tod}.bin" > rpointer.ocn
    for file in "${case_name}.docn.r.${ymd}_${tod}.nc" "${case_name}.docn.rs1.${ymd}_${tod}.bin"; do
      if [[ -s "${ref_dir}/${file}" ]]; then stage_file_atomic "${ref_dir}/${file}" "${file}"; fi
    done
  else
    echo "${ocn_time}" > rpointer.ocn
    stage_file_atomic "${ref_dir}/${case_name}.mpaso.rst.${ymd}_${tod}.nc" "${case_name}.mpaso.rst.${ymd}_${tod}.nc"
  fi
  # EAM starts from ncdata in hybrid mode; the eam.r placeholder is not staged.
  for file in "${case_name}.elm.r.${ymd}-${tod}.nc" "${case_name}.mosart.r.${ymd}-${tod}.nc" \
              "${case_name}.cpl.r.${ymd}-${tod}.nc" "${case_name}.mpassi.rst.${ymd}_${tod}.nc"; do
    stage_file_atomic "${ref_dir}/${file}" "${file}"
  done

  cd "${case_dir}"
  configure_member_case "${archive}" "${run_hours}"
  user_eam_nl "${ref_dir}/${case_name}.eam.i.${ymd}-${tod}.nc"
  user_elm_nl "${ref_dir}/${case_name}.elm.r.${ymd}-${tod}.nc"
  user_mosart_nl "${ref_dir}/${case_name}.mosart.r.${ymd}-${tod}.nc"
  user_mpassi_nl "${ymd}" "${hour}"
  if [[ "${my_runtype}" == "Full-CPL" ]]; then user_mpaso_nl "${ymd}" "${hour}" "${ocn_time}"; fi
  ./preview_namelists
  patch_mpas_streams "${run_dir}" streams.seaice src.mpassi "${ref_dir}/${case_name}.mpassi.rst.${ymd}_${tod}.nc" "${case_dir}" 11 38
  if [[ "${my_runtype}" == "Full-CPL" ]]; then
    patch_mpas_streams "${run_dir}" streams.ocean src.mpaso "${ref_dir}/${case_name}.mpaso.rst.${ymd}_${tod}.nc" "${case_dir}" 12
  fi
  echo "Prepared hybrid start for ${enstr} at ${ymd}-${tod}"
)

# Resume an interrupted member from its latest archived restart.
prepare_member_continue() (
  local enstr="$1" run_hours="$2" state="$3"
  local case_name="${my_casename}.${enstr}" case_dir run_dir archive rest_dir file atm_rst
  case_dir=$(member_case_dir "${enstr}"); run_dir=$(member_run_dir "${enstr}"); archive=$(member_archive "${enstr}")
  rest_dir="${archive}/rest/${state}"
  atm_rst="${case_name}.eam.r.${state}.nc"
  set -e

  # st_archive leaves the latest restart set and rpointers in the run
  # directory; restage from the archive only when that state is not intact.
  if [[ "$(head -n 1 "${run_dir}/rpointer.atm" 2>/dev/null)" != "${atm_rst}" || ! -s "${run_dir}/${atm_rst}" ]]; then
    echo "Restaging ${enstr} restart set from ${rest_dir}"
    cd "${run_dir}"
    rm -f rpointer*
    for file in "${rest_dir}"/*; do
      case "$(basename "${file}")" in
        *.h[0-9].*|*.h[0-9][0-9].*) continue ;;   # skip in-progress history copies
      esac
      stage_file_atomic "${file}" "$(basename "${file}")"
    done
    echo "${atm_rst}" > rpointer.atm
    echo "./${case_name}.elm.r.${state}.nc" > rpointer.lnd
    echo "./${case_name}.mosart.r.${state}.nc" > rpointer.rof
    echo "${case_name}.cpl.r.${state}.nc" > rpointer.drv
    echo "${state:0:10}_$(printf '%02d' $((10#${state:11:5} / 3600))):00:00" > rpointer.ice
    if [[ "${my_runtype}" == "AMIP" ]]; then
      printf '%s\n%s\n' "${case_name}.docn.r.${state:0:10}_${state:11:5}.nc" "${case_name}.docn.rs1.${state:0:10}_${state:11:5}.bin" > rpointer.ocn
    else
      cp rpointer.ice rpointer.ocn
    fi
  fi

  # MPAS-SI and MOSART buildnml require the hybrid reference restarts in the
  # run directory even with CONTINUE_RUN=TRUE; st_archive removed them.
  cd "${run_dir}"
  for file in "${case_name}.mpassi.rst.${my_e3sm_start_date}_${my_e3sm_start_tod}.nc" \
              "${case_name}.mosart.r.${START_STAMP}.nc"; do
    [[ -s "${file}" ]] || stage_file_atomic "${archive}/rest/${START_STAMP}/${file}" "${file}"
  done

  cd "${case_dir}"
  configure_member_case "${archive}" "${run_hours}"
  ./xmlchange CONTINUE_RUN=TRUE
  echo "Prepared continuation for ${enstr} from ${state}"
)

# Configure one member to run from its current state to the end time.
prepare_member() {
  local enstr="$1" state="$2" run_hours
  run_hours=$(( (END_EPOCH - $(stamp_to_epoch "${state}")) / 3600 ))
  if [[ "${state}" == "${START_STAMP}" ]]; then
    prepare_member_start "${enstr}" "${run_hours}"
  else
    prepare_member_continue "${enstr}" "${run_hours}" "${state}"
  fi
}

# Members that have not reached the end time, as "ENxx:state" words. Callers
# assign the output first (PENDING=$(pending_members) || exit 1) so a failure
# here stops the driver.
pending_members() {
  local i enstr state
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}")
    [[ -x "$(member_case_dir "${enstr}")/xmlchange" ]] || fail "missing member case from Step 1: $(member_case_dir "${enstr}")"
    state=$(member_state_stamp "${enstr}")
    if [[ "${state}" != "${END_STAMP}" ]]; then echo "${enstr}:${state}"; fi
  done
}

capture_member_failure() {
  local enstr="$1" run_dir log diag
  run_dir=$(member_run_dir "${enstr}")
  diag="${LOG_DIR}/failure.${enstr}.${SLURM_JOB_ID:-$$}.free.log"
  log=$(find "${run_dir}" -maxdepth 1 -type f -name 'e3sm.log.*' -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk 'NR==1 {print $2}')
  {
    echo "Failure diagnostics for ${enstr} at $(date '+%Y-%m-%d %H:%M:%S')"
    tail -n 40 "$(member_case_dir "${enstr}")/CaseStatus" 2>/dev/null || true
    if [[ -n "${log}" ]]; then
      echo "===== ${log} ====="
      grep -nEi 'error|fatal|abort|segmentation|killed|bad termination' "${log}" | tail -n 60 || true
      tail -n 100 "${log}"
    fi
  } > "${diag}" 2>&1
  echo "Captured ${enstr} failure diagnostics: ${diag}"
}

#------------------------------------------------------------------------------
# Modes
#------------------------------------------------------------------------------
show_status() {
  local i enstr state jobs done_count=0
  printf '%-5s %-17s %-6s %s\n' member state done jobs
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}")
    state=$(member_state_stamp "${enstr}")
    jobs=$(member_jobs "${enstr}")
    if [[ "${state}" == "${END_STAMP}" ]]; then done_count=$((done_count + 1)); fi
    printf '%-5s %-17s %-6s %s\n' "${enstr}" "${state}" "$([[ "${state}" == "${END_STAMP}" ]] && echo yes || echo no)" "${jobs:--}"
  done
  echo "${done_count}/${my_ensnum} members at ${END_STAMP}; layout ${LAYOUT}, ${MEMBER_NODES} nodes per member (${JOB_QUEUE}, ${JOB_WALLTIME})"
  jobs=$(driver_jobs)
  [[ -z "${jobs}" ]] || echo "single_job allocation: ${jobs}"
  [[ -s "${COMPLETE_MARKER}" ]] && echo "Completion record: ${COMPLETE_MARKER}"
  return 0
}

# member_jobs: one CIME batch job per pending member.
submit_member_jobs() {
  local entry enstr state jobs setup_log submitted=0
  local -a failed=()
  for entry in ${PENDING}; do
    enstr=${entry%%:*}; state=${entry#*:}
    jobs=$(member_jobs "${enstr}")
    if [[ -n "${jobs}" ]]; then
      echo "${enstr}: jobs active (${jobs}); not resubmitting"
      continue
    fi
    setup_log="${LOG_DIR}/free_member_setup.${enstr}.$(date +%Y%m%d%H%M%S).log"
    echo "${enstr}: ${state} -> ${END_STAMP}; setup log ${setup_log}"
    # A failing member is reported and skipped so the others still submit.
    prepare_member "${enstr}" "${state}" > "${setup_log}" 2>&1 || { echo "ERROR: setup failed for ${enstr}; see ${setup_log}" >&2; failed+=("${enstr}"); continue; }
    if ! (cd "$(member_case_dir "${enstr}")" && ./case.submit) >> "${setup_log}" 2>&1; then
      echo "ERROR: case.submit failed for ${enstr}; see ${setup_log}" >&2
      failed+=("${enstr}")
      continue
    fi
    grep -E "Submitted job" "${setup_log}" | tail -n 2 || true
    submitted=$((submitted + 1))
  done
  echo "Submitted ${submitted} member job(s). Check progress with: $0 status"
  (( ${#failed[@]} == 0 )) || fail "${#failed[@]} member(s) not submitted: ${failed[*]}; fix and rerun $0 (submitted members are skipped)"
}

# single_job: submit one allocation that runs every pending member.
submit_single_job() {
  local i jobs job_id
  jobs=$(driver_jobs)
  [[ -z "${jobs}" ]] || { echo "single_job allocation already queued or running: ${jobs}"; return 0; }
  for ((i=1; i<=my_ensnum; i++)); do
    jobs=$(member_jobs "$(member_name "${i}")")
    [[ -z "${jobs}" ]] || fail "$(member_name "${i}") still has member jobs (${jobs}); wait for them or cancel them first"
  done
  job_id=$(sbatch --parsable --job-name="${DRIVER_JOB_NAME}" --account="${my_project}" \
    --partition="${JOB_QUEUE}" --time="${JOB_WALLTIME}" --nodes="${my_job_nnodes}" --exclusive \
    --output="${my_wkdir}/e3sm_ctrl_free.%j" "${SCRIPT_PATH}" run "${my_wkdir}")
  echo "Submitted single_job allocation ${job_id%%;*} (${my_job_nnodes} nodes, ${JOB_QUEUE}, ${JOB_WALLTIME}); log ${my_wkdir}/e3sm_ctrl_free.${job_id%%;*}"
}

submit_run() {
  [[ ! -s "${COMPLETE_MARKER}" ]] || { echo "Free run already complete: ${COMPLETE_MARKER}"; return 0; }
  PENDING=$(pending_members) || exit 1
  if [[ -z "${PENDING}" ]]; then
    echo "All members reached ${END_STAMP}; run: $0 finalize"
    return 0
  fi
  if [[ "${LAYOUT}" == "member_jobs" ]]; then submit_member_jobs; else submit_single_job; fi
}

# Inside the single_job allocation: run all pending members in parallel waves.
run_single_job() {
  local entry enstr idx case_name archived
  local -a pending=() pids=() names=() failed=()
  local -A member_state=()
  trap 'echo "ERROR: received a termination signal"; jobs -pr | xargs -r kill -TERM 2>/dev/null; wait; exit 143' INT TERM
  echo "== Free run ${my_casename}: ${START_STAMP} -> ${END_STAMP}, job ${SLURM_JOB_ID} =="
  date
  [[ ! -s "${COMPLETE_MARKER}" ]] || { echo "Free run already complete: ${COMPLETE_MARKER}"; return 0; }
  [[ "${SLURM_JOB_NUM_NODES:-${my_job_nnodes}}" == "${my_job_nnodes}" ]] || fail "allocation has ${SLURM_JOB_NUM_NODES} nodes; configuration requires ${my_job_nnodes}"
  PENDING=$(pending_members) || exit 1
  for entry in ${PENDING}; do
    enstr=${entry%%:*}; member_state[${enstr}]=${entry#*:}; pending+=("${enstr}")
  done

  for enstr in "${pending[@]}"; do
    echo "${enstr}: ${member_state[${enstr}]} -> ${END_STAMP}"
    prepare_member "${enstr}" "${member_state[${enstr}]}" > "${LOG_DIR}/free_setup.${enstr}.${SLURM_JOB_ID}.log" 2>&1 &
    pids+=("$!"); names+=("${enstr}")
    if (( ${#pids[@]} >= MAX_PARALLEL_SETUP )) || [[ "${enstr}" == "${pending[-1]}" ]]; then
      for idx in "${!pids[@]}"; do
        wait "${pids[$idx]}" || fail "setup failed for ${names[$idx]}; see ${LOG_DIR}/free_setup.${names[$idx]}.${SLURM_JOB_ID}.log"
      done
      pids=(); names=()
    fi
  done

  for enstr in "${pending[@]}"; do
    (cd "$(member_case_dir "${enstr}")" && ./case.submit --no-batch) \
      > "${LOG_DIR}/e3sm.${enstr}.o${SLURM_JOB_ID}.free.log" 2>&1 &
    pids+=("$!"); names+=("${enstr}")
    if (( ${#pids[@]} >= BATCH_SIZE )) || [[ "${enstr}" == "${pending[-1]}" ]]; then
      echo "Waiting for members: ${names[*]}"
      for idx in "${!pids[@]}"; do
        wait "${pids[$idx]}" || echo "WARNING: case.submit returned non-zero for ${names[$idx]}"
      done
      pids=(); names=()
    fi
  done

  # A member succeeded only if its end-time restart reached the archive.
  for enstr in "${pending[@]}"; do
    case_name="${my_casename}.${enstr}"
    archived="$(member_archive "${enstr}")/rest/${END_STAMP}/${case_name}.eam.r.${END_STAMP}.nc"
    if [[ ! -s "${archived}" && -s "$(member_run_dir "${enstr}")/${case_name}.eam.r.${END_STAMP}.nc" ]]; then
      echo "Running short-term archive for ${enstr}"
      (cd "$(member_case_dir "${enstr}")" && srun -N 1 -n 1 ./case.st_archive) || true
    fi
    if [[ -s "${archived}" ]]; then
      echo "${enstr} reached ${END_STAMP}"
    else
      capture_member_failure "${enstr}"
      failed+=("${enstr}")
    fi
  done
  (( ${#failed[@]} == 0 )) || fail "run failed for: ${failed[*]}; fix and rerun '$0 submit' to resume from the last archived restart"
  finalize_run
  date
}

# With the template's 6-hourly streams (mfilt=4), EAM/ELM open a new history
# file at 00Z each day, plus one last file at the end time.
FREE_HIST_STREAMS=("atm:eam.h0" "atm:eam.h1" "atm:eam.h2" "lnd:elm.h0" "lnd:elm.h1")
finalize_run() {
  local i enstr case_name archive epoch stamp missing=0 dir dir_epoch stream file jobs
  [[ ! -s "${COMPLETE_MARKER}" ]] || { echo "Free run already complete: ${COMPLETE_MARKER}"; return 0; }
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}")
    [[ "$(member_state_stamp "${enstr}")" == "${END_STAMP}" ]] || fail "${enstr} has not reached ${END_STAMP}; run '$0 status'"
    jobs=$(member_jobs "${enstr}")
    [[ -z "${jobs}" ]] || fail "${enstr} still has active jobs (${jobs}); wait for short-term archiving to finish"
  done
  echo "All members reached ${END_STAMP}; verifying archived history"
  for ((i=1; i<=my_ensnum; i++)); do
    enstr=$(member_name "${i}"); case_name="${my_casename}.${enstr}"; archive=$(member_archive "${enstr}")
    for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=86400)); do
      stamp=$(epoch_to_stamp "${epoch}")
      for stream in "${FREE_HIST_STREAMS[@]}"; do
        file="${archive}/${stream%%:*}/hist/${case_name}.${stream#*:}.${stamp}.nc"
        [[ -s "${file}" ]] || { echo "ERROR: missing history file: ${file}" >&2; missing=$((missing + 1)); }
      done
    done
  done
  (( missing == 0 )) || fail "${missing} history file(s) missing; restart sets were kept for recovery"

  if [[ "${KEEP_INTERMEDIATE}" == "FALSE" ]]; then
    for ((i=1; i<=my_ensnum; i++)); do
      archive=$(member_archive "$(member_name "${i}")")
      for dir in "${archive}"/rest/????-??-??-?????; do
        [[ -d "${dir}" ]] || continue
        dir_epoch=$(stamp_to_epoch "$(basename "${dir}")") || continue
        if (( dir_epoch > START_EPOCH && dir_epoch < END_EPOCH )); then
          echo "Removing intermediate restart set: ${dir}"
          rm -rf -- "${dir}"
        fi
      done
    done
  fi

  printf 'valid_time=%s\nstart_time=%s\ncase=%s\nensemble_size=%s\nslurm_job_id=%s\ncompleted_at=%s\n' \
    "${END_STAMP}" "${START_STAMP}" "${my_casename}" "${my_ensnum}" "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${COMPLETE_MARKER}.tmp"
  mv -f "${COMPLETE_MARKER}.tmp" "${COMPLETE_MARKER}"
  echo "Wrote completion record: ${COMPLETE_MARKER}"
}

case "${MODE}" in
  submit) submit_run ;;
  status) show_status ;;
  finalize) finalize_run ;;
  run)
    [[ "${LAYOUT}" == "single_job" ]] || fail "run mode requires my_free_layout=single_job"
    run_single_job
    ;;
esac
