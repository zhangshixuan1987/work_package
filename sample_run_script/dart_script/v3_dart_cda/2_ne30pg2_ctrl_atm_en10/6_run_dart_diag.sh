#!/bin/bash
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_diag
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --output=e3sm_dart_diag.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail

# Step 6 modes:
#   seq      - preprocess cycle observation sequences
#   obs      - run standard observation-space diagnostics
#   common   - run diagnostics for observations common across cycles
#   closest  - independently select the ensemble member closest to the mean
#   all      - run seq, obs, and common; also run closest when enabled below
# Each mode has its own range-specific completion record and can be restarted
# independently. A valid record is skipped unless DIAG_FORCE=TRUE.
# submit like sbatch --export=ALL,DIAG_MODE=common 6_run_dart_diag.sh

DIAG_MODE="${DIAG_MODE:-all}"
# Leave empty to use my_eam_dart_diag_start and my_eam_dart_diag_end.
DIAG_START="${DIAG_START:-}"
DIAG_END="${DIAG_END:-}"
# TRUE bypasses valid completion records and reruns every selected mode.
DIAG_FORCE="${DIAG_FORCE:-TRUE}"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"
CONFIG_FILE="${WORK_DIR}/create_and_setup_case.sh"
[[ -r "${CONFIG_FILE}" ]] || fail "missing workflow configuration: ${CONFIG_FILE}"
source "${CONFIG_FILE}"
mkdir -p "${my_log_dir}" "${my_status_dir}" "${my_lock_dir}"

for cmd in awk bash date flock mkdir mv readlink; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done

DIAG_START="${DIAG_START:-${my_eam_dart_diag_start:-}}"
DIAG_END="${DIAG_END:-${my_eam_dart_diag_end:-}}"
DIAG_TASKS="${DIAG_TASKS:-24}"
CLOSEST_STEP_HOURS="${my_eam_dart_diag_closest_hours:-24}"
[[ "${DIAG_MODE}" =~ ^(seq|obs|common|closest|all)$ ]] || fail "DIAG_MODE must be seq, obs, common, closest, or all; got: ${DIAG_MODE}"
[[ "${DIAG_TASKS}" =~ ^[1-9][0-9]*$ ]] || fail "DIAG_TASKS must be a positive integer"
[[ "${CLOSEST_STEP_HOURS}" =~ ^[1-9][0-9]*$ ]] || fail "my_eam_dart_diag_closest_hours must be positive"
DIAG_FORCE="${DIAG_FORCE^^}"
[[ "${DIAG_FORCE}" == "TRUE" || "${DIAG_FORCE}" == "FALSE" ]] || fail "DIAG_FORCE must be TRUE or FALSE; got: ${DIAG_FORCE}"
if [[ "${DIAG_FORCE}" == "TRUE" ]]; then
  echo "Force mode enabled: requested diagnostics will be reprocessed"
fi

stamp_to_epoch() {
  local stamp="$1" ymd tod hour minute second
  [[ "${stamp}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-[0-9]{5}$ ]] || return 1
  ymd="${stamp%-*}"
  tod="${stamp##*-}"
  (( 10#${tod} < 86400 )) || return 1
  hour=$((10#${tod} / 3600))
  minute=$(((10#${tod} % 3600) / 60))
  second=$((10#${tod} % 60))
  date -u -d "${ymd} $(printf '%02d:%02d:%02d' "${hour}" "${minute}" "${second}")" +%s
}

epoch_to_stamp() {
  local epoch="$1" ymd hour minute second tod
  read -r ymd hour minute second < <(date -u -d "@${epoch}" '+%Y-%m-%d %H %M %S')
  tod=$((10#${hour} * 3600 + 10#${minute} * 60 + 10#${second}))
  printf '%s-%05d\n' "${ymd}" "${tod}"
}

START_EPOCH=$(stamp_to_epoch "${DIAG_START}") || fail "invalid DIAG_START: ${DIAG_START}"
END_EPOCH=$(stamp_to_epoch "${DIAG_END}") || fail "invalid DIAG_END: ${DIAG_END}"
EXPERIMENT_START_EPOCH=$(stamp_to_epoch "${my_e3sm_start_date}-${my_e3sm_start_tod}") || fail "invalid configured E3SM start time"
(( START_EPOCH <= END_EPOCH )) || fail "DIAG_START is after DIAG_END"
(( START_EPOCH >= EXPERIMENT_START_EPOCH )) || fail "DIAG_START precedes the configured E3SM start time"
[[ "${my_eam_dart_cycle_hours:-}" =~ ^[1-9][0-9]*$ ]] || fail "my_eam_dart_cycle_hours must be positive"
STEP_SECONDS=$((my_eam_dart_cycle_hours * 3600))
CLOSEST_STEP_SECONDS=$((CLOSEST_STEP_HOURS * 3600))
(( (END_EPOCH - START_EPOCH) % STEP_SECONDS == 0 )) || fail "diagnostic range is not aligned to the ${my_eam_dart_cycle_hours}-hour EAM DA interval"
(( (END_EPOCH - EXPERIMENT_START_EPOCH) % STEP_SECONDS == 0 )) || fail "DIAG_END is not aligned to the EAM DA cadence from the experiment start"
EAM_DA_COMPLETED_CYCLES=$(((END_EPOCH - EXPERIMENT_START_EPOCH) / STEP_SECONDS))

STATUS_DIR="${my_status_dir}"
mkdir -p "${STATUS_DIR}"
exec 8>"${my_lock_dir}/step6_diagnostics.lock"
flock -n 8 || fail "another Step 6 diagnostic driver is active"

validate_obs_input_range() {
  local epoch stamp cycle_dir obs_seq
  for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=STEP_SECONDS)); do
    stamp=$(epoch_to_stamp "${epoch}")
    cycle_dir="${my_eam_dart_run_dir}/${stamp}"
    obs_seq="${cycle_dir}/${my_casename}.dart.e.eam_obs_seq_final.${stamp}"
    [[ -s "${obs_seq}" ]] || fail "missing or empty diagnostic input: ${obs_seq}"
  done
}

validate_closest_input_range() {
  local epoch stamp member member_name case_name restart
  local -a restart_matches
  [[ "${my_eam_dart_diag_run_closest_member}" == ".true." ]] || fail \
    "DIAG_MODE=closest requires my_eam_dart_diag_run_closest_member=.true."
  (( (END_EPOCH - START_EPOCH) % CLOSEST_STEP_SECONDS == 0 )) || fail \
    "closest diagnostic range is not aligned to the ${CLOSEST_STEP_HOURS}-hour cadence"
  for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=CLOSEST_STEP_SECONDS)); do
    stamp=$(epoch_to_stamp "${epoch}")
    for ((member=1; member<=my_ensnum; member++)); do
      member_name=$(printf 'EN%02d' "${member}")
      case_name="${my_casename}.${member_name}"
      restart="${my_modeldir}/${member_name}/archive/rest/${stamp}/${case_name}.*eam.i*${stamp}.nc"
      restart_matches=( ${restart} )
      (( ${#restart_matches[@]} == 1 )) || fail "expected one closest-member input matching: ${restart}"
      [[ -s "${restart_matches[0]}" ]] || fail "missing or empty closest-member input: ${restart_matches[0]}"
    done
  done
}

echo "Preflighting Step 6 input range"
case "${DIAG_MODE}" in
  seq|obs|common) validate_obs_input_range ;;
  closest) validate_closest_input_range ;;
  all)
    validate_obs_input_range
    if [[ "${my_eam_dart_diag_run_closest_member}" == ".true." ]]; then
      validate_closest_input_range
    fi
    ;;
esac
echo "Step 6 preflight completed successfully"

run_mode() {
  local mode="$1" worker status_file progress_file tmp_file
  worker="${my_workflow_lib}/diagnostics/dart_diag_${mode}_worker.sh"
  [[ -x "${worker}" ]] || fail "missing diagnostic worker: ${worker}"
  status_file="${STATUS_DIR}/diagnostic_complete.${mode}.${DIAG_START}-${DIAG_END}"
  progress_file="${STATUS_DIR}/diagnostic_in_progress.${mode}.${DIAG_START}-${DIAG_END}"

  if [[ "${DIAG_FORCE}" != "TRUE" ]] &&
     [[ -s "${status_file}" ]] &&
     [[ "$(awk -F= '$1 == "case" {print $2}' "${status_file}")" == "${my_casename}" ]] &&
     [[ "$(awk -F= '$1 == "ensemble_size" {print $2}' "${status_file}")" == "${my_ensnum}" ]]; then
    echo "Diagnostic mode ${mode} already complete: ${status_file}"
    return 0
  fi

  rm -f -- "${status_file}"
  printf 'mode=%s\nstart=%s\nend=%s\ncase=%s\nensemble_size=%s\njob_id=%s\nstarted_at=%s\n' \
    "${mode}" "${DIAG_START}" "${DIAG_END}" "${my_casename}" "${my_ensnum}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${progress_file}"

  echo "== Running Step 6 mode ${mode}: ${DIAG_START} through ${DIAG_END} =="
  STEP6_DRIVER_ACTIVE=TRUE MY_DART_DIAG_START="${DIAG_START}" MY_DART_DIAG_END="${DIAG_END}" EAM_DA_COMPLETED_CYCLES="${EAM_DA_COMPLETED_CYCLES}" DIAG_TASKS="${DIAG_TASKS}" "${worker}"

  tmp_file="${status_file}.tmp.${SLURM_JOB_ID:-$$}"
  printf 'mode=%s\nstart=%s\nend=%s\ncase=%s\nensemble_size=%s\njob_id=%s\ncompleted_at=%s\n' \
    "${mode}" "${DIAG_START}" "${DIAG_END}" "${my_casename}" "${my_ensnum}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${tmp_file}"
  mv -f "${tmp_file}" "${status_file}"
  rm -f -- "${progress_file}"
  echo "Completed Step 6 mode ${mode}: ${status_file}"
}

case "${DIAG_MODE}" in
  seq) run_mode seq ;;
  obs) run_mode obs ;;
  common) run_mode common ;;
  closest) run_mode closest ;;
  all)
    run_mode seq
    run_mode obs
    run_mode common
    if [[ "${my_eam_dart_diag_run_closest_member}" == ".true." ]]; then
      run_mode closest
    else
      echo "Skipping Step 6 mode closest: my_eam_dart_diag_run_closest_member=${my_eam_dart_diag_run_closest_member}"
    fi
    ;;
esac
