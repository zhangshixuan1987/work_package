#!/bin/bash
#SBATCH --account=esmd
#SBATCH --time=24:00:00
#SBATCH --partition=slurm
#SBATCH --job-name=e3sm_dart_interp_init
#SBATCH --nodes=1
#SBATCH --output=e3sm_dart_interp_init.%j
#SBATCH --exclusive
#SBATCH --no-kill
#SBATCH --requeue

set -Eeuo pipefail

# User settings. Step 8 processes midnight restart files, one calendar day at a time.
INTERP_START="2011-12-15-00000"
INTERP_END="2011-12-28-00000"
MAX_CONCURRENT_WORKERS=4
OVERWRITE_EXISTING="FALSE"
PROCESS_EAM="TRUE"
PROCESS_ELM="TRUE"
REQUIRE_UPSTREAM_MARKERS="TRUE"

readonly EAM_VARIABLES="lat_d,lon_d,U,V,T,Q,PS"
readonly ELM_VARIABLES="cols1d_ityp,cols1d_active,cols1d_ityplun,cols1d_landunit_index,cols1d_topounit_index,cols1d_gridcell_index,cols1d_jxy,cols1d_ixy,cols1d_lat,cols1d_lon,DZSNO,ZSNO,ZISNO,T_SOISNO,H2OSOI_LIQ,H2OSOI_ICE"

fail() { echo "ERROR: $*" >&2; exit 1; }
normalize_bool() {
  local name="$1" value="${2^^}"
  [[ "${value}" == "TRUE" || "${value}" == "FALSE" ]] || fail "${name} must be TRUE or FALSE, got: $2"
  printf '%s' "${value}"
}
stamp_to_epoch() {
  local stamp="$1"
  [[ "${stamp}" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}-00000$ ]] || return 1
  date -u -d "${stamp%-*} 00:00:00" +%s
}
epoch_to_stamp() { date -u -d "@$1" '+%Y-%m-%d-00000'; }
record_value() { awk -F= -v key="$2" '$1 == key {sub(/^[^=]*=/, ""); print; exit}' "$1"; }

if [[ -n "${SLURM_SUBMIT_DIR:-}" ]]; then
  WORK_DIR=$(readlink -f "${SLURM_SUBMIT_DIR}") || fail "cannot resolve SLURM_SUBMIT_DIR"
else
  SCRIPT_PATH=$(readlink -f "${BASH_SOURCE[0]}") || fail "cannot resolve script path"
  WORK_DIR=$(dirname "${SCRIPT_PATH}")
fi
cd "${WORK_DIR}"
[[ -r create_and_setup_case.sh ]] || fail "missing create_and_setup_case.sh"
source ./create_and_setup_case.sh
[[ -n "${my_conda_setup_file:-}" ]] || fail "my_conda_setup_file is unset"
[[ -r "${my_conda_setup_file}" ]] || fail "configured Conda setup is not readable: ${my_conda_setup_file}"
[[ -n "${my_analysis_conda_env:-}" ]] || fail "my_analysis_conda_env is unset"
echo "Activating configured analysis environment: ${my_analysis_conda_env}"
source "${my_conda_setup_file}"
conda activate "${my_analysis_conda_env}" || fail "could not activate Conda environment: ${my_analysis_conda_env}"

for cmd in awk date flock mktemp mv ncdump ncks ncrename readlink rm; do
  command -v "${cmd}" >/dev/null 2>&1 || fail "required command not found: ${cmd}"
done
for setting in OVERWRITE_EXISTING PROCESS_EAM PROCESS_ELM REQUIRE_UPSTREAM_MARKERS; do
  printf -v "${setting}" '%s' "$(normalize_bool "${setting}" "${!setting}")"
done
[[ "${MAX_CONCURRENT_WORKERS}" =~ ^[1-9][0-9]*$ ]] || fail "MAX_CONCURRENT_WORKERS must be a positive integer"
[[ "${PROCESS_EAM}" == "TRUE" || "${PROCESS_ELM}" == "TRUE" ]] || fail "at least one component must be enabled"
[[ "${my_ensnum}" =~ ^[1-9][0-9]*$ ]] || fail "invalid ensemble size: ${my_ensnum}"

START_EPOCH=$(stamp_to_epoch "${INTERP_START}") || fail "INTERP_START must be a valid midnight timestamp"
END_EPOCH=$(stamp_to_epoch "${INTERP_END}") || fail "INTERP_END must be a valid midnight timestamp"
(( START_EPOCH <= END_EPOCH )) || fail "INTERP_START is after INTERP_END"

OUT_DIR="${my_dart_root}/post/init_rgd"
STATUS_DIR="${my_status_dir}"
mkdir -p "${OUT_DIR}" "${STATUS_DIR}" "${my_log_dir}" "${my_lock_dir}"

exec 8>"${my_lock_dir}/interp_init.lock"
flock -n 8 || fail "another Step 8 interpolation driver is active"

validate_upstream_marker() {
  local stamp="$1" marker marker_time
  if [[ "${stamp}" == "${my_e3sm_start_date}-${my_e3sm_start_tod}" ]]; then
    marker="${STATUS_DIR}/perturb_complete.${stamp}"
  else
    marker="${STATUS_DIR}/cycle_complete.${stamp}"
  fi
  if [[ ! -s "${marker}" ]]; then
    [[ "${REQUIRE_UPSTREAM_MARKERS}" == "FALSE" ]] || fail "missing upstream completion record: ${marker}"
    echo "WARNING: accepting legacy restart without completion record: ${stamp}" >&2
    return 0
  fi
  marker_time=$(record_value "${marker}" valid_time)
  [[ "${marker_time}" == "${stamp}" ]] || fail "invalid valid_time in upstream record: ${marker}"
}
preflight_date() {
  local stamp="$1" member case_name member_archive restart_path eam_in elm_in
  validate_upstream_marker "${stamp}"
  for ((preflight_i=1; preflight_i<=my_ensnum; preflight_i++)); do
    member=$(printf 'EN%02d' "${preflight_i}")
    case_name="${my_casename}.${member}"
    member_archive="${my_modeldir}/${member}/archive"
    restart_path="${member_archive}/rest/${stamp}"
    [[ -d "${restart_path}" ]] || fail "missing restart directory for ${member}: ${restart_path}"
    eam_in="${restart_path}/${case_name}.eam.i.${stamp}.nc"
    elm_in="${restart_path}/${case_name}.elm.r.${stamp}.nc"
    [[ "${PROCESS_EAM}" != "TRUE" || -s "${eam_in}" ]] || fail "missing EAM input: ${eam_in}"
    [[ "${PROCESS_ELM}" != "TRUE" || -s "${elm_in}" ]] || fail "missing ELM input: ${elm_in}"
    [[ "${PROCESS_EAM}" != "TRUE" ]] || ncks -m -v "${EAM_VARIABLES}" "${eam_in}" >/dev/null 2>&1 || fail "invalid EAM input or variables: ${eam_in}"
    [[ "${PROCESS_ELM}" != "TRUE" ]] || ncks -m -v "${ELM_VARIABLES}" "${elm_in}" >/dev/null 2>&1 || fail "invalid ELM input or variables: ${elm_in}"
  done
}


validate_output() {
  local component="$1" file="$2"
  [[ -s "${file}" ]] || return 1
  ncdump -h "${file}" >/dev/null 2>&1 || return 1
  case "${component}" in
    eam)
      ncks -m -v lat,lon,U,V,T,Q,PS "${file}" >/dev/null 2>&1
      ;;
    elm)
      ncks -m -v "${ELM_VARIABLES}" "${file}" >/dev/null 2>&1
      ;;
    *) return 1 ;;
  esac
}

validate_completion_record() {
  local record="$1" member="$2" stamp="$3" eam_out="$4" elm_out="$5"
  [[ -s "${record}" ]] || return 1
  [[ "$(record_value "${record}" valid_time)" == "${stamp}" ]] || return 1
  [[ "$(record_value "${record}" case)" == "${my_casename}" ]] || return 1
  [[ "$(record_value "${record}" member)" == "${member}" ]] || return 1
  [[ "$(record_value "${record}" ensemble_size)" == "${my_ensnum}" ]] || return 1
  [[ "$(record_value "${record}" process_eam)" == "${PROCESS_EAM}" ]] || return 1
  [[ "$(record_value "${record}" process_elm)" == "${PROCESS_ELM}" ]] || return 1
  [[ "${PROCESS_EAM}" != "TRUE" ]] || validate_output eam "${eam_out}" || return 1
  [[ "${PROCESS_ELM}" != "TRUE" ]] || validate_output elm "${elm_out}" || return 1
}

process_member() {
  local stamp="$1" member="$2" case_name member_archive restart_path date_out
  local eam_in elm_in eam_out elm_out record progress tmp_record tmpdir
  case_name="${my_casename}.${member}"
  member_archive="${my_modeldir}/${member}/archive"
  restart_path="${member_archive}/rest/${stamp}"
  date_out="${OUT_DIR}/${stamp}"
  eam_in="${restart_path}/${case_name}.eam.i.${stamp}.nc"
  elm_in="${restart_path}/${case_name}.elm.r.${stamp}.nc"
  eam_out="${date_out}/${case_name}.eam.i.${stamp}.nc"
  elm_out="${date_out}/${case_name}.elm.r.${stamp}.nc"
  record="${STATUS_DIR}/interp_init_complete.${member}.${stamp}"
  progress="${STATUS_DIR}/interp_init_in_progress.${member}.${stamp}"

  exec 9>"${my_lock_dir}/interp_init.${member}.${stamp}.lock"
  flock -n 9 || { echo "ERROR: interpolation already active for ${member} at ${stamp}" >&2; return 1; }

  if [[ "${OVERWRITE_EXISTING}" == "FALSE" ]] && validate_completion_record "${record}" "${member}" "${stamp}" "${eam_out}" "${elm_out}"; then
    echo "Skipping verified ${member} at ${stamp}"
    return 0
  fi

  [[ "${PROCESS_EAM}" != "TRUE" || -s "${eam_in}" ]] || { echo "ERROR: missing EAM input: ${eam_in}" >&2; return 1; }
  [[ "${PROCESS_ELM}" != "TRUE" || -s "${elm_in}" ]] || { echo "ERROR: missing ELM input: ${elm_in}" >&2; return 1; }
  [[ "${PROCESS_EAM}" != "TRUE" ]] || ncdump -h "${eam_in}" >/dev/null 2>&1 || { echo "ERROR: invalid EAM input: ${eam_in}" >&2; return 1; }
  [[ "${PROCESS_ELM}" != "TRUE" ]] || ncdump -h "${elm_in}" >/dev/null 2>&1 || { echo "ERROR: invalid ELM input: ${elm_in}" >&2; return 1; }

  mkdir -p "${date_out}"
  rm -f -- "${record}"
  printf 'valid_time=%s\ncase=%s\nmember=%s\nensemble_size=%s\nprocess_eam=%s\nprocess_elm=%s\njob_id=%s\nstarted_at=%s\n' \
    "${stamp}" "${my_casename}" "${member}" "${my_ensnum}" "${PROCESS_EAM}" "${PROCESS_ELM}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${progress}"

  tmpdir=$(mktemp -d "${date_out}/.interp.${member}.${stamp}.XXXXXX") || return 1
  trap '[[ -z "${tmpdir:-}" ]] || rm -rf -- "${tmpdir}"' EXIT
  if [[ "${PROCESS_ELM}" == "TRUE" ]]; then
    if ! ncks -O -v "${ELM_VARIABLES}" "${elm_in}" "${tmpdir}/elm.nc" || ! validate_output elm "${tmpdir}/elm.nc"; then
      rm -rf -- "${tmpdir}"
      echo "ERROR: ELM extraction failed for ${member} at ${stamp}" >&2
      return 1
    fi
  fi
  if [[ "${PROCESS_EAM}" == "TRUE" ]]; then
    if ! ncks -O -v "${EAM_VARIABLES}" "${eam_in}" "${tmpdir}/eam.nc" ||
       ! ncrename -O -v lat_d,lat -v lon_d,lon -d ncol_d,ncol "${tmpdir}/eam.nc" ||
       ! validate_output eam "${tmpdir}/eam.nc"; then
      rm -rf -- "${tmpdir}"
      echo "ERROR: EAM extraction failed for ${member} at ${stamp}" >&2
      return 1
    fi
  fi

  [[ "${PROCESS_ELM}" != "TRUE" ]] || mv -f "${tmpdir}/elm.nc" "${elm_out}"
  [[ "${PROCESS_EAM}" != "TRUE" ]] || mv -f "${tmpdir}/eam.nc" "${eam_out}"
  rmdir "${tmpdir}"
  tmpdir=""

  tmp_record="${record}.tmp.${SLURM_JOB_ID:-$$}"
  printf 'valid_time=%s\ncase=%s\nmember=%s\nensemble_size=%s\nprocess_eam=%s\nprocess_elm=%s\njob_id=%s\ncompleted_at=%s\n' \
    "${stamp}" "${my_casename}" "${member}" "${my_ensnum}" "${PROCESS_EAM}" "${PROCESS_ELM}" \
    "${SLURM_JOB_ID:-none}" "$(date '+%Y-%m-%d %H:%M:%S')" > "${tmp_record}"
  mv -f "${tmp_record}" "${record}"
  rm -f -- "${progress}"
  echo "Completed ${member} at ${stamp}"
}

run_date() {
  local stamp="$1" member index status failed=0
  local -a pids=() members=()
  validate_upstream_marker "${stamp}"
  echo "Processing ${stamp}"

  for ((i=1; i<=my_ensnum; i++)); do
    member=$(printf 'EN%02d' "${i}")
    process_member "${stamp}" "${member}" &
    pids+=("$!")
    members+=("${member}")
    if (( ${#pids[@]} >= MAX_CONCURRENT_WORKERS )); then
      for index in "${!pids[@]}"; do
        wait "${pids[$index]}" || { status=$?; echo "ERROR: ${members[$index]} failed at ${stamp} (${status})" >&2; failed=1; }
      done
      pids=()
      members=()
    fi
  done
  for index in "${!pids[@]}"; do
    wait "${pids[$index]}" || { status=$?; echo "ERROR: ${members[$index]} failed at ${stamp} (${status})" >&2; failed=1; }
  done
  (( failed == 0 )) || fail "one or more members failed at ${stamp}"
}

echo "Preflighting all requested restart inputs"
for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=86400)); do
  preflight_date "$(epoch_to_stamp "${epoch}")"
done
echo "Preflight completed successfully"

for ((epoch=START_EPOCH; epoch<=END_EPOCH; epoch+=86400)); do
  run_date "$(epoch_to_stamp "${epoch}")"
done

echo "Step 8 completed: ${INTERP_START} through ${INTERP_END}"
echo "Outputs: ${OUT_DIR}"
