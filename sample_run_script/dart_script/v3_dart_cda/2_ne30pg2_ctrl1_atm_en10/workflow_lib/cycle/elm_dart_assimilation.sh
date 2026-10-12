#!/bin/bash
# ELM-DART assimilation for one coupled Step 4 cycle.
# Sourced inside a dedicated subshell by e3sm_dart_single_cycle.sh.

#source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh
#source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-haswell.sh
[[ -n "${my_dart_env_file:-}" ]] || { echo "ERROR: my_dart_env_file is unset" >&2; exit 1; }
[[ -r "${my_dart_env_file}" ]] || { echo "ERROR: configured DART environment is not readable: ${my_dart_env_file}" >&2; exit 1; }
echo "Using configured DART machine environment: ${my_dart_env_file}"
source "${my_dart_env_file}"

elm_fail() { echo "ERROR: ELM DART: $*" >&2; return 1; }
elm_require_file() { [[ -s "$1" ]] || elm_fail "missing or empty file: $1"; }
elm_require_executable() { [[ -x "$1" ]] || elm_fail "missing executable: $1"; }
elm_positive_int() { [[ "${2:-}" =~ ^[1-9][0-9]*$ ]] || elm_fail "$1 must be a positive integer, got: ${2:-unset}"; }

elm_configure_namelist() {
  local nml="$1" logical_value period_days period_seconds
  period_days=$((my_elm_dart_cycle_hours / 24))
  period_seconds=$(((my_elm_dart_cycle_hours % 24) * 3600))

  sed -i -E \
    -e "s@^[[:space:]]*ens_size[[:space:]]*=.*@   ens_size = ${my_ensnum}@" \
    -e "s@^[[:space:]]*num_output_state_members[[:space:]]*=.*@   num_output_state_members = ${my_ensnum}@" \
    -e "s@^[[:space:]]*num_output_obs_members[[:space:]]*=.*@   num_output_obs_members = ${my_ensnum}@" \
    -e "s@^[[:space:]]*assimilation_period_days[[:space:]]*=.*@   assimilation_period_days = ${period_days}@" \
    -e "s@^[[:space:]]*assimilation_period_seconds[[:space:]]*=.*@   assimilation_period_seconds = ${period_seconds}@" \
    "${nml}" || return 1

  if [[ "${ELM_APPLY_LND_PROFILE}" == "TRUE" ]]; then
    for logical_value in "${lnd_da_output_sequential_prior_post}" "${lnd_da_use_sequential_prior_post}" "${lnd_da_perturb_from_single_instance}" "${lnd_da_spread_restoration}" "${lnd_da_sampling_error_correction}" "${lnd_da_horiz_dist_only}" "${lnd_da_strongly_coupled}"; do
      [[ "${logical_value}" =~ ^\.(true|false)\.$ ]] || { elm_fail "invalid land-pass Fortran logical: ${logical_value}"; return 1; }
    done
    [[ "${lnd_da_perturbation_amplitude}" =~ ^[0-9]+([.][0-9]+)?$ ]] || { elm_fail "invalid lnd_da_perturbation_amplitude"; return 1; }
    [[ "${lnd_da_perturbation_method}" == "uniform" || "${lnd_da_perturbation_method}" == "model" ]] || { elm_fail "lnd_da_perturbation_method must be uniform or model"; return 1; }
    [[ "${lnd_da_inf_flavor_prior}" =~ ^[0-9]+$ && "${lnd_da_inf_flavor_posterior}" =~ ^[0-9]+$ ]] || { elm_fail "land inflation flavors must be non-negative integers"; return 1; }
    [[ "${lnd_da_cutoff}" =~ ^[0-9]+([.][0-9]+)?$ ]] || { elm_fail "invalid lnd_da_cutoff"; return 1; }
    if [[ "${lnd_da_output_sequential_prior_post,,}" == ".true." && "${lnd_da_use_sequential_prior_post,,}" == ".true." ]]; then
      elm_fail "land filter cannot both output and consume sequential priors"
      return 1
    fi
    case "${lnd_da_state_model}:${lnd_da_obs_model}" in
      Atmosphere:Atmosphere|Atmosphere:Land|Land:Atmosphere|Land:Land) ;;
      *) elm_fail "land state/obs models must each be Atmosphere or Land"; return 1 ;;
    esac

    sed -i -E \
      -e "s@^[[:space:]]*perturb_from_single_instance[[:space:]]*=.*@   perturb_from_single_instance = ${lnd_da_perturb_from_single_instance}@" \
      -e "s@^[[:space:]]*perturbation_amplitude[[:space:]]*=.*@   perturbation_amplitude = ${lnd_da_perturbation_amplitude}@" \
      -e "s@^[[:space:]]*obs_sequence_in_name[[:space:]]*=.*@   obs_sequence_in_name = '${lnd_da_obs_sequence_in_name}'@" \
      -e "s@^[[:space:]]*inf_flavor[[:space:]]*=.*@   inf_flavor = ${lnd_da_inf_flavor_prior}, ${lnd_da_inf_flavor_posterior}@" \
      -e "s@^[[:space:]]*cutoff[[:space:]]*=.*@   cutoff = ${lnd_da_cutoff}@" \
      -e "s@^[[:space:]]*spread_restoration[[:space:]]*=.*@   spread_restoration = ${lnd_da_spread_restoration}@" \
      -e "s@^[[:space:]]*sampling_error_correction[[:space:]]*=.*@   sampling_error_correction = ${lnd_da_sampling_error_correction}@" \
      -e "s@^[[:space:]]*horiz_dist_only[[:space:]]*=.*@   horiz_dist_only = ${lnd_da_horiz_dist_only}@" \
      "${nml}" || return 1

    sed -i "/^[[:space:]]*perturbation_amplitude[[:space:]]*=/a\\   perturbation_method = '${lnd_da_perturbation_method}'\\n   output_sequential_prior_post = ${lnd_da_output_sequential_prior_post}\\n   use_sequential_prior_post = ${lnd_da_use_sequential_prior_post}" "${nml}" || return 1
    if grep -Eiq '^[[:space:]]*&strongly_coupled_localization_nml' "${nml}"; then
      elm_fail "input.nml already defines strongly_coupled_localization_nml; refusing to create a duplicate group"
      return 1
    fi
    printf '\n&strongly_coupled_localization_nml\n   strongly_coupled = %s\n   state_model = '\''%s'\''\n   obs_model = '\''%s'\''\n   /\n' "${lnd_da_strongly_coupled}" "${lnd_da_state_model}" "${lnd_da_obs_model}" >> "${nml}" || return 1

    grep -Eq "^[[:space:]]*ens_size[[:space:]]*=[[:space:]]*${my_ensnum}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*assimilation_period_days[[:space:]]*=[[:space:]]*${period_days}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*assimilation_period_seconds[[:space:]]*=[[:space:]]*${period_seconds}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*perturb_from_single_instance[[:space:]]*=[[:space:]]*${lnd_da_perturb_from_single_instance}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*perturbation_amplitude[[:space:]]*=[[:space:]]*${lnd_da_perturbation_amplitude}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*output_sequential_prior_post[[:space:]]*=[[:space:]]*${lnd_da_output_sequential_prior_post}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*use_sequential_prior_post[[:space:]]*=[[:space:]]*${lnd_da_use_sequential_prior_post}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*perturbation_method[[:space:]]*=[[:space:]]*'${lnd_da_perturbation_method}'" "${nml}" || return 1
    grep -Eq "^[[:space:]]*obs_sequence_in_name[[:space:]]*=[[:space:]]*'${lnd_da_obs_sequence_in_name}'" "${nml}" || return 1
    grep -Eq "^[[:space:]]*inf_flavor[[:space:]]*=[[:space:]]*${lnd_da_inf_flavor_prior},[[:space:]]*${lnd_da_inf_flavor_posterior}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*cutoff[[:space:]]*=[[:space:]]*${lnd_da_cutoff}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*spread_restoration[[:space:]]*=[[:space:]]*${lnd_da_spread_restoration}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*sampling_error_correction[[:space:]]*=[[:space:]]*${lnd_da_sampling_error_correction}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*horiz_dist_only[[:space:]]*=[[:space:]]*${lnd_da_horiz_dist_only}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*strongly_coupled[[:space:]]*=[[:space:]]*${lnd_da_strongly_coupled}" "${nml}" || return 1
    grep -Eq "^[[:space:]]*state_model[[:space:]]*=[[:space:]]*'${lnd_da_state_model}'" "${nml}" || return 1
    grep -Eq "^[[:space:]]*obs_model[[:space:]]*=[[:space:]]*'${lnd_da_obs_model}'" "${nml}" || return 1
  fi

  if [[ -z "${my_elm_vector_history_stream}" ]]; then
    sed -i \
      -e "/^[[:space:]]*'vector_files\.txt'/d" \
      -e "s/,'elm_vector_history\.nc'//g" \
      -e "/'vector'[[:space:]]*,[[:space:]]*'NO_COPY_BACK'/d" \
      "${nml}" || return 1
  fi
}

elm_find_observation() {
  local stamp="$1" yyyymm="${stamp:0:4}${stamp:5:2}" candidate
  for candidate in \
    "${my_elm_dart_obsdir}/${yyyymm}_6H/obs_seq.${stamp}" \
    "${my_elm_dart_obsdir}/${yyyymm}/obs_seq.${stamp}" \
    "${my_elm_dart_obsdir}/obs_seq.${stamp}"; do
    [[ -s "${candidate}" ]] && { printf '%s\n' "${candidate}"; return 0; }
  done
  elm_fail "no land observation sequence found for ${stamp} under ${my_elm_dart_obsdir}"
}


elm_target_iso_time() {
  local seconds hour minute second
  seconds=$((10#${DART_SECONDS}))
  hour=$((seconds / 3600))
  minute=$(((seconds % 3600) / 60))
  second=$((seconds % 60))
  printf '%04d-%02d-%02dT%02d:%02d:%02d\n' \
    "$((10#${DART_YEAR}))" "$((10#${DART_MONTH}))" "$((10#${DART_DAY}))" \
    "${hour}" "${minute}" "${second}"
}

elm_iso_to_stamp() {
  local iso="$1" ymd hour minute second seconds
  [[ "${iso}" =~ ^([0-9]{4}-[0-9]{2}-[0-9]{2})T([0-9]{2}):([0-9]{2}):([0-9]{2})$ ]] || {
    elm_fail "invalid ELM history time: ${iso}"
    return 1
  }
  ymd=${BASH_REMATCH[1]}
  hour=$((10#${BASH_REMATCH[2]}))
  minute=$((10#${BASH_REMATCH[3]}))
  second=$((10#${BASH_REMATCH[4]}))
  (( hour < 24 && minute < 60 && second < 60 )) || {
    elm_fail "ELM history time is out of range: ${iso}"
    return 1
  }
  seconds=$((hour * 3600 + minute * 60 + second))
  printf '%s-%05d\n' "${ymd}" "${seconds}"
}

elm_find_history_record() {
  local enstr="$1" stream="$2" target="$3" stamp filename candidate member_archive
  member_archive="${my_modeldir}/${enstr}/archive"
  stamp=$(elm_iso_to_stamp "${target}") || return 1
  filename="${my_casename}.${enstr}.elm.${stream}.${stamp}.nc"
  for candidate in \
    "${member_archive}/lnd/hist/${filename}" \
    "${my_modeldir}/${enstr}/run/${filename}"; do
    if [[ -s "${candidate}" ]]; then
      printf '%s|0\n' "${candidate}"
      return 0
    fi
  done
  elm_fail "no ${stream} history file for ${enstr} at ${target} (expected ${filename})"
}

elm_link_single_history_record() {
  local record="$1" output="$2" expected="$3" source index stamp
  source=${record%|*}
  index=${record##*|}
  stamp=$(elm_iso_to_stamp "${expected}") || return 1
  [[ "${index}" == "0" && -s "${source}" && "${source}" == *".${stamp}.nc" ]] || {
    elm_fail "history source does not match expected timestamp ${stamp}: ${source}"
    return 1
  }
  [[ -z "${output}" ]] || ln -s "${source}" "${output}" || return 1
}


elm_validate_vector_record() {
  local record="$1" source header variable
  local -a required_variables=(NEE H2OSNO TLAI TWS SMP TV RH2M_R PBOT TBOT)
  source=${record%|*}
  header=$(ncdump -h "${source}" 2>/dev/null) || {
    elm_fail "invalid vector history file: ${source}"
    return 1
  }
  for variable in "${required_variables[@]}"; do
    grep -Eq "^[[:space:]]*(byte|char|short|int|int64|float|double|ubyte|ushort|uint|uint64|string)[[:space:]]+${variable}\\(" <<< "${header}" || {
      elm_fail "vector history file is missing ${variable}: ${source}"
      return 1
    }
  done
}

elm_preflight() {
  local name executable i enstr restart history_record vector_record target_iso
  for name in my_ensnum my_elm_dart_nnodes my_task_per_node; do
    elm_positive_int "${name}" "${!name:-}" || return 1
  done
  if [[ "${ELM_USES_SEQUENTIAL_PRIOR}" != "TRUE" ]]; then
    [[ -n "${my_elm_dart_obsdir:-}" ]] || elm_fail "my_elm_dart_obsdir is unset" || return 1
  fi
  [[ -n "${my_elm_history_stream:-}" ]] || elm_fail "my_elm_history_stream is unset" || return 1
  for executable in filter elm_to_dart; do
    elm_require_executable "${ELM_WORK}/${executable}" || return 1
  done
  if [[ "${ELM_APPLY_LND_PROFILE}" != "TRUE" ]] || (( lnd_da_inf_flavor_prior != 0 || lnd_da_inf_flavor_posterior != 0 )); then
    elm_require_executable "${ELM_WORK}/fill_inflation_restart" || return 1
  fi
  elm_require_file "${my_elm_filter_nml}" || return 1
  target_iso=$(elm_target_iso_time) || return 1
  if [[ "${ELM_DART_PREFLIGHT_STATIC_ONLY:-FALSE}" != "TRUE" ]]; then
    if [[ "${ELM_USES_SEQUENTIAL_PRIOR}" == "TRUE" ]]; then
      ELM_OBS="${my_eam_dart_run_dir}/${ELM_STAMP}/${my_casename}.dart.e.eam_obs_seq_final.${ELM_STAMP}"
      elm_require_file "${ELM_OBS}" || return 1
    else
      ELM_OBS=$(elm_find_observation "${ELM_STAMP}") || return 1
    fi
  fi
  for i in $(seq 1 "${my_ensnum}"); do
    enstr=$(printf 'EN%02d' "${i}")
    member_archive="${my_modeldir}/${enstr}/archive"
    restart="${member_archive}/rest/${ELM_STAMP}/${my_casename}.${enstr}.elm.r.${ELM_STAMP}.nc"
    elm_require_file "${restart}" || return 1
    ncdump -h "${restart}" >/dev/null 2>&1 || elm_fail "invalid restart: ${restart}" || return 1
    history_record=$(elm_find_history_record "${enstr}" "${my_elm_history_stream}" "${target_iso}") || return 1
    elm_link_single_history_record "${history_record}" "" "${target_iso}" || return 1
    if [[ -n "${my_elm_vector_history_stream}" ]]; then
      vector_record=$(elm_find_history_record "${enstr}" "${my_elm_vector_history_stream}" "${target_iso}") || return 1
      elm_validate_vector_record "${vector_record}" || return 1
      elm_link_single_history_record "${vector_record}" "" "${target_iso}" || return 1
    fi
  done
}

ELM_WORK="${my_elm_dart_code}/models/${my_elm_dart_model}/work"
ELM_STAMP=$(printf '%04d-%02d-%02d-%05d' "$((10#${DART_YEAR}))" "$((10#${DART_MONTH}))" "$((10#${DART_DAY}))" "$((10#${DART_SECONDS}))")
ELM_RUNROOT="${my_elm_dart_run_dir}"
ELM_DADIR="${ELM_RUNROOT}/${ELM_STAMP}"
ELM_MARKER="${DA_TRANSACTION_DIR}/.dart_elm_filter_in_progress"
ELM_COMPLETE="${my_status_dir}/elm_assim_complete.${ELM_STAMP}"
ELM_APPLY_LND_PROFILE="FALSE"
ELM_USES_SEQUENTIAL_PRIOR="FALSE"
case "${strongly_coupled_on,,}" in
  on)
    ELM_APPLY_LND_PROFILE="TRUE"
    [[ "${lnd_da_use_sequential_prior_post,,}" != ".true." ]] || ELM_USES_SEQUENTIAL_PRIOR="TRUE"
    ;;
  off) ;;
  *) elm_fail "strongly_coupled_on must be on or off, got: ${strongly_coupled_on:-unset}"; return 1 ;;
esac
ELM_OBS_LINK_NAME="obs_seq.out"
[[ "${ELM_APPLY_LND_PROFILE}" != "TRUE" ]] || ELM_OBS_LINK_NAME="${lnd_da_obs_sequence_in_name}"

elm_preflight || return 1
if [[ "${ELM_DART_PREFLIGHT_ONLY:-FALSE}" == TRUE ]]; then
  echo "ELM DART ${ELM_DART_PREFLIGHT_STATIC_ONLY:+static }preflight passed for ${ELM_STAMP}"
  return 0
fi
if [[ "${ELM_DART_TRANSACTION_STARTED:-FALSE}" == "TRUE" ]]; then
  [[ -s "${ELM_MARKER}" ]] || elm_fail "sequential transaction marker is missing: ${ELM_MARKER}" || return 1
else
  [[ ! -e "${ELM_MARKER}" ]] || elm_fail "stale marker requires a full forecast rebuild: ${ELM_MARKER}" || return 1
fi
mkdir -p "${ELM_DADIR}" "${my_status_dir}" "${DA_TRANSACTION_DIR}" || return 1
find "${ELM_DADIR}" -mindepth 1 -maxdepth 1 -delete || return 1
rm -f -- "${ELM_COMPLETE}"
printf 'cycle=%s\nvalid_time=%s\nslurm_job_id=%s\nstarted_at=%s\nphase=elm_filter\n' \
  "${DATA_ASSIMILATION_CYCLES}" "${ELM_STAMP}" "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${ELM_MARKER}" || return 1

cp -p "${my_elm_filter_nml}" "${ELM_DADIR}/input.nml" || return 1
elm_configure_namelist "${ELM_DADIR}/input.nml" || elm_fail "could not configure input.nml" || return 1
ln -s "${ELM_OBS}" "${ELM_DADIR}/${ELM_OBS_LINK_NAME}" || return 1
for executable in filter elm_to_dart; do
  cp -p "${ELM_WORK}/${executable}" "${ELM_DADIR}/${executable}" || return 1
done
if [[ "${ELM_APPLY_LND_PROFILE}" != "TRUE" ]] || (( lnd_da_inf_flavor_prior != 0 || lnd_da_inf_flavor_posterior != 0 )); then
  cp -p "${ELM_WORK}/fill_inflation_restart" "${ELM_DADIR}/fill_inflation_restart" || return 1
fi

cd "${ELM_DADIR}" || return 1
: > restart_files.txt
: > history_files.txt
[[ -z "${my_elm_vector_history_stream}" ]] || : > vector_files.txt
ELM_TARGET_ISO=$(elm_target_iso_time) || return 1
for i in $(seq 1 "${my_ensnum}"); do
  enstr=$(printf 'EN%02d' "${i}"); inst=$(printf '%04d' "${i}")
  member_archive="${my_modeldir}/${enstr}/archive"
  restart_source="${member_archive}/rest/${ELM_STAMP}/${my_casename}.${enstr}.elm.r.${ELM_STAMP}.nc"
  ln -s "${restart_source}" elm.nc || return 1
  ./elm_to_dart > "elm_to_dart.${enstr}.log" 2>&1 || elm_fail "elm_to_dart failed for ${enstr}" || return 1
  mv -f elm.nc "elm_restart_${inst}.nc" || return 1
  history_record=$(elm_find_history_record "${enstr}" "${my_elm_history_stream}" "${ELM_TARGET_ISO}") || return 1
  elm_link_single_history_record "${history_record}" "elm_history_${inst}.nc" "${ELM_TARGET_ISO}" || return 1
  printf 'elm_restart_%s.nc\n' "${inst}" >> restart_files.txt
  printf 'elm_history_%s.nc\n' "${inst}" >> history_files.txt
  if [[ -n "${my_elm_vector_history_stream}" ]]; then
    vector_record=$(elm_find_history_record "${enstr}" "${my_elm_vector_history_stream}" "${ELM_TARGET_ISO}") || return 1
    elm_validate_vector_record "${vector_record}" || return 1
    elm_link_single_history_record "${vector_record}" "elm_vector_history_${inst}.nc" "${ELM_TARGET_ISO}" || return 1
    printf 'elm_vector_history_%s.nc\n' "${inst}" >> vector_files.txt
  fi
done

# Initialize inflation only when the configured land pass uses it.
if [[ "${ELM_APPLY_LND_PROFILE}" != "TRUE" ]] || (( lnd_da_inf_flavor_prior != 0 || lnd_da_inf_flavor_posterior != 0 )); then
  need_prior=FALSE
  need_posterior=FALSE
  [[ "${ELM_APPLY_LND_PROFILE}" == "TRUE" && "${lnd_da_inf_flavor_prior}" == "0" ]] || need_prior=TRUE
  [[ "${ELM_APPLY_LND_PROFILE}" != "TRUE" || "${lnd_da_inf_flavor_posterior}" == "0" ]] || need_posterior=TRUE
  ln -s elm_restart_0001.nc elm_restart.nc
  ln -s elm_history_0001.nc elm_history.nc
  [[ -z "${my_elm_vector_history_stream}" ]] || ln -s elm_vector_history_0001.nc elm_vector_history.nc
  previous_dir=""
  while IFS= read -r candidate; do
    [[ "${need_prior}" != "TRUE" ]] || compgen -G "${candidate}/output_priorinf_*.nc" >/dev/null || continue
    [[ "${need_posterior}" != "TRUE" ]] || compgen -G "${candidate}/output_postinf_*.nc" >/dev/null || continue
    previous_dir="${candidate}"
    break
  done < <(find "${ELM_RUNROOT}" -mindepth 1 -maxdepth 1 -type d ! -path "${ELM_DADIR}" -print | sort -r)
  if [[ -n "${previous_dir}" ]]; then
    if [[ "${need_prior}" == "TRUE" ]]; then
      for source in "${previous_dir}"/output_priorinf_*.nc; do
        target=${source##*/}; target=${target/output_/input_}; cp -p "${source}" "${target}" || return 1
      done
    fi
    if [[ "${need_posterior}" == "TRUE" ]]; then
      for source in "${previous_dir}"/output_postinf_*.nc; do
        target=${source##*/}; target=${target/output_/input_}; cp -p "${source}" "${target}" || return 1
      done
    fi
  else
    ./fill_inflation_restart > fill_inflation_restart.log 2>&1 || elm_fail "could not initialize inflation" || return 1
  fi
  rm -f elm_restart.nc elm_history.nc elm_vector_history.nc
fi

# The ELM model interface reads the ensemble time from the conventional
# single-member filenames during filter initialization, even though the
# ensemble itself is supplied through restart_files.txt/history_files.txt.
# Keep these aliases present for the filter run.
ln -s elm_restart_0001.nc elm_restart.nc || return 1
ln -s elm_history_0001.nc elm_history.nc || return 1
[[ -z "${my_elm_vector_history_stream}" ]] || ln -s elm_vector_history_0001.nc elm_vector_history.nc || return 1

ELM_DART_NTASKS=$((my_elm_dart_nnodes * my_task_per_node))
echo "$(date) -- BEGIN ELM FILTER (${my_elm_dart_nnodes} nodes, ${ELM_DART_NTASKS} tasks)"
srun --exclusive --nodes="${my_elm_dart_nnodes}" --ntasks="${ELM_DART_NTASKS}" \
  --mpi=pmi2 --kill-on-bad-exit -l --cpu_bind=cores -c 1 -m plane="${my_task_per_node}" ./filter \
  || elm_fail "filter failed for ${ELM_STAMP}" || return 1
echo "$(date) -- END ELM FILTER"
rm -f elm_restart.nc elm_history.nc elm_vector_history.nc || return 1

# Filter writes each posterior through elm_restart_NNNN.nc, a symlink to
# the archived target-time restart. A stale marker forces full reforecasting.
for i in $(seq 1 "${my_ensnum}"); do
  enstr=$(printf 'EN%02d' "${i}"); inst=$(printf '%04d' "${i}")
  member_archive="${my_modeldir}/${enstr}/archive"
  restart_target="${member_archive}/rest/${ELM_STAMP}/${my_casename}.${enstr}.elm.r.${ELM_STAMP}.nc"
  [[ -L "elm_restart_${inst}.nc" ]] || elm_fail "working restart is not a symlink for ${enstr}" || return 1
  [[ "$(readlink -f "elm_restart_${inst}.nc")" == "$(readlink -f "${restart_target}")" ]] || elm_fail "working restart target mismatch for ${enstr}" || return 1
  ncdump -h "${restart_target}" >/dev/null 2>&1 || elm_fail "invalid analyzed restart for ${enstr}" || return 1
done

[[ -s obs_seq.final ]] && mv -f obs_seq.final "elm_obs_seq.${ELM_STAMP}.final"
[[ -s dart_log.out ]] && mv -f dart_log.out "elm_dart_log.${ELM_STAMP}.out"
tmp_record="${ELM_COMPLETE}.tmp.${SLURM_JOB_ID:-$$}"
printf 'valid_time=%s\ncase=%s\nensemble_size=%s\ncycle=%s\nda_mode=on\nslurm_job_id=%s\ncompleted_at=%s\n' \
  "${ELM_STAMP}" "${my_casename}" "${my_ensnum}" "${DATA_ASSIMILATION_CYCLES}" \
  "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${tmp_record}" || return 1
mv -f "${tmp_record}" "${ELM_COMPLETE}" || return 1
rm -f "${ELM_MARKER}" || return 1
echo "Validated ${my_ensnum} ELM analyses and wrote ${ELM_COMPLETE}"
