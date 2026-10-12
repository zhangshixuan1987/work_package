#!/bin/bash
# Two-stage, four-pass bidirectional EAM/ELM assimilation.
# Stage 1 native passes run concurrently; stage 2 cross passes run concurrently.

bidirectional_fail() { echo "ERROR: bidirectional DART: $*" >&2; return 1; }

BD_STAMP=$(printf '%04d-%02d-%02d-%05d' "$((10#${DART_YEAR}))" "$((10#${DART_MONTH}))" "$((10#${DART_DAY}))" "$((10#${DART_SECONDS}))")
BD_EAM_NATIVE_DIR="${my_eam_dart_run_dir}/${BD_STAMP}/eam_atmosphere"
BD_ELM_NATIVE_DIR="${my_elm_dart_run_dir}/${BD_STAMP}/elm_land"
BD_EAM_PRIOR="${BD_EAM_NATIVE_DIR}/${my_casename}.dart.e.eam_obs_seq_final.${BD_STAMP}"
BD_ELM_PRIOR="${BD_ELM_NATIVE_DIR}/elm_obs_seq.${BD_STAMP}.final"
BD_LOG_ROOT="${LOG_DIR}/bidirectional.${BD_STAMP}"
BD_MARKER="${DA_TRANSACTION_DIR}/.dart_bidirectional_in_progress"
mkdir -p "${BD_LOG_ROOT}" "${DA_TRANSACTION_DIR}" \
  || bidirectional_fail "could not create bidirectional working directories" \
  || return 1
[[ ! -e "${BD_MARKER}" ]] || bidirectional_fail "stale bidirectional marker requires a full forecast rebuild: ${BD_MARKER}" || return 1

run_bd_pass() {
  local pass_name="$1" component_script="$2"
  shift 2
  ( set -Eeuo pipefail; cd "${my_wkdir}"; export "$@"; . "${component_script}" ) \
    > "${BD_LOG_ROOT}/${pass_name}.log" 2>&1
}

wait_bd_pass() {
  local pass_name="$1" pid="$2" status=0
  wait "${pid}" || status=$?
  if (( status != 0 )); then
    echo "ERROR: ${pass_name} failed with status ${status}; log: ${BD_LOG_ROOT}/${pass_name}.log" >&2
    tail -n 60 "${BD_LOG_ROOT}/${pass_name}.log" >&2 || true
    return "${status}"
  fi
  echo "$(date) -- completed bidirectional pass ${pass_name}"
}

printf 'cycle=%s\nvalid_time=%s\nslurm_job_id=%s\nstarted_at=%s\nphase=stage1\n' \
  "${DATA_ASSIMILATION_CYCLES}" "${BD_STAMP}" "${SLURM_JOB_ID:-none}" "$(date '+%F %T')" > "${BD_MARKER}" || return 1
echo "$(date) -- BEGIN BIDIRECTIONAL STAGE 1: native component passes"

# EAM and ELM are independent during the native stage. Start EAM immediately
# while ELM performs its component-specific preflight; an ELM preflight failure
# must not prevent EAM from completing. Both native outputs are still required
# before either cross-component pass can start.
run_bd_pass eam_atmosphere "${my_workflow_lib}/cycle/eam_dart_assimilation.sh" \
  DART_PASS_TAG=eam_atmosphere \
  atm_da_compute_posterior=.false. \
  atm_da_output_sequential_prior_post=.true. atm_da_use_sequential_prior_post=.false. \
  atm_da_output_mean=.true. atm_da_output_sd=.true. atm_da_output_members=.true. \
  atm_da_inf_flavor_prior=0 atm_da_inf_flavor_posterior=0 \
  atm_da_cutoff="${atm_da_cutoff}" atm_da_spread_restoration=.false. \
  atm_da_sampling_error_correction=.false. atm_da_horiz_dist_only=.true. \
  atm_da_no_normalization_of_scale_heights=.true. \
  atm_da_obs_sequence_in_name=obs_seq.out \
  atm_da_strongly_coupled=.false. atm_da_state_model=Atmosphere atm_da_obs_model=Atmosphere &
bd_eam_native_pid=$!

run_bd_pass elm_land_preflight "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" \
  DART_PASS_TAG=elm_land ELM_DART_PREFLIGHT_ONLY=TRUE ELM_DART_TRANSACTION_STARTED=FALSE \
  lnd_da_output_sequential_prior_post=.true. lnd_da_use_sequential_prior_post=.false. \
  lnd_da_obs_sequence_in_name=obs_seq.out lnd_da_sequential_obs_file= \
  lnd_da_strongly_coupled=.false. lnd_da_state_model=Land lnd_da_obs_model=Land &
bd_elm_preflight_pid=$!

bd_stage1_status=0
bd_elm_native_pid=
if wait_bd_pass elm_land_preflight "${bd_elm_preflight_pid}"; then
  run_bd_pass elm_land "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" \
    DART_PASS_TAG=elm_land ELM_DART_TRANSACTION_STARTED=FALSE \
    lnd_da_output_sequential_prior_post=.true. lnd_da_use_sequential_prior_post=.false. \
    lnd_da_obs_sequence_in_name=obs_seq.out lnd_da_sequential_obs_file= \
    lnd_da_strongly_coupled=.false. lnd_da_state_model=Land lnd_da_obs_model=Land &
  bd_elm_native_pid=$!
else
  echo "ERROR: native ELM preflight failed; EAM native assimilation will be allowed to finish" >&2
  bd_stage1_status=1
fi

wait_bd_pass eam_atmosphere "${bd_eam_native_pid}" || bd_stage1_status=1
if [[ -n "${bd_elm_native_pid}" ]]; then
  wait_bd_pass elm_land "${bd_elm_native_pid}" || bd_stage1_status=1
fi
(( bd_stage1_status == 0 )) || bidirectional_fail "stage 1 failed; cross passes were not started" || return 1
[[ -s "${BD_EAM_PRIOR}" ]] || bidirectional_fail "missing EAM sequential-prior output: ${BD_EAM_PRIOR}" || return 1
[[ -s "${BD_ELM_PRIOR}" ]] || bidirectional_fail "missing ELM sequential-prior output: ${BD_ELM_PRIOR}" || return 1

# Validate the ELM cross input and updated land ensemble before either cross
# pass begins. The EAM pass validates the same ELM-prior file independently.
run_bd_pass elm_atmosphere_cross_preflight "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" \
  DART_PASS_TAG=elm_atmosphere_cross ELM_DART_PREFLIGHT_ONLY=TRUE ELM_DART_TRANSACTION_STARTED=FALSE \
  lnd_da_output_sequential_prior_post=.false. lnd_da_use_sequential_prior_post=.true. \
  lnd_da_obs_sequence_in_name=obs_seq.sequential_prior_eam \
  lnd_da_sequential_obs_file="${BD_EAM_PRIOR}" \
  lnd_da_strongly_coupled=.true. lnd_da_state_model=Land lnd_da_obs_model=Atmosphere \
  || bidirectional_fail "cross ELM preflight failed" || return 1

sed -i 's/^phase=.*/phase=stage2/' "${BD_MARKER}" || return 1
echo "$(date) -- BEGIN BIDIRECTIONAL STAGE 2: cross-component passes"
run_bd_pass elm_atmosphere_cross "${my_workflow_lib}/cycle/elm_dart_assimilation.sh" \
  DART_PASS_TAG=elm_atmosphere_cross ELM_DART_TRANSACTION_STARTED=FALSE \
  lnd_da_output_sequential_prior_post=.false. lnd_da_use_sequential_prior_post=.true. \
  lnd_da_obs_sequence_in_name=obs_seq.sequential_prior_eam \
  lnd_da_sequential_obs_file="${BD_EAM_PRIOR}" \
  lnd_da_strongly_coupled=.true. lnd_da_state_model=Land lnd_da_obs_model=Atmosphere &
bd_elm_cross_pid=$!

run_bd_pass eam_land_cross "${my_workflow_lib}/cycle/eam_dart_assimilation.sh" \
  DART_PASS_TAG=eam_land_cross \
  atm_da_compute_posterior=.false. \
  atm_da_output_sequential_prior_post=.false. atm_da_use_sequential_prior_post=.true. \
  atm_da_output_mean=.true. atm_da_output_sd=.true. atm_da_output_members=.true. \
  atm_da_inf_flavor_prior=0 atm_da_inf_flavor_posterior=0 \
  atm_da_cutoff=0.4 atm_da_spread_restoration=.false. \
  atm_da_sampling_error_correction=.false. atm_da_horiz_dist_only=.true. \
  atm_da_no_normalization_of_scale_heights=.false. \
  atm_da_obs_sequence_in_name=obs_seq.sequential_prior_elm \
  atm_da_sequential_obs_file="${BD_ELM_PRIOR}" \
  atm_da_strongly_coupled=.true. atm_da_state_model=Atmosphere atm_da_obs_model=Land &
bd_eam_cross_pid=$!

bd_stage2_status=0
wait_bd_pass elm_atmosphere_cross "${bd_elm_cross_pid}" || bd_stage2_status=1
wait_bd_pass eam_land_cross "${bd_eam_cross_pid}" || bd_stage2_status=1
(( bd_stage2_status == 0 )) || bidirectional_fail "stage 2 failed; cycle handoff is blocked" || return 1

rm -f "${BD_MARKER}" || return 1
echo "Completed all bidirectional DART passes for ${BD_STAMP}"
