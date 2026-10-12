# Shared configuration for all numbered workflow stages.
# This file is sourced; it should define settings without running workflow work.

# --- Workflow and runtime paths ---------------------------------------------
# Update this absolute root if the workflow directory is moved.
export my_workflow_root="/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/4_ne30pg2_dart_scp_en10"
export my_workflow_lib="${my_workflow_root}/workflow_lib"
export my_runtime_dir="${my_workflow_root}/runtmp"
export my_log_dir="${my_runtime_dir}/logs"
export my_status_dir="${my_runtime_dir}/status"
export my_lock_dir="${my_runtime_dir}/locks"
export my_handoff_dir="${my_runtime_dir}/handoff"
export my_run_script_dir="${my_runtime_dir}/run_scripts"

################################################################################
# --- Machine and Slurm defaults ---------------------------------------------
# Explicit analysis environment used by Steps 2 and 3 for NCO and related tools.
################################################################################
export my_machine=compy
export my_project="esmd"
export my_jobqueue="short"
export my_walltime="02:00:00"
export my_task_per_node=40
export my_job_nnodes=20
export my_layout="custom-2_1x6_nhours"

################################################################################
# --- Machine and Slurm defaults ---------------------------------------------
# Explicit analysis environment used by Steps 2 and 3 for NCO and related tools.
################################################################################
export my_conda_setup_file="/people/zhan391/miniforge3/etc/profile.d/conda.sh"
export my_analysis_conda_env="e3sm_analysis"
export my_dart_env_file="${my_workflow_lib}/env/env_${my_machine}_specific.sh"
export my_eam_filter_nml="${my_workflow_lib}/namelists/eam/filter.nml"
export my_eam_perturb_nml="${my_workflow_lib}/namelists/eam/perturb.nml"
export my_eam_diag_nml="${my_workflow_lib}/namelists/eam/diagnostics.nml"
export my_elm_filter_nml="${my_workflow_lib}/namelists/elm/filter.nml"

################################################################################
# --- Ensemble and cycle execution -------------------------------------------
# Initial model time, ensemble size, concurrency, retries, and cycle batching.
################################################################################
export my_ensnum=10
export my_casedate="2011-12-01"
export my_casetod="00000"
export my_nodes_per_member=2
export my_retry_forecast_timeout_sec=12600
export my_wait_poll_interval_sec=20
export my_max_parallel_setup=10
export my_max_parallel_handoff=10
# Cycles attempted in one Step 4 allocation; use 1 for one cycle per job.
export my_cycles_per_job=48
# Start another cycle only when this runtime plus the shutdown margin remains.
export my_min_cycle_time_sec=12600
export my_cycle_shutdown_margin_sec=900

################################################################################
# --- E3SM experiment identity and output paths -------------------------------
# my_runtype must be Full-CPL or AMIP; keep compset and reference files consistent.
################################################################################
export my_e3sm_code="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
export my_runtype="AMIP"
export my_compset="F20TR"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/compyfs/zhan391/v3_dart_cda_scratch"
export my_casename="DARTSCPEN${my_ensnum}_${my_compset}_${my_resolution}_${my_machine}"

################################################################################
# Model experiment directory and shared executable. Step 1 builds E3SM once
# and reuses that executable for the ensemble members.
################################################################################
export my_modeldir="${my_runpath}/${my_casename}"
export my_modelcase="${my_modeldir}/case_scripts"
export my_modelexe="${my_modeldir}/build/e3sm.exe"
export my_eam_topography_file="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_eam_se_mapping_file="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_eam_cs_grid_file="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_elm_hist_mapping_file="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_traave.20231110.nc"
export my_eam_hist_mapping_file="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_eam_init_mapping_file="/compyfs/zhan391/acme_init/map_file/map_ne30np4_to_cmip6_180x360_aave.c231227.nc"

################################################################################
# --- Reference initial conditions -------------------------------------------
# Full-CPL uses the coupled component restart set below. AMIP does not require
# MPAS-O, but the active initialization checks still require MPAS-I and coupler.
################################################################################
export my_refcase="20241109.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.pm-cpu"
export my_refdate=${my_casedate}
export my_reftod=${my_casetod}
export my_refdir="/compyfs/zhan391/acme_init/E3SMv3_INT/${my_refdate}-${my_reftod}"
export my_refeam_in="/compyfs/zhan391/acme_init/E3SMv3_INT/0101-01-01-00000/v3.LR.piControl.eam.i.0101-01-01-00000.nc"
export my_refeam_ic="${my_refdir}/v3.LR.ne30pg2.ERA5.eam.i.${my_refdate}-${my_reftod}.nc"
export my_refelm_in="${my_refdir}/${my_refcase}.elm.r.${my_refdate}-${my_reftod}.nc"
export my_refrof_in="${my_refdir}/${my_refcase}.mosart.r.${my_refdate}-${my_reftod}.nc"
export my_refocn_in="${my_refdir}/${my_refcase}.mpaso.rst.${my_refdate}_${my_reftod}.nc"
export my_refice_in="${my_refdir}/${my_refcase}.mpassi.rst.${my_refdate}_${my_reftod}.nc"
export my_refcpl_in="${my_refdir}/${my_refcase}.cpl.r.${my_refdate}-${my_reftod}.nc"

################################################################################
# --- DART installation and observations -------------------------------------
# Live completed-cycle counter. Step 4 updates this value transactionally.
################################################################################
export my_e3sm_completed_cycles=1
export my_e3sm_cycle_hours=6 # shared E3SM forecast/handoff interval 
# One-time recovery override. When TRUE, Step 4 may advance uniformly older
# case XML clocks to the configured completed-cycle state after validating all
# member restarts and archived DART evidence. A successful cycle resets it.
export my_allow_forward_xml_recovery="FALSE"
# Forecast reuse controls:
#   my_forecast_ready_for_da=FALSE keeps the readiness checks enabled.
#   my_skip_completed_members=TRUE reuses each valid archived member forecast
#   and reruns only members whose target forecast is missing or invalid.
#   FALSE forces every member forecast to be rebuilt. A stale DART transaction
#   marker also forces a full rebuild, regardless of this configured value.
export my_forecast_ready_for_da="FALSE"
export my_skip_completed_members="FALSE"
# Raw CIME output is owned by each ENxx case. The layout value is recorded in
# completion markers so downstream stages can reject incompatible old state.
export my_dart_root="${my_modeldir}/dart_en$(printf '%02d' "${my_ensnum}")"
export my_e3sm_start_date=${my_casedate}
export my_e3sm_start_tod=${my_casetod}
export my_e3sm_end_date="2011-12-01"
export my_e3sm_end_tod="21600"

################################################################################
# --- EAM DART assimilation --------------------------------------------------
# Step 4 runs this component only when its cadence is due; simultaneous due components split the nodes.
################################################################################
export my_eam_dart_da="on"
export my_eam_dart_cycle_hours=6
export my_eam_dart_end_date="${my_e3sm_end_date}"
export my_eam_dart_end_tod="${my_e3sm_end_tod}"
export my_eam_dart_run_dir="${my_dart_root}/eam"
export my_eam_dart_model="eam-se"
export my_eam_dart_pgrid=".true."
export my_eam_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART_SCP"
# EAM model-top observation screening and vertical-localization controls.
export my_eam_no_obs_assim_above_level=5
export my_eam_use_log_vertical_scale=".true."
export my_eam_vert_normalization_scale_height="1.5"
export my_eam_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"
# Optional per-cycle EAM DA settings. Keys combine the exact valid time and
# parameter name. Omitted parameters retain their normal defaults.
declare -Ag my_eam_cycle_overrides=(
  # ["2011-12-02-00000:localization_cutoff"]="0.10"
  # ["2011-12-02-00000:inflation_damping"]="0.9"
  # ["2011-12-02-00000:no_obs_assim_above_level"]="7"
)


################################################################################
# --- ELM DART assimilation --------------------------------------------------
# Step 4 runs this component only when its cadence is due; simultaneous due components split the nodes.
################################################################################
# Disabled for this reverse-coupling quick test: the pre-existing ELM
# sequential-prior file below stands in for a live ELM DART pass.
export my_elm_dart_da="on"
export my_elm_dart_cycle_hours=6
export my_elm_dart_end_date="${my_e3sm_end_date}"
export my_elm_dart_end_tod="${my_e3sm_end_tod}"
export my_elm_dart_run_dir="${my_dart_root}/elm"
export my_elm_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART_SCP"
export my_elm_dart_obsdir="/compyfs/zhan391/acme_init/Observations/SMAP"
export my_elm_dart_model="elm"
export my_elm_history_stream="h1"
export my_elm_vector_history_stream="h2"

################################################################################
# --- Strongly coupled DA setup ----------------------------------------------
################################################################################
export strongly_coupled_on="on"
# legacy uses the individual ATM/LND profile below. Bidirectional executes
# native EAM+ELM passes first, followed by both cross-component passes.
export dart_coupling_mode="bidirectional"

################################################################################
# Stage 1 native EAM profile. The bidirectional orchestrator changes these
# settings only inside the Stage 2 EAM land-cross subprocess.
################################################################################
# Shared EAM settings used by both native and cross passes.
export atm_da_compute_posterior=".false."
export atm_da_output_mean=".true."
export atm_da_output_sd=".true."
export atm_da_output_members=".true."
export atm_da_inf_flavor_prior="0"
export atm_da_inf_flavor_posterior="0"
export atm_da_cutoff="0.4"
export atm_da_spread_restoration=".false."
export atm_da_sampling_error_correction=".false."
export atm_da_horiz_dist_only=".true."
# Stage 1 native EAM input; replaced with the ELM sequential-prior filename
# for the Stage 2 EAM land-to-atmosphere cross pass.
# Stage 1: obs_sequence_in_name = 'obs_seq.out'
# Stage 2: obs_sequence_in_name = 'obs_seq.sequential_prior_elm'
export atm_da_obs_sequence_in_name="obs_seq.out"
# Stage 1 native EAM setup. For the Stage 2 strongly coupled EAM cross pass,
# output/use sequential-prior settings are swapped and strongly_coupled is enabled.
export atm_da_output_sequential_prior_post=".true."
export atm_da_use_sequential_prior_post=".false."
export atm_da_strongly_coupled=".false."
# Strongly coupled Land-to-Atmosphere DA:
# update the atmospheric state using land-observation sequential priors.
export atm_da_state_model="Atmosphere"
export atm_da_obs_model="Land"

################################################################################
# Stage 1 native ELM profile. The bidirectional orchestrator changes these
# settings only inside the Stage 2 ELM atmosphere-cross subprocess.
################################################################################
export lnd_da_compute_posterior=".false."
export lnd_da_output_mean=".true."
export lnd_da_output_sd=".true."
export lnd_da_output_members=".true."
export lnd_da_inf_flavor_prior="0"
export lnd_da_inf_flavor_posterior="0"
export lnd_da_cutoff="0.4"
export lnd_da_spread_restoration=".false."
export lnd_da_sampling_error_correction=".false."
export lnd_da_horiz_dist_only=".true."
# Stage 1 native ELM input. For the Stage 2 Atmosphere-to-Land cross pass,
# this is replaced with the EAM atmospheric sequential-prior filename.
# Stage 1: obs_sequence_in_name = 'obs_seq.out'
# Stage 2: obs_sequence_in_name = 'obs_seq.sequential_prior_eam'
export lnd_da_obs_sequence_in_name="obs_seq.out"
# Stage 1 native ELM setup. For the Stage 2 strongly coupled ELM cross pass,
# output/use sequential-prior settings are swapped and coupling is enabled.
export lnd_da_output_sequential_prior_post=".true."
export lnd_da_use_sequential_prior_post=".false."
export lnd_da_strongly_coupled=".false."
# Stage 1 updates Land using land observations.
# Stage 2 overrides obs_model to Atmosphere, updating Land using
# atmospheric-observation sequential priors.
export lnd_da_state_model="Land"
export lnd_da_obs_model="Atmosphere"

################################################################################
# --- DART diagnostic defaults ------------------------------------------------
################################################################################
# Step 6 uses the configured experiment range when DIAG_START and DIAG_END are empty.
export my_eam_dart_diag_start="2011-12-01-00000"
export my_eam_dart_diag_end="2012-01-01-00000"
# The following values are Fortran logical strings consumed by DART tools.
export my_eam_dart_diag_use_custom_range=".true."
# Enable the independent Step 6 `closest` mode. With DIAG_MODE=all, closest is
# added after seq, obs, and common. Its summaries replace older results for the
# same quantity only after the new result completes successfully.
export my_eam_dart_diag_run_closest_member=".false."
export my_eam_dart_diag_run_obs2netcdf=".false."
export my_eam_dart_diag_run_obs_diag=".true."
