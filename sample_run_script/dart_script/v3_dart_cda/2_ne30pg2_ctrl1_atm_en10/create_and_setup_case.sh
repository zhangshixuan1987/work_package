# Shared configuration for all numbered workflow stages.
# This file is sourced; it should define settings without running workflow work.

################################################################################
# --- Workflow and runtime paths ---------------------------------------------
# Update this absolute root if the workflow directory is moved.
################################################################################
export my_workflow_root="/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/2_ne30pg2_ctrl1_atm_en10"
export my_workflow_lib="${my_workflow_root}/workflow_lib"
export my_runtime_dir="${my_workflow_root}/runtmp"
export my_log_dir="${my_runtime_dir}/logs"
export my_status_dir="${my_runtime_dir}/status"
export my_lock_dir="${my_runtime_dir}/locks"
export my_handoff_dir="${my_runtime_dir}/handoff"
export my_run_script_dir="${my_runtime_dir}/run_scripts"

################################################################################
# Step 4 --nodes must match my_job_nnodes; other stages use their own allocations.
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
export my_casedate="2012-05-01"
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
export my_casename="CTRFEN${my_ensnum}S1_${my_compset}_${my_resolution}_${my_machine}"

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
# --- Shared E3SM cycling controls -------------------------------------------
# Live completed-cycle counter. Step 4 updates this value transactionally.
################################################################################
export my_e3sm_completed_cycles=0
export my_e3sm_cycle_hours=6
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
# Fixed workflow invariant; all archive readers and writers use ENxx/archive.
export my_dart_root="${my_modeldir}/dart_en$(printf '%02d' "${my_ensnum}")"
export my_e3sm_start_date=${my_casedate}
export my_e3sm_start_tod=${my_casetod}
export my_e3sm_end_date="2012-06-01"
export my_e3sm_end_tod="00000"

################################################################################
# --- Continuous free run (Step 4: 4_run_e3sm_continuous.sh) ------------------
# This case has no DA cycling. After Steps 1-3, every member runs one continuous
# forecast from the perturbed IC at my_e3sm_start_date to my_e3sm_end_date in a
# single model run (restart written only at the end time). No eam.i (inithist)
# files are written. my_free_layout selects how members are run:
#   member_jobs  each member is its own CIME batch job on my_free_member_nodes
#                nodes in my_free_job_queue (~40 min per month on 20 nodes)
#   single_job   one my_job_nnodes-node allocation in my_free_single_job_queue
#                runs all members together, my_nodes_per_member nodes each
#                (several hours per month on 2 nodes per member)
################################################################################
export my_free_layout="member_jobs"
export my_free_member_nodes=20
export my_free_job_queue="short"
export my_free_job_walltime="02:00:00"
export my_free_single_job_queue="slurm"
export my_free_single_job_walltime="24:00:00"
# Delete restart sets between start and end at finalize (TRUE keeps them).
export my_free_keep_intermediate_restarts="FALSE"

################################################################################
# --- EAM DART assimilation --------------------------------------------------
# Step 4 derives component nodes from my_job_nnodes and both DA switches:
# both enabled = equal halves; one enabled = all nodes; both disabled = no DA launch.
################################################################################
export my_eam_dart_da="on"
export my_eam_dart_cycle_hours=6
export my_eam_dart_end_date="${my_e3sm_end_date}"
export my_eam_dart_end_tod="${my_e3sm_end_tod}"
export my_eam_dart_run_dir="${my_dart_root}/eam"
export my_eam_dart_model="eam-se"
export my_eam_dart_pgrid=".true."
export my_eam_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART"
# Lowest EAM model level at which observations may be assimilated; levels
# above it (smaller indices) are excluded near the diffusive model top.
export my_eam_no_obs_assim_above_level=5
export my_eam_use_log_vertical_scale=".true."
export my_eam_vert_normalization_scale_height="1.5"
# Use a workflow-owned loader rather than a generated file in DARTs work tree.
export my_eam_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"
# Optional per-cycle EAM DA settings. Keys combine the exact valid time and
# parameter name. Omitted parameters retain their normal defaults.
declare -Ag my_eam_cycle_overrides=(
  # Cycle: 2011-12-28-21600
  ["2011-12-28-21600:localization_cutoff"]="0.1"
  ["2011-12-28-21600:inflation_damping"]="0.5"
  ["2011-12-28-21600:no_obs_assim_above_level"]="7"
  # Cycle: 2011-12-28-43200
  ["2011-12-28-43200:localization_cutoff"]="0.1"
  ["2011-12-28-43200:inflation_damping"]="0.5"
  ["2011-12-28-43200:no_obs_assim_above_level"]="7"
  # Cycle: 2011-12-28-64800
  ["2011-12-28-64800:localization_cutoff"]="0.1"
  ["2011-12-28-64800:inflation_damping"]="0.5"
  ["2011-12-28-64800:no_obs_assim_above_level"]="7"
)

################################################################################
# --- ELM DART assimilation --------------------------------------------------
# ELM DA links single-record h1 history and h2 vector files from each member archive.
################################################################################
export my_elm_dart_da="off"
export my_elm_dart_cycle_hours=6
export my_elm_dart_end_date="${my_e3sm_end_date}"
export my_elm_dart_end_tod="${my_e3sm_end_tod}"
export my_elm_dart_run_dir="${my_dart_root}/elm"
export my_elm_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART_SCP"
export my_elm_dart_obsdir="/compyfs/zhan391/acme_init/Observations/SMAP"
export my_elm_dart_model="elm"
export my_elm_history_stream="h1"
export my_elm_vector_history_stream="h2"
declare -Ag my_elm_dart_cycle_overrides=(
)

################################################################################
# --- Strongly coupled DA setup ----------------------------------------------
################################################################################
export strongly_coupled_on="off"

export atm_da_compute_posterior=".false."
export atm_da_output_sequential_prior_post=".false."
export atm_da_use_sequential_prior_post=".false."
export atm_da_output_mean=".true."
export atm_da_output_sd=".true."
export atm_da_output_members=".true."
export atm_da_strongly_coupled=".false."
export atm_da_state_model="Atmosphere"
export atm_da_obs_model="Atmosphere"

export lnd_da_output_sequential_prior_post=".false."
export lnd_da_use_sequential_prior_post=".false."
export lnd_da_perturb_from_single_instance=".true."
export lnd_da_perturbation_amplitude="0.2"
export lnd_da_perturbation_method="uniform"
export lnd_da_obs_sequence_in_name="obs_seq.out"
export lnd_da_inf_flavor_prior="0"
export lnd_da_inf_flavor_posterior="0"
export lnd_da_cutoff="0.4"
export lnd_da_spread_restoration=".false."
export lnd_da_sampling_error_correction=".false."
export lnd_da_horiz_dist_only=".true."
export lnd_da_strongly_coupled=".false."
export lnd_da_state_model="Land"
export lnd_da_obs_model="Land"

################################################################################
# --- EAM DART diagnostic defaults ------------------------------------------
# Step 6 uses this range when DIAG_START and DIAG_END are empty.
################################################################################
#export my_eam_dart_diag_start="2012-05-01-00000"
export my_eam_dart_diag_start="2012-05-16-00000"
export my_eam_dart_diag_end="2012-06-01-00000"
# Diagnostic switches use Fortran logical strings because the workers pass them to DART tools.
export my_eam_dart_diag_use_custom_range=".true."
# Enable the independent Step 6 `closest` mode. With DIAG_MODE=all, closest is
# added after seq, obs, and common. Its summaries replace older results for the
# same quantity only after the new result completes successfully.
export my_eam_dart_diag_run_closest_member=".true."
export my_eam_dart_diag_closest_hours=24
export my_eam_dart_diag_run_obs2netcdf=".false."
export my_eam_dart_diag_run_obs_diag=".true."
