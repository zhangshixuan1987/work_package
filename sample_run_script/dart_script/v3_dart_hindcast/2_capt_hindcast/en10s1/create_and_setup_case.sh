# Shared configuration for all numbered hindcast workflow stages.
# This file is sourced; it should define settings without running workflow work.
#
# New DA experiment: copy this directory (without runtmp/), then edit at least
# my_refcase, my_ensnum and my_casedate (and my_runtype/my_compset/my_job_nnodes).

################################################################################
# --- Workflow and runtime paths ---------------------------------------------
# The root is the directory holding this file, so a copied workflow needs no edit.
################################################################################
export my_workflow_root="$(cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")" && pwd)"
export my_workflow_lib="${my_workflow_root}/workflow_lib"
export my_runtime_dir="${my_workflow_root}/runtmp"
export my_log_dir="${my_runtime_dir}/logs"
export my_status_dir="${my_runtime_dir}/status"
export my_lock_dir="${my_runtime_dir}/locks"
export my_run_script_dir="${my_runtime_dir}/run_scripts"

################################################################################
# --- Machine and Slurm defaults (model run of each member) -------------------
################################################################################
export my_machine=compy
export my_project="esmd"
export my_jobqueue="short"
export my_walltime="02:00:00"
export my_task_per_node=40
export my_job_nnodes=10            # nodes per ensemble member
export my_max_parallel_setup=4     # members created/configured concurrently by Steps 1-2
export my_e3sm_unified="/share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh"  # Step 4 (NCO)
export my_machine_env_file="${my_workflow_lib}/env/env_${my_machine}_specific.sh"       # Step 5 (nccopy)

################################################################################
# --- Ensemble, initial time and forecast length ------------------------------
################################################################################
export my_ensnum=10
export my_enstart=1
export my_casedate="2012-06-01"
export my_casetod="00000"
export my_leadtime="15day"
export my_stop_option="nmonths"    # ndays, nmonths or nyears
export my_stop_n=2
export my_rest_option="nmonths"
export my_rest_n=1

################################################################################
# --- E3SM experiment identity -------------------------------------------------
# my_runtype must be AMIP (F20TR) or Full-CPL (WCYCL20TR); the case template is
# workflow_lib/run_template/compile_and_setup_e3sm.${my_compset}.sh
################################################################################
export my_e3sm_code="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
export my_runtype="AMIP"
export my_compset="F20TR"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/compyfs/zhan391/v3_dart_cda_scratch"

################################################################################
# --- Initial conditions -------------------------------------------------------
# Perturbed ICs (Steps 0a/0b): see the block at the end of this file.
################################################################################
export my_ic_source="perturbed"       # dart_analysis | perturbed (Steps 0a/0b, see README)

################################################################################
# --- Hindcast case -----------------------------------------------------------
# Named after the DA case (DARTEN10_F20TR_... -> DARTEN10_15day_F20TR_...), so a
# copied workflow pointing at a new DA run gets its own case. Set my_casetag
# for a second hindcast from the same analysis.
################################################################################
export my_casetag="CAPTEN10S1"
export my_casename="${my_casetag}_${my_leadtime}_${my_compset}_${my_resolution}_${my_machine}"
export my_modeldir="${my_runpath}/${my_casename}"   # hindcast case root (members under it)

################################################################################
# --- Ocean and sea-ice initial condition (Full-CPL only) ---------------------
################################################################################
export my_refmpas="20241109.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.pm-cpu"
export my_refdir1="/compyfs/zhan391/acme_init/E3SMv3_INT/${my_casedate}-${my_casetod}"
export my_refocn_in="${my_refdir1}/${my_refmpas}.mpaso.rst.${my_casedate}_${my_casetod}.nc"
export my_refice_in="${my_refdir1}/${my_refmpas}.mpassi.rst.${my_casedate}_${my_casetod}.nc"

################################################################################
# --- Perturbed initial conditions (my_ic_source=perturbed: Steps 0a/0b) -------
# Step 0a stages, for every member, the EAM state my_refeam_in with
# PS,U,V,T,Q,CLDLIQ,CLDICE taken from my_refeam_ic, plus the ELM/MOSART/MPAS/cpl
# restarts below, under ${my_modeldir}/ENxx/archive/rest/<date>; Step 0b perturbs
# the EAM state with DART. The hindcast then starts from its own case.
################################################################################
export my_refcase="${my_casename}"
export my_refdir="${my_modeldir}"
export my_refdate="${my_casedate}"
export my_reftod="${my_casetod}"
export my_refeam_in="/compyfs/zhan391/v3_dart_cda_scratch/CTRL_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/init/2012-01-01-00000/v3.LR.amip_0101.eam.i.2012-01-01-00000.nc"
export my_refeam_ic="${my_refdir1}/v3.LR.ne30pg2.ERA5.eam.i.${my_casedate}-${my_casetod}.nc"
export my_refelm_in="${my_refdir1}/${my_refmpas}.elm.r.${my_casedate}-${my_casetod}.nc"
export my_refrof_in="${my_refdir1}/${my_refmpas}.mosart.r.${my_casedate}-${my_casetod}.nc"
export my_refcpl_in="${my_refdir1}/${my_refmpas}.cpl.r.${my_casedate}-${my_casetod}.nc"
# DART perturbation (Step 0b)
export my_conda_setup_file="/people/zhan391/miniforge3/etc/profile.d/conda.sh"
export my_analysis_conda_env="e3sm_analysis"
export my_dart_root="${my_modeldir}/dart_en$(printf '%02d' "${my_ensnum}")"
export my_eam_perturb_nml="${my_workflow_lib}/namelists/eam/perturb.nml"
export my_eam_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART"
export my_eam_dart_model="eam-se"
export my_eam_dart_run_dir="${my_dart_root}/eam"
export my_eam_dart_cycle_hours=6
export my_eam_dart_pgrid=".true."
export my_eam_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"
export my_eam_topography_file="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_eam_se_mapping_file="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_eam_cs_grid_file="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"

################################################################################
# --- Forecast validation (Step 3) ---------------------------------------------
# History streams that must exist for every day / month of the forecast before
# a member is marked complete ("component:stream", component = archive subdir).
################################################################################
export my_expected_daily_streams="atm:eam.h1 atm:eam.h2"
export my_expected_monthly_streams="atm:eam.h0 lnd:elm.h0"

################################################################################
# --- Post-processing defaults (Step 4) ---------------------------------------
################################################################################
export my_post_start="${my_casedate}"                                    # first day
export my_post_end=$(date -d "${my_casedate} + ${my_stop_n} ${my_stop_option#n}" +%F)   # exclusive
export my_e3sm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20231110.nc"
