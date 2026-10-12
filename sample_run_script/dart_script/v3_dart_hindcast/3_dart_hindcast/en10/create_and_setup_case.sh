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
export my_jobqueue="short" #"slurm"
export my_walltime="02:00:00" #"04:00:00"
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
export my_casedate="2012-01-01"
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
# --- DA analysis the hindcast starts from ------------------------------------
# Restarts: ${my_refdir}/ENxx/archive/rest/${my_casedate}-${my_casetod}
# (the older shared layout ${my_refdir}/archive/rest/... is also accepted).
################################################################################
export my_ic_source="dart_analysis"   # dart_analysis | perturbed (Steps 0a/0b, see README)
export my_refcase="DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy"
export my_refdir="${my_runpath}/${my_refcase}"

################################################################################
# --- Hindcast case -----------------------------------------------------------
# Named after the DA case (DARTEN10_F20TR_... -> DARTEN10_15day_F20TR_...), so a
# copied workflow pointing at a new DA run gets its own case. Set my_casetag
# for a second hindcast from the same analysis.
################################################################################
export my_casetag="DARTEN${my_ensnum}"
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
