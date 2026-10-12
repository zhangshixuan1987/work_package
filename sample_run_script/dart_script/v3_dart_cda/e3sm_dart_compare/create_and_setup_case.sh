# Configuration for the production E3SM-DART observation comparison.
# This file is sourced by 1_run_dart_diag_compare.sh.

# -----------------------------------------------------------------------------
# Machine and launch settings
# -----------------------------------------------------------------------------
export my_machine="compy"
export my_project="esmd"
export my_jobqueue="short"
export my_walltime="12:00:00"
export my_task_per_node=40
export my_job_nnodes=1
export my_job_ntasks=24

# -----------------------------------------------------------------------------
# Production experiments
# -----------------------------------------------------------------------------
# Each workflow experiment supplies its authoritative case name, ensemble size,
# DART run directory, machine environment, and diagnostics namelist through its
# own create_and_setup_case.sh. The production driver prepends the archived 2011
# CAM6 assimilation as an additional comparison input. Order is significant and
# is preserved in DART's obs_common_subset input lists.
comparison_experiment_roots=(
  "/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/1_ne30pg2_ctrl_atm_en10"
  "/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/1_ne30pg2_dart_atm_en10"
  "/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/1_ne30pg2_dart_atm_en20"
  "/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/1_ne30pg2_dart_atm_en40"
  "/qfs/people/zhan391/e3sm_dart_work/run_script/v3_dart_cda/3_ne30pg2_dart_cpl_en40"
)

# -----------------------------------------------------------------------------
# Shared comparison range and stages
# -----------------------------------------------------------------------------
# Both endpoints are inclusive and must be aligned to my_dart_window.
export my_dart_diag_ymds="2011-12-01-00000"
export my_dart_diag_ymde="2011-12-31-64800"
export my_dart_window=6

# Generate *.common sequences, optionally convert them to NetCDF, then run
# observation-space diagnostics.
export my_dartdiag_common=".true."
export my_dartdiag_obs2nc=".false."
export my_dartdiag_obdiag=".true."

# -----------------------------------------------------------------------------
# DART tools and shared EAM grid assets
# -----------------------------------------------------------------------------
export my_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART"
export my_dart_eam="eam-se"
export my_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"
export my_e3sm_topo="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_e3sm_semap="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_e3sm_csgrid="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_dart_pgrid=".true."
export my_dart_ntask=24

# Timeline values retained for the legacy date-reporting code in the driver.
export my_dartymds="2011-12-01"
export my_darttods="00000"
export my_dart_cycle=124
