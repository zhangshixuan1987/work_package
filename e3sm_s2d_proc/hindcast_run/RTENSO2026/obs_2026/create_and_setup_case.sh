#initial time and ensemble size for the model 
export my_ensnum=10
export my_runtype="hybrid"
export my_workdir=${PWD}

# run configuration
export my_machine=pm-cpu
export my_project="e3sm"
export my_jobqueue="regular"
export my_task_per_node=128 # pm-cpu
export my_job_nnodes=270    # total nodes (nens x node per ens): 3 dates x 10 members x 9 nodes
export my_layout="custom-9_1x2_nyears"
export my_min_nodes_per_sim=9
#export my_layout="custom-28_1x2_nyears"
#export my_min_nodes_per_sim=28

export my_runopt="nyears" # simulation setup
export my_runn="2" # simulation length
export my_restopt="nyears" # restart setup
export my_restn="1" # restart setup
export my_initopt="YEARLY" #"MONTHLY"
export my_walltime="24:00:00"

export my_short_archive=true #short archive
export do_compile_setup_only=false

# Step 1 of restart: set true to reconfigure cases, then resubmit with do_ensemble_run=true
# See also block below for restart run configuration setup 
export do_continue_setup=false #true

# For short spin-up setup, when turned on the phase 1 will do a short-spinup with shorter ocean time step  
# see also block below for the other setup 
export do_short_spinup=true # must be false for continuation setup/run
export do_skip_completed_spinup=true # skip members with complete spinup restart files

# For final simulation stage, when short spinup is turned on this is phase 2  
export do_ensemble_run=true #run simulation

######################################################
# restart run configuration
# only used when do_continue_setup = true 
# note: restart ymd is corresponding to leadymd
######################################################
export my_restart_ymd=("2026-03-01" "2026-06-01" "2026-09-01")
export my_restart_tod="00000"
export my_continue_recompile=false #true # recompile because the layout changed
export my_continue_runopt="nyears"
export my_continue_runn="1"  
export my_continue_restopt="nyears" 
export my_continue_restn="1" 
export my_continue_walltime="12:00:00"
############################################

#model configuration (compset, resolution etc.)
export my_e3sm_code="/global/homes/z/zhan391/code/E3SMv3"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_compset="WCYCL20TR"
#export my_runpath="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
export my_runpath="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"
#change time step size if needed
export my_mpaso_dt="00:30:00" # default: "00:30:00"
export my_mpassi_dt=1800      # default: 1800

export my_leadymd=("2026-02-01" "2026-05-01" "2026-08-01")
export my_leadtod="00000"
export my_leadhh=`echo "${my_leadtod}" / 3600 | bc`

# RTENSO2026 initial conditions, one directory per component (../../ensemble_init/RTENSO2026):
#   <my_refdir>/YYYYMMDD/atm/, ocn/, lnd/<method>/, rof/<method>/, ice/<method>/, cpl/<method>/
# my_ic_method: land, river, and coupler method (clim, obs, hybrid_ic); in the directory name
# my_ice_method: sea-ice method, the same by default or nsidc (ice/nsidc/, make_ice_ic.bash)
export my_leadref="RTENSO2026"
export my_ic_method="obs"
export my_ice_method="${my_ic_method}"
export my_leadcase="${my_compset}_${my_resolution}_${my_leadref}_${my_ic_method}"
[[ ${my_ice_method} == "${my_ic_method}" ]] || my_leadcase="${my_leadcase}_ice_${my_ice_method}"
export my_refdir="/global/cfs/cdirs/e3sm/S2S2D/RTENSO2026/initial_conditions"
export my_refcase1="v3.LR.amip_0101.HICCUP.atm_era5" # for atmosphere
export my_refcase2="20251024_s2d_spinup" # for other

################################################
# spinup run configuration
################################################
export my_spinup_mpaso_dt="00:15:00" # default: "00:30:00"
export my_spinup_mpassi_dt=1800      # default: 1800
export my_spinup_runopt="nmonths"      # spinup stop option
export my_spinup_runn="1"            # spinup length
export my_spinup_restopt="nmonths"     # spinup restart option
export my_spinup_restn="1"           # spinup restart N
###################################################################

#diagnotic configuration
export my_eami_rgdmap="/global/homes/z/zhan391/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
export my_e3sm_rgdmap="/global/homes/z/zhan391/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/global/homes/z/zhan391/maps/map_r05_to_cmip6_180x360_aave.20200901.nc"
