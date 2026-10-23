#initial time and ensemble size for the model 
export my_ensnum=10
export my_runtype="hybrid"
export my_workdir=${PWD}

# run configuration
export my_machine=pm-cpu
export my_project="e3sm"
export my_jobqueue="regular"
export my_task_per_node=128 # pm-cpu
export my_job_nnodes=280    # total nodes (nens x node per ens)
export my_layout="custom-7_1x2_nyears"
export my_min_nodes_per_sim=7

export my_runopt="nyears" # simulation setup
export my_runn="2" # simulation length
export my_restopt="nyears" # restart setup
export my_restn="1" # restart setup
export my_initopt="YEARLY" #"MONTHLY"
export my_walltime="24:00:00"

export my_short_archive=true #short archive
export do_compile_setup_only=false #true
export do_continue_setup=false #see also restart setup
export do_ensemble_run=true  #run simulation 

######################################################
# restart run configuration
# only used when do_continue_run = True 
# note: restart ymd is corresponding to leadymd
######################################################
export my_restart_ymd=("" "" "" "")
export my_restart_tod="00000"
export my_continue_recompile=false # recompile if use a different layout
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
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
#export my_runpath="/global/cfs/cdirs/e3smdata/simulations/S2S2D"
export my_runpath="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_S2D"

export my_leadymd=("2017-05-01" "2017-11-01" "2018-05-01" "2018-11-01")
export my_leadtod="00000"
export my_leadhh=`echo "${my_leadtod}" / 3600 | bc`

# option: BruteForce or Coupled-DA or JRA55-FOSIRL or ERA5-FOSIRL 
export my_leadref="JRA55-FOSIRL" 
export my_leadcase="${my_compset}_${my_resolution}_${my_leadref//-/_}"
export my_refdir="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/${my_leadref}"
export my_refcase1="v3.LR.amip_0101.HICCUP.atm_era5" # for atmosphere
export my_refcase2="20251024_s2d_spinup" # for other 

#diagnotic configuration
export my_eami_rgdmap="/global/homes/z/zhan391/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
export my_e3sm_rgdmap="/global/homes/z/zhan391/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/global/homes/z/zhan391/maps/map_r05_to_cmip6_180x360_aave.20200901.nc"
