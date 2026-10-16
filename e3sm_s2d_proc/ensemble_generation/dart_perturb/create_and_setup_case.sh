#initial time and ensemble size for the model 
export my_ensnum=10
export my_casedate="1980-05-01"
export my_casetod="00000"

#model configuration (compset, resolution etc.)
export my_e3sm_code="/pscratch/sd/z/zhan391/e3sm_project/code/E3SMv3"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D" 

export my_e3sm_topo="/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_e3sm_semap="/pscratch/sd/z/zhan391/e3sm_project/code/HOMME/SEMapping.nc"
export my_e3sm_csgrid="/pscratch/sd/z/zhan391/e3sm_project/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_e3sm_rgdmap="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D/maps/map_r05_to_cmip6_180x360_aave.20200901.nc"

#single initial condition files to start the simulation 
export my_refcase1="v3.LR.ne30pg2.ERA5"
export my_refdate1=${my_casedate}
export my_reftod1=${my_casetod}
export my_refdir1="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.S2D/${my_refdate1}-${my_reftod1}"
export my_refeam_ic="${my_refdir1}/${my_refcase1}.eam.i.${my_refdate1}-${my_reftod1}.nc"
#export my_refeam_in="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.amip_0101/1980-01-01-00000/v3.LR.amip_0101.eam.i.1980-01-01-00000.nc"

#the spin-up below start from 0 for 1958-2023
export my_refcase2="20251024_s2d_spinup"
export my_refdir2="/global/cfs/cdirs/e3sm/lvroekel/biweekly_restarts/archive/rest/${my_refdate2}-${my_reftod2}"
export my_refelm_in="${my_refdir2}/${my_refcase2}.elm.r.${my_refdate2}-${my_reftod2}.nc"
export my_refrof_in="${my_refdir2}/${my_refcase2}.mosart.r.${my_refdate2}-${my_reftod2}.nc"
export my_refocn_in="${my_refdir2}/${my_refcase2}.mpaso.rst.${my_refdate2}_${my_reftod2}.nc"
export my_refice_in="${my_refdir2}/${my_refcase2}.mpassi.rst.${my_refdate2}_${my_reftod2}.nc"
export my_refcpl_in="${my_refdir2}/${my_refcase2}.cpl.r.${my_refdate2}-${my_reftod2}.nc"

export my_dart_code="/pscratch/sd/z/zhan391/e3sm_project/code/DART"
export my_dart_case="${my_refcase1}"
export my_dart_runpath="${my_refdir1}"
export my_dart_eam="eam-se"
export my_dartymds="${my_casedate}" #start ymd of eam-dart da
export my_darttods="${my_casetod}"  #start ymd of eam-dart da
export my_dartymde="${my_casedate}" #end ymd of dart da
export my_darttode="${my_casetod}"  #end tod of dart da
export my_dart_env="${my_dart_code}/models/eam-se/work/env_mach_specific.sh"

#dart run configuration
export my_task_per_node=64
export my_job_nnodes=4
export my_dart_ntask=$((my_task_per_node * my_job_nnodes))
export my_dart_pgrid=".true."
export my_dart_queue=${my_jobqueue}
export my_dart_project=${my_project}
export my_dart_machine=${my_machine}
