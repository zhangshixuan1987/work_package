#initial time and ensemble size for the model 
export my_ensnum=9
export my_casedate="${my_casedate:-1980-02-01}"
export my_casetod="00000"

#model configuration (compset, resolution etc.)
export my_e3sm_code="/global/homes/z/zhan391/code/E3SMv3"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="${my_runpath:-/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT}"
########################################################
# my_fosopt: must be the name below:
#            BruteForce or Coupled-DA or 
#            JRA55-FOSIRL or ERA5-FOSIRL
#my_fosopt: ERA5 or JRA55 
#           needed to determine initial conditions for
#           sea-ice land and river for BruteForce case
#my_eamopt: ERA5-Tang or CUSTOMIZE
#           ERA5-Tang is default, and CUSTOMIZE option 
#           needs extra setup in later section    
########################################################
export my_runcase="${my_runcase:-JRA55-FOSIRL}"
export my_fosopt="${my_fosopt:-JRA55}"
export my_eamopt="${my_eamopt:-ERA5-Tang2}"
export my_ocnopt="${my_ocnopt:-}"
export my_overwrite="${my_overwrite:-false}"

export my_e3sm_topo="/global/cfs/cdirs/e3sm/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_e3sm_semap="/global/homes/z/zhan391/code/HOMME/SEMapping.nc"
export my_e3sm_csgrid="/global/homes/z/zhan391/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_eami_rgdmap="/global/homes/z/zhan391/maps/map_ne30np4_to_cmip6_180x360_aave.20181001.nc"
export my_e3sm_rgdmap="/global/homes/z/zhan391/maps/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/global/homes/z/zhan391/maps/map_r05_to_cmip6_180x360_aave.20200901.nc"

#atmosphere initial condition files (ERA5-Tang1)
export my_atmcase1="v3.LR.amip_0101.HICCUP.atm_era5"
export my_atmdate1=${my_casedate}
export my_atmtod1="00"
export my_atmdir1="/global/cfs/cdirs/e3sm/tang30/S2D/atmIC_ERA5" 
export my_eam_ic1="${my_atmdir1}/${my_atmcase1}.${my_atmdate1}.${my_atmtod1}.ne30np4.L80.nc"

#atmosphere initial condition files (ERA5-Tang2)
export my_atmcase2="v3.LR.amip_0101.HICCUP.atm_era5"
export my_atmdate2=${my_casedate}
export my_atmtod2="00"
export my_atmdir2="/pscratch/sd/t/tang30/share/ForShixuan/ERA5_IC"
export my_eam_ic2="${my_atmdir2}/${my_atmcase2}.${my_atmdate1}.${my_atmtod1}.ne30np4.L80.nc"

##atmosphere initial condition files (CUSTOMIZE)
#export my_atmcasex="v3.LR.amip_0101.HICCUP.atm_era5"
#export my_atmdatex=${my_casedate}
#export my_atmtodx="00"
#export my_atmdirx="/pscratch/sd/t/tang30/share/ForShixuan/ERA5_IC"
#export my_eam_icx="${my_atmdir2}/${my_atmcase2}.${my_atmdate1}.${my_atmtod1}.ne30np4.L80.nc"

# template file for eam.i (from existing simulations)
export my_eami_template="/pscratch/sd/z/zhan391/e3sm_project/input/v3.LR.amip_0101/1980-01-01-00000/v3.LR.amip_0101.eam.i.1980-01-01-00000.nc"

#fosirl ocean/land initial condition files (JRA55-forced)
export my_foscase1="20251024_s2d_spinup"
export my_fosrefy1="1980"
export my_fossimy1="0353"
#export my_refdir1="/global/cfs/cdirs/e3sm/lvroekel/biweekly_restarts/archive/rest
export my_fosdir1="/global/cfs/cdirs/e3sm/lvroekel/archive/rest"

#fosirl ocean/land initial condition files (ERA5-forced,plance holder)
export my_foscase2=""
export my_fosrefy2=""
export my_fossimy2=""
export my_fosdir2=""

#bruteforce ocean initial condition files (Ben Kirtman)
export my_ocndir1="/pscratch/sd/l/lvroekel/kirtman_restarts"
export my_ocncase1="OceanAnal_Restart_%s_E3SM_Grid.nc"

#bruteforce ocean initial condition files (Luke Van Roekel, selected years)
export my_ocndir2="/pscratch/sd/l/lvroekel/kirtman_restarts/new_test"
export my_ocncase2="20251024_s2d_spinup.mpaso.rst.%s_00000.nc"

#bruteforce CFSRv2 ocean initial condition files (Luke Van Roekel, 1980-2018)
export my_ocndir3="/pscratch/sd/l/lvroekel/ics"
export my_ocncase3="20251024_s2d_spinup.mpaso.rst.%s_00000.cfsrv2.nc"

#4DEnVar ocean/land/atmosphere initial-condition files
export my_daroot="/global/cfs/cdirs/e3sm/zhan391/v3.LR.S2D.ENSINT/4DEnVar_ICs_SHIPF"
export my_dacase="assm_anom_ocn_4DEnVar_loc100_v3"
export my_daocn="${my_dacase}.mpaso.rst.%s_00000_new.nc"
export my_dalnd="${my_dacase}.elm.r.%s-00000.nc"
export my_daeam_1980s="WCYCL20TR_ne30pg2_r05_IcoswISC30E3r5_4DEnVar.eam.i.%s-01800.nc"
export my_daeam="${my_dacase}.eam.i.%s-00000.nc"

export my_dart_code="/global/homes/z/zhan391/code/DART"
export my_dart_case="${my_atmcase1}"
export my_dart_runpath="${my_runpath}/${my_runcase}/${my_atmdate1}-00000"
export my_dart_eam="eam-se"
export my_dartymds="${my_casedate}" #start ymd of eam-dart da
export my_darttods="${my_casetod}"  #start ymd of eam-dart da
export my_dartymde="${my_casedate}" #end ymd of dart da
export my_darttode="${my_casetod}"  #end tod of dart da
export my_dart_env="${my_dart_code}/models/eam-se/mach_env/env_pm-cpu_specific.sh"
export my_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"

#dart run configuration
export my_task_per_node=64
export my_job_nnodes=1
export my_dart_ntask=$((my_task_per_node * my_job_nnodes))
export my_dart_pgrid=".true."
export my_dart_queue="${my_jobqueue:-}"
export my_dart_project="${my_project:-}"
export my_dart_machine="${my_machine:-}"
