export my_task_per_node=40
export my_job_nnodes=40
export my_job_ntasks=40
export my_layout="custom-2_1x6_nhours"
export my_machine=compy
export my_project="esmd"
export my_jobqueue="short"
export my_walltime="02:00:00"
export my_date=`date +"%m-%d-%y"`

#initial time and ensemble size for the model 
export my_ensnum=40
export my_casedate="2011-12-01"
export my_casetod="00000"

#model configuration (compset, resolution etc.)
export my_e3sm_code="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
export my_runtype="AMIP" # or "Full-CPL"
export my_compset="F20TR"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/compyfs/zhan391/v3_dart_cda_scratch"
export my_casename="DARTEN${my_ensnum}_${my_compset}_${my_resolution}_${my_machine}"

#my_modeldir is the directory for dart simulation. 
#my_modelexe is the pre-existing e3sm.exe file. We only need to compile e3sm once for ensemble #1
export my_modeldir="${my_runpath}/${my_casename}"
export my_modelcase="${my_modeldir}/case_scripts"
export my_modelexe="${my_modeldir}/build/e3sm.exe"
export my_e3sm_topo="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_e3sm_semap="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_e3sm_csgrid="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_e3sm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"

#single initial condition files to start the simulation 
#need a sample eam.i file from existing model simulation
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

#dart configuration 
export my_dart_cycle=8
export my_dart_arcdir="archive"
export my_dart_window=6 # 6hourly cycling da 
export my_dart_code="/qfs/people/zhan391/e3sm_dart_work/code/DART"
export my_ensnum="${my_ensnum}"
export my_dart_runpath="${my_modeldir}/archive"
export my_dart_eam="eam-se"
export my_dart_obsdir="/compyfs/zhan391/acme_init/Observations/NCEP+ACARS+GPS+AIRS"

#configuration for DART cycling data assimilation 
export my_dartymds=${my_casedate} #start ymd of eam-dart da
export my_darttods=${my_casetod}  #start ymd of eam-dart da
export my_dartymde="2012-01-01"   #end ymd of dart da
export my_darttode="00000"        #end tod of dart da

#configuration for DART diagnostics 
#export my_dart_diag_ymds="2011-12-16-00000"
#export my_dart_diag_ymde="2012-01-01-00000"
#export my_dart_diag_ymds="2011-12-01-00000"
#export my_dart_diag_ymde="2012-01-01-00000"
export my_dartdiag_custom=".true."

#dart run configuration
export my_dart_ntask=$((my_task_per_node * my_job_nnodes))
export my_dart_pgrid=".true."
export my_dart_queue=${my_jobqueue}
export my_dart_project=${my_project}
export my_dart_machine=${my_machine}
