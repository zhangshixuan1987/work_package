export my_task_per_node=40
export my_job_nnodes=20
export my_job_ntasks=40
export my_layout="custom-10_1x1_nmonths"
export my_machine=compy
export my_project="esmd"
export my_jobqueue="short"
export my_walltime="02:00:00"
export my_date=`date +"%m-%d-%y"`

#initial time and ensemble size for the model 
export my_ensnum=20
export my_enstart=1
export my_casedate="2012-01-01"
export my_casetod="00000"
export my_leadtime="15day" 

#model configuration (compset, resolution etc.)
export my_e3sm_code="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
export my_runtype="AMIP" # or "Full-CPL"
export my_compset="F20TR"
export my_resolution="ne30pg2_r05_IcoswISC30E3r5"
export my_runpath="/compyfs/zhan391/v3_dart_cda_scratch"
export my_casename="DARTEN${my_ensnum}_${my_leadtime}_${my_compset}_${my_resolution}_${my_machine}"

#my_modeldir is the directory for dart simulation. 
#my_modelexe is the pre-existing e3sm.exe file. 
#We only need to compile e3sm once for ensemble #1
export my_modeldir="${my_runpath}/${my_casename}"
export my_modelcase="${my_modeldir}/case_scripts"
export my_modelexe="${my_modeldir}/build/e3sm.exe"
export my_e3sm_topo="/compyfs/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4pg2_x6t-SGH.c20210614.nc"
export my_e3sm_semap="/qfs/people/zhan391/e3sm_dart_work/code/HOMME/SEMapping.nc"
export my_e3sm_csgrid="/qfs/people/zhan391/e3sm_dart_work/code/DART/models/eam-se/work/SEMapping_cs_grid_NE30.nc"
export my_e3sm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"
export my_elm_rgdmap="/compyfs/zhan391/acme_init/map_file/map_r05_to_cmip6_180x360_aave.20200901.nc"

#initial condition files to start the simulation 
export my_refcase="DARTEN20_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy"
export my_refdate=${my_casedate}
export my_reftod=${my_casetod}
export my_refdir="/compyfs/zhan391/v3_dart_cda_scratch/${my_refcase}/archive"
