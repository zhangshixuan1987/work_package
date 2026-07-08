#!/bin/csh
date
setenv PROJECT         'ESMD'
setenv CESM_PROJ       'ESMD'
setenv CESM_EMAIL      'shixuan.zhang@pnnl.gov'

setenv MACH            'compy'
setenv COMPILER        'intel'
setenv RESOLUTION      'ne30_ne30'
setenv NTASKS_PER_INST '1600'
setenv NINST           '1'
setenv NTHRDS         '1'
         
setenv CCSMTAG  'E3SM_MAINT1.0'
setenv COMPSET  'FC5AV1C-04P2'
setenv CCSMROOT '/compyfs/zhan391/code/E3SM_MAINT1.0'
setenv PTMP     '/compyfs/zhan391/e3sm_simulation/e3sm_mam_cnvg'
setenv EXEDIR   '/compyfs/zhan391/e3sm_simulation/e3sm_mam_cnvg/exe//compile_E3SM_MAINT1.0_FC5AV1C-04P2_ne30_ne30_compy_intel_1600proc'
setenv CSMDATA  '/compyfs/inputdata/'
         
setenv CASE      'e3smv1_tsc_ne30_ne30_E3SM_MAINT1.0_FC5AV1C-04P2_intel_compy_DT0030_06_fullmodel_rad'
setenv CASEROOT  '/compyfs/zhan391/e3sm_simulation/e3sm_mam_cnvg/cases/e3smv1_tsc_ne30_ne30_E3SM_MAINT1.0_FC5AV1C-04P2_intel_compy_DT0030_06_fullmodel_rad'
setenv RUNDIR    '/compyfs/zhan391/e3sm_simulation/e3sm_mam_cnvg/run/e3smv1_tsc_ne30_ne30_E3SM_MAINT1.0_FC5AV1C-04P2_intel_compy_DT0030_06_fullmodel_rad'

setenv postCIME 5

set nbundles = -1         
set case_setup_script = '/compyfs/zhan391/run_script/e3sm_cnvg_202007/create_and_setup_bundled_case.csh' 

set dtime    = 30
set atm_init = /compyfs/zhan391/acme_init/FC5AV1C-04P2_init_201907/FC5AV1C-04P2_ne30_ne30_intel_cori-knl.cam.i.0001-06-01-00000.nc
set lnd_init = /compyfs/zhan391/acme_init/FC5AV1C-04P2_init_201907/FC5AV1C-04P2_ne30_ne30_intel_cori-knl.clm2.r.0001-06-01-00000.nc

set ncycle   = 1          # 1 cycle in total, i.e., no restart

@ nlen = 2 * 3600 / 30   # run model for 60 minutes.
set stop_o = 'nsteps'

@ nhtfrq  = 1800 / 30   # output every 30 minutes
set mfilt = 1               # 1 time step per file 

@ nhtfrq1 = 60  / 30     # output every 24s 
@ mfilt1  = 360 / 60         # 1 time step per file 
#@ nhtfrq1  = 1
#@ mfilt1   = 360 / 30     # 1 time step per file

set nl_file = /compyfs/zhan391/run_script/e3sm_cnvg_202007/namelist_files/cam_nl_e3smv1_tsc
#====================================================================
# Create and set up new case. No need to build the model. 
#====================================================================
source ${case_setup_script}

cd $CASEROOT

#cp $MODSROOT/*90 SourceMods/src.cam/
#cp $MODSROOT/namelist_definition.xml   $CCSMROOT/models/atm/cam/bld/namelist_files/namelist_definition.xml
##cp $MODSROOT/namelist_defaults_cam.xml $CCSMROOT/models/atm/cam/bld/namelist_files/namelist_defaults_cam.xml

./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

#-----------------------------------------
# Runtime options: edit env_run.xml
#-----------------------------------------
cd $CASEROOT

#./xmlchange  -file env_run.xml -id RUN_STARTDATE -val 'startdate_dummy'

@ nresub = $ncycle - 1

./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o 
./xmlchange  -file env_run.xml -id  REST_N       -val 1
./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'ndays' 
./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

@ ncpl = 86400 / $dtime

./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "04:00:00"
#./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val  "short"  #"slurm"

#--------------------
# Namelist variables
#--------------------
cat > user_nl_cam <<EOF
 empty_htapes       = .true.
 avgflag_pertape    = 'I',
 nhtfrq             = $nhtfrq,
 mfilt              = $mfilt,
 ndens              = 1,
 iradsw             = 2,
 iradlw             = 2,
 ncdata             = '$atm_init'
 history_clubb      = .true.
!clubb_history      = .true.
!clubb_rad_history  = .true.
 fincl1             = 'PS','U','V','T','Q','RELHUM','RHREFHT','QFLX','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','num_a4','LANDFRAC',
!tstep_type = 5, 
!qsplit     = 1, 
!rsplit     = 3, 
!se_nsplit  = 2,
!hypervis_subcycle = 3
EOF
cat $nl_file >> user_nl_cam

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF

#============
# Run model 
#============
cd $CASEROOT

if ( $nbundles > 0 ) then
   echo CASEROOT is $CASEROOT
else

   if ( $postCIME <= 2 ) then
      ./$CASE.submit
   else
      ./case.submit
   endif

endif

date
