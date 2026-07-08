#!/bin/csh
date
setenv PROJECT         'ESMD'
setenv CESM_PROJ       'ESMD'
setenv CESM_EMAIL      'shixuan.zhang@pnnl.gov'

setenv MACH            'compy'
setenv COMPILER        'intel'
setenv RESOLUTION      '1.9x2.5_gx1v6'
setenv NTASKS_PER_INST '96'
setenv NINST           '1'
setenv NTHRDS         '4'
         
setenv CCSMTAG  'E3SM_CAM4_SGR'
setenv COMPSET  'F'
setenv CCSMROOT '/compyfs/zhan391/code/E3SM_CAM4_SGR'
setenv PTMP     '/compyfs/zhan391/dart_e3sm_cam4'
setenv EXEDIR   '/compyfs/zhan391/dart_e3sm_cam4/exe/compile_E3SM_CAM4_SGR_F_1.9x2.5_gx1v6_compy_intel_96proc'
setenv CSMDATA  '/compyfs/zhan391/acme_init/csmdata'
         
setenv CASE      'dart_e3sm_assm_1.9x2.5_gx1v6_DT1800'
setenv CASEROOT  '/compyfs/zhan391/dart_e3sm_cam4/cases/dart_e3sm_assm_1.9x2.5_gx1v6_DT1800'
setenv RUNDIR    '/compyfs/zhan391/dart_e3sm_cam4/run/dart_e3sm_assm_1.9x2.5_gx1v6_DT1800'

setenv postCIME 5

set case_setup_script = '/compyfs/zhan391/run_script/dart_e3sm_test/create_and_setup_bundled_dart.csh' 

set run_refcase = 'dart_e3sm_cntl_1.9x2.5_gx1v6_DT1800'
set run_refdate = '2009-01-01'
set run_reftod  = '21600'

set dtime    = 1800
set atm_init = /compyfs/zhan391/dart_e3sm_cam4/run/dart_e3sm_assm_1.9x2.5_gx1v6_DT1800/dart_e3sm_cntl_1.9x2.5_gx1v6_DT1800.cam.i.2009-01-01-21600.nc
set lnd_init = /compyfs/zhan391/dart_e3sm_cam4/run/dart_e3sm_assm_1.9x2.5_gx1v6_DT1800/dart_e3sm_cntl_1.9x2.5_gx1v6_DT1800.clm2.r.2009-01-01-21600.nc

set ncycle   = 1          # 1 cycle in total, i.e., no restart

#@ nlen = 3600 / 1800   # run model for 60 minutes.
#set stop_o   = 'nsteps'
set nlen      = 6 # run model for 5 years
set stop_o    = 'nhours'

#@ nhtfrq     = 1800 / 1800    # output every 30 minutes
set nhtfrq    = -6               # monthly average
set mfilt     = 1                # 1 time step per file 

set nl_file = /compyfs/zhan391/run_script/dart_e3sm_test/namelist_files/cam_nl_dart_e3sm_assm
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
./xmlchange  -file env_run.xml -id  REST_N       -val 6
./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'nhours'
./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

@ ncpl = 86400 / $dtime

./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "1:00:00"
#./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val "slurm"
#./xmlchange  -file env_batch.xml -id PROJECT            -val "m3089"
#./xmlchange  -file env_batch.xml -id CHARGE_ACCOUNT     -val "m3089"

./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ATM  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_CPL  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_OCN  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_WAV  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_GLC  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ICE  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ROF  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_LND  -val netcdf
./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ESP  -val netcdf

#--------------------
# Namelist variables
#--------------------
# fincl1             = 'PS','U','V','T','Q','CLDLIQ','CLDICE','NUMLIQ','NUMICE','num_a1','num_a2','num_a3','LANDFRAC'
cat > user_nl_cam <<EOF
 ncdata             = '$atm_init'
 avgflag_pertape    = 'A','I'
 nhtfrq             = $nhtfrq,$nhtfrq
 mfilt              = $mfilt,$mfilt
 ndens              = 1,1
 fincl2             = 'PHIS'
EOF
cat $nl_file >> user_nl_cam

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
 hist_empty_htapes = .true.
 hist_fincl1 = 'TSA'
 hist_nhtfrq = $nhtfrq
 hist_mfilt  = $mfilt
 hist_avgflag_pertape = 'I'
EOF

#============
# Run model 
#============
cd $CASEROOT

if ( $postCIME <= 2 ) then
   ./$CASE.submit
else
   ./case.submit
endif

date
