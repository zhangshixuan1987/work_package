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
#!ncdata             = '$atm_init'
cat > user_nl_cam <<EOF
 avgflag_pertape    = 'A','I'
 nhtfrq             = $nhtfrq,-6
 mfilt              = $mfilt,1
 ndens              = 1,1
 fincl2             = 'PHIS'
EOF
cat $nl_file >> user_nl_cam

#!finidat = '$lnd_init'
cat > user_nl_clm <<EOF
 hist_empty_htapes = .true.
 hist_fincl1 = 'TSA'
 hist_nhtfrq = -6
 hist_mfilt  = 1
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
