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

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "4:00:00"

#--------------------
# Namelist variables
#--------------------
cat > user_nl_cam <<EOF
 empty_htapes       = .true.
 avgflag_pertape    = 'I',
 nhtfrq             = $nhtfrq,
 mfilt              = $mfilt,
 ndens              = 1,
 rsplit               = $rsplit
 se_nsplit            = $se_nsplit
 cld_macmic_num_steps = $cld_macmic_num_steps
 iradsw             = 2,
 iradlw             = 2,
 ncdata             = '$atm_init'
 history_clubb      = .true.
!clubb_history      = .true.
!clubb_rad_history  = .true.
 fincl1             = 'WP2_CLUBB','WP3_CLUBB','WPRCP_CLUBB','WPRTP_CLUBB','WPTHLP_CLUBB','WPTHVP_CLUBB','VP2_CLUBB','VPWP_CLUBB','UP2_CLUBB','UPWP_CLUBB','THLP2_CLUBB','STEND_CLUBB','RTPTHLP_CLUBB','RTP2_CLUBB','RCMTEND_CLUBB','RCM_CLUBB','RHO_CLUBB','RIMTEND_CLUBB','RVMTEND_CLUBB','CLOUDFRAC_CLUBB','CLOUDCOVER_CLUBB','PS','U','V','T','Q','RELHUM','RHREFHT','QFLX','CLDLIQ','CLDICE','NUMLIQ','NUMICE','CDNUMC','CCN3','IWC','LWC','TS','num_a1','num_a2','num_a3','LANDFRAC','ICWMR','ICIMR','ICEFRAC','FREQI','FREQL','FREQR','FREQS','DTENDTQ','DTENDTH','DTCOND','FICE','PRECC','PRECL','PRECSH','THETAL','PBLH','QT','SL','CLDST','ZMDLF','TTENDICE','QVTENDICE', 'QITENDICE','NITENDICE','DPDLFLIQ','DPDLFICE','DPDLFT','DP_CLD','DQP','RELVAR','CONCLD','CMELIQ','CLOUD','CLDTOT','CLDMED','CLDLOW','CLDST','LHFLX','OMEGA','OMEGA500','OMEGAT'
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
