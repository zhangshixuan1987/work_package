#====================================================================
# Create and set up new case. No need to build the model. 
#====================================================================
source ${case_setup_script}

cd $CASEROOT

./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

#-----------------------------------------
# Runtime options: edit env_run.xml
#-----------------------------------------
cd $CASEROOT

@ nresub = $ncycle - 1

./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o 
./xmlchange  -file env_run.xml -id  REST_N       -val 1 
./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'nhours' 
./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

@ ncpl = 86400 / $dtime

./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 

./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "12:00:00"

#--------------------
# Namelist variables
#--------------------
cat > user_nl_cam <<EOF
 avgflag_pertape    = 'I','A'
 nhtfrq             = $nhtfrq,$nhtfrq
 mfilt              = $mfilt, $mfilt
 ndens              = 1, 1
 empty_htapes       = .true.,
 fincl1             = 'PS', 'U', 'V', 'T', 'Q', 'CLDLIQ', 'LANDFRAC', 'DTCORE', 'PTTEND', 'RKZ_ql','RKZ_qv', 'RKZ_qsat', 'RKZ_RH', 'RKZ_f', 'RKZ_dfdRH', 'RKZ_term_A', 'RKZ_term_B', 'RKZ_term_C', 'RKZ_dfdt', 'RKZ_ql_incld','RKZ_Al', 'RKZ_Av', 'RKZ_AT', 'RKZ_zqme', 'RKZ_qme', 'RKZ_qme_lm4_qv', 'RKZ_qme_lm4_ql', 'RKZ_qme_lm5_ng', 'RKZ_qme_lm5_ps', 'RKZ_Taf', 'RKZ_Tbf'
 fincl2             = 'PS', 'U', 'V', 'T', 'Q', 'CLDLIQ', 'LANDFRAC', 'DTCORE', 'PTTEND', 'RKZ_ql','RKZ_qv', 'RKZ_qsat', 'RKZ_RH', 'RKZ_f', 'RKZ_dfdRH', 'RKZ_term_A', 'RKZ_term_B', 'RKZ_term_C', 'RKZ_dfdt', 'RKZ_ql_incld','RKZ_Al', 'RKZ_Av', 'RKZ_AT', 'RKZ_zqme', 'RKZ_qme', 'RKZ_qme_lm4_qv', 'RKZ_qme_lm4_ql', 'RKZ_qme_lm5_ng', 'RKZ_qme_lm5_ps', 'RKZ_Taf', 'RKZ_Tbf'
 ncdata             = '$atm_init'
 deep_scheme        = 'off',
 shallow_scheme     = 'off',
 l_tracer_aero      = .false.
 l_vdiff            = .false.
 l_rayleigh         = .false.
 l_gw_drag          = .false.
 l_ac_energy_chk    = .true.
 l_bc_energy_fix    = .true.
 l_dry_adj          = .false.
 l_st_mac           = .true.
 l_st_mic           = .false.
 l_rad              = .false.
EOF
cat $nl_file >> user_nl_cam

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
