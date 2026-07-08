#!/bin/csh
date
setenv PROJECT         'ESMD'
setenv CESM_PROJ       'ESMD'
setenv CESM_EMAIL      'hui.wan@pnnl.gov'

setenv MACH            'compy'
setenv COMPILER        'intel'
setenv RESOLUTION      'ne30_ne30'
setenv NTASKS_PER_INST '576'
setenv NINST           '1'
setenv NTHRDS         '4'
         
setenv CCSMTAG  'ACME_condensation_cjvogl'
setenv COMPSET  'FC5AQUAP'
setenv CCSMROOT '/compyfs/zhan391/code/ACME_condensation_cjvogl'
setenv PTMP     '/compyfs/zhan391/cnvg_condens_master_shixuan'
setenv EXEDIR   '/compyfs/zhan391/cnvg_condens_master_shixuan/exe/compile_ACME_condensation_cjvogl_FC5AQUAP_ne30_ne30_compy_intel_576proc'
setenv CSMDATA  '/compyfs/inputdata'
         
setenv CASE      'RKZ_A88_B88_C88_ql17_lmt4_ne30_ne30_FC5AQUAP_intel_compy_DT0001_01_dycore_mac'
setenv CASEROOT  '/compyfs/zhan391/cnvg_condens_master_shixuan/cases/RKZ_A88_B88_C88_ql17_lmt4_ne30_ne30_FC5AQUAP_intel_compy_DT0001_01_dycore_mac'
setenv RUNDIR    '/compyfs/zhan391/cnvg_condens_master_shixuan/run/RKZ_A88_B88_C88_ql17_lmt4_ne30_ne30_FC5AQUAP_intel_compy_DT0001_01_dycore_mac'

setenv postCIME 5

set case_setup_script = '/compyfs/zhan391/run_script/condens_master_shixuan/create_and_setup_bundled_case.csh' 

set dtime    = 1
set atm_init = /compyfs/zhan391/acme_init/ne30_FC5_init/init_gen_clim_FC5_default.cam.i.2008-01-01-00000.nc

set ncycle   = 1       # 1 cycle in total, i.e., no restart

set nlen     = 12
set stop_o   = 'nhours'

set nhtfrq   = -1
set mfilt    = 1                # 1 time step per file 

#@ nlen = 3600 / 1   # run model for 60 minutes.
#set stop_o = 'nsteps'

#@ nhtfrq  = 1800 / 1    # output every 30 minutes
#set mfilt = 1                # 1 time step per file 

set nl_file = /compyfs/zhan391/run_script/condens_master_shixuan/namelist_files/cam_nl_RKZ_A88_B88_C88_ql17_lmt4
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
