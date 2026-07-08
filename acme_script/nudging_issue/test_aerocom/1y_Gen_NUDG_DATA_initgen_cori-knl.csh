#!/bin/csh 
date

set my_branch_shortname = "master"
#=====================================
# HOW TO USE THIS SCRIPT:
#
# 1. If the model source code has already been checked out 
#    to ~/codes/$my_branch_shortname, set fetch_code to 0 (off);
#    otherwise, set it to 1 (on) to obtain code.
#
set fetch_code    = 0   # 0 = No, >0 = Yes

# 2. If the executable has already been created under
#    $PTMP/exe/, set compile_model to 0 (off);
#    otherwise set it to 1 to build the model.
#
set compile_model = 0  # 0 = No, >0 = Yes

# 3. 
#
set run_model  = 1   # 0 = No, >0 = Yes

#-----------------------------
set debug = 'TRUE'
set debug = 'FALSE'

set git_dir = `pwd`

set dtime = 1800

####################################################################
# Fetch code
####################################################################
setenv CCSMTAG E3SM_MAINT1.0_NDG
setenv CCSMROOT /compyfs/zhan391/code/$CCSMTAG

if ($fetch_code > 0) then
   mkdir -p $CCSMROOT:h
   cd $CCSMROOT:h
   git clone -b huiwanpnnl/atm/$my_branch_shortname git@github.com:ACME-Climate/ACME.git $CCSMTAG

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   ESMD
setenv CESM_PROJ $PROJECT

setenv postCIME 1
setenv COMPSET F20TRC5-CMIP6
setenv RESOLUTION ne30_ne30 #ne30_oECv3
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  1024 #2016
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  24

setenv WORKDIR /compyfs/zhan391

setenv CSMDATA /compyfs/inputdata/

set taskname = F20TRC5-CMIP6_NUDG_ISSUE

setenv initDir /compyfs/zhan391/acme_init

set casename = AMIP_NDG_DATA_${COMPSET}_DT`printf "%04d" ${dtime}`

#-----------------
#setup run information
set dtime1      = 1800 
set run_refdate = "2009-10-01"
set run_reftod  = "00000"
set atm_file    = "${initDir}/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800.cam.i.2010-10-01-00000.nc"
set lnd_file    = "${initDir}/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800.clm2.r.2010-10-01-00000.nc"

setenv PTMP     $WORKDIR/$taskname
setenv EXELOC   $PTMP/exe/

#-----------------
mkdir -p $PTMP

if ($NINST > 1) then
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NINST}x${NTASKS_PER_INST}x${NTHRDS}bundle
else
   set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_${MACH}_${COMPILER}_${NTASKS_PER_INST}x${NTHRDS}threads
endif

if ($debug == 'TRUE') then
   set execase = ${execase}_debug
endif

setenv EXEDIR ${EXELOC}/$execase

set script_output_dir = $PTMP/scripts
if ( -e $script_output_dir )  then
   mv $script_output_dir ${script_output_dir}_bak_`date +%F-%H%M%S-%N`
endif
mkdir -p $script_output_dir

set case_setup_script = "create_and_setup_bundled_gen.csh"
cp $git_dir/$case_setup_script $script_output_dir/$case_setup_script

####################################################################
# Compile model
####################################################################
if ($compile_model > 0) then

   setenv CASE     $execase
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case

   source ${case_setup_script}

   # Build the model

   cd $CASEROOT
   ./xmlchange -file env_build.xml -id GMAKE_J -val '16'
   ./xmlchange -file env_build.xml -id DEBUG   -val $debug
   ./case.build

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then

   cd $git_dir

   setenv CASE     $casename
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   set atm_init   = ${atm_file}
   set lnd_init   = ${lnd_file}

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'
   #./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'

   # Runtime options: edit env_run.xml
   #
   # Note: with 2048 cores (128 nodes) the L72 model integrates at 
   # the speed of roughly 1 hour wall clock time per model month.
   # Set restart cycle to 6 months to make full use of the max.
   # run time for the job size bin.
   #./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "$CSMDATA/ocn/docn7/SSTDATA/sst_ice_CMIP6_DECK_E3SM_1x1_c20180213.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_DATA_FILENAME -val "/compyfs/zhan391/acme_init/csmdata/ocn/docn7/SSTDATA/sst_weekly_cdcunits_1x1_1990-01-01-2013-08-11_n3.nc"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_ALIGN -val "1990"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_START -val "1990"
   ./xmlchange  -file env_run.xml -id SSTICE_YEAR_END   -val "2013"

  # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N           -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION      -val 'nmonths'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '1'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'nmonths'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "3:00:00"
   # Properly change time step for both atm and lnd.
   
   @ ncpl = 86400 / $dtime

   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime

   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ATM  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_CPL  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_OCN  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_WAV  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_GLC  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ICE  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ROF  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_LND  -val netcdf
   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME_ESP  -val netcdf

   # Namelist variables


cat <<EOF >> user_nl_cam
 ncdata               = '$atm_init'
 history_amwg         = .true.
 docosp               = .false.
 cosp_amwg            = .false.
 history_verbose      = .false.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .false.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = 2,
 iradlw               = 2,
 inithist             = 'MONTHLY'
 inithist_all         = .true.
 empty_htapes         = .false.
 nhtfrq               = 0,
 mfilt                = 1,
 ndens                = 2,
EOF

cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif

