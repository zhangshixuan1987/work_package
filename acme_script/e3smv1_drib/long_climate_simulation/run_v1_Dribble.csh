#!/bin/csh 
date

#=====================================
# HOW TO USE THIS SCRIPT:
#
# 1. Assume that the model source code has already been checked out 
#    set compile_model to 0 (off);
#    otherwise set it to 1 to build the model.
set compile_model = 1   # 0 = No, >0 = Yes

# 2. run the model 
#
set run_model  = 1   # 0 = No, >0 = Yes

#-----------------------------
set debug = 'TRUE'
set debug = 'FALSE'

set git_dir = `pwd`

set dtime = 1800

####################################################################
set branch_hash = 684ade
set branch_date = 20201008

set taskname = FC5AV1C-L_201806_10YR_FRE
set expname  = FC5AV1C-L_NEW_TQDRB_DT1800

##work directory for code and model output;;; 
set HOME    = /compyfs/zhan391

setenv CSMDATA /compyfs/inputdata/

####################################################################
# code
####################################################################
setenv CCSMTAG "E3SMv1_commit_${branch_hash}"
setenv CCSMROOT $HOME/code/$CCSMTAG 

####################################################################
# Machine, compset, PE layout etc.
####################################################################

setenv CESM_EMAIL shixuan.zhang@pnnl.gov
setenv PROJECT   ESMD
setenv CESM_PROJ $PROJECT

setenv postCIME 1
setenv COMPSET FC5AV1C-L 
setenv RESOLUTION ne30_ne30
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST 1200
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  24

setenv WORKDIR $HOME

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

set case_setup_script = "create_and_setup_bundled_case.csh"
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

   setenv CASE     $expname
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

   # Runtime options: edit env_run.xml
   #
   # Note: with 2048 cores (128 nodes) the L72 model integrates at 
   # the speed of roughly 1 hour wall clock time per model month.
   # Set restart cycle to 6 months to make full use of the max.
   # run time for the job size bin.

   ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val '0004-10-01'
   ./xmlchange  -file env_run.xml -id  REST_N       -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'nmonths'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '124'   #'288'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'nmonths'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val  '0' #'5'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val  "48:00:00" #"144:00:00"
   # Properly change time step for both atm and lnd.
   
   @ ncpl = 86400 / $dtime

   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime

   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME      -val netcdf

   # Namelist variables

cat > user_nl_cam <<EOF
 l_dribble_tend_into_macmic_loop = .true.
 dribble_start_step   = 1 
 history_amwg         = .true.
 history_verbose      = .true.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = -1,
 iradlw               = -1,
 inithist             = 'MONTHLY'
 inithist_all         = .true.
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


