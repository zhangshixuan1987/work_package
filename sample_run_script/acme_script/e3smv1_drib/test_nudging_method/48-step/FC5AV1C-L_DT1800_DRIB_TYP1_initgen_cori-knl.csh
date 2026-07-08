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
setenv CCSMTAG  E3SMv1_201806
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
setenv COMPSET FC5AV1C-L
setenv RESOLUTION ne30_ne30 #ne30_oECv3
setenv MACH compy
setenv COMPILER intel
setenv scheduler 'PBS'

setenv NTASKS_PER_INST  384
setenv NINST  1
setenv NTHRDS 1
setenv NCORES_PER_NODE  24

setenv WORKDIR /compyfs/zhan391

setenv CSMDATA /compyfs/inputdata/

setenv INPUT_NUDGING /compyfs/zhan391/FC5AV1C-L_201806_initgen_84nodes/run/NDGDATA_FC5AV1C-L_org_DT1800


set taskname = E3SMv1_DRIBBLING_48STEP_TEST
set casename = DRIB_TYP1_${COMPSET}_NDG_CLIM_DT`printf "%04d" ${dtime}`

#-----------------
#setup branch run information
set dtime1      = 1800 
set run_refcase = FC5AV1C-04P2_ne30_ne30_intel_cori-knl
set run_refdate = "0001-08-01"
set run_reftod  = "00000"
set initDir     = /compyfs/zhan391/acme_init/FC5AV1C-04P2_init_201907/

setenv PTMP     $WORKDIR/$taskname
setenv EXELOC   $PTMP/exe

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
   
   ./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'
   #./xmlchange -file env_build.xml -id GMAKE_J -val '16'
   #./xmlchange -file env_build.xml -id DEBUG   -val $debug
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

###copy the required restart files to the directory for branch run#####
   if ( ! -d $RUNDIR ) then
      mkdir -p $RUNDIR
   endif
   set camiin  = "$run_refcase.cam.i.$run_refdate-$run_reftod.nc"
   set camrin  = "$run_refcase.cam.r.$run_refdate-$run_reftod.nc"
   set camrsin = "$run_refcase.cam.rs.$run_refdate-$run_reftod.nc"
  #set camrhin = "$run_refcase.cam.rh0.$run_refdate-$run_reftod.nc"
   set lndiin  = "$run_refcase.clm2.r.$run_refdate-$run_reftod.nc"
  #set lndrin  = "$run_refcase.clm2.rh0.$run_refdate-$run_reftod.nc"
   set cicerin = "$run_refcase.cice.r.$run_refdate-$run_reftod.nc"
   set cplrin  = "$run_refcase.cpl.r.$run_refdate-$run_reftod.nc"
   set docnrin  = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.nc"
   set docnrsin = "$run_refcase.docn.rs1.$run_refdate-$run_reftod.bin"
   ln -sf ${initDir}${camiin}  $RUNDIR/${camiin}
   ln -sf ${initDir}${camrin}  $RUNDIR/${camrin}
   ln -sf ${initDir}${camrsin} $RUNDIR/${camrsin}
  #ln -sf ${initDir}${camrhin} $RUNDIR/${camrhin}
   ln -sf ${initDir}${lndiin}  $RUNDIR/${lndiin}
  #ln -sf ${initDir}${lndrin}  $RUNDIR/${lndrin}
   ln -sf ${initDir}${cicerin} $RUNDIR/${cicerin}
   ln -sf ${initDir}${cplrin}  $RUNDIR/${cplrin}
   ln -sf ${initDir}${docnrsin} $RUNDIR/${docnrsin}

   echo "${camrin}"     >! $RUNDIR/rpointer.atm
   echo "${lndiin}"     >! $RUNDIR/rpointer.lnd
   echo "${cicerin}"    >! $RUNDIR/rpointer.ice
   echo "${cplrin}"     >! $RUNDIR/rpointer.drv
   echo "${docnrin}"    >! $RUNDIR/rpointer.ocn
   echo "${docnrsin}"   >> $RUNDIR/rpointer.ocn

   set atm_init   = "$RUNDIR/${camiin}"
   set lnd_init   = "$RUNDIR/${lndiin}"
   set ice_init   = "${cicerin}"

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'
   ./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'

   # Runtime options: edit env_run.xml
   #
   # Note: with 2048 cores (128 nodes) the L72 model integrates at 
   # the speed of roughly 1 hour wall clock time per model month.
   # Set restart cycle to 6 months to make full use of the max.
   # run time for the job size bin.

 # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val 'nmonths'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '48'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'nsteps'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_QUEUE          -val 'short'
   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "2:00:00"
   # Properly change time step for both atm and lnd.
   
   @ ncpl = 86400 / $dtime

   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime

   ./xmlchange  -file env_run.xml -id  PIO_TYPENAME      -val netcdf

   # Namelist variables


cat <<EOF >> user_nl_cam
 ncdata               = '$atm_init'
 history_amwg         = .true.
 history_verbose      = .false.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .true.
 history_budget       = .true.
!history_microphysics = .true.
 cosp_lite            = .true.
 docosp               = .true.
 iradsw               = 2,
 iradlw               = 2,
 inithist             = 'MONTHLY'
 inithist_all         = .true.
 empty_htapes         = .false.
 clubb_history        = .true.
 clubb_rad_history    = .true. 
! do_extra_macmic_diag = .true.
 dribble_tend_into_macmic_loop = 1
 dribble_start_step            = 24
 nhtfrq = 1
 mfilt  = 1

EOF

cat <<EOF >> user_nl_cam
 !.......................................................
 ! nudging
 !.......................................................
  Nudge_Model  = .True.
  Nudge_Method = 'Linear'
  Nudge_Tau    = 6.0  !nudging tau in h 
  Nudge_Loc    = .True.
  Nudge_Curr   = .True.
  Num_Slice    = 1
  Nudge_Path   = '${INPUT_NUDGING}/'
  Nudge_File_Template = 'NDGDATA_FC5AV1C-L_org_DT1800.cam.h1.%y-%m-%d-%s.nc'
  Nudge_Times_Per_Day =  8 !! nudging input data frequency
  Model_Times_Per_Day = 48 !! should not be larger than 48 if dtime = 1800s
  Nudge_Uprof = 1
  Nudge_Ucoef = 1.
  Nudge_Vprof = 1
  Nudge_Vcoef = 1.
  Nudge_Tprof = 0
  Nudge_Tcoef = 0.
  Nudge_Qprof = 0
  Nudge_Qcoef = 0.
  Nudge_PSprof = 0
  Nudge_PScoef = 0.
  Nudge_Beg_Year = 0
  Nudge_Beg_Month = 1
  Nudge_Beg_Day = 1
  Nudge_End_Year = 9999
  Nudge_End_Month = 12
  Nudge_End_Day = 31
EOF

# ==============
# CLM Namelist 
# With a RUN_TYPE=hybrid the finidat is automatically specified
# ==============
cat > user_nl_clm <<EOF
 finidat = '$lnd_init'
EOF

# ==============
#  CICE Namelist 
# ==============
cat > user_nl_cice <<EOF
ice_ic = '$ice_init'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


