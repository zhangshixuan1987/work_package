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
setenv CCSMTAG ACME_master_20180611
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

setenv INPUT_NUDGING /compyfs/zhan391/acme_init/

set taskname = FC5AV1C-L_201806_initgen_84nodes
set casename = CPDRIB_202010_${COMPSET}_rev_DT`printf "%04d" ${dtime}`

#-----------------
#setup branch run information
set dtime1      = 1800 
set run_refcase = FC5AV1C-04P2_ne30_ne30_intel_cori-knl
set run_refdate = "0001-05-01"
set run_reftod  = "00000"
set initDir     = ${INPUT_NUDGING}/FC5AV1C-04P2_init_201907/

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
   cp -r ${initDir}${camiin}  $RUNDIR/${camiin}
   cp -r ${initDir}${camrin}  $RUNDIR/${camrin}
   cp -r ${initDir}${camrsin} $RUNDIR/${camrsin}
  #cp -r ${initDir}${camrhin} $RUNDIR/${camrhin}
   cp -r ${initDir}${lndiin}  $RUNDIR/${lndiin}
  #cp -r ${initDir}${lndrin}  $RUNDIR/${lndrin}
   cp -r ${initDir}${cicerin} $RUNDIR/${cicerin}
   cp -r ${initDir}${cplrin}  $RUNDIR/${cplrin}
   cp -r ${initDir}${docnrsin} $RUNDIR/${docnrsin}

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
   #./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'

   # Runtime options: edit env_run.xml
   #
   # Note: with 2048 cores (128 nodes) the L72 model integrates at 
   # the speed of roughly 1 hour wall clock time per model month.
   # Set restart cycle to 6 months to make full use of the max.
   # run time for the job size bin.

 # ./xmlchange  -file env_run.xml -id  RUN_STARTDATE  -val $start_date
   ./xmlchange  -file env_run.xml -id  REST_N         -val '1'
   ./xmlchange  -file env_run.xml -id  REST_OPTION    -val 'ndays'

  # Set sim length so that a cycle can be finished in 8 hours
  #
   @ nmon = $dtime * 16 / 2000
   ./xmlchange  -file env_run.xml -id  STOP_N       -val '1'      #'240'
   ./xmlchange  -file env_run.xml -id  STOP_OPTION  -val 'ndays'  #'nsteps'

   @ nresub = 68 / $nmon
   ./xmlchange  -file env_run.xml -id  RESUBMIT     -val '0'

   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "2:00:00"
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
 l_dribling_tend      = .true.
 l_dribling_uv        = .false.
 l_dribling_w         = .false.
 history_amwg         = .true.
 docosp               = .false.
 cosp_amwg            = .false.
 history_verbose      = .true.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .false.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = -1,
 iradlw               = -1,
 inithist             = 'DAILY'
 inithist_all         = .true.
 empty_htapes         = .false.
 pergro_test_active   = .false.
 l_aerosol_cldgrow    = .true.
 l_aerosol_cldshnk    = .true.
 l_aerosol_oldcld     = .true.
 l_aerosol_mixing     = .true.
 clubb_history        = .true.
 clubb_rad_history    = .true. 
 macmic_extra_diag    = .true.
 macmic_clubb_diag    = .false.
 macmic_mg2_diag      = .false.
 nhtfrq               =    1,   
 mfilt                =    1,  
 avgflag_pertape      =  'I', 
EOF

cat <<EOF >> user_nl_cam
 !.......................................................
 ! nudging
 !.......................................................
  Nudge_Model  = .False.
  Nudge_Method = 'Linear'
  Nudge_Path   = '${INPUT_NUDGING}/NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800/'
  Nudge_File_Template = 'NDRB_NUDG_DATA_F20TRC5-CMIP6_DT1800.cam.h1.%y-%m-%d-%s.nc'
  Nudge_Times_Per_Day = 4  !! nudging input data frequency
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
  Nudge_Beg_Year = 2009
  Nudge_Beg_Month = 12
  Nudge_Beg_Day = 1
  Nudge_End_Year = 2010
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
cat > $user_cice <<EOF
ice_ic = '$ice_init'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

date
endif


