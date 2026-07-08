#!/bin/csh
date

##############################################################################################
# This is the driver scripts for the SCM test simulations which
# - compiles the model (just once);
# - generates one script for each simulation; 
# - runs those scripts to create cases and get ready for production run;
# Author: shixuan.zhang@pnnl.gov
##############################################################################################
set compile_model = 0
set run_model     = 1

setenv debug  'TRUE'
setenv debug  'FALSE'

#####determine which queue for the job###########
setenv queue  'debug'
setenv queue  'regular'

#setenv submit_to_queue  'FALSE'
setenv MACH              compy
setenv submit_to_queue  'TRUE'
setenv wall_time        '03:00:00'

# $taskname is used to organize the exe/run/case directories. Not used for other purposes.
setenv taskname  "F_SCAM5_201806"
# $runcase is the user specified case and configurations for SCM simulation 
# A configuration file ${runcase}_default_config_info.csh is created by 
# the user to save the test case information
setenv runcase   "DYCOMSrf01_default"

#job name 
set jobname = ${runcase}_no_gw_convect

# Number of vertical levels 
setenv NVLEV    72

# model time step size 
#set dtime = 1800 # default for SE dycore 
set dtime = 1200  # default for FV dycore

# simulation length 
setenv stop_option   nhours
setenv stop_n        96 # the forcing data tsec = 0, 345600 ;

#---------------------------------
# user specified configuration 
#---------------------------------

# User enter any needed modules to load or use below
#  EXAMPLE:
#module load python

# What version of E3SM? (v1 or v2)
#  Select v1 if you are running with E3SMv1 RELEASE code
#  Select v2 if you are running with up-to-date E3SM master or v2 release code
#  (If you are running with non-up-to-date master code then you may need to modify
#    aspects of this script to get it to compile.)
# NOTE: v2 tunings provided are subject to change as model release is finalized
#  (and if so will be updated here)
# This script does NOT clone the model code from the E3SM repo.
# # Instead, it assumes that the code is located at USER's code directory 
setenv e3sm_version v1
setenv CCSMTAG maint-1.0_radheat_cpl
setenv CCSMROOT /compyfs/zhan391/code/${CCSMTAG}

# Set the dynamical core
#   1) Select "Eulerian" ONLY IF you are running E3SMv1 release code 
#   2) Select "SE" IF you are running code from recent E3SM master or v2
setenv dycore Eulerian
#  WARNING:  EULERIAN DYCORE SCM IS NO LONGER SUPPORTED. You are only safe
#  to use Eulerian dycore SCM if you are using E3SMv1 release code.  Else,
#  user be(very)ware

# COSP, set to false unless user really wants it
setenv do_cosp  false

#====================================================================
# Paths to  model input/output
#====================================================================
if ($MACH == "cori-knl") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS 1

   setenv CSMDATA  /project/projectdirs/acme/inputdata
   setenv WORKDIR /compyfs/zhan391
   setenv PTMP   $WORKDIR/$wkdrnam
   setenv EXELOC $PTMP/exe
   set initDir = /global/cscratch1/sd/zhan391/acme_input/FC5AV1C-L_init_201712

else if ($MACH == "cori-haswell") then

   setenv CESM_PROJ  m3089
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS 1

   setenv CSMDATA  /project/projectdirs/acme/inputdata
   setenv WORKDIR /compyfs/zhan391
   setenv PTMP   $WORKDIR/$wkdrnam
   setenv EXELOC $PTMP/exe
   set initDir = /global/cscratch1/sd/zhan391/acme_input/FC5AV1C-L_init_201712

else if ($MACH == "constance") then

   setenv CESM_PROJ uq_climate
   setenv PROJECT $CESM_PROJ
   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS  1
   setenv NCORES_PER_NODE  1
   set scheduler = "SLURM"

   setenv CSMDATA  /pic/projects/climate/csmdata
   setenv WORKDIR /compyfs/zhan391
   setenv PTMP     $WORKDIR/$taskname
   setenv EXELOC   $PTMP/exe
   set initDir = /pic/projects/uq_climate/wanh895/acme_input/FC5AV1C-L_init_201712

else if ($MACH == "compy") then

   setenv CESM_EMAIL shixuan.zhang@pnnl.gov
   setenv PROJECT   ESMD
   setenv CESM_PROJ $PROJECT

   setenv COMPILER   intel
   setenv NTHRDS 1
   setenv NTASKS  1
   setenv NCORES_PER_NODE  1
   set scheduler = "PBS"

   setenv CSMDATA  /compyfs/inputdata
   setenv WORKDIR  /compyfs/zhan391
   setenv PTMP     $WORKDIR/$taskname
   setenv EXELOC   $PTMP/exe
   set initDir = /compyfs/zhan391/acme_init/ne30_FC5_init

else

   echo Specify paths for MACH $MACH ! Abort.
   exit

endif

if ($e3sm_version == v1) then
  setenv COMPSET  F_SCAM5
  setenv atm_mod  cam
  setenv physset  cam5
else
  setenv COMPSET  F_SCAM
  setenv atm_mod  eam
  setenv physset  default
endif

#-----------------------------
if ($dycore == Eulerian) then
  setenv RESOLUTION  T42_T42
endif

if ($dycore == SE) then
  setenv RESOLUTION ne4_ne4
endif

#-----------------
set git_dir = `pwd`

mkdir -p $PTMP 

set casename = ${jobname}_${RESOLUTION}_L${NVLEV}_DT`printf "%04d" ${dtime}`

set execase = compile_${CCSMTAG}_${COMPSET}_${RESOLUTION}_L${NVLEV}_${MACH}_${COMPILER}

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

#------------------
#load case info and basic configuration 
set case_config_script = "${runcase}_config_info.csh"
source ${case_config_script}

####################################################################
# Compile model (just once)
####################################################################
if ($compile_model > 0) then

   setenv CASE     $execase
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case
   source ${case_setup_script}

   # Build the model
   cd $CASEROOT

   #./xmlchange -file env_build.xml -id GMAKE_J -val '8'
   ./xmlchange -file env_build.xml -id DEBUG   -val $debug

   echo
   echo Start to build model
   echo

   ./case.build

    echo
    echo Finished building the model.
    echo

endif

#####################################################################
# Conduct simulation(s)
#####################################################################
if ($run_model > 0) then

   #---------------------------------
   cd $git_dir

   setenv CASE     $casename
   setenv CASEROOT ${PTMP}/cases/$CASE
   setenv RUNDIR   ${PTMP}/run/$CASE

   # Create and set up new case; no need to build model

   source ${case_setup_script}

   cd $CASEROOT

   ./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'
   #./xmlchange -file env_build.xml -id CAM_CONFIG_OPTS -append --val='-cosp'
   ./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val $wall_time

  # Set sim length so that a cycle can be finished in 8 hours
   ./xmlchange  -file env_run.xml -id  STOP_N          -val $stop_n      
   ./xmlchange  -file env_run.xml -id  STOP_OPTION     -val $stop_option
   ./xmlchange  -file env_run.xml -id  SAVE_TIMING_DIR -val $PTMP
   ./xmlchange  -file env_run.xml -id  DOUT_S          -val 'FALSE'
   ./xmlchange  -file env_run.xml -id  DOUT_L_MS       -val 'FALSE'

  # Set coupling frequency properly based on the time step size 
   @ ncpl = 86400 / $dtime
   ./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
   ./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime
   ./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime

   # User specified namelist variables
cat <<EOF >> user_nl_cam
   history_verbose      = .true.
   history_aero_optics  = .true.
   history_aerosol      = .true.
   history_clubb        = .true.
   history_budget       = .true.
   clubb_history        = .true.
   clubb_rad_history    = .true. 
   use_gw_convect       = .false.
   gw_drag_file         = ''
   nhtfrq               =  1, 
   mfilt                = 100000,
   avgflag_pertape      = 'I', 
   fincl1               = 'FUL', 'FDL'
EOF

   # Run model 
   cd $CASEROOT
   ./case.submit

endif
