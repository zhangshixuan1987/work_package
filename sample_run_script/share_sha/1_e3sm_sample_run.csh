#!/bin/csh
date

set echo verbose

set fetch_code    = 0   # 0 = No, >0 = Yes

set setup_model   = 1   # 0 = No, >0 = Yes
set compile_model = 1   # 0 = No, >0 = Yes
set run_model     = 1   # 0 = No, >0 = Yes

#refer the following directory to simulation setup 
#https://acme-climate.atlassian.net/wiki/spaces/NGDNA/pages/3565846567/CONUS+ne32x32+simulation+script

####################################################################
# Fetch code
####################################################################
setenv CODE_ROOT /pscratch/sd/z/zhan391/DARPA_project/e3sm_model/code
setenv CCSM_TAG  E3SM
setenv CCSM_ROOT ${CODE_ROOT}/${CCSM_TAG}

setenv BRANCH $CCSM_TAG 

if ($fetch_code > 0) then

   cd ${CODE_ROOT}
   git git@github.com:E3SM-Project/E3SM.git $CCSM_TAG

   # Setup git hooks
   rm -rf .git/hooks
   git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
   git config commit.template .git/hooks/commit.template
   git checkout ${BRANCH}
   # Bring in all submodule components
   git submodule update --init --recursive

endif

####################################################################
# Machine, compset, PE layout etc.
####################################################################
setenv COMPSET    "F20TR"
setenv RESOLUTION "ne30pg2_EC30to60E2r2" 
setenv MACH       pm-cpu 
setenv COMPILER   gnu
setenv PROJECT    "m1867"

setenv PTMP             /pscratch/sd/z/zhan391/darpa_scratch
setenv CASE_NAME        ${COMPSET}_${RESOLUTION}
setenv CASE_ROOT        $PTMP/$CASE_NAME/cases
setenv CASE_RUN_DIR     $PTMP/$CASE_NAME/run
setenv CASE_BUILD_DIR   $PTMP/$CASE_NAME/build

####################################################################
# Run options setup 
####################################################################
set MODEL_START_TYPE = "hybrid" #change to "startup" if use random ic from E3SM data directory 
set MODEL_START_YEAR = "2007"
set MODEL_START_DATE = "${MODEL_START_YEAR}-01-01"

#Additional options for 'branch' and 'hybrid'
set GET_REFCASE      = TRUE
# reference case name
set RUN_REFCASE      = "v2.LR.historical_0101"
# same as MODEL_START_DATE for 'branch', can be different for 'hybrid'
set RUN_REFDATE      = "2007-01-01"
# reference case run directory  
set RUN_REFDIR       = "/pscratch/sd/z/zhan391/DARPA_project/SCREAM_INIT/E3SMv2/2007-01-01"

set lnd_data         = "/pscratch/sd/z/zhan391/DARPA_project/SCREAM_INIT/LAND/betacast/SCREAM-IELM-ne30pg2_20100101_0120.elm.r.2010-01-01-00000.nc"
set lnd_surdat       = "/global/cfs/cdirs/e3sm/inputdata/lnd/clm2/surfdata_map/surfdata_ne30pg2_simyr2000_c210402.nc"

set DEBUG            = 'FALSE'
set JOB_QUEUE        = "debug" #"regular"
set WALLTIME         = "00:30:00"
set STOP_OPTION      = "nmonths"
set STOP_N           = "1"
set RESUBMIT         = "0"
set REST_OPTION      = "nmonths"
set REST_N           = "12" 

set NTASKS           = 120
set NTHRDS           = 2
set NNODES_ATM       = 1
set NNODES_OCN       = 1
set NTASKS_PER_NODE  = 64

#specific setups for model 
set RUN_START_DATE   = ${MODEL_START_DATE}
set RUN_START_YEAR   = ${MODEL_START_YEAR}
set RUN_START_TOD    = 0

set DTIME             = 1800 #unit:s
set NCPL_ATM          = `expr 24 \* 60 \* 60 \/ $DTIME`

set NUTOP             = 2.5e5
set SEDTIME           = `expr $DTIME \/ 6`
set EPS_AGRID         = 1e-9
set DT_TRACER_FACTOR  = 6
set HYPERVIS_SUBCLE_Q = 6
set IRADSW            = 2 #radiation coupling every 2 steps (i.e. 1hr)
set IRADLW            = 2 #radiation coupling every 2 steps (i.e. 1hr)

#derive total tasks
set NTASKS_ATM       = `expr ${NNODES_ATM} \* ${NTASKS_PER_NODE}`
set NTASKS_OCN       = `expr ${NNODES_OCN} \* ${NTASKS_PER_NODE}`
set NTASKS_PD        = `expr ${NTASKS_PER_NODE} \* ${NTHRDS}`
#echo $NTASKS_ATM $NNODES_OCN $NTASKS_PD
if ( ${NTASKS_OCN} != ${NTASKS_ATM} ) then
    set NNODES=`expr ${NNODES_ATM} + ${NNODES_OCN}`
else
    set NNODES=${NNODES_ATM}
endif 

set PELAYOUT = "${NNODES}x${NTASKS_PER_NODE}x${NTHRDS}"
set EMAIL    = "shixuan.zhang@pnnl.gov"


####################################################################
# Create model
####################################################################
if ( -d $CASE_ROOT ) then 
  rm -rvf $CASE_ROOT
endif

${CCSM_ROOT}/cime/scripts/create_newcase \
                --case $CASE_NAME \
                --output-root ${CASE_ROOT} \
                --script-root ${CASE_ROOT} \
                --compset ${COMPSET} \
                --res ${RESOLUTION} \
                --machine ${MACH} \
                --compiler ${COMPILER} \
                --queue ${JOB_QUEUE} \
                --project ${PROJECT} \
                --walltime ${WALLTIME}

# Copy this script to case directory
#=============================================
cp -v `basename $0` ${CASE_ROOT}/

####################################################################
# Setup model 
####################################################################
if ($setup_model > 0 ) then

   cd $CASE_ROOT

   ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'

   ./xmlchange   JOB_WALLCLOCK_TIME=$WALLTIME
   ./xmlchange   JOB_QUEUE=$JOB_QUEUE

   ./xmlchange   RUNDIR=$CASE_RUN_DIR
   ./xmlchange   EXEROOT=$CASE_BUILD_DIR
   ./xmlchange   RUN_STARTDATE=$RUN_START_DATE
   ./xmlchange   START_TOD=$RUN_START_TOD
   ./xmlchange   STOP_N=$STOP_N
   ./xmlchange   STOP_OPTION=$STOP_OPTION
   ./xmlchange   REST_N=$REST_N
   ./xmlchange   REST_OPTION=$REST_OPTION
   ./xmlchange   DOUT_S='FALSE'
   ./xmlchange   BUDGETS='TRUE'

    if (${NTASKS_OCN} != ${NTASKS_ATM}) then
        ./xmlchange NTASKS=${NTASKS_OCN}
        ./xmlchange NTASKS_ATM=${NTASKS_ATM}
        ./xmlchange ROOTPE_ATM=${NNODES_OCN}
    else
        ./xmlchange NTASKS=${NTASKS_ATM}
    endif

    ./xmlchange NTHRDS_ATM=${NTHRDS}
    ./xmlchange MAX_MPITASKS_PER_NODE=${NTASKS_PER_NODE}
    ./xmlchange MAX_TASKS_PER_NODE=${NTASKS_PD}

#  ./xmlchange --id DEBUG --val ${debug_compile}
   # Set PIO format, use PIO version 2, and increase PIO buffer size 
   ./xmlchange PIO_NETCDF_FORMAT="64bit_data"
   ./xmlchange PIO_VERSION="2"
   #./xmlchange PIO_BUFFER_SIZE_LIMIT=64200000
   #./xmlchange PIO_REARR_COMM_MAX_PEND_REQ_COMP2IO=64
   ./xmlchange ATM_NCPL=${NCPL_ATM}
   #./xmlchange SCREAM_CMAKE_OPTIONS="SCREAM_NP 4 SCREAM_NUM_VERTICAL_LEV 128 SCREAM_NUM_TRACERS 10"
   #./xmlchange DEBUG=TRUE #debug rather than optimized build.

   ./xmlchange CAM_DYCORE=se
   ./xmlchange EPS_AGRID=${EPS_AGRID}

   ./xmlchange CAM_TARGET=theta-l
   #./xmlchange SSTICE_DATA_FILENAME=${sst_data}
   #./xmlchange SSTICE_YEAR_START=${RUN_START_YEAR},SSTICE_YEAR_END=${RUN_START_YEAR},SSTICE_YEAR_ALIGN=${RUN_START_YEAR}
   ./xmlchange GLC_AVG_PERIOD=glc_coupling_period

   # use following line if use elm init from a previous run 
   #./xmlchange ELM_NAMELIST_OPTS="use_init_interp=.true. init_interp_method='general'"
   

    # Edit CAM namelist to set dycore options for new grid
    cat <<EOF >> user_nl_eam

history_amwg        = .true.
history_aero_optics = .false.
history_aerosol     = .false.
history_budget      = .true.
history_verbose     = .true.
inithist            = 'MONTHLY'
inithist_all        = .true.
nhtfrq          =   0,  -6,  
mfilt           =   1, 120, 
avgflag_pertape = 'A', 'A', 

fexcl1 = 'CFAD_SR532_CAL', 'LINOZ_DO3', 'LINOZ_DO3_PSC', 'LINOZ_O3CLIM', 'LINOZ_O3COL', 'LINOZ_SSO3', 'hstobie_linoz',
fincl1 = 'extinct_sw_inp','extinct_lw_bnd7','extinct_lw_inp','CLD_CAL', 'TREFMNAV', 'TREFMXAV',
fincl2 = 'PRECT','PRECC','TMQ','TUQ','TVQ','TCO', 'SCO', 'SHFLX', 'LHFLX', 'QFLX', 
         'TGCLDCWP', 'TGCLDLWP', 'TGCLDIWP', 'TH7001000','THE7001000','FSNT', 'FLNT','FSDS', 'FLUT','PBLH',
         'QRL','QRS', 'FLDS', 'FLNS', 'FSDS', 'FSNS', 'FSNTOA','FSNTOAC','FSUTOA','FSUTOAC','FLNTC',
         'FSNTC','FSNSC','FLNSC','FLDS','SWCF','LWCF','CLDTOT', 'CLDHGH','CLDLOW','CLDMED','CLDTOT', 'AODVIS',
         'Nudge_U', 'Nudge_V','Nudge_T', 'Nudge_Q',

EOF

cat <<EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
 hist_mfilt  = 1, !365
 hist_nhtfrq = 0, !-24
 hist_avgflag_pertape = 'A', !'A'
 ! Override
 fsurdat="${lnd_surdat}"
 finidat="${lnd_data}"
 check_finidat_fsurdat_consistency = .false.
EOF

cat << EOF >> user_nl_mosart
 rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
 rtmhist_mfilt = 1, !365
 rtmhist_ndens = 2
 rtmhist_nhtfrq = 0, !-24
EOF

   # Finally, run setup
   ./case.setup 

   # The run location is determined in the bowels of CIME
   # Symlink to that location from user-chosen $case_root (=current dir)
   ln -s `./xmlquery -value RUNDIR` run

   # This disables the logic that sets tprof_n and tprof_options internally.
   ./xmlchange --file env_run.xml TPROF_TOTAL=-1
   echo "tprof_n = 1" >> user_nl_cpl
   echo "tprof_option = 'ndays'" >> user_nl_cpl

endif

#====================================================================
# Compile 
#====================================================================
if ($compile_model > 0) then # Build the model

   rm -rvf $CASE_BUILD_DIR
   cd $CASE_ROOT
   ./case.build

    # Set file striping on run dir for writing large files
    ls -l run > /dev/null #just using this command as a check that run dir exists
    lfs setstripe -S 1m -c 64 run
    lfs getstripe run >& lfs_run.txt

endif

#####################################################################
# Conduct simulation
#####################################################################
if ($run_model > 0) then
  cd $CASE_ROOT
  if ( $RESUBMIT > 0 ) then
    ./xmlchange  RESUBMIT=$RESUBMIT
  endif 
  if ( $MODEL_START_TYPE  == "initial" ) then
    ./xmlchange RUN_TYPE='startup'
    ./xmlchange CONTINUE_RUN='FALSE'
  else if ( $MODEL_START_TYPE  == "continue" ) then 
    ./xmlchange  CONTINUE_RUN='TRUE'
  else if ( $MODEL_START_TYPE  == "branch" || $MODEL_START_TYPE  == "hybrid" ) then 
    ./xmlchange  RUN_TYPE=$MODEL_START_TYPE
    ./xmlchange  GET_REFCASE=$GET_REFCASE
    ./xmlchange  RUN_REFDIR=$RUN_REFDIR
    ./xmlchange  RUN_REFCASE=$RUN_REFCASE
    ./xmlchange  RUN_REFDATE=$RUN_REFDATE
    ./xmlchange  CONTINUE_RUN='FALSE'
  else 
    echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
    exit 380
  endif 

  #Run the model 
  ./case.submit --batch-args="--mail-type=ALL --mail-user=${EMAIL}"
endif

echo "Done working in ${CASE_ROOT}"

