#!/bin/bash -fe
ulimit -S -c 0

# E3SM Water Cycle v2 run_e3sm script template.
#
# Inspired by v1 run_e3sm script as well as SCREAM group simplified run script.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=compy
readonly PROJECT="e3sm"

# Simulation
readonly COMPSET="F2010"
readonly RESOLUTION="SGP_ne32x16pg2"
readonly CASE_NAME="v2.SGP.ne512.F2010.p30"
readonly CASE_GROUP="v2_SGP"

# Code and compilation
readonly CHECKOUT="E3SM_maint-2.0"
readonly BRANCH="bcde284dc9ded281edf5d398265b236b3e3c4bbf" # master as of 20211022
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# PE layout
readonly NTASKS=1200
readonly NTHRDS=1

# Run options
readonly MODEL_START_TYPE="initial"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="0001-01-01"

# Additional options for 'branch' and 'hybrid'
#readonly GET_REFCASE=TRUE
#readonly RUN_REFDIR="/lcrc/group/e3sm/ac.golaz/E3SMv2/v2.LR.piClim-histaer_0021/init"
#readonly RUN_REFCASE="v2.LR.piClim-control"
#readonly RUN_REFDATE="0021-01-01"   # same as MODEL_START_DATE for 'branch', can be different for 'hybrid'

# Set paths
readonly CODE_ROOT="${HOME}/compy/model/sp/${CHECKOUT}"
readonly CASE_ROOT="${HOME}/compy/model/sp/${CHECKOUT}/cases/${CASE_NAME}"
readonly CRUN_ROOT="/compyfs/zhan524/e3sm_scratch/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CRUN_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CRUN_ROOT}/archive

# Define type of run
#  short tests: 'XS_2x5_ndays', 'XS_1x10_ndays', 'S_1x10_ndays', 
#               'M_1x10_ndays', 'L_1x10_ndays', 'XL_1x10_ndays'
#  or 'production' for full simulation
readonly run='F2010'

if [ "${run}" != "production" ]; then

  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}
  readonly CASE_RUN_DIR=${CRUN_ROOT}/run
  readonly PELAYOUT="M"
  readonly WALLTIME="01:59:00"
  readonly STOP_OPTION="ndays"
  readonly STOP_N=1
  readonly REST_OPTION="ndays"
  readonly REST_N=1
  readonly RESUBMIT=0
  readonly DO_SHORT_TERM_ARCHIVING=false

else

  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="L"
  readonly WALLTIME="48:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="55"
  readonly REST_OPTION="nyears"
  readonly REST_N="5"
  readonly RESUBMIT="2"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="nyears"
readonly HIST_N="5"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=false

# --- Now, do the work ---

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Setup
case_setup

# Build
case_build

# Configure runtime options
runtime_options

# Copy script into case_script directory for provenance
copy_script

# Submit
case_submit

# All done
echo $'\n----- All done -----\n'

}

# =======================
# Custom user_nl settings
# =======================

user_nl() {

cat << EOF >> user_nl_eam
nu_top=4e4 
se_tstep=33.33333333333
drydep_srf_file = '/compyfs/inputdata/atm/cam/chem/trop_mam/atmsrf_SGP_ne32x8pg2_20230123.nc'
inithist = 'MONTHLY'
inithist_all = .true.
fincl1  = 'PS',
          'AODVIS',  !! AeroCom indirect diagnostics 
          'angstrm',
          'cod',
          'cdr',
          'cdnc',
          'cdnum',
          'icnum',
          'clt',
          'lcc',
          'lwp',
          'iwp',
          'icc',
          'icnc',
          'icr',
          'LHFLX',
          'SHFLX',
          'OMEGA500',
          'rh700',
          'colrv',
          'ccn',
          'ccn.1bl',
          'ccn.3bl',
          'ptop',
          'ttop',
          'rwp',
          'lwp2',
          'iwp2',
          'autoconv',
          'accretn',
          'FSUTOA_d1',
          'FSUTOAC_d1',
          'FSUTOA',
          'FSUTOAC',
          'FLUTC',
          'FLUT', 
          'PRECC', 
          'PRECL', 
          'PRECT', 
          'TH7001000',
docosp    = .false.,
cosp_lite = .true.,
cosp_ncolumns        = 10
cosp_nradsteps       = 3
cosp_lmisr_sim       = .true.
cosp_lisccp_sim      = .true.
cosp_lmodis_sim      = .true.
cosp_llidar_sim      = .true.
history_amwg         = .true.
history_aero_optics  = .true.
history_aerosol      = .true.
history_clubb        = .true.
history_budget       = .true.
history_verbose      = .true.
hist_hetfrz_classnuc = .true. 
do_aerocom_ind3      = .true.
rad_diag_1           = 'A:Q:H2O', 'N:O2:O2', 'N:CO2:CO2', 'A:O3:O3', 'N:N2O:N2O', 'N:CH4:CH4', 'N:CFC11:CFC11', 'N:CFC12:CFC12',

!!!!.......................................................
!!!! nudging
!!!!
!!!! Default setup in E3SMv1 
!!!!    Nudge_Tau         = -999
!!!!    Nudge_Loc_PhysOut = .False.
!!!!    Nudge_CurrentStep = .False.
!!!!    Nudge_File_Ntime  = 1
!!!!    Nudge_Method      = ‘Step’
!!!!
!!!! Setup for MERRA2 
!!!!    Nudge_Path           = /qfs/projects/eagles/pma/merra2/ne30np4/'
!!!!    Nudge_File_Template  = 'MERRA2_ne30np4_%y-%m-%d-%s.nc'
!!!!    Nudge_Times_Per_Day  = 8  !! nudging input data frequency
!!!!    Nudge_File_Ntime     = 1 
!!!!.......................................................
!!!Nudge_Model          = .True.
!!!Nudge_Path           = '/compyfs/zhan524/myinput/ndata/eraint_ne30_pg2L72/'
!!!Nudge_File_Template  = 'eraint_ne30_pg2L72_%y%m%d00.nc'
!!!Nudge_Times_Per_Day  = 4  !! nudging input data frequency
!!!Model_Times_Per_Day  = 48 !! should not be larger than 48 if dtime = 1800s
!!!Nudge_Uprof          = 2
!!!Nudge_Ucoef          = 1.
!!!Nudge_Vprof          = 2
!!!Nudge_Vcoef          = 1.
!!!Nudge_Tprof          = 0
!!!Nudge_Tcoef          = 0.
!!!Nudge_Qprof          = 0
!!!Nudge_Qcoef          = 0.
!!!Nudge_PSprof         = 0
!!!Nudge_PScoef         = 0.
!!!Nudge_Beg_Year       = 0001
!!!Nudge_Beg_Month      = 1
!!!Nudge_Beg_Day        = 1
!!!Nudge_End_Year       = 9999
!!!Nudge_End_Month      = 1
!!!Nudge_End_Day        = 1
!!!Nudge_Vwin_Lindex    = 0.
!!!Nudge_Vwin_Hindex    = 70.
!!!Nudge_Vwin_Ldelta    = 0.1
!!!Nudge_Vwin_Hdelta    = 0.1
!!!Nudge_Vwin_lo        = 0.
!!!Nudge_Vwin_hi        = 1.
!!!Nudge_Method         = 'Linear'
!!!Nudge_Loc_PhysOut    = .True.
!!!Nudge_Tau            = 6.        !! relaxation time scale, unit: 6h
!!!Nudge_CurrentStep    = .False.
!!!Nudge_File_Ntime     = 4

! Historical, vs single forcing configurations

! | Configuration      | GHGs      | Aerosols and | Ozone     | Solar     | Volcanoes | Land use
! |                    |           | precursors   |           |           |           |         
! -----------------------------------------------------------------------------------------------
! | historical         | varying   | varying      | varying   | varying   | varying   | varying
! | hist-GHG           | varying   | 1850         | 1850      | 1850      | 1850      | 1850
! | hist-aer           | 1850      | varying      | 1850      | 1850      | 1850      | 1850
! | hist-all-xGHG-xaer | 1850      | 1850         | varying   | varying   | varying   | varying

!!!! (1) GHGs settings
!!!
!!! bndtvghg		= ' '
!!! ch4vmr		= 808.249e-9
!!! co2vmr		= 284.317000e-6
!!! f11vmr		= 32.1102e-12
!!! f12vmr		= 0.0
!!! flbc_list		= ' '
!!! n2ovmr		= 273.0211e-9
!!! scenario_ghg		= 'FIXED'

! (2) aeorosols and precursors

!!!! (3) ozone
!!!
!!! chlorine_loading_fixed_ymd		= 18500101
!!! chlorine_loading_type		= 'FIXED'
!!!
!!! linoz_data_cycle_yr		= 1850
!!! linoz_data_type		= 'CYCLICAL'
!!!
!!!! (4) solar
!!!
!!! solar_data_file		= '${input_data_dir}/atm/cam/solar/Solar_1850control_input4MIPS_c20181106.nc'
!!! solar_data_type		= 'FIXED'
!!! solar_data_ymd		= 18500101
!!!
!!!! (5) volcanoes
!!!
!!! prescribed_volcaero_cycle_yr		= 1
!!! prescribed_volcaero_file		= 'CMIP_DOE-ACME_radiation_average_1850-2014_v3_c20171204.nc'
!!! prescribed_volcaero_type		= 'CYCLICAL'

EOF

cat << EOF >> user_nl_elm
! hist_dov2xy = .true.,.true.
! hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
! hist_mfilt = 1,365
! hist_nhtfrq = 0,-24
! hist_avgflag_pertape = 'A','A'

! Override
check_finidat_fsurdat_consistency = .false.
!!!finidat = '/compyfs/liji711/RRM_IC_generation/SGP_ne32x8_elmini.nc'

EOF

###cat << EOF >> user_nl_mosart
### rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
### rtmhist_mfilt = 1,365
### rtmhist_ndens = 2
### rtmhist_nhtfrq = 0,-24
###EOF

# Override SST and sea-ice datasets
#./xmlchange SSTICE_DATA_FILENAME=$input_data_dir/ocn/docn7/SSTDATA/sst_ice_v2.LR.piControl_0.5x0.5_climo_0001-0500.nc
#./xmlchange SSTICE_GRID_FILENAME=$input_data_dir/ocn/docn7/domain.ocn.0.5x0.5.c211007.nc
#./xmlchange SSTICE_YEAR_ALIGN=1
#./xmlchange SSTICE_YEAR_START=0
#./xmlchange SSTICE_YEAR_END=0

}

patch_mpas_streams() {

echo

}

######################################################
### Most users won't need to change anything below ###
######################################################

#-----------------------------------------------------
fetch_code() {

    if [ "${do_fetch_code,,}" != "true" ]; then
        echo $'\n----- Skipping fetch_code -----\n'
        return
    fi

    echo $'\n----- Starting fetch_code -----\n'
    local path=${CODE_ROOT}
    local repo=e3sm

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .
    
    # Setup git hooks
    rm -rf .git/hooks
    git clone git@github.com:E3SM-Project/E3SM-Hooks.git .git/hooks
    git config commit.template .git/hooks/commit.template

    # Check out desired branch
    git checkout ${BRANCH}

    # Custom addition
    if [ "${CHERRY}" != "" ]; then
        echo ----- WARNING: adding git cherry-pick -----
        for commit in "${CHERRY[@]}"
        do
            echo ${commit}
            git cherry-pick ${commit}
        done
        echo -------------------------------------------
    fi

    # Bring in all submodule components
    git submodule update --init --recursive

    popd
}

#-----------------------------------------------------
create_newcase() {

    if [ "${do_create_newcase,,}" != "true" ]; then
        echo $'\n----- Skipping create_newcase -----\n'
        return
    fi

    echo $'\n----- Starting create_newcase -----\n'

    ${CODE_ROOT}/cime/scripts/create_newcase \
        --case ${CASE_NAME} \
        --case-group ${CASE_GROUP} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} \
        --project ${PROJECT} \
        --walltime ${WALLTIME} #\
#       --pecount ${PELAYOUT}

    if [ $? != 0 ]; then
      echo $'\nNote: if create_newcase failed because sub-directory already exists:'
      echo $'  * delete old case_script sub-directory'
      echo $'  * or set do_newcase=false\n'
      exit 35
    fi

}

#-----------------------------------------------------
case_setup() {

    if [ "${do_case_setup,,}" != "true" ]; then
        echo $'\n----- Skipping case_setup -----\n'
        return
    fi

    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # Setup some CIME directories
    ./xmlchange EXEROOT=${CASE_BUILD_DIR}
    ./xmlchange RUNDIR=${CASE_RUN_DIR}

    # Short term archiving
    ./xmlchange DOUT_S=${DO_SHORT_TERM_ARCHIVING^^}
    ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

   ./xmlchange NTASKS_ATM=$NTASKS
   ./xmlchange NTHRDS_ATM=$NTHRDS
   ./xmlchange ROOTPE_ATM='0'

   ./xmlchange NTASKS_LND=$NTASKS
   ./xmlchange NTHRDS_LND=$NTHRDS
   ./xmlchange ROOTPE_LND='0'

   ./xmlchange NTASKS_ROF=$NTASKS
   ./xmlchange NTHRDS_ROF=$NTHRDS
   ./xmlchange ROOTPE_ROF='0'

   ./xmlchange NTASKS_ICE=$NTASKS
   ./xmlchange NTHRDS_ICE=$NTHRDS
   ./xmlchange ROOTPE_ICE='0'

   ./xmlchange NTASKS_OCN=$NTASKS
   ./xmlchange NTHRDS_OCN=$NTHRDS
   ./xmlchange ROOTPE_OCN='0'

   ./xmlchange NTASKS_GLC=$NTASKS
   ./xmlchange NTHRDS_GLC=$NTHRDS
   ./xmlchange ROOTPE_GLC='0'

   ./xmlchange NTASKS_WAV=$NTASKS
   ./xmlchange NTHRDS_WAV=$NTHRDS
   ./xmlchange ROOTPE_WAV='0'

   ./xmlchange NTASKS_CPL=$NTASKS
   ./xmlchange NTHRDS_CPL=$NTHRDS
   ./xmlchange ROOTPE_CPL='0'

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    # Custom user_nl
    user_nl

    # Finally, run CIME case.setup
    ./case.setup --reset

    popd
}

#-----------------------------------------------------
case_build() {

    pushd ${CASE_SCRIPTS_DIR}

    # do_case_build = false
    if [ "${do_case_build,,}" != "true" ]; then

        echo $'\n----- case_build -----\n'

        if [ "${OLD_EXECUTABLE}" == "" ]; then
            # Ues previously built executable, make sure it exists
            if [ -x ${CASE_BUILD_DIR}/e3sm.exe ]; then
                echo 'Skipping build because $do_case_build = '${do_case_build}
            else
                echo 'ERROR: $do_case_build = '${do_case_build}' but no executable exists for this case.'
                exit 297
            fi
        else
            # If absolute pathname exists and is executable, reuse pre-exiting executable
            if [ -x ${OLD_EXECUTABLE} ]; then
                echo 'Using $OLD_EXECUTABLE = '${OLD_EXECUTABLE}
                cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
            else
                echo 'ERROR: $OLD_EXECUTABLE = '$OLD_EXECUTABLE' does not exist or is not an executable file.'
                exit 297
            fi
        fi
        echo 'WARNING: Setting BUILD_COMPLETE = TRUE.  This is a little risky, but trusting the user.'
        ./xmlchange BUILD_COMPLETE=TRUE

    # do_case_build = true
    else

        echo $'\n----- Starting case_build -----\n'

        # Turn on debug compilation option if requested
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi

        # Run CIME case.build
        ./case.build

        # Some user_nl settings won't be updated to *_in files under the run directory
        # Call preview_namelists to make sure *_in and user_nl files are consistent.
        ./preview_namelists

    fi

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

    # set larger threshold
    ./xmlchange EPS_AGRID=1e-9

    # Set simulation start date
    ./xmlchange RUN_STARTDATE=${START_DATE}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Set resubmissions
    if (( RESUBMIT > 0 )); then
        ./xmlchange RESUBMIT=${RESUBMIT}
    fi

    # Run type
    # Start from default of user-specified initial conditions
    if [ "${MODEL_START_TYPE,,}" == "initial" ]; then
        ./xmlchange RUN_TYPE="startup"
        ./xmlchange CONTINUE_RUN="FALSE"

    # Continue existing run
    elif [ "${MODEL_START_TYPE,,}" == "continue" ]; then
        ./xmlchange CONTINUE_RUN="TRUE"

    elif [ "${MODEL_START_TYPE,,}" == "branch" ] || [ "${MODEL_START_TYPE,,}" == "hybrid" ]; then
        ./xmlchange RUN_TYPE=${MODEL_START_TYPE,,}
        ./xmlchange GET_REFCASE=${GET_REFCASE}
	./xmlchange RUN_REFDIR=${RUN_REFDIR}
        ./xmlchange RUN_REFCASE=${RUN_REFCASE}
        ./xmlchange RUN_REFDATE=${RUN_REFDATE}
        echo 'Warning: $MODEL_START_TYPE = '${MODEL_START_TYPE} 
	echo '$RUN_REFDIR = '${RUN_REFDIR}
	echo '$RUN_REFCASE = '${RUN_REFCASE}
	echo '$RUN_REFDATE = '${START_DATE}
 
    else
        echo 'ERROR: $MODEL_START_TYPE = '${MODEL_START_TYPE}' is unrecognized. Exiting.'
        exit 380
    fi

    # Patch mpas streams files
    patch_mpas_streams

    popd
}

#-----------------------------------------------------
case_submit() {

    if [ "${do_case_submit,,}" != "true" ]; then
        echo $'\n----- Skipping case_submit -----\n'
        return
    fi

    echo $'\n----- Starting case_submit -----\n'
    pushd ${CASE_SCRIPTS_DIR}
    
    # Run CIME case.submit
    ./case.submit

    popd
}

#-----------------------------------------------------
copy_script() {

    echo $'\n----- Saving run script for provenance -----\n'

    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    local this_script_name=`basename $0`
    local script_provenance_name=${this_script_name}.`date +%Y%m%d-%H%M%S`
    cp -vp ${this_script_name} ${script_provenance_dir}/${script_provenance_name}

}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Now, actually run the script
#-----------------------------------------------------
main



