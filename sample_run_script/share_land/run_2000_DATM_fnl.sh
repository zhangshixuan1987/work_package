#!/bin/bash -fe

# E3SM Coupled Model Group run_e3sm script template.
#
# Bash coding style inspired by:
# http://kfirlavi.herokuapp.com/blog/2012/11/14/defensive-bash-programming

main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=pm-cpu
readonly PROJECT="m3525" #"e3sm"

########################################################################################
readonly COMPSET="20TR_DATM%ERA56HR_ELM%CNPRDCTCBCTOP_SICE_SOCN_MOSART_SGLC_SWAV_SIAC_SESP"
readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"
readonly CASE_NAME="test.v3-LR.DATM.ne30pg2_r05_IcoswISC30E3r5.${MACHINE}"

# If this is part of a simulation campaign, ask your group lead about using a case_group label
# otherwise, comment out
#readonly CASE_GROUP="v3.LR"

# Code and compilation
readonly CHECKOUT="20240924"
readonly BRANCH="v3.0.0"  # master as of 2024-03-04 399d4301138617088dd93214123d6c025e061302
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
readonly MODEL_START_TYPE='hybrid' #, 'continue', 'branch', 'hybrid'
readonly MODEL_START_YEAR=1979
readonly START_DATE=`printf "%04d" ${MODEL_START_YEAR}`"-01-01"
readonly RUN_ICYCLE=0

# Additional options for 'branch' and 'hybrid'
if [[ $RUN_ICYCLE == 0 ]];then 
  readonly RUN_REFDATE1="1979-01-01"
  readonly RUN_REFCASE1="v3.LR.historical_0101"
  readonly RUN_REFDIR1="/pscratch/sd/z/zhan391/e3sm_project/input/${RUN_REFCASE1}/${RUN_REFDATE1}-00000"
  readonly LND_INIT_FILE1="${RUN_REFDIR1}/${RUN_REFCASE1}.elm.r.${RUN_REFDATE1}-00000.nc"
  readonly RTM_INIT_FILE1="${RUN_REFDIR1}/${RUN_REFCASE1}.mosart.r.${RUN_REFDATE1}-00000.nc"
else 
  readonly RUN_REFDATE1="1980-01-01" 
  readonly RUN_REFCASE1=${CASE_NAME} 
  readonly RUN_REFDIR1="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_testings/${RUN_REFCASE1}/archive/rest/${RUN_REFDATE1}-00000"
  readonly LND_INIT_FILE1="${RUN_REFDIR1}/${RUN_REFCASE1}.elm.r.${RUN_REFDATE1}-00000.nc"
  readonly RTM_INIT_FILE1="${RUN_REFDIR1}/${RUN_REFCASE1}.mosart.r.${RUN_REFDATE1}-00000.nc"
fi 
readonly GET_REFCASE=TRUE
readonly RUN_REFDATE=${START_DATE}
readonly RUN_REFCASE=${RUN_REFCASE1}
readonly RUN_REFDIR="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_testings/${RUN_REFCASE}/archive/rest/${RUN_REFDATE}-00000"
if [ ! -d ${RUN_REFDIR} ];then
  mkdir -p ${RUN_REFDIR}
fi
#now copy initial condition file to targe location 
readonly LND_INIT_FILE="${RUN_REFDIR}/${RUN_REFCASE1}.elm.r.${RUN_REFDATE}-00000.nc"
readonly RTM_INIT_FILE="${RUN_REFDIR}/${RUN_REFCASE1}.mosart.r.${RUN_REFDATE}-00000.nc"
echo $LND_INIT_FILE
echo $RTM_INIT_FILE
echo ${RUN_REFDATE}_00000 > ${RUN_REFDIR}/rpointer.ocn
echo ${RUN_REFDATE}_00000 > ${RUN_REFDIR}/rpointer.ice
ff=`basename ${LND_INIT_FILE}`; echo ./${ff} > ${RUN_REFDIR}/rpointer.lnd
ff=`basename ${RTM_INIT_FILE}`; echo ${ff} > ${RUN_REFDIR}/rpointer.rof
#if [ -f ${LND_INIT_FILE} ]; then
#  mv ${LND_INIT_FILE} ${LND_INIT_FILE}.cyc${RUN_ICYCLE}
#fi
#if [ -f ${RTM_INIT_FILE} ]; then
#  mv ${RTM_INIT_FILE} ${RTM_INIT_FILE}.cyc${RUN_ICYCLE}
#fi
#cp -rp ${LND_INIT_FILE1} ${LND_INIT_FILE}
#cp -rp ${RTM_INIT_FILE1} ${RTM_INIT_FILE}

# Addition setup for data atmosphere 
readonly DATM_MODE="ERA56HR"
readonly DATM_HIST_YR_ALIGN="${MODEL_START_YEAR}"
readonly DATM_HIST_YR_START="1980"
readonly DATM_HIST_YR_END="1980"
readonly DATM_CO2_TSERIES="20tr"
readonly DATM_PRESAERO="clim_2000"

# Set paths
readonly CODE_ROOT="/pscratch/sd/z/zhan391/e3sm_project/code/${CHECKOUT}"
readonly CASE_ROOT="/pscratch/sd/z/zhan391/e3sm_project/E3SMv3_testings/${CASE_NAME}"

# Sub-directories
readonly CASE_BUILD_DIR=${CASE_ROOT}/build
readonly CASE_ARCHIVE_DIR=${CASE_ROOT}/archive

# Define type of run
#  short tests: 'XS_1x10_ndays', 'XS_2x5_ndays', 'S_1x10_ndays', 'M_1x10_ndays', 'L_1x10_ndays'
#  or 'production' for full simulation
#readonly run='L_1x10_ndays'  # build with this to ensure non-threading
#readonly run='S_1x10_ndays'
#readonly run='S_2x5_ndays'
#readonly run='M_1x10_ndays'
#readonly run='custom-40_1x10_ndays'
readonly run='production'

if [[ "${run}" != "production" ]]; then
  echo "setting up Short test simulations: ${run}"
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/tests/${run}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/tests/${run}/run
  readonly PELAYOUT=${layout}
  readonly WALLTIME="00:30:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=false
else
  echo "setting up ${run}"
  # Production simulation
  readonly CASE_SCRIPTS_DIR=${CASE_ROOT}/case_scripts
  readonly CASE_RUN_DIR=${CASE_ROOT}/run
  readonly PELAYOUT="custom-20" #"L"
  readonly WALLTIME="30:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="20"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="0"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="nyears"
readonly HIST_N="1"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=""

# --- Toggle flags for what to do ----
do_fetch_code=false
do_create_newcase=true
do_case_setup=true
do_case_build=true
do_case_submit=true

# --- Now, do the work ---
if [ "${do_create_newcase,,}" == "true" ]; then
  if [ -d ${CASE_SCRIPTS_DIR} ];then
    if [ ! -d ${CASE_SCRIPTS_DIR}.`date +%Y%m%d` ]; then 	  
      mv ${CASE_SCRIPTS_DIR} ${CASE_SCRIPTS_DIR}.`date +%Y%m%d`
    else 
      rm -rvf ${CASE_SCRIPTS_DIR}
    fi 
  fi
fi 

if [ "${do_case_build,,}" == "true" ]; then
  if [ -d ${CASE_BUILD_DIR} ];then
    rm -rvf ${CASE_BUILD_DIR}
  fi
fi  

# Make directories created by this script world-readable
umask 022

# Fetch code from Github
fetch_code

# Create case
create_newcase

# Custom PE layout
custom_pelayout

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

cat << EOF >> user_nl_datm
 !add setup for datm here if needed
EOF

cat << EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.
 hist_fexcl1 = 'AGWDNPP','ALTMAX_LASTYEAR','AVAIL_RETRANSP','AVAILC','BAF_CROP',
               'BAF_PEATF','BIOCHEM_PMIN_TO_PLANT','CH4_SURF_AERE_SAT','CH4_SURF_AERE_UNSAT','CH4_SURF_DIFF_SAT',
               'CH4_SURF_DIFF_UNSAT','CH4_SURF_EBUL_SAT','CH4_SURF_EBUL_UNSAT','CMASS_BALANCE_ERROR','cn_scalar',
               'COL_PTRUNC','CONC_CH4_SAT','CONC_CH4_UNSAT','CONC_O2_SAT','CONC_O2_UNSAT',
               'cp_scalar','CWDC_HR','CWDC_LOSS','CWDC_TO_LITR2C','CWDC_TO_LITR3C',
               'CWDC_vr','CWDN_TO_LITR2N','CWDN_TO_LITR3N','CWDN_vr','CWDP_TO_LITR2P',
               'CWDP_TO_LITR3P','CWDP_vr','DWT_CONV_CFLUX_DRIBBLED','F_CO2_SOIL','F_CO2_SOIL_vr',
               'F_DENIT_vr','F_N2O_DENIT','F_N2O_NIT','F_NIT_vr','FCH4_DFSAT',
               'FINUNDATED_LAG','FPI_P_vr','FPI_vr','FROOTC_LOSS','HR_vr',
               'LABILEP_TO_SECONDP','LABILEP_vr','LAND_UPTAKE','LEAF_MR','leaf_npimbalance',
               'LEAFC_LOSS','LEAFC_TO_LITTER','LFC2','LITR1_HR','LITR1C_TO_SOIL1C',
               'LITR1C_vr','LITR1N_TNDNCY_VERT_TRANS','LITR1N_TO_SOIL1N','LITR1N_vr','LITR1P_TNDNCY_VERT_TRANS',
               'LITR1P_TO_SOIL1P','LITR1P_vr','LITR2_HR','LITR2C_TO_SOIL2C','LITR2C_vr',
               'LITR2N_TNDNCY_VERT_TRANS','LITR2N_TO_SOIL2N','LITR2N_vr','LITR2P_TNDNCY_VERT_TRANS','LITR2P_TO_SOIL2P',
               'LITR2P_vr','LITR3_HR','LITR3C_TO_SOIL3C','LITR3C_vr','LITR3N_TNDNCY_VERT_TRANS',
               'LITR3N_TO_SOIL3N','LITR3N_vr','LITR3P_TNDNCY_VERT_TRANS','LITR3P_TO_SOIL3P','LITR3P_vr',
               'M_LITR1C_TO_LEACHING','M_LITR2C_TO_LEACHING','M_LITR3C_TO_LEACHING','M_SOIL1C_TO_LEACHING','M_SOIL2C_TO_LEACHING',
               'M_SOIL3C_TO_LEACHING','M_SOIL4C_TO_LEACHING','NDEPLOY','NEM','nlim_m',
               'o2_decomp_depth_unsat','OCCLP_vr','PDEPLOY','PLANT_CALLOC','PLANT_NDEMAND',
               'PLANT_NDEMAND_COL','PLANT_PALLOC','PLANT_PDEMAND','PLANT_PDEMAND_COL','plim_m',
               'POT_F_DENIT','POT_F_NIT','POTENTIAL_IMMOB','POTENTIAL_IMMOB_P','PRIMP_TO_LABILEP',
               'PRIMP_vr','PROD1P_LOSS','QOVER_LAG','RETRANSN_TO_NPOOL','RETRANSP_TO_PPOOL',
               'SCALARAVG_vr','SECONDP_TO_LABILEP','SECONDP_TO_OCCLP','SECONDP_vr','SMIN_NH4_vr',
               'SMIN_NO3_vr','SMINN_TO_SOIL1N_L1','SMINN_TO_SOIL2N_L2','SMINN_TO_SOIL2N_S1','SMINN_TO_SOIL3N_L3',
               'SMINN_TO_SOIL3N_S2','SMINN_TO_SOIL4N_S3','SMINP_TO_SOIL1P_L1','SMINP_TO_SOIL2P_L2','SMINP_TO_SOIL2P_S1',
               'SMINP_TO_SOIL3P_L3','SMINP_TO_SOIL3P_S2','SMINP_TO_SOIL4P_S3','SMINP_vr','SOIL1_HR','SOIL1C_TO_SOIL2C','SOIL1C_vr','SOIL1N_TNDNCY_VERT_TRANS','SOIL1N_TO_SOIL2N','SOIL1N_vr',
               'SOIL1P_TNDNCY_VERT_TRANS','SOIL1P_TO_SOIL2P','SOIL1P_vr','SOIL2_HR','SOIL2C_TO_SOIL3C',
               'SOIL2C_vr','SOIL2N_TNDNCY_VERT_TRANS','SOIL2N_TO_SOIL3N','SOIL2N_vr','SOIL2P_TNDNCY_VERT_TRANS',
               'SOIL2P_TO_SOIL3P','SOIL2P_vr','SOIL3_HR','SOIL3C_TO_SOIL4C','SOIL3C_vr',
               'SOIL3N_TNDNCY_VERT_TRANS','SOIL3N_TO_SOIL4N','SOIL3N_vr','SOIL3P_TNDNCY_VERT_TRANS','SOIL3P_TO_SOIL4P',
               'SOIL3P_vr','SOIL4_HR','SOIL4C_vr','SOIL4N_TNDNCY_VERT_TRANS','SOIL4N_TO_SMINN',
               'SOIL4N_vr','SOIL4P_TNDNCY_VERT_TRANS','SOIL4P_TO_SMINP','SOIL4P_vr','SOLUTIONP_vr',
               'TCS_MONTH_BEGIN','TCS_MONTH_END','TOTCOLCH4','water_scalar','WF',
               'wlim_m','WOODC_LOSS','WTGQ'
 hist_fincl1 = 'SNOWDP','COL_FIRE_CLOSS','NPOOL','PPOOL','TOTPRODC'
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA'
 hist_mfilt = 1,365
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 create_crop_landunit = .false.
 check_finidat_pct_consistency = .false.
 finidat = "${LND_INIT_FILE}"
 
! fatmlndfrc = '\${DIN_LOC_ROOT}/share/domains/domain.lnd.r05_IcoswISC30E3r5.231121.nc'
 fsurdat = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr1850_c200609_with_TOP.nc'
 fsnowaging = '\${DIN_LOC_ROOT}/lnd/clm2/snicardata/snicar_drdt_bst_fit_60_c070416.nc'
 fsnowoptics = '\${DIN_LOC_ROOT}/lnd/clm2/snicardata/snicar_optics_5bnd_mam_c160322.nc'
 fsoilordercon = '\${DIN_LOC_ROOT}/lnd/clm2/paramdata/CNP_parameters_c180529.nc'

 ! default landuse.timeseries file below probably was created using CMIP5 raw land data.
 flanduse_timeseries = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_hist_simyr1850-2015_c240308.nc'

 !below used by v3.LR fully coupled model
 flndtopo = '\$DIN_LOC_ROOT/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr2000_c200609_with_TOP.nc'
 stream_fldfilename_ndep = '\$DIN_LOC_ROOT/lnd/clm2/ndepdata/fndep_elm_cbgc_exp_simyr1849-2101_1.9x2.5_c190103.nc'
 stream_fldfilename_popdens = '\$DIN_LOC_ROOT/lnd/clm2/firedata/elmforc.ssp5_hdm_0.5x0.5_simyr1850-2100_c190109.nc'
 tw_irr = .false.
 use_atm_downscaling_to_topunit = .false.
EOF

cat << EOF >> user_nl_mosart 
 coupling_period = 10800
 finidat_rtm = "${RTM_INIT_FILE}"
EOF

cat << EOF >> user_nl_cpl
  flds_co2a = .true.
EOF

}

# =====================================
# Customize MPAS stream files if needed
# =====================================

patch_mpas_streams() {

echo
echo 'Modifying MPAS streams files'

}

# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout(){

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
    if [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
        hthrd=2  # hyper-threading
    elif [ "${MACHINE}" == "pm-cpu" ]; then
        ncore=128
        hthrd=1  # hyper-threading
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for current custom PE layout setting.'
        exit 400
    fi

    # Layout-specific customization
    # Applicable to all custom layouts
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=1
    ./xmlchange NTHRDS=1
    ./xmlchange ROOTPE=0
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$(( $ncore * $hthrd))

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}
    ./xmlchange NTASKS=$(( $nnodes * $ncore ))

    popd

fi

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
    local repo=E3SM

    echo "Cloning $repo repository branch $BRANCH under $path"
    if [ -d "${path}" ]; then
        echo "ERROR: Directory already exists. Not overwriting"
        exit 20
    fi
    mkdir -p ${path}
    pushd ${path}

    # This will put repository, with all code
    git clone git@github.com:E3SM-Project/${repo}.git .

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

    if [[ ${PELAYOUT} == custom-* ]];
    then
        layout="M" # temporary placeholder for create_newcase
    else
        layout=${PELAYOUT}

    fi

    # Base arguments
    args=" --case ${CASE_NAME} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} \
        --walltime ${WALLTIME} \
        --pecount ${PELAYOUT}"

    # Oprional arguments
    if [ ! -z "${PROJECT}" ]; then
      args="${args} --project ${PROJECT}"
    fi
    if [ ! -z "${CASE_GROUP}" ]; then
      args="${args} --case-group ${CASE_GROUP}"
    fi
    if [ ! -z "${QUEUE}" ]; then
      args="${args} --queue ${QUEUE}"
    fi

    ${CODE_ROOT}/cime/scripts/create_newcase ${args}

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

    # Build with COSP, except for a data atmosphere (datm)
    if [ `./xmlquery --value COMP_ATM` == "datm"  ]; then 
      echo $'\nThe specified configuration uses a data atmosphere, so cannot activate COSP simulator\n'
    else
      echo $'\nConfiguring E3SM to use the COSP simulator\n'
      ./xmlchange --id CAM_CONFIG_OPTS --append --val='-cosp'
    fi

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    #DATM setup   
    ./xmlchange DATM_MODE=${DATM_MODE} 
    ./xmlchange DATM_CLMNCEP_YR_ALIGN=${DATM_HIST_YR_ALIGN}
    ./xmlchange DATM_CLMNCEP_YR_START=${DATM_HIST_YR_START}
    ./xmlchange DATM_CLMNCEP_YR_END=${DATM_HIST_YR_END}
    ./xmlchange DATM_CO2_TSERIES=${DATM_CO2_TSERIES}
    ./xmlchange DATM_PRESAERO=${DATM_PRESAERO}

    #initialization of elm model setup 
    ./xmlchange ELM_FORCE_COLDSTART="off"

    #configuration followed by v3.LR fully coupled model 
    ./xmlchange ROF_NCPL=8
    ./xmlchange CCSM_BGC="CO2A"
    ./xmlchange ELM_CO2_TYPE="diagnostic"
    ./xmlchange BARRIER_OPTION="never"
    ./xmlchange --id ELM_BLDNML_OPTS --append --val='-bgc bgc -nutrient cnp -nutrient_comp_pathway rd -soil_decomp ctc -methane -solar_rad_scheme top'

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

    fi

    # Some user_nl settings won't be updated to *_in files under the run directory
    # Call preview_namelists to make sure *_in and user_nl files are consistent.
    echo $'\n----- Preview namelists -----\n'
    ./preview_namelists

    popd
}

#-----------------------------------------------------
runtime_options() {

    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}

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
	echo "Prepare the restart files - copy restart-point files over to ../run for the relocated case"

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

