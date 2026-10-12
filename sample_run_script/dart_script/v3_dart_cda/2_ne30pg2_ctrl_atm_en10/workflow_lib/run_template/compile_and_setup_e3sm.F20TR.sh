#!/bin/bash -fe
main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE=my_machine #compy
readonly PROJECT=my_project #"esmd" 
readonly JOB_SLURM=my_jobqueue #"short"
readonly JOB_NTASKS=my_task_per_node
#readonly COMPILER="--compiler intel" # spec must include --compiler, otherwise empty it or comment out the whole line

# Simulation
readonly DATE=`date +"%m-%d-%y"`
readonly COMPSET=my_compset #"F20TR"
readonly RESOLUTION=my_resolution #"ne11np4_EC30to60E2r2"
readonly CASE_NAME=my_casename #"${COMPSET}_${RESOLUTION}_${ENSNAM}_${MACHINE}"

# Code and compilation
readonly CHECKOUT="20230715" # 
readonly BRANCH=""           #  4361928bd23968e00e3d6c04212a0ceec21e6bf6
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
#readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE=my_casedate
readonly START_TOD=my_casetod
readonly START_HOUR=`echo "${START_TOD}" / 3600 | bc`

readonly ATM_TIME=${RUN_REFDATE}-${RUN_REFTOD}
readonly OCN_TIME=${RUN_REFDATE}_${RUN_REFTOD}

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE #FALSE
readonly RUN_REFDATE=my_refdate
readonly RUN_REFTOD=my_reftod
readonly RUN_REFCASE=my_refcase
readonly RUN_REFDIR=my_refdir

readonly MODEL_START_TYPE="branch"
readonly rof_init="${RUN_REFDIR}/${RUN_REFCASE}.mosart.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
readonly atm_init="${RUN_REFDIR}/${RUN_REFCASE}.eam.i.${RUN_REFDATE}-${RUN_REFTOD}.nc"
readonly ice_init="${RUN_REFDIR}/${RUN_REFCASE}.mpassi.rst.${RUN_REFDATE}_${RUN_REFTOD}.nc"
readonly lnd_init="${RUN_REFDIR}/${RUN_REFCASE}.elm.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
echo $atm_init 
echo $lnd_init
echo $rof_init
echo $ice_init

#use below if current day 
readonly lnd_surdat="surfdata_0.5x0.5_simyr2010_c230922_with_TOP.nc"
#use below if pre-industrial
#readonly lnd_surdat="surfdata_0.5x0.5_simyr1850_c200609_with_TOP.nc"

#readonly sst_data="/compyfs/zhan391/acme_init/SST_forcing/sst_weekly_cdcunits_1x1_1990-01-01-2013-08-11_n3.nc"
#readonly sst_grid="/compyfs/zhan391/acme_init/SST_forcing/domain.ocn.1x1.111007.nc"
#readonly sst_syear="1990"
#readonly sst_eyear="2013"

readonly sst_data="/compyfs/zhan391/acme_init/SST_forcing/sst_ice_NOAA_AVHRR_E3SM_1x1_c20231225.nc"
readonly sst_grid="/compyfs/zhan391/acme_init/SST_forcing/domain.ocn.1x1.111007.nc"
readonly sst_syear="1981"
readonly sst_eyear="2020"

# Set paths
#"/qfs/people/zhan391/e3sm_dart_work/code/E3SM-maint-2.1"
readonly CODE_ROOT=my_e3sm_code
readonly CASE_ROOT=my_case_root
readonly CASE_BUILD_DIR=my_case_build
readonly CASE_RUN_DIR=my_case_run 
readonly CASE_ARCHIVE_DIR=my_case_archive
readonly CASE_SCRIPTS_DIR=my_case_script

readonly run=my_layout #'custom-24_1x6_nhours'

if [ "${run}" != "production" ]; then
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}
  readonly PELAYOUT=${layout}
  readonly WALLTIME="06:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=true
else
  # Production simulation
  readonly PELAYOUT="L"
  readonly WALLTIME="06:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="1"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="10"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="ndays" 
readonly HIST_N="5"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE=old_modelexe

# --- Now, do the work ---

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

cat << EOF >> user_nl_eam
 ncdata          = '$atm_init'
 inithist        = '6-HOURLY'
 inithist_all    = .true.

 cosp_lite = .true.

 empty_htapes = .true.

 ! All EAM history is 6-hourly: each DA cycle is a separate 6-h run, so daily
 ! or monthly averages would only cover the last cycle. Daily/monthly means are
 ! built in post-processing. T/U/V/Q/PS/CLDLIQ/CLDICE come from the DART output.
 avgflag_pertape = 'A','A','I'
 nhtfrq = -6,-6,-6
 mfilt  = 4,4,4

 ! 6hourly (h0) A: core aerosol/chemistry diagnostics (2-D)
 fincl1 = 'AODALL','AODVIS','AODABS','AODBC','AODDUST','AODPOM','AODSO4','AODSOA','AODSS',
          'ABURDENBC','ABURDENDUST','ABURDENMOM','ABURDENPOM','ABURDENSEASALT','ABURDENSO4','ABURDENSOA',
          'O3_SRF','CO_SRF','TCO','SCO','CDNUMC',
          'SFDMS','SFSO2','TROP_P','SOLIN','PHIS',
          'LANDFRAC','OCNFRAC','ICEFRAC'

 ! 6hourly (h1): fields not in DART or h2; 6-h min/max for daily extremes
 fincl2 = 'OMEGA1000:I','OMEGA975:I','OMEGA950:I','OMEGA925:I','OMEGA900:I','OMEGA850:I',
          'OMEGA800:I','OMEGA700:I','OMEGA600:I','OMEGA500:I','OMEGA400:I',
          'OMEGA300:I','OMEGA200:I','OMEGA100:I',
          'TTQ:I','TREFHTMN:M','TREFHTMX:X','PRECTMX:X'

 ! 6hourly (h2) A
 fincl3 = 'RHREFHT','TTOP','TOZ','PS','PSL','OMEGA500:A','PRECT:A','PRECL:A','PRECC:A','PRECSL:A',
          'TS:A','CLDLOW:A','CLDMED:A','CLDHGH:A','CLDTOT:A','TGCLDLWP:A','TGCLDIWP:A','TGCLDCWP:A',
          'TMQ:A','CAPE:A','CIN:A','TS:A','TREFHT:A','QREFHT:A','U10:A','TAUX:A','TAUY:A',
          'SHFLX:A','LHFLX:A','TUQ:A','TVQ:A',
          'FSNT:A','FLNT:A','FSNTOA:A','FLUT:A','FSUTOA:A',
          'FSNTC:A','FLNTC:A','FSNTOAC:A','FLUTC:A',
          'SWCF:A','LWCF:A', 
          'FSNS:A','FLNS:A','FSDS:A','FLDS:A',
          'FSNSC:A','FLNSC:A','FSDSC:A','FLDSC:A',  
          'QFLX','U90M','V90M',
          'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANTAU_ISCCP','MEANPTOP_ISCCP','MEANTB_ISCCP',
          'CLDTOT_CAL','CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN','CLDHGH_CAL',
          'CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN','CLDMED_CAL','CLDMED_CAL_LIQ',
          'CLDMED_CAL_ICE','CLDMED_CAL_UN','CLDLOW_CAL','CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE',
          'CLDLOW_CAL_UN'

 ! -- chemUCI settings ------------------
 history_chemdyg_summary = .true.
 history_gaschmbudget_2D = .false.
 history_gaschmbudget_2D_levels = .false.
 history_gaschmbudget_num = 6 !! no impact if  history_gaschmbudget_2D = .false.

 ! -- MAM5 settings ------------------    
 is_output_interactive_volc = .true.       

EOF

cat << EOF >> user_nl_elm
 hist_dov2xy = .true.,.true.,.false.,
 ! All ELM history is 6-hourly averages: each DA cycle is a separate 6-h run,
 ! so longer averages (daily/monthly) never span more than the last cycle.
 ! Daily/monthly means are built in post-processing from these 6-h means.
 hist_empty_htapes    = .true.
 hist_fincl1 = 'SNOWDP','COL_FIRE_CLOSS','NPOOL','PPOOL','TOTPRODC',
               'GPP','NPP','NEE','NEP','NBP','ER','HR','AR',
               'TOTVEGC','TOTSOMC','TOTLITC','TOTECOSYSC','LEAFC','FROOTC',
               'TOTVEGN','SMINN','TOTVEGP','SMINP',
               'FPG','FPG_P','FPI','FPI_P',
               'ELAI','BTRAN','TV','FPSN',
               'FCEV','FCTR','FGEV','FSH_V','FSH_G',
               'QOVER','QDRAI','QINFL','FSAT','H2OCAN','WA',
               'FCH4','NFIRE','FAREA_BURNED'
 hist_fincl2 = 'H2OSNO','FSNO','QRUNOFF','QSNOMELT','FSNO_EFF','SNORDSL','SNOW','FSA','FSDS','FSR',
               'FLDS','FIRE','FIRA','SOILWATER_10CM','SOILLIQ','SOILICE','QSOIL','U10','U10WITHGUSTS',
               'TWS','TSOI_10CM','TSA','TLAI','THBOT','TSOI','TSOI_ICE','SOILLIQ_ICE','SOILICE_ICE',
               'TAUX','TAUY','FSH','HC','HCSOI','EFLX_LH_TOT','SNOW_DEPTH','RH2M','RAIN','QVEGE','QVEGT',
               'QBOT','Q2M','H2OSOI','H2OSFC','ZWT','ZBOT','TBOT','TG','PBOT'

 hist_mfilt           = 4,4
 hist_nhtfrq          = -6,-6
 hist_avgflag_pertape = 'A','A'
 
 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 
 ! Unless using the above finidat for 1850, also set the following, esp. the 2nd one
 finidat="${lnd_init}"
 check_finidat_fsurdat_consistency = .false.
 check_finidat_pct_consistency   = .false.
 flanduse_timeseries = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_hist_simyr1850-2015_c240308.nc'
 fsurdat = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/${lnd_surdat}'
EOF

cat << EOF >> user_nl_mosart 
 finidat_rtm = "${rof_init}" 
 rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
 rtmhist_mfilt  = 4,4
 rtmhist_ndens  = 2
 rtmhist_nhtfrq = -6,-6   ! 6-h averages; DA cycles are separate 6-h runs
EOF

cat << EOF >> user_nl_mpassi
 config_calendar_type = 'gregorian'
 config_initial_condition_type = 'restart'
 config_start_time = '${START_DATE}_${START_HOUR}'
 config_do_restart = .true.
 config_do_restart_bgc = .true.
 config_do_restart_hbrine = .true.
 config_do_restart_snow_density = .true.
 config_do_restart_snow_grain_radius = .true.
 config_do_restart_zsalinity = .true.
 config_restart_timestamp_name = 'rpointer.ice'
EOF

}



# =====================================
# Customize MPAS stream files if needed
# =====================================

patch_mpas_streams() {
echo
echo 'Modifying MPAS streams files'
pushd ${CASE_RUN_DIR}
# change streams.seaice file
patch streams.seaice << EOF
--- streams.seaice    
+++ streams.seaice    
@@ -11,1 +11,1 @@
-                  filename_template="/compyfs/inputdata/ice/mpas-seaice/IcoswISC30E3r5/mpassi.IcoswISC30E3r5.20231120.nc"
+                  filename_template="${ice_init}"
@@ -34,1 +34,8 @@
-
+ 
+<immutable_stream name="restart_ic"
+                  type="input"
+                  io_type="pnetcdf"
+                  filename_template="${ice_init}"
+                  filename_interval="none"
+                  input_interval="initial_only" />
+
EOF

# copy to SourceMods
cp streams.seaice ${CASE_SCRIPTS_DIR}/SourceMods/src.mpassi/

popd

}

# =====================================================
# Custom PE layout: custom-N where N is number of nodes
# =====================================================

custom_pelayout() {

if [[ ${PELAYOUT} == custom-* ]];
then
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'

    # Number of cores per node (machine specific)
    if [ "${MACHINE}" == "chrysalis" ]; then
        ncore=64
    elif [ "${MACHINE}" == "compy" ]; then
        ncore=40
    elif [ "${MACHINE}" == "pm-cpu" ]; then
        ncore=64
    else
        echo 'ERROR: MACHINE = '${MACHINE}' is not supported for custom PE layout.' 
        exit 400
    fi

    # Extract number of nodes
    tmp=($(echo ${PELAYOUT} | tr "-" " "))
    nnodes=${tmp[1]}

    # Customize
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( $nnodes * ${JOB_NTASKS} ))
    ./xmlchange NTHRDS=1
    ./xmlchange MAX_MPITASKS_PER_NODE=$ncore
    ./xmlchange MAX_TASKS_PER_NODE=$ncore
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

    ${CODE_ROOT}/cime/scripts/create_newcase \
        --case ${CASE_NAME} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} ${COMPILER}\
        --project ${PROJECT} \
        --walltime ${WALLTIME} \
        --pecount ${layout}

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

    # Source Mods copy .F90 files to src.*
    # cp ${HOME}/source_files/E3SMv2_UCI-MZT-MSC_20220629/mo_gas_phase_chemdr.F90   ${CASE_SCRIPTS_DIR}/SourceMods/src.eam/mo_gas_phase_chemdr.F90

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

    # Change the calendar 
    ./xmlchange --id CALENDAR  --val 'GREGORIAN'

    # Extracts input_data_dir in case it is needed for user edits to the namelist later
    local input_data_dir=`./xmlquery DIN_LOC_ROOT --value`

    #changing chemistry mechanism
    #local usr_mech_infile="$CODE_ROOT/components/eam/chem_proc/inputs/pp_chemUCI_linozv3_mam5_vbs.in"
    #echo 'Changing chemistry to :'${usr_mech_infile}
    #./xmlchange --id CAM_CONFIG_OPTS --append --val='-microphys p3 -chem superfast_mam5_resus_mom_soag -vbs -usr_mech_infile '${usr_mech_infile}

   # change SST files 
   ./xmlchange SSTICE_DATA_FILENAME=${sst_data}
   ./xmlchange SSTICE_GRID_FILENAME=${sst_grid}
   ./xmlchange SSTICE_YEAR_START=${sst_syear},SSTICE_YEAR_END=${sst_eyear},SSTICE_YEAR_ALIGN=${sst_syear}

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

        # Increase parallelism on compute nodes
        if [ -z "${SLURMD_NODENAME}" ]
        then 
            echo $'\nCompiling on login node: GMAKE_J=8\n' 
            ./xmlchange GMAKE_J=8
        else 
            echo $'\nCompiling on compute node: GMAKE_J=64\n' 
            ./xmlchange GMAKE_J=64
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
    ./xmlchange START_TOD=${START_TOD}

    # Segment length
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}

    # Restart frequency
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}

    # Coupler history
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}

    # Coupler budgets (always on)
    ./xmlchange BUDGETS=TRUE

    # Job Queue 
    ./xmlchange JOB_QUEUE="${JOB_SLURM}"

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
        ./xmlchange RUN_REFTOD=${RUN_REFTOD}
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
