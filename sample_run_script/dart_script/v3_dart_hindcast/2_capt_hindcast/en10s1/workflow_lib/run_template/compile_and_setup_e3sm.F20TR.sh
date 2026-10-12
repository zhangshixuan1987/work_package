#!/bin/bash -fe
#===============================================================================
# E3SM case template for the hindcast workflow, my_compset=F20TR.
#
# Step 1 (1_run_e3sm_ensemble_setup.sh) copies this file to
# runtmp/run_scripts/compile_and_setup_e3sm.sh and runs it once for the base
# case (compiles e3sm.exe) and once per member (reuses that executable).
# All settings arrive in HC_* environment variables from
# create_and_setup_case.sh; nothing in this file is rewritten by sed.
#
# Edit here: user_nl() history output / namelists, prescribed SST files.
# Member initial conditions (ncdata, finidat, finidat_rtm, MPAS restart
# streams, rpointers, hybrid start) are applied and verified per member by
# workflow_lib/common.sh:configure_member, so they are not set here.
#===============================================================================
main() {

# --- Configuration (from the workflow driver) ----
readonly MACHINE=${HC_MACHINE:?}
readonly PROJECT=${HC_PROJECT:?}
readonly JOB_QUEUE=${HC_JOB_QUEUE:?}
readonly WALLTIME=${HC_WALLTIME:?}
readonly NODES=${HC_NODES:?}
readonly TASK_PER_NODE=${HC_TASK_PER_NODE:?}

readonly COMPSET=${HC_COMPSET:?}
readonly RESOLUTION=${HC_RESOLUTION:?}
readonly CASE_NAME=${HC_CASE_NAME:?}
readonly DEBUG_COMPILE=false

readonly START_DATE=${HC_START_DATE:?}
readonly START_TOD=${HC_START_TOD:?}
readonly START_HOUR=$(( 10#${START_TOD} / 3600 ))
readonly STOP_OPTION=${HC_STOP_OPTION:?}
readonly STOP_N=${HC_STOP_N:?}
readonly REST_OPTION=${HC_REST_OPTION:?}
readonly REST_N=${HC_REST_N:?}
readonly DO_SHORT_TERM_ARCHIVING=true

# Coupler history
readonly HIST_OPTION="nhours"
readonly HIST_N="6"

# Initial-condition files written into user_nl (reset per member later)
readonly atm_init=${HC_ATM_INIT:?}
readonly lnd_init=${HC_LND_INIT:?}
readonly rof_init=${HC_ROF_INIT:?}
readonly OCN_DATE=${HC_OCN_DATE:?}

# Prescribed SST / sea ice (AMIP)
readonly sst_data="/compyfs/zhan391/acme_init/SST_forcing/sst_ice_NOAA_AVHRR_E3SM_1x1_c20231225.nc"
readonly sst_grid="/compyfs/zhan391/acme_init/SST_forcing/domain.ocn.1x1.111007.nc"
readonly sst_syear="1981"
readonly sst_eyear="2020"

# Paths
readonly CODE_ROOT=${HC_CODE_ROOT:?}
readonly CASE_ROOT=${HC_CASE_ROOT:?}
readonly CASE_SCRIPTS_DIR=${HC_CASE_SCRIPTS_DIR:?}
readonly CASE_BUILD_DIR=${HC_CASE_BUILD_DIR:?}
readonly CASE_RUN_DIR=${HC_CASE_RUN_DIR:?}
readonly CASE_ARCHIVE_DIR=${HC_CASE_ARCHIVE_DIR:?}

# Empty: compile e3sm.exe (base case). Otherwise: reuse this executable (members).
readonly OLD_EXECUTABLE=${HC_OLD_EXECUTABLE:-}

# --- Now, do the work ---
umask 022
create_newcase
custom_pelayout
case_setup
case_build
runtime_options
copy_script
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

 avgflag_pertape = 'A','A','I'
 nhtfrq = 0,-24,-6
 mfilt  = 1,1,4

 fincl1 = 'AODALL','AODBC','AODDUST','AODPOM','AODSO4','AODSOA','AODSS','AODVIS',
          'CLDLOW','CLDMED','CLDHGH','CLDTOT',
          'CLDHGH_CAL','CLDLOW_CAL','CLDMED_CAL','CLD_MISR','CLDTOT_CAL',
          'CLMODIS','FISCCP1_COSP','FLDS','FLNS','FLNSC','FLNT','FLUT',
          'FLUTC','FSDS','FSDSC','FSNS','FSNSC','FSNT','FSNTOA','FSNTOAC','FSNTC',
          'ICEFRAC','LANDFRAC','LWCF','OCNFRAC','OMEGA','PRECC','PRECL','PRECSC','PRECSL','PS','PSL','Q',
          'QFLX','QREFHT','RELHUM','SCO','SHFLX','SOLIN','SWCF','T','TAUX','TAUY','TCO',
          'TGCLDLWP','TMQ','TREFHT','TREFMNAV','TREFMXAV','TS','U','U10','V','Z3',
          'dst_a1DDF','dst_a3DDF','dst_c1DDF','dst_c3DDF','dst_a1SFWET','dst_a3SFWET','dst_c1SFWET','dst_c3SFWET',
          'O3','LHFLX',
          'O3_2DTDA_trop','O3_2DTDB_trop','O3_2DTDD_trop','O3_2DTDE_trop','O3_2DTDI_trop','O3_2DTDL_trop',
          'O3_2DTDN_trop','O3_2DTDO_trop','O3_2DTDS_trop','O3_2DTDU_trop','O3_2DTRE_trop','O3_2DTRI_trop',
          'O3_SRF','NO_2DTDS','NO_TDLgt','NO2_2DTDD','NO2_2DTDS','NO2_TDAcf','CO_SRF','TROPE3D_P','TROP_P',
          'CDNUMC','SFDMS','so4_a1_sfgaex1','so4_a2_sfgaex1','so4_a3_sfgaex1','so4_a5_sfgaex1','soa_a1_sfgaex1',
          'soa_a2_sfgaex1','soa_a3_sfgaex1','GS_soa_a1','GS_soa_a2','GS_soa_a3','AQSO4_H2O2','AQSO4_O3',
          'SFSO2','SO2_CLXF','SO2','DF_SO2','AQ_SO2','GS_SO2','WD_SO2','ABURDENSO4_STR','ABURDENSO4_TRO',
          'ABURDENSO4','ABURDENBC','ABURDENDUST','ABURDENMOM','ABURDENPOM','ABURDENSEASALT',
          'ABURDENSOA','AODSO4_STR','AODSO4_TRO',
          'EXTINCT','AODABS','AODABSBC','CLDICE','CLDLIQ','CLD_CAL_TMPLIQ','CLD_CAL_TMPICE','Mass_bc_srf',
          'Mass_dst_srf','Mass_mom_srf','Mass_ncl_srf','Mass_pom_srf','Mass_so4_srf','Mass_soa_srf','Mass_bc_850',
          'Mass_dst_850','Mass_mom_850','Mass_ncl_850','Mass_pom_850','Mass_so4_850','Mass_soa_850','Mass_bc_500',
          'Mass_dst_500','Mass_mom_500','Mass_ncl_500','Mass_pom_500','Mass_so4_500','Mass_soa_500','Mass_bc_330',
          'Mass_dst_330','Mass_mom_330','Mass_ncl_330','Mass_pom_330','Mass_so4_330','Mass_soa_330','Mass_bc_200',
          'Mass_dst_200','Mass_mom_200','Mass_ncl_200','Mass_pom_200','Mass_so4_200','Mass_soa_200',
          'O3_2DTDD','O3_2DCIP','O3_2DCIL','CO_2DTDS','CO_2DTDD','CO_2DCEP','CO_2DCEL','NO_2DTDD',
          'FLNTC','SAODVIS',
          'H2OLNZ',
          'dst_a1SF','dst_a3SF',
          'PHIS','CLOUD','TGCLDIWP','TGCLDCWP','AREL',
          'CLDTOT_ISCCP','MEANCLDALB_ISCCP','MEANPTOP_ISCCP','CLD_CAL',
          'CLDTOT_CAL_LIQ','CLDTOT_CAL_ICE','CLDTOT_CAL_UN',
          'CLDHGH_CAL_LIQ','CLDHGH_CAL_ICE','CLDHGH_CAL_UN',
          'CLDMED_CAL_LIQ','CLDMED_CAL_ICE','CLDMED_CAL_UN',
          'CLDLOW_CAL_LIQ','CLDLOW_CAL_ICE','CLDLOW_CAL_UN',
          'CLWMODIS','CLIMODIS'

 !daily (h1) I
 fincl2 = 'TCO:A','SCO:A','FLUT:A','PRECT:A',
          'PS:I','PSL:I','TS:I','TBOT:I','TUQ:I','TVQ:I','UBOT:I','VBOT:I','ZBOT:I','U10:I',
          'TAUX:I','TAUY:I','TREFHT:I','QREFHT:I','TMQ:I','CAPE:I','CIN:I','TTQ:I',
          'U1000:I','U975:I','U950:I','U925:I','U900:I','U850:I','U800:I','U700:I','U600:I',
          'U500:I','U400:I','U300:I','U200:I','U100:I','U050:I','U010:I',
          'V1000:I','V975:I','V950:I','V925:I','V900:I','V850:I','V800:I','V700:I','V600:I',
          'V500:I','V400:I','V300:I','V200:I','V100:I','V050:I','V010:I',
          'Z1000:I','Z975:I','Z950:I','Z925:I','Z900:I','Z850:I','Z800:I','Z700:I','Z600:I',
          'Z500:I','Z400:I','Z300:I','Z200:I','Z100:I','Z050:I','Z010:I'
          'T1000:I','T975:I','T950:I','T925:I','T900:I','T850:I','T800:I','T700:I','T600:I',
          'T500:I','T400:I','T300:I','T200:I','T100:I','T050:I','T010:I',
          'Q1000:I','Q975:I','Q950:I','Q925:I','Q900:I','Q850:I','Q800:I','Q700:I','Q600:I',
          'Q500:I','Q400:I','Q300:I','Q200:I','Q100:I','Q050:I','Q010:I',
          'OMEGA1000:I','OMEGA975:I','OMEGA950:I','OMEGA925:I','OMEGA900:I','OMEGA850:I',
          'OMEGA800:I','OMEGA700:I','OMEGA600:I','OMEGA500:I','OMEGA400:I',
          'OMEGA300:I','OMEGA200:I','OMEGA100:I',
          'TREFHTMN:M','TREFHTMX:X','PRECTMX:X',

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
 hist_fincl2 = 'H2OSNO', 'FSNO', 'QRUNOFF', 'QSNOMELT', 'FSNO_EFF', 'SNORDSL', 'SNOW', 'FSA', 'FSDS', 'FSR', 'FLDS', 'FIRE', 'FIRA', 'SOILWATER_10CM','SOILLIQ','SOILICE','QSOIL','U10','U10WITHGUSTS','TWS','TSOI_10CM','TSA','TLAI','THBOT','TSOI','TSOI_ICE','SOILLIQ_ICE','SOILICE_ICE','TAUX','TAUY','FSH','HC','HCSOI','EFLX_LH_TOT','SNOW_DEPTH','RH2M','RAIN','QVEGE','QVEGT','QBOT','Q2M','H2OSOI','H2OSFC','ZWT','ZBOT','TBOT','TG','THBOT','PBOT'

 hist_mfilt = 1,365
 hist_nhtfrq = 0,-24
 hist_avgflag_pertape = 'A','A'

 check_finidat_year_consistency = .false.
 check_dynpft_consistency = .false.
 ! Unless using the above finidat for 1850, also set the following, esp. the 2nd one
 finidat="${lnd_init}"
 check_finidat_fsurdat_consistency = .false.
 check_finidat_pct_consistency   = .false.
 flanduse_timeseries = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/landuse.timeseries_0.5x0.5_hist_simyr1850-2015_c240308.nc'
 fsurdat = '\${DIN_LOC_ROOT}/lnd/clm2/surfdata_map/surfdata_0.5x0.5_simyr2010_c230922_with_TOP.nc' 
EOF

cat << EOF >> user_nl_mosart 
 finidat_rtm = "${rof_init}" 
 rtmhist_fincl2 = 'RIVER_DISCHARGE_OVER_LAND_LIQ'
 rtmhist_mfilt  = 1,1
 rtmhist_ndens  = 2
 rtmhist_nhtfrq = 0,-24
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

# =====================================================
# Custom PE layout: NODES x TASK_PER_NODE MPI tasks
# =====================================================
custom_pelayout() {
    echo $'\n CUSTOMIZE PROCESSOR CONFIGURATION:'
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange NTASKS=$(( NODES * TASK_PER_NODE ))
    ./xmlchange NTHRDS=1
    ./xmlchange MAX_MPITASKS_PER_NODE=${TASK_PER_NODE}
    ./xmlchange MAX_TASKS_PER_NODE=${TASK_PER_NODE}
    popd
}

#-----------------------------------------------------
create_newcase() {
    echo $'\n----- Starting create_newcase -----\n'
    ${CODE_ROOT}/cime/scripts/create_newcase \
        --case ${CASE_NAME} \
        --output-root ${CASE_ROOT} \
        --script-root ${CASE_SCRIPTS_DIR} \
        --handle-preexisting-dirs u \
        --compset ${COMPSET} \
        --res ${RESOLUTION} \
        --machine ${MACHINE} \
        --project ${PROJECT} \
        --walltime ${WALLTIME} \
        --pecount M
}

#-----------------------------------------------------
case_setup() {
    echo $'\n----- Starting case_setup -----\n'
    pushd ${CASE_SCRIPTS_DIR}

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

    ./xmlchange --id CALENDAR --val 'GREGORIAN'

    # Prescribed SST / sea ice (AMIP)
    ./xmlchange SSTICE_DATA_FILENAME=${sst_data}
    ./xmlchange SSTICE_GRID_FILENAME=${sst_grid}
    ./xmlchange SSTICE_YEAR_START=${sst_syear},SSTICE_YEAR_END=${sst_eyear},SSTICE_YEAR_ALIGN=${sst_syear}

    # Custom user_nl
    user_nl

    ./case.setup --reset
    popd
}

#-----------------------------------------------------
case_build() {
    pushd ${CASE_SCRIPTS_DIR}
    if [ -n "${OLD_EXECUTABLE}" ]; then
        echo $'\n----- Reusing executable -----\n'
        if [ ! -x "${OLD_EXECUTABLE}" ]; then
            echo "ERROR: OLD_EXECUTABLE=${OLD_EXECUTABLE} does not exist or is not executable."
            exit 297
        fi
        mkdir -p ${CASE_BUILD_DIR}
        cp -fp ${OLD_EXECUTABLE} ${CASE_BUILD_DIR}/
        ./xmlchange BUILD_COMPLETE=TRUE
    else
        echo $'\n----- Starting case_build -----\n'
        if [ "${DEBUG_COMPILE^^}" == "TRUE" ]; then
            ./xmlchange DEBUG=${DEBUG_COMPILE^^}
        fi
        if [ -z "${SLURMD_NODENAME}" ]; then
            echo $'\nCompiling on login node: GMAKE_J=8\n'
            ./xmlchange GMAKE_J=8
        else
            echo $'\nCompiling on compute node: GMAKE_J=64\n'
            ./xmlchange GMAKE_J=64
        fi
        ./case.build
    fi
    # Call preview_namelists to make sure *_in and user_nl files are consistent.
    echo $'\n----- Preview namelists -----\n'
    ./preview_namelists
    popd
}

#-----------------------------------------------------
runtime_options() {
    echo $'\n----- Starting runtime_options -----\n'
    pushd ${CASE_SCRIPTS_DIR}
    ./xmlchange RUN_STARTDATE=${START_DATE}
    ./xmlchange START_TOD=${START_TOD}
    ./xmlchange STOP_OPTION=${STOP_OPTION,,},STOP_N=${STOP_N}
    ./xmlchange REST_OPTION=${REST_OPTION,,},REST_N=${REST_N}
    ./xmlchange HIST_OPTION=${HIST_OPTION,,},HIST_N=${HIST_N}
    ./xmlchange BUDGETS=TRUE
    ./xmlchange JOB_QUEUE="${JOB_QUEUE}"
    popd
}

#-----------------------------------------------------
copy_script() {
    echo $'\n----- Saving run script for provenance -----\n'
    local script_provenance_dir=${CASE_SCRIPTS_DIR}/run_script_provenance
    mkdir -p ${script_provenance_dir}
    cp -vp "$0" ${script_provenance_dir}/$(basename "$0").$(date +%Y%m%d-%H%M%S)
}

#-----------------------------------------------------
# Silent versions of popd and pushd
pushd() {
    command pushd "$@" > /dev/null
}
popd() {
    command popd "$@" > /dev/null
}

# Run the script
main
