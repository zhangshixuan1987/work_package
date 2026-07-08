#!/bin/bash -fe
main() {

# For debugging, uncomment libe below
#set -x

# --- Configuration flags ----

# Machine and project
readonly MACHINE="compy"
readonly PROJECT="esmd"
readonly JOB_SLURM="short"
readonly JOB_NTASKS="40"
#readonly COMPILER="--compiler intel" # spec must include --compiler, otherwise empty it or comment out the whole line

# Simulation
readonly DATE=`date +"%m-%d-%y"`
readonly COMPSET="F20TR"
readonly RESOLUTION="ne30pg2_r05_IcoswISC30E3r5"
readonly CASE_NAME="DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy.EN10"

# Code and compilation
readonly CHECKOUT="20230715" # 
readonly BRANCH=""           #  4361928bd23968e00e3d6c04212a0ceec21e6bf6
readonly CHERRY=( )
readonly DEBUG_COMPILE=false

# Run options
#readonly MODEL_START_TYPE="hybrid"  # 'initial', 'continue', 'branch', 'hybrid'
readonly START_DATE="2011-12-01"
readonly START_TOD="00000"

readonly ATM_TIME=${RUN_REFDATE}-${RUN_REFTOD}
readonly OCN_TIME=${RUN_REFDATE}_${RUN_REFTOD}
readonly OCN_DATE=${START_DATE}_`printf "%02d" ${START_HOUR}`:00:00

# Additional options for 'branch' and 'hybrid'
readonly GET_REFCASE=TRUE
readonly RUN_REFDATE="2011-12-01"
readonly RUN_REFTOD="00000"
readonly RUN_REFCASE="DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy.EN10"
readonly RUN_REFDIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/rest/2011-12-01-00000"

readonly MODEL_START_TYPE="branch"
readonly rof_init="${RUN_REFDIR}/${RUN_REFCASE}.mosart.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
readonly atm_init="${RUN_REFDIR}/${RUN_REFCASE}.eam.i.${RUN_REFDATE}-${RUN_REFTOD}.nc"
readonly ice_init="${RUN_REFDIR}/${RUN_REFCASE}.mpassi.rst.${RUN_REFDATE}_${RUN_REFTOD}.nc"
readonly lnd_init="${RUN_REFDIR}/${RUN_REFCASE}.elm.r.${RUN_REFDATE}-${RUN_REFTOD}.nc"
echo $atm_init 
echo $lnd_init
echo $rof_init
echo $ice_init

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
readonly CODE_ROOT="/qfs/people/zhan391/e3sm_dart_work/code/E3SMv3"
readonly CASE_ROOT="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy"
readonly CASE_BUILD_DIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/EN10/build"
readonly CASE_RUN_DIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/EN10/run"
readonly CASE_ARCHIVE_DIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive"
readonly CASE_SCRIPTS_DIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/EN10/case_scripts"

readonly run="custom-2_1x6_nhours"

if [ "${run}" != "production" ]; then
  # Short test simulations
  tmp=($(echo $run | tr "_" " "))
  layout=${tmp[0]}
  units=${tmp[2]}
  resubmit=$(( ${tmp[1]%%x*} -1 ))
  length=${tmp[1]##*x}
  readonly PELAYOUT=${layout}
  readonly WALLTIME="02:00:00"
  readonly STOP_OPTION=${units}
  readonly STOP_N=${length}
  readonly REST_OPTION=${STOP_OPTION}
  readonly REST_N=${STOP_N}
  readonly RESUBMIT=${resubmit}
  readonly DO_SHORT_TERM_ARCHIVING=true
else
  # Production simulation
  readonly PELAYOUT="L"
  readonly WALLTIME="02:00:00"
  readonly STOP_OPTION="nyears"
  readonly STOP_N="1"
  readonly REST_OPTION="nyears"
  readonly REST_N="1"
  readonly RESUBMIT="10"
  readonly DO_SHORT_TERM_ARCHIVING=false
fi

# Coupler history 
readonly HIST_OPTION="nhours" 
readonly HIST_N="6"

# Leave empty (unless you understand what it does)
readonly OLD_EXECUTABLE="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/build/e3sm.exe"

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

 avgflag_pertape = 'A','A','I'
 nhtfrq = 0,-24,-6
 mfilt  = 1,1,1

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
    ./xmlchange START_TOD="00000"

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
        ./xmlchange GET_REFCASE=TRUE
        ./xmlchange RUN_REFDIR="/compyfs/zhan391/v3_dart_cda_scratch/DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy/archive/rest/2011-12-01-00000"
        ./xmlchange RUN_REFCASE="DARTEN10_F20TR_ne30pg2_r05_IcoswISC30E3r5_compy.EN10"
        ./xmlchange RUN_REFDATE="2011-12-01"
        ./xmlchange RUN_REFTOD="00000"
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
