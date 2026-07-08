#====================================================================
# create new case
#====================================================================
rm -rf $CASEROOT

cd  $CCSMROOT/cime/scripts

./create_newcase -case $CASEROOT -mach $MACH -project $PROJECT \
                 -res $RESOLUTION -compset $COMPSET -compiler $COMPILER --walltime $wall_time

#====================================================================
# set up case
#====================================================================
cd $CASEROOT

# SCM must run in serial mode
if ($dycore == Eulerian) then
  ./xmlchange --id MPILIB --val mpi-serial
endif

# Define executable and run directories
./xmlchange --id RUNDIR  --val "${RUNDIR}"
./xmlchange --id EXEROOT --val "${EXEDIR}"

# Set to debug, only on certain MACHs 
if ($MACH =~ 'cori*') then
  ./xmlchange --id JOB_QUEUE --val 'debug'
endif

if ($MACH == 'quartz' || $MACH == 'syrah') then
  ./xmlchange --id JOB_QUEUE --val 'pdebug'
endif

# Get local input data directory path
./xmlchange -file env_run.xml   -id DIN_LOC_ROOT          -val $CSMDATA
./xmlchange -file env_run.xml   -id DIN_LOC_ROOT_CLMFORC  -val $CSMDATA

# need to use single thread
@ nproc = $NTASKS
foreach component ( ATM LND ICE OCN CPL GLC ROF WAV )
  ./xmlchange  NTASKS_$component=$nproc,NTHRDS_$component=$NTHRDS
end

# CAM configure options.  By default set up with settings the same as E3SMv1
set CAM_CONFIG_OPTS="-phys cam5 -scam -nlev $NVLEV -clubb_sgs"
if ( $dycore == Eulerian ) then
  set CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS -nospmd -nosmp"
endif

if ( $do_cosp == true ) then
  set  CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS -cosp -verbose"
endif

# This option ONLY to be used for the REPLAY mode
if ($init_aero_type == none) then
  set CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS -chem linoz_mam4_resus_mom_soag -rain_evap_to_coarse_aero -bc_dep_to_snow_updates"
endif

if ($init_aero_type == cons_droplet || $init_aero_type == prescribed || $init_aero_type == observed) then
  set CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS -chem none"
endif

./xmlchange CAM_CONFIG_OPTS="$CAM_CONFIG_OPTS"

if ($e3sm_version == v2) then
  ./xmlchange CAM_TARGET=theta-l
endif

#-----------------------------
if ($dycore == Eulerian) then
  setenv clubb_micro_steps  8
endif

if ($dycore == SE) then 
  setenv clubb_micro_steps  6
endif

# User enter CAM namelist options
# Add additional output here for example
cat <<EOF >> user_nl_${atm_mod}
 cld_macmic_num_steps = $clubb_micro_steps
 cosp_lite = .true.
 use_gw_front = .true.
 iopfile = '$CSMDATA/$iop_path/$iop_file'
 scm_iop_srf_prop = $do_iop_srf_prop 
 scm_relaxation = $do_scm_relaxation
 iradlw = 1
 iradsw = 1
 swrad_off = $do_turnoff_swrad 
 lwrad_off = $do_turnoff_lwrad
 precip_off = $do_turnoff_precip
 scmlat = $lat 
 scmlon = $lon
EOF

#  Settings shared by v1 and v2
cat <<EOF >> user_nl_${atm_mod}
 use_hetfrz_classnuc = .true.
 micro_mg_dcs_tdep = .true.
 microp_aero_wsub_scheme = 1
 sscav_tuning = .true.
 convproc_do_aer = .true.
 demott_ice_nuc = .true.
 liqcf_fix = .true.
 regen_fix = .true.
 resus_fix = .false.
 mam_amicphys_optaa = 1
 fix_g1_err_ndrop = .true.
 ssalt_tuning = .true.
 relvar_fix = .true.
 mg_prc_coeff_fix = .true.
 rrtmg_temp_fix = .true.
 mam_amicphys_optaa = 1
 fix_g1_err_ndrop = .true.
 ssalt_tuning = .true.
 use_rad_dt_cosz = .true.
 ice_sed_ai = 500.0
 do_tms = .false.
 n_so4_monolayers_pcage = 8.0D0
 se_ftype = 2
EOF

# Set v1 parameters if running EAMv1
if ($e3sm_version == v1) then
cat <<EOF >> user_nl_${atm_mod}
 cldfrc_dp1 = 0.045D0
 clubb_ice_deep = 16.e-6
 clubb_ice_sh = 50.e-6
 clubb_liq_deep = 8.e-6
 clubb_liq_sh = 10.e-6
 clubb_C2rt = 1.75D0
 zmconv_c0_lnd = 0.007
 zmconv_c0_ocn = 0.007
 zmconv_dmpdz = -0.7e-3
 zmconv_ke = 1.5E-6
 effgw_oro = 0.25
 seasalt_emis_scale = 0.85
 dust_emis_fact = 2.05D0
 clubb_gamma_coef = 0.32
 clubb_C8 = 4.3
 cldfrc2m_rhmaxi = 1.05D0
 clubb_c_K10 = 0.3
 effgw_beres = 0.4
 so4_sz_thresh_icenuc = 0.075e-6
 micro_mg_accre_enhan_fac = 1.5D0
 zmconv_tiedke_add = 0.8D0
 zmconv_cape_cin = 1
 zmconv_mx_bot_lyr_adj = 2
 taubgnd = 2.5D-3
 clubb_C1 = 1.335
 raytau0 = 5.0D0
 prc_coef1 = 30500.0D0
 prc_exp = 3.19D0
 prc_exp1 = -1.2D0
 se_ftype = 2
 clubb_C14 = 1.3D0
EOF

else

# V2 tunings.  WARNING: may not be final yet
cat <<EOF >> user_nl_${atm_mod}
 zmconv_trigdcape_ull = .true.
 cld_sed   = 1.0D0
 effgw_beres  = 0.35
 gw_convect_hcf  = 12.5
 effgw_oro  = 0.375
 dust_emis_fact =  1.50D0
 linoz_psc_T = 197.5
 clubb_tk1 = 253.15D0
 clubb_c1               = 2.4
 clubb_c11              = 0.70
 clubb_c11b             = 0.20
 clubb_c11c             = 0.85
 clubb_c14              = 2.5D0
 clubb_c1b              = 2.8
 clubb_c1c              = 0.75
 clubb_c6rtb            = 7.50
 clubb_c6rtc            = 0.50
 clubb_c6thlb           = 7.50
 clubb_c6thlc           = 0.50
 clubb_c8               = 5.2
 clubb_c_k10            = 0.35
 clubb_c_k10h           = 0.35
 clubb_gamma_coef       = 0.12D0
 clubb_gamma_coefb      = 0.28D0
 clubb_gamma_coefc      = 1.2
 clubb_mu               = 0.0005
 clubb_wpxp_l_thresh    = 100.0D0
 clubb_ice_deep         = 14.e-6
 clubb_ice_sh           = 50.e-6
 clubb_liq_deep         = 8.e-6
 clubb_liq_sh           = 10.e-6
 clubb_C2rt             = 1.75D0
 clubb_use_sgv          = .true.
 seasalt_emis_scale     = 0.6
 zmconv_dmpdz  = -0.7e-3
 zmconv_c0_lnd          = 0.0020
 zmconv_c0_ocn          = 0.0020
 zmconv_ke              = 5.0E-6
 zmconv_alfa            = 0.14D0
 zmconv_tp_fac          = 2.0D0
 zmconv_tiedke_add      = 0.8D0
 zmconv_cape_cin        = 1
 zmconv_mx_bot_lyr_adj  = 1
 prc_coef1               = 30500.0D0
 prc_exp                 = 3.19D0
 prc_exp1                = -1.40D0
 micro_mg_accre_enhan_fac = 1.75D0
 microp_aero_wsubmin     = 0.001D0
 so4_sz_thresh_icenuc    = 0.080e-6
 micro_mg_berg_eff_factor = 0.7D0
 cldfrc_dp1              = 0.018D0
 cldfrc2m_rhmaxi         = 1.05D0
 taubgnd                 = 2.5D-3
 raytau0                 = 5.0D0
EOF

endif


# if constant droplet was selected then modify name list to reflect this
if ($init_aero_type == cons_droplet) then

cat <<EOF >> user_nl_${atm_mod}
  micro_do_nccons = .true.
  micro_do_nicons = .true.
  micro_nccons = $micro_nccons_val
  micro_nicons = $micro_nicons_val
EOF

endif

# if prescribed or observed aerosols set then need to put in settings for prescribed aerosol model
if ($init_aero_type == cons_droplet || $init_aero_type == prescribed ||$init_aero_type == observed) then

cat <<EOF >> user_nl_${atm_mod}
  use_hetfrz_classnuc = .false.
  aerodep_flx_type = 'CYCLICAL'
  aerodep_flx_datapath = '$CSMDATA/$presc_aero_path'
  aerodep_flx_file = '$presc_aero_file'
  aerodep_flx_cycle_yr = 01
  prescribed_aero_type = 'CYCLICAL'
  prescribed_aero_datapath='$CSMDATA/$presc_aero_path'
  prescribed_aero_file='$presc_aero_file'
  prescribed_aero_cycle_yr = 01
EOF

endif

# if observed aerosols then set flag
if ($init_aero_type == observed) then

cat <<EOF >> user_nl_${atm_mod}
  scm_observed_aero = .true.
EOF

endif

# avoid the monthly cice file from writing as this 
#   appears to be currently broken for SCM
cat <<EOF >> user_nl_cice
  histfreq='y','x','x','x','x'
EOF

# Set correct version of CLM or ELM (depending on version of code base)
if ($e3sm_version == v1) then
  set CLM_CONFIG_OPTS="-phys clm4_5"
  ./xmlchange CLM_CONFIG_OPTS="$CLM_CONFIG_OPTS"
else
  set ELM_CONFIG_OPTS="-phys elm"
  ./xmlchange ELM_CONFIG_OPTS="$ELM_CONFIG_OPTS"
endif

# Modify the run start and duration parameters for the desired case
./xmlchange RUN_STARTDATE="$startdate",START_TOD="$start_in_sec",STOP_OPTION="$stop_option",STOP_N="$stop_n"

# Modify the latitude and longitude for the particular case
./xmlchange PTS_MODE="TRUE",PTS_LAT="$lat",PTS_LON="$lon"
./xmlchange MASK_GRID="USGS"
./xmlchange CALENDAR="GREGORIAN"

./case.setup

# Don't want to write restarts as this appears to be broken for 
#  CICE model in SCM.  For now set this to a high value to avoid
./xmlchange PIO_TYPENAME="netcdf"
./xmlchange REST_N=300000

# Modify some parameters for CICE to make it SCM compatible 
./xmlchange CICE_AUTO_DECOMP="FALSE"
./xmlchange CICE_DECOMPTYPE="blkrobin"
./xmlchange --id CICE_BLCKX --val 1
./xmlchange --id CICE_BLCKY --val 1
./xmlchange --id CICE_MXBLCKS --val 1
./xmlchange CICE_CONFIG_OPTS="-nodecomp -maxblocks 1 -nx 1 -ny 1"


