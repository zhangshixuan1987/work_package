#!/bin/csh
date
setenv PROJECT         'ESMD'
setenv CESM_PROJ       'ESMD'
setenv CESM_EMAIL      'shixuan.zhang@pnnl.gov'

setenv do_cosp         'false'
setenv dycore          'Eulerian'
setenv MACH            'compy'
setenv COMPILER        'intel'
setenv RESOLUTION      'T42_T42'
setenv NTASKS_PER_INST '1'
setenv NINST           '1'
setenv NTHRDS         '1'
         
setenv CCSMTAG  'ACME_master_20180611'
setenv COMPSET  'F_SCAM5'
setenv CCSMROOT '/compyfs/zhan391/code/ACME_master_20180611'
setenv PTMP     '/compyfs/zhan391/F_SCAM5_201806/DYCOMS-RF01_RERUN'
setenv EXEDIR   '/compyfs/zhan391/F_SCAM5_201806/DYCOMS-RF01_RERUN/exe//compile_ACME_master_20180611_F_SCAM5_T42_T42_compy_intel_1proc'
setenv CSMDATA  '/compyfs/inputdata/'
         
setenv CASE      'drib_model_DYCOMS-RF01_RERUN_T42_T42_ACME_master_20180611_F_SCAM5_intel_compy_DT1200_01'
setenv CASEROOT  '/compyfs/zhan391/F_SCAM5_201806/DYCOMS-RF01_RERUN/cases/drib_model_DYCOMS-RF01_RERUN_T42_T42_ACME_master_20180611_F_SCAM5_intel_compy_DT1200_01'
setenv RUNDIR    '/compyfs/zhan391/F_SCAM5_201806/DYCOMS-RF01_RERUN/run/drib_model_DYCOMS-RF01_RERUN_T42_T42_ACME_master_20180611_F_SCAM5_intel_compy_DT1200_01'

setenv postCIME 5

set case_setup_script = '/compyfs/zhan391/run_script/e3smv1_nudging/scm_run/DYCOMS_RF01_RERUN/create_and_setup_bundled_case.csh' 

####setup the ouptput frequency###
set dtime    = 1200
set ncycle   = 1                         # 1 cycle in total, i.e., no restart
@ nlen       = 96 * 3600 / 1200   # run model for 60 minutes.
set stop_o   = 'nsteps'
set nhtfrq   = 1
@ mfilt      = 96 * 3600 / 1200 + 1           # 1 time step per file

#====================================================================
# Create and set up new case. No need to build the model. 
#====================================================================
source ${case_setup_script}

cd $CASEROOT

#cp $MODSROOT/*90 SourceMods/src.cam/
#cp $MODSROOT/namelist_definition.xml   $CCSMROOT/models/atm/cam/bld/namelist_files/namelist_definition.xml
##cp $MODSROOT/namelist_defaults_cam.xml $CCSMROOT/models/atm/cam/bld/namelist_files/namelist_defaults_cam.xml

./xmlchange -file env_build.xml -id BUILD_COMPLETE  -val 'TRUE'

#-----------------------------------------
# Runtime options: edit env_run.xml
#-----------------------------------------
cd $CASEROOT

#./xmlchange  -file env_run.xml -id RUN_STARTDATE -val 'startdate_dummy'

@ nresub = $ncycle - 1

./xmlchange  -file env_run.xml -id  STOP_N       -val $nlen
./xmlchange  -file env_run.xml -id  STOP_OPTION  -val $stop_o 
#./xmlchange  -file env_run.xml -id  REST_N       -val $nlen
#./xmlchange  -file env_run.xml -id  REST_OPTION  -val 'ndays' 
#./xmlchange  -file env_run.xml -id  RESUBMIT     -val $nresub 
./xmlchange  -file env_run.xml -id  DOUT_S       -val 'FALSE'
./xmlchange  -file env_run.xml -id  DOUT_L_MS    -val 'FALSE'

@ ncpl = 86400 / $dtime

./xmlchange  -file env_run.xml -id  ATM_NCPL          -val $ncpl
./xmlchange  -file env_run.xml -id  CAM_NAMELIST_OPTS -val dtime=$dtime 
./xmlchange  -file env_run.xml -id  CLM_NAMELIST_OPTS -val dtime=$dtime 

#./xmlchange  -file env_batch.xml -id JOB_WALLCLOCK_TIME -val "6:00:00"
./xmlchange  -file env_run.xml -id  USER_REQUESTED_WALLTIME -val "6:00:00"

#--------------------
# Namelist variables
#--------------------
cat <<EOF >> user_nl_cam
 !eul_nsplit           = 2
 !rsplit               = 3
 !qsplit               = 1
 !inithist           = 'ENDOFRUN'
 !inithist_all       = .false.
 l_dribling_tend      = .true.
 l_dribling_uv        = .false.
 l_dribling_w         = .false.
 history_amwg         = .true.
 docosp               = .false.
 cosp_amwg            = .false.
 history_verbose      = .true.
 history_aero_optics  = .false.
 history_aerosol      = .false.
 history_clubb        = .true.
 history_budget       = .true.
!history_microphysics = .true.
 iradsw               = 2,
 iradlw               = 2,
 inithist             = 'MONTHLY'
 inithist_all         = .true.
 pergro_test_active   = .false.
 l_aerosol_cldgrow    = .true.
 l_aerosol_cldshnk    = .true.
 l_aerosol_oldcld     = .true.
 l_aerosol_mixing     = .true.
 clubb_history        = .true.
 clubb_rad_history    = .true.
 macmic_extra_diag    = .true.
 macmic_clubb_diag    = .true.
 macmic_mg2_diag      = .false.
 empty_htapes       = .false.
 avgflag_pertape    = 'I','I','I'
 nhtfrq             = $nhtfrq,$nhtfrq,$nhtfrq
 mfilt              = $mfilt,$mfilt,$mfilt
 ndens              = 1,1,1
 fincl2               = 'thlm', 'thvm', 'rtm', 'rcm', 'rvm', 'um', 'vm', 'cloud_frac', 'cloud_cover',
                        'rcm_in_layer', 'rcm_in_cloud', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt',
                        'Lscale', 'Lscale_pert_1', 'Lscale_pert_2', 'T_in_K', 'rel_humidity',
                        'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down',
                        'wp2thvp', 'wp2rcp', 'wprtpthlp', 'sigma_sqd_w_zt', 'rho', 'radht', 'radht_LW', 'radht_SW',
                        'rsat', 'rsati', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm',
                        'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_forcing', 'rtm_sdmp','rtm_mc',
                        'rtm_pd', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_bt', 'thlm_ma', 'thlm_ta', 'thlm_mfl',
                        'thlm_tacl', 'thlm_cl', 'thlm_forcing', 'thlm_sdmp','thlm_mc', 'thlm_old', 'thlm_without_ta',
                        'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl',
                        'rtm_old', 'rtm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl',
                        'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'vm_bt', 'vm_ma',
                        'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp',
                        'wp3_ac', 'wp3_bp1', 'wp3_bp2', 'wp3_pr1', 'wp3_pr2', 'wp3_dp1', 'wp3_cl', 'mixt_frac', 'w_1',
                        'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2',
                        'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2',
                        'cloud_frac_1', 'cloud_frac_2', 'a3_coef_zt', 'wp3_on_wp2_zt', 'chi_1', 'chi_2', 'stdev_chi_1',
                        'stdev_chi_2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2',
                        'corr_chi_eta_1', 'corr_chi_eta_2', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'precip_frac',
                        'precip_frac_1', 'precip_frac_2', 'Ncnm',  'C11_Skw_fnc', 'wp2', 'rtp2', 'thlp2', 'rtpthlp',
                        'wprtp', 'wpthlp', 'wp4', 'up2', 'vp2', 'wpthvp', 'rtpthvp', 'thlpthvp', 'wprcp', 'rc_coef',
                        'wm_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'rho_zm', 'sigma_sqd_w', 'Skw_velocity',
                        'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C7_Skw_fnc', 'C1_Skw_fnc', 'a3_coef',
                        'wp3_on_wp2','em', 'mean_w_up', 'mean_w_down','shear','wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac',
                        'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_dp1', 'wp2_dp2', 'wp2_cl', 'wp2_pd', 'wp2_sf',
                        'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_cl',
                        'vp2_pd', 'vp2_sf', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1',
                        'up2_pr2', 'up2_cl', 'up2_pd', 'up2_sf', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp',
                        'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl',
                        'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma',
                        'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3',
                        'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'rtp2_bt',
                        'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf',
                        'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1',
                        'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt',
                        'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2',
                        'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'wpthlp_entermfl',
                        'wpthlp_exit_mfl', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max',
                        'wprtp_mfl_min', 'wprtp_mfl_max'

 fincl3               = 'wp3_01', 'wp3_02', 'wp3_03', 'wp3_04', 'wp3_05', 'wp3_06',
                        'wp2_01', 'wp2_02', 'wp2_03', 'wp2_04', 'wp2_05', 'wp2_06',
                        'rc_1_01', 'rc_1_02', 'rc_1_03', 'rc_1_04', 'rc_1_05', 'rc_1_06',
                        'rc_2_01', 'rc_2_02', 'rc_2_03', 'rc_2_04', 'rc_2_05', 'rc_2_06',
                        'w_1_01',   'w_1_02', 'w_1_03',   'w_1_04',  'w_1_05',  'w_1_06',
                        'w_2_01',   'w_2_02', 'w_2_03',   'w_2_04',  'w_2_05',  'w_2_06',
                        'wp2rcp_01', 'wp2rcp_02', 'wp2rcp_03', 'wp2rcp_04', 'wp2rcp_05', 'wp2rcp_06',
                        'wprcp_01', 'wprcp_02', 'wprcp_03', 'wprcp_04', 'wprcp_05', 'wprcp_06',
                        'wp2rtp_01', 'wp2rtp_02', 'wp2rtp_03', 'wp2rtp_04', 'wp2rtp_05', 'wp2rtp_06',
                        'wprtp_01', 'wprtp_02', 'wprtp_03', 'wprtp_04', 'wprtp_05', 'wprtp_06',
                        'wpthlp_01', 'wpthlp_02', 'wpthlp_03', 'wpthlp_04', 'wpthlp_05', 'wpthlp_06',
                        'wp2thvp_01', 'wp2thvp_02', 'wp2thvp_03', 'wp2thvp_04', 'wp2thvp_05', 'wp2thvp_06',
                        'wpthvp_01', 'wpthvp_02', 'wpthvp_03', 'wpthvp_04', 'wpthvp_05', 'wpthvp_06',
                        'wp3_bt_01', 'wp3_ma_01', 'wp3_ta_01', 'wp3_tp_01', 'wp3_ac_01', 'wp3_bp1_01', 'wp3_bp2_01','wp3_pr1_01', 'wp3_pr2_01', 'wp3_dp1_01', 'wp3_cl_01',
                        'wp3_bt_02', 'wp3_ma_02', 'wp3_ta_02', 'wp3_tp_02', 'wp3_ac_02', 'wp3_bp1_02', 'wp3_bp2_02','wp3_pr1_02', 'wp3_pr2_02', 'wp3_dp1_02', 'wp3_cl_02',
                        'wp3_bt_03', 'wp3_ma_03', 'wp3_ta_03', 'wp3_tp_03', 'wp3_ac_03', 'wp3_bp1_03', 'wp3_bp2_03','wp3_pr1_03', 'wp3_pr2_03', 'wp3_dp1_03', 'wp3_cl_03',
                        'wp3_bt_04', 'wp3_ma_04', 'wp3_ta_04', 'wp3_tp_04', 'wp3_ac_04', 'wp3_bp1_04', 'wp3_bp2_04','wp3_pr1_04', 'wp3_pr2_04', 'wp3_dp1_04', 'wp3_cl_04',
                        'wp3_bt_05', 'wp3_ma_05', 'wp3_ta_05', 'wp3_tp_05', 'wp3_ac_05', 'wp3_bp1_05', 'wp3_bp2_05','wp3_pr1_05', 'wp3_pr2_05', 'wp3_dp1_05', 'wp3_cl_05',
                        'wp3_bt_06', 'wp3_ma_06', 'wp3_ta_06', 'wp3_tp_06', 'wp3_ac_06', 'wp3_bp1_06', 'wp3_bp2_06','wp3_pr1_06', 'wp3_pr2_06', 'wp3_dp1_06', 'wp3_cl_06',
                        'wp2_bt_01', 'wp2_ma_01', 'wp2_ta_01', 'wp2_ac_01', 'wp2_bp_01', 'wp2_pr1_01', 'wp2_pr2_01','wp2_pr3_01', 'wp2_dp1_01', 'wp2_dp2_01', 'wp2_cl_01',
                        'wp2_pd_01', 'wp2_sf_01', 'wp2_bt_02', 'wp2_ma_02', 'wp2_ta_02', 'wp2_ac_02', 'wp2_bp_02', 'wp2_pr1_02', 'wp2_pr2_02','wp2_pr3_02', 'wp2_dp1_02',
                        'wp2_dp2_02', 'wp2_cl_02','wp2_pd_02', 'wp2_sf_02', 'wp2_bt_03', 'wp2_ma_03', 'wp2_ta_03', 'wp2_ac_03', 'wp2_bp_03', 'wp2_pr1_03', 'wp2_pr2_03',
                        'wp2_pr3_03', 'wp2_dp1_03', 'wp2_dp2_03', 'wp2_cl_03','wp2_bt_04', 'wp2_ma_04', 'wp2_ta_04', 'wp2_ac_04', 'wp2_bp_04', 'wp2_pr1_04', 'wp2_pr2_04',
                        'wp2_pr3_04', 'wp2_dp1_04', 'wp2_dp2_04', 'wp2_cl_04','wp2_pd_04', 'wp2_sf_04', 'wp2_bt_05', 'wp2_ma_05', 'wp2_ta_05', 'wp2_ac_05', 'wp2_bp_05',
                        'wp2_pr1_05', 'wp2_pr2_05','wp2_pr3_05', 'wp2_dp1_05', 'wp2_dp2_05', 'wp2_cl_05', 'wp2_bt_06', 'wp2_ma_06', 'wp2_ta_06', 'wp2_ac_06', 'wp2_bp_06',
                        'wp2_pr1_06', 'wp2_pr2_06','wp2_pr3_06', 'wp2_dp1_06', 'wp2_dp2_06', 'wp2_cl_06',

clubb_vars_zt    = 'thlm', 'thvm', 'rtm', 'rcm', 'rvm', 'um', 'vm', 'um_ref','vm_ref','ug', 'vg', 'cloud_frac', 'cloud_cover', 'rcm_in_layer', 'rcm_in_cloud', 'p_in_Pa', 'exner', 'rho_ds_zt', 'thv_ds_zt', 'Lscale', 'Lscale_pert_1', 'Lscale_pert_2', 'T_in_K', 'rel_humidity', 'wp3', 'wpthlp2', 'wp2thlp', 'wprtp2', 'wp2rtp', 'Lscale_up', 'Lscale_down', 'tau_zt', 'Kh_zt', 'wp2thvp', 'wp2rcp', 'wprtpthlp', 'sigma_sqd_w_zt', 'rho', 'radht', 'radht_LW', 'radht_SW', 'Ncm', 'Nc_in_cloud', 'Nc_activated', 'snowslope', 'sed_rcm', 'rsat', 'rsati', 'diam', 'mass_ice_cryst', 'rcm_icedfs', 'u_T_cm', 'rtm_bt', 'rtm_ma', 'rtm_ta', 'rtm_mfl', 'rtm_tacl', 'rtm_cl', 'rtm_forcing', 'rtm_sdmp','rtm_mc', 'rtm_pd', 'rvm_mc', 'rcm_mc', 'rcm_sd_mg_morr', 'thlm_bt', 'thlm_ma', 'thlm_ta', 'thlm_mfl', 'thlm_tacl', 'thlm_cl', 'thlm_forcing', 'thlm_sdmp','thlm_mc', 'thlm_old', 'thlm_without_ta', 'thlm_mfl_min', 'thlm_mfl_max', 'thlm_enter_mfl', 'thlm_exit_mfl', 'rtm_old', 'rtm_without_ta', 'rtm_mfl_min', 'rtm_mfl_max', 'rtm_enter_mfl', 'rtm_exit_mfl', 'um_bt', 'um_ma', 'um_gf', 'um_cf', 'um_ta', 'um_f', 'um_sdmp', 'um_ndg', 'vm_bt', 'vm_ma', 'vm_gf', 'vm_cf', 'vm_ta', 'vm_f', 'vm_sdmp', 'vm_ndg', 'wp3_bt', 'wp3_ma', 'wp3_ta', 'wp3_tp', 'wp3_ac', 'wp3_bp1', 'wp3_bp2', 'wp3_pr1', 'wp3_pr2', 'wp3_dp1', 'wp3_cl', 'mixt_frac', 'w_1', 'w_2', 'varnce_w_1', 'varnce_w_2', 'thl_1', 'thl_2', 'varnce_thl_1', 'varnce_thl_2', 'rt_1', 'rt_2', 'varnce_rt_1', 'varnce_rt_2', 'rc_1', 'rc_2', 'rsatl_1', 'rsatl_2', 'cloud_frac_1', 'cloud_frac_2', 'a3_coef_zt', 'wp3_on_wp2_zt', 'chi_1', 'chi_2', 'stdev_chi_1', 'stdev_chi_2', 'stdev_eta_1', 'stdev_eta_2', 'covar_chi_eta_1', 'covar_chi_eta_2', 'corr_chi_eta_1', 'corr_chi_eta_2', 'crt_1', 'crt_2', 'cthl_1', 'cthl_2', 'precip_frac', 'precip_frac_1', 'precip_frac_2', 'Ncnm', 'wp2_zt', 'thlp2_zt', 'wpthlp_zt', 'wprtp_zt', 'rtp2_zt', 'rtpthlp_zt', 'up2_zt', 'vp2_zt', 'upwp_zt', 'vpwp_zt', 'C11_Skw_fnc'
clubb_vars_zm    = 'wp2', 'rtp2', 'thlp2', 'rtpthlp', 'wprtp', 'wpthlp', 'wp4', 'up2', 'vp2', 'wpthvp', 'rtpthvp', 'thlpthvp', 'tau_zm', 'Kh_zm', 'wprcp', 'rc_coef', 'wm_zm', 'thlprcp', 'rtprcp', 'rcp2', 'upwp', 'vpwp', 'rho_zm', 'sigma_sqd_w', 'Skw_velocity', 'gamma_Skw_fnc', 'C6rt_Skw_fnc', 'C6thl_Skw_fnc', 'C7_Skw_fnc', 'C1_Skw_fnc', 'a3_coef', 'wp3_on_wp2', 'rcm_zm', 'rtm_zm', 'thlm_zm', 'cloud_frac_zm', 'rho_ds_zm', 'thv_ds_zm', 'em', 'mean_w_up', 'mean_w_down', 'shear', 'wp3_zm', 'Frad', 'Frad_LW', 'Frad_SW', 'Frad_LW_up', 'Frad_SW_up', 'Frad_LW_down', 'Frad_SW_down', 'Fprec', 'Fcsed', 'wp2_bt', 'wp2_ma', 'wp2_ta', 'wp2_ac', 'wp2_bp', 'wp2_pr1', 'wp2_pr2', 'wp2_pr3', 'wp2_dp1', 'wp2_dp2', 'wp2_cl', 'wp2_pd', 'wp2_sf', 'vp2_bt', 'vp2_ma', 'vp2_ta', 'vp2_tp', 'vp2_dp1', 'vp2_dp2', 'vp2_pr1', 'vp2_pr2', 'vp2_cl', 'vp2_pd', 'vp2_sf', 'up2_bt', 'up2_ma', 'up2_ta', 'up2_tp', 'up2_dp1', 'up2_dp2', 'up2_pr1', 'up2_pr2', 'up2_cl', 'up2_pd', 'up2_sf', 'wprtp_bt', 'wprtp_ma', 'wprtp_ta', 'wprtp_tp', 'wprtp_ac', 'wprtp_bp', 'wprtp_pr1', 'wprtp_pr2', 'wprtp_pr3', 'wprtp_dp1', 'wprtp_mfl', 'wprtp_cl', 'wprtp_sicl', 'wprtp_pd', 'wprtp_forcing', 'wprtp_mc', 'wpthlp_bt', 'wpthlp_ma', 'wpthlp_ta', 'wpthlp_tp', 'wpthlp_ac', 'wpthlp_bp', 'wpthlp_pr1', 'wpthlp_pr2', 'wpthlp_pr3', 'wpthlp_dp1', 'wpthlp_mfl', 'wpthlp_cl', 'wpthlp_sicl', 'wpthlp_forcing', 'wpthlp_mc', 'rtp2_bt', 'rtp2_ma', 'rtp2_ta', 'rtp2_tp', 'rtp2_dp1', 'rtp2_dp2', 'rtp2_cl', 'rtp2_pd', 'rtp2_sf', 'rtp2_forcing', 'rtp2_mc', 'thlp2_bt', 'thlp2_ma', 'thlp2_ta', 'thlp2_tp', 'thlp2_dp1', 'thlp2_dp2', 'thlp2_cl', 'thlp2_pd', 'thlp2_sf', 'thlp2_forcing', 'thlp2_mc', 'rtpthlp_bt', 'rtpthlp_ma', 'rtpthlp_ta', 'rtpthlp_tp1', 'rtpthlp_tp2', 'rtpthlp_dp1', 'rtpthlp_dp2', 'rtpthlp_cl', 'rtpthlp_sf', 'rtpthlp_forcing', 'rtpthlp_mc', 'wpthlp_entermfl', 'wpthlp_exit_mfl', 'wprtp_enter_mfl', 'wprtp_exit_mfl', 'wpthlp_mfl_min', 'wpthlp_mfl_max', 'wprtp_mfl_min', 'wprtp_mfl_max'
!tstep_type = 5, 
!qsplit     = 1, 
!rsplit     = 3, 
!se_nsplit  = 2,
!hypervis_subcycle = 3
EOF

#============
# Run model 
#============
cd $CASEROOT

# Submit case to queue if set, else submit
#   via the case.run script
  if ($submit_to_queue == TRUE) then
    ./case.submit
  else
    ./case.submit --no-batch
  endif

date
