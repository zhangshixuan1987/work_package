;;;;;3year(2010-2012) ERA5 reanalysis information;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ERA5_DIR = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/ERA5/2010-2012/regrid_data/"
ERA5_SFC_LEV     = ERA5_DIR + "era5_analysis_land_monthly_2010-2012.nc"
ERA5_SIG_LEV     = ERA5_DIR + "era5_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV     = ERA5_DIR + "era5_analysis_pressure_level_2010-2012.nc"
ERA5_SIG_LEV_ENS = ERA5_DIR + "era5_ensemble_spread_single_level_3hourly_2010-2012.nc"
                  ;"era5_ensemble_member_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_pressure_level_monthly_2010-2012.nc"
if(l_era5_single_level)then
 ERA5_HRES   = ERA5_SIG_LEV
 ERA5_ENS    = ERA5_SIG_LEV_ENS
end if 

if(l_era5_pressure_level)then
 ERA5_HRES   = ERA5_PRE_LEV
 ERA5_ENS    = ERA5_PRE_LEV_ENS
end if

if(l_era5_surface_level)then
 ERA5_HRES   = ERA5_SIG_LEV
 ERA5_ENS    = "" 
end if

;;;;;;3-year monthly data;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ModelOutDir = "/global/cscratch1/sd/zhan391/run_e3sm_crysphere/"
E3SM_NE30        = ModelOutDir + "E3SMV1_NE30_to_180x360.cam.h0.2010-2012.nc" 
E3SM_RRM         = ModelOutDir + "E3SMV1_RRM_to_180x360.cam.h0.2010-2012.nc"
E3SM_NE120       = ModelOutDir + "E3SMV1_NE120_to_180x360.cam.h0.2010-2012.nc"

;;;;Climo data;;;
ERA5_CLIM      = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/ERA5"
CERES_EBAF_TOA = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/ceres_ebaf_toa_v4.1/"
CERES_EBAF_SFC = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/ceres_ebaf_surface_v4.1/"
GPCP_v2p3      = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/GPCP_v2.3"
ISCCPCOSP      = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/ISCCPCOSP"
OAFlux         = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/OAFlux"
HadISST        = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/HadISST_PD"
SSMI           = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/SSMI"
CRU            = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/CRU"
MERRA2         = "/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/MERRA2"

;;;;;;SAM index data;;;
SAM_ERA5       = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/local_diagnostics/sam_index_analysis/normalized_with_30year_ERA5/SAM_Norm_AOI_Monthly_2010-2012_ERA5.txt"
SAM_E3SM_NE30  = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/local_diagnostics/sam_index_analysis/normalized_with_30year_ERA5/SAM_Norm_AOI_Monthly_2010-2012_E3SM_NE30.txt"
SAM_E3SM_RRM   = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/local_diagnostics/sam_index_analysis/normalized_with_30year_ERA5/SAM_Norm_AOI_Monthly_2010-2012_E3SM_RRM.txt"
SAM_E3SM_NE120 = "/global/cfs/cdirs/m3089/shixuan/run_e3sm_crysphere/local_diagnostics/sam_index_analysis/normalized_with_30year_ERA5/SAM_Norm_AOI_Monthly_2010-2012_E3SM_NE120.txt"

;;;;;;E3SM experiments;;;;;;;;;;;
;confg_string     = (/"ERA5","E3SM_NE30", "E3SM_RRM", "E3SM_NE120"/)
confg_string     = (/"ERA5","ne30pg2","Antarctica RRM", "ne120pg2"/)
nexps            = dimsizes (confg_string)
Exp_data_files   = new (nexps, string)
Exp_short_names  = new (nexps, string)
Exp_long_names   = new (nexps, string)
SAM_indx_files   = new (nexps, string)

iexp = 0
Exp_data_files(iexp)   = ERA5_HRES
Exp_short_names(iexp)  = confg_string(iexp)
Exp_long_names(iexp)   = "ERA5 reanalysis"
SAM_indx_files(iexp)   = SAM_ERA5

iexp = iexp+1
Exp_data_files(iexp)   = E3SM_NE30
Exp_short_names(iexp)  = confg_string(iexp)
Exp_long_names(iexp)   = "E3SM (Low  resolution)"
SAM_indx_files(iexp)   = SAM_E3SM_NE30

iexp = iexp +1
Exp_data_files(iexp)   = E3SM_RRM
Exp_short_names(iexp)  = confg_string(iexp)
Exp_long_names(iexp)   = "E3SM (Regional refinement)"
SAM_indx_files(iexp)   = SAM_E3SM_RRM

iexp = iexp +1
Exp_data_files(iexp)   = E3SM_NE120
Exp_short_names(iexp)  = confg_string(iexp)
Exp_long_names(iexp)   = "E3SM (High resolution)"
SAM_indx_files(iexp)   = SAM_E3SM_NE120

