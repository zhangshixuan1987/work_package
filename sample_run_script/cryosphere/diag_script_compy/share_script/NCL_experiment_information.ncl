;;;;;ERA5 reanalysis information;;;;;;
ERA5_DIR = "/compyfs/zhan391/E3SM_PROJECT/ERA5/2010-2012/regrid_data/"
ERA5_SFC_LEV     = ERA5_DIR + "era5_analysis_land_monthly_2010-2012.nc"
ERA5_SIG_LEV     = ERA5_DIR + "era5_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV     = ERA5_DIR + "era5_analysis_pressure_level_2010-2012.nc"
ERA5_SIG_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_pressure_level_monthly_2010-2012.nc"

ModelOutDir = "/compyfs/zhan391/E3SM_PROJECT/E3SM_Wuyin/post_process/"
E3SM_NE120  = ModelOutDir + "E3SMV1_NE30_to_180x360.cam.h0.2010-2012.nc" 
E3SM_NE30   = ModelOutDir + "E3SMV1_NE120_to_180x360.cam.h0.2010-2012.nc"

;;;;;;E3SM experiments;;;;;;;;;;;
confg_string   = (/"ERA5","E3SM_NE30", "E3SM_NE120"/)
nexps          = dimsizes (confg_string)

Exp_short_Names  = new (nexps, string)
Exp_long_names   = new (nexps, string)

Exp_short_Names  = confg_string

Exp_long_names(0)  = "ERA5 (reanalysis)"
Exp_long_names(1)  = "E3SM (Low  resolution)" 
Exp_long_names(2)  = "E3SM (High resolution)" 

