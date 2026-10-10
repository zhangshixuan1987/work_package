;;;;;ERA5 reanalysis information;;;;;;
ERA5_DIR = "/compyfs/zhan391/E3SM_PROJECT/ERA5/2010-2012/regrid_data/"
ERA5_SFC_LEV     = ERA5_DIR + "era5_analysis_land_monthly_2010-2012.nc"
ERA5_SIG_LEV     = ERA5_DIR + "era5_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV     = ERA5_DIR + "era5_analysis_pressure_level_2010-2012.nc"
ERA5_SIG_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_pressure_level_monthly_2010-2012.nc"

;;;;;;E3SM experiments;;;;;;;;;;;
confg_string = (/"ENSNDG_ERA5", "EDAUVT"/)
nconfg       = dimsizes (confg_string)

hincast_string = (/"Hind_d0","Hind_d1","Hind_d2","Hind_d5"/)
nhind = dimsizes(hincast_string)

nexps = 1 * nhind + 1
OUT_Exp_Names = new (nexps, string)
Exp_long_names= new (nexps, string)

OUT_Exp_Names(0)  = confg_string(0) + "_"+ hincast_string(0)
OUT_Exp_Names(1:) = confg_string(1) + "_"+ hincast_string(:)

Exp_long_names(0) = confg_string(0)
Exp_long_names(1:)= confg_string(1) + "_"+ hincast_string(:)  

