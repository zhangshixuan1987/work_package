;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "../share_script/NCL_experiment_information.ncl"
load "../share_script/Table_Var_info.ncl"
begin

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Regions      = (/"Polar (SH)"/)
  LatN         = -55.0
  LatS         = -90.0
  nregs        = dimsizes(Regions)

;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups       = Exp_short_Names
  ngrps        = dimsizes(Groups)

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Varname      = VarList
  nvars        = dimsizes(Varname)

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
ERA5_DIR = "/compyfs/zhan391/E3SM_PROJECT/ERA5/2010-2012/regrid_data/"
ERA5_SFC_LEV     = ERA5_DIR + "era5_analysis_land_monthly_2010-2012.nc"
ERA5_SIG_LEV     = ERA5_DIR + "era5_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV     = ERA5_DIR + "era5_analysis_pressure_level_2010-2012.nc"
ERA5_SIG_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_sig_level_monthly_2010-2012.nc"
ERA5_PRE_LEV_ENS = ERA5_DIR + "era5_ensemble_member_analysis_pressure_level_monthly_2010-2012.nc"

ModelOutDir = "/compyfs/zhan391/E3SM_PROJECT/E3SM_Wuyin/post_process/"
E3SM_NE120  = ModelOutDir + "E3SMV1_NE30_to_180x360.cam.h0.2010-2012.nc"
E3SM_NE30   = ModelOutDir + "E3SMV1_NE120_to_180x360.cam.h0.2010-2012.nc"


 
 FERA5_SFC_ENM  = addfile(ERA5_SFC_LEV,"r")
 FERA5_SIG_ENM  = addfile(ERA5_SIG_LEV,"r")
 FERA5_PRE_ENM  = addfile(ERA5_PRE_LEV,"r")
 FERA5_SIG_E10  = addfile(ERA5_SIG_LEV_ENS,"r")
 FERA5_PRE_E10  = addfile(ERA5_PRE_LEV_ENS,"r")

 FE3SM_NE30     = addfile(E3SM_NE30,"r")
 FE3SM_NE120    = addfile(E3SM_NE120,"r")

 fname = "global_mean.txt"
 system("/bin/rm -f "+fname)   ; remove any pre-existing file

 print(nvars)

  do iv = 0,nvars-1

   varnam    = Varname(iv)
   print(varnam)
   gw        = File0->gw
   lat       = File0->lat
   lon       = File0->lon
   var1      = File0->$varnam$(:,:,:)
   var2      = File1->$varnam$(:,:,:)
   var3      = File2->$varnam$(:,:,:)

   avar1     = dim_avg_n_Wrap(var1,0)
   avar2     = dim_avg_n_Wrap(var2,0)
   avar3     = dim_avg_n_Wrap(var3,0)
  
   if(varnam.eq."OMEGA500")then
    avar1    = avar1 *3600.0*24.0/100.0 
    avar2    = avar2 *3600.0*24.0/100.0
    avar3    = avar3 *3600.0*24.0/100.0
   end if 

   if(varnam.eq."PRECC".or.varnam.eq."PRECL")then
    avar1    = avar1 *3600.0*24.0 *1000.0
    avar2    = avar2 *3600.0*24.0 *1000.0
    avar3    = avar3 *3600.0*24.0 *1000.0
   end if

   if(varnam.eq."QREFHT") then
    avar1    = avar1 *1000.0
    avar2    = avar2 *1000.0
    avar3    = avar3 *1000.0
   end if 

   if(varnam.eq."TGCLDCWP".or.varnam.eq."TGCLDLWP".or.varnam.eq."TGCLDIWP") then
    avar1    = avar1 *1000.0
    avar2    = avar2 *1000.0
    avar3    = avar3 *1000.0
   end if

   glbmean1  = wgt_areaave_Wrap(avar1,  gw, 1.0, 1) 
   glbmean2  = wgt_areaave_Wrap(avar2,  gw, 1.0, 1)
   glbmean3  = wgt_areaave_Wrap(avar3,  gw, 1.0, 1)

   alist = [/varnam,glbmean1,glbmean2,glbmean3/]
  
   if(iv.eq.0)then
    write_table(fname, "w", alist, "%s,%16.2f,%16.2f,%16.2f")
   else
    write_table(fname, "a", alist, "%s,%16.2f,%16.2f,%16.2f")
   end if 

  end do 

end
