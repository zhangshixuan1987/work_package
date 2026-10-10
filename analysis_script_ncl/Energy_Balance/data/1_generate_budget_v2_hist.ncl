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
begin

  Exp_Nam   = "v2.LR.historical"
  DATA_DIR  = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  E3SM_FILE = systemfunc ("ls "+DATA_DIR+Exp_Nam+"*.nc") ; file paths

  ;constant 
  lvap = 2.501e6
  lice = 3.337e5
  lsub = lvap + lice
  rhow = 1e3 ; water density kg/m3
  
  f   = addfiles(E3SM_FILE,"r")
  lat = f[0]->lat
  lon = f[0]->lon
  ListSetType (f, "join")

 ;;;;;;;;;;;;;;;;;;;;;;;;; 
 ;;latent heat fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;  
  ;;direct output by EAMv2;;;; 
  lhf1 = f[:]->LHFLX ; upward 
  lhf1 = lhf1
  ;;recalculate with evaporation fluxes, precipitation and snow 
  prect = f[:]->PRECC 
  prect = f[:]->PRECC + f[:]->PRECL
  prest = f[:]->PRECSC
  prest = f[:]->PRECSC + f[:]->PRECSL
  qflx  = f[:]->QFLX ; evaporation 
  lhf2  = lhf1
  lhf2  = lsub * qflx - lice * rhow *(prect - prest) 
 ;print(min(lhf1) + " " + max(lhf1))
 ;print(min(lhf2) + " " + max(lhf2))
 ;print (avg(qflx) + "  " +avg(prect)* rhow +  " " +avg(prest)* rhow)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;sensible heat fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  shf = f[:]->SHFLX ; upward
  print(min(shf) + " " + max(shf))

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;Longwave radiation fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lwsfc = f[:]->FLNS  ; All sky (upward)
  lwsfc_clr = f[:]->FLNSC  ; Clear sky; upward 

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;shortwave radiation fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  swsfc = f[:]->FSNS  ; downward; All sky 
  swsfc_clr = f[:]->FSNSC   ; Clear sky

  netrad = lwsfc 
  netrad = swsfc - lwsfc 

  netrad_clr = lwsfc_clr 
  netrad_clr = swsfc_clr - lwsfc_clr 

  netflx1 = lwsfc
  netflx2 = lwsfc
  netflx1 = netrad - lhf1 - shf ; note: lhf1 is used here
  netflx2 = netrad - lhf2 - shf ; note: lhf2 is used here

  swsfc@units      = "W/m2"
  swsfc_clr@units  = "W/m2"
  lwsfc@units      = "W/m2"
  lwsfc_clr@units  = "W/m2"
  shf@units        = "W/m2"
  lhf1@units       = "W/m2"
  lhf2@units       = "W/m2"
  netrad@units     = "W/m2"
  netrad_clr@units = "W/m2"
  netflx1@units    = "W/m2"
  netflx2@units    = "W/m2"

  swsfc@description = "Net shortwave radiation fluxes at the surface"
  lwsfc@description = "Net longwave radiation fluxes at the surface"
  swsfc_clr@description = "Net shortwave radiation fluxes (clear-sky) at the surface"
  lwsfc_clr@description = "Net shortwave radiation fluxes (clear-sky) at the surface"
  netrad@description = "Net radiation fluxes at the surface"
  netrad_clr@description = "Net radiation fluxes (clear-sky) at the surface"
  shf@description  = "Sensible heat fluxes (land model) at the surface"
  lhf2@description = "Latent heat fluxes (Lv+Lf) QFLX - Lf 1e3 (PRECC + PRECL - PRECSC - PRECSL)"
  lhf1@description = "Latent heat fluxes (land model)"
  
  netflx1@description = "Net surfalce flux: FSNS - FLNS - SHFLX - LHFLX (Land surface)"
  netflx2@description = "Net surfalce flux: FSNS - FLNS - SHFLX - LHFLX (New formula)"

  prect = prect * 1000 * 86400 
  prest = prest * 1000 * 86400

  prect@units = "mm/day"
  prest@units = "mm/day" 
  qflx@units  = "kg/m2/s"

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;open netcdf file to save the data for plots;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_name = "Surface_energy_fluxes_"+Exp_Nam+".nc"
  system("rm " + out_file_name)

  fout = addfile(out_file_name,"c")
  fout->PRECT      = prect
  fout->PREST      = prest
  fout->QFLX       = qflx
  fout->LHF        = lhf2
  fout->LHF_LAND   = lhf1
  fout->SHF        = shf
  fout->SWSFC      = swsfc
  fout->LWSFC      = lwsfc
  fout->SWSFC_C    = swsfc_clr
  fout->LWSFC_C    = lwsfc_clr
  fout->NETCFSFC   = netrad
  fout->NETCFSFC_C = netrad_clr
  fout->NETFLX1    = netflx1
  fout->NETFLX2    = netflx2

end 
