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

  DATA_DIR  = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  HRES_FL   = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  LENS_FL   = DATA_DIR + "ERA5_ensemble_360x180_monthly_1981-2010.nc" 

  ;constant 
  lvap = 2.501e6
  lice = 3.337e5
  lsub = lvap + lice ; sublimation 
  rhow = 1e3 ; water density kg/m3
  
  f1   = addfile(HRES_FL,"r")
  f2   = addfile(LENS_FL,"r")

  lat  = f1->lat
  lon  = f1->lon
  time = f1->time 

  nens = 11
  nlat = dimsizes(lat)
  nlon = dimsizes(lon)
  ntim = dimsizes(time) 

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; 
 ;;latent heat fluxes directly output by ERA5          ;;
 ;;notes: the mslhf in ERA5 is quite similar to slhf   ;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  lhf1      = new((/ntim,nens,nlat,nlon/),typeof(f1->slhf))
  lhf1!0    = "time"
  lhf1!1    = "nens"
  lhf1!2    = "lat"
  lhf1!3    = "lon"
  lhf1&time = time 
  lhf1&lat  = lat
  lhf1&lon  = lon 

  lhf1(:,0,:,:)  = (/f1->slhf/) ; upward 
  lhf1(:,1:,:,:) = (/f2->slhf/) ; upward
  lhf1           = - lhf1 / 86400.0
  lhf1@units     = "W/m2"
  printMinMax(lhf1,False)

  ;;recalculate with evaporation fluxes, precipitation and snow 
  lhf2           = lhf1
  lhf2(:,0,:,:)  = (/-lsub * f1->mer - lice * f1->mcpr - lice * f1->mlspr/)
  lhf2(:,1:,:,:) = (/-lsub * f2->mer - lice * f2->mcpr - lice * f2->mlspr/)
  lhf2@units     = "W/m2"
  printMinMax(lhf2,False) 

  ;;precipitation 
  prect           = lhf1
  prect(:,0,:,:)  = (/f1->mcpr + f1->mlspr/)
  prect(:,1:,:,:) = (/f2->mcpr + f2->mlspr/)
  prect           = prect * 86400
  prect@units     = "mm/day"
  printMinMax(prect,False)

  prest           = lhf1
  prest(:,0,:,:)  = (/f1->mlssr + f1->mcsr/)
  prest(:,1:,:,:) = (/f2->mlssr + f2->mcsr/)
  prest           = prest * 86400
  prest@units     = "mm/day"
  printMinMax(prest,False)

  qflx            = lhf1
  qflx(:,0,:,:)   = (/f1->mer/)
  qflx(:,1:,:,:)  = (/f2->mer/)
  qflx@units      = "kg/m2/s"
  printMinMax(qflx,False)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;sensible heat fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;; 
  shf           = lhf1
  shf(:,0,:,:)  = (/f1->sshf/)
  shf(:,1:,:,:) = (/f2->sshf/)
  shf           = - shf / 86400.0
  shf@units     = "W/m2"
  printMinMax(shf,False)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;Longwave radiation fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  lwsfc           = lhf1
  lwsfc(:,0,:,:)  = (/f1->str/)  ; All sky (surface net upward)
  lwsfc(:,1:,:,:) = (/f2->str/) 
  lwsfc           = -lwsfc/86400.0         ; change to upward
  lwsfc@units     = "W/m2"
  printMinMax(lwsfc,False) 
 
  lwsfc_clr           = lhf1
  lwsfc_clr(:,0,:,:)  = (/f1->strc/)    ; Clear sky; (surface net upward)
  lwsfc_clr(:,1:,:,:) = (/f2->strc/)
  lwsfc_clr           = - lwsfc_clr/86400.0 ; change to upward
  lwsfc_clr@units     = "W/m2"
  printMinMax(lwsfc_clr,False)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;shortwave radiation fluxes;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  swsfc           = lhf1
  swsfc(:,0,:,:)  = (/f1->ssr/)  ; All sky (surface net downward)
  swsfc(:,1:,:,:) = (/f2->ssr/) 
  swsfc           = swsfc/86400.0 ; change to upward 
  swsfc@units     = "W/m2"
  printMinMax(swsfc,False)

  swsfc_clr           = lhf1
  swsfc_clr(:,0,:,:)  = (/f1->ssrc/) ; Clear sky
  swsfc_clr(:,1:,:,:) = (/f2->ssrc/)
  swsfc_clr           = swsfc_clr/86400.0 ; change to upward 
  swsfc_clr@units     = "W/m2"
  printMinMax(swsfc_clr,False)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;net radiation fluxes;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;
  netrad           = lwsfc 
  netrad           = swsfc - lwsfc 
  netrad@units     = "W/m2"

  netrad_clr       = lwsfc_clr 
  netrad_clr       = swsfc_clr - lwsfc_clr 
  netrad_clr@units = "W/m2"

 ;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;net surface fluxes;;;;;
 ;;positive upward;;;;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;
  netflx1       = lwsfc
  netflx2       = lwsfc
  netflx1       = netrad - lhf1 - shf ; note: lhf1 is used here 
  netflx2       = netrad - lhf2 - shf ; note: lhf2 is used here
  netflx1@units = "W/m2"
  netflx2@units = "W/m2"
 
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
  
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;open netcdf file to save the data for plots;;;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_name = "Surface_energy_fluxes_ERA5.nc"
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
