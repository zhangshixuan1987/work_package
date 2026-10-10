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

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.
  grav      = 9.80665

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DATA_DIR    = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  ERA5_File   = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  EAMv2_File   = DATA_DIR + "20210727.v2.LR.F20TR_eam.h0_monthly_1981-2010.nc"
  intyp       = 2

  f1          = addfile(ERA5_File,"r")
  pfix        = f1->lev

  f2          = addfile(EAMv2_File,"r")
  lat         = f2->lat
  phis        = f2->PHIS(:,:,:)
  psrm        = f2->PS(:,:,:)
  tmod        = f2->T(:,:,:,:) 
  zmod        = f2->Z3(:,:,:,:)
  
  zsmd        = phis
  zsmd        = zsmd /grav
  zsmd@units  = "m"
  pmod        = pfix * 100.0 ; convert to Pa

  dimx = dimsizes(tmod)
  ntime = dimx(0)
  nlev  = dimx(1)
  nlat  = dimx(2)
  nlon  = dimx(3)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;interpolate to ERA5 pressure levels
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  tbot   = tmod(:,nlev-1,:,:)
  intyp  = 2                             ; 1=linear, 2=log, 3=log-log
  kxtrp  = True                          ; True=extrapolate
  varflg = 1                             ; temperature is variable
  tpm    = vinth2p_ecmwf(tmod,f2->hyam,f2->hybm,pfix,psrm,intyp,f2->P0/100.0, \
               1,kxtrp,varflg,tbot,zsmd)
  varflg = -1
  zpm    = vinth2p_ecmwf(zmod,f2->hyam,f2->hybm,pfix,psrm,intyp,f2->P0/100.0, \
               1,kxtrp,varflg,tbot,zsmd)
  printVarSummary(tpm)
  printVarSummary(zpm)

  ptop   = min(pfix)
  p2tdz  = stdatmus_p2tdz(ptop)
  ztop   = p2tdz(2)

  dpmod  = zpm 
  dzmod  = zpm
  dpmod  = dpres_plevel(pmod, psrm, ptop, 0) 
  dzmod  = tofloat(dz_height(zpm,  conform(zpm,zsmd,(/0,2,3/)), ztop, 0))
  
  dpmod@units = "Pa"
  dzmod@units = "m"

  ;;open netcdf file to save the data for plots;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_name = "Hydrostatic_analysis_out_EAMv2.nc"

  system("rm " + out_file_name)
  fout = addfile(out_file_name,"c")
  fout->T    = tpm
  fout->Z    = zpm
  fout->plev = pfix
  fout->PHIS = phis
  fout->ZS   = zsmd
  fout->PS   = psrm
  fout->dz   = dzmod
  fout->dp   = dpmod

  delete([/tpm,zpm,pfix,zsmd,psrm,dzmod,dpmod/])

exit

end 

