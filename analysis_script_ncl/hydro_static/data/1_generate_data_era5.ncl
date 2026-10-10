;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
external MAKEIC "./MAKEIC.so"
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
  f1          = addfile(ERA5_File,"r")
  pfix        = f1->lev
  phis        = f1->phis(:,:,:)
  ter5        = f1->t(:,:,:,:)
  psr5        = f1->sp(:,:,:)
  zer5        = f1->z(:,:,:,:)
  zer5        = zer5/grav 
  zsr5        = phis
  zsr5        = phis/grav
  per5        = pfix * 100.0 ; convert to Pa

  psr5@units  = "Pa"
  per5@units  = "Pa"
  zsr5@units  = "m"
  zer5@units  = "m"

  dimx = dimsizes(ter5)
  ntime = dimx(0)
  nlev  = dimx(1)
  nlat  = dimx(2)
  nlon  = dimx(3)

  ptop   = min(pfix)
  p2tdz  = stdatmus_p2tdz(ptop)
  ztop   = p2tdz(2)
  dper5  = zer5
  dzer5  = zer5
  dper5  = dpres_plevel(pfix, psr5, ptop, 0)  ; dp(30)
  dzer5  = tofloat(dz_height(zer5, conform(zer5,zsr5,(/0,2,3/)), ztop, 0))

  ;;open netcdf file to save the data for plots;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_name = "Hydrostatic_analysis_out_ERA5.nc"
  system("rm " + out_file_name)
  fout = addfile(out_file_name,"c")
  fout->T    = ter5
  fout->Z    = zer5
  fout->plev = pfix
  fout->PHIS = phis
  fout->ZS   = zsr5
  fout->PS   = psr5
  fout->dz   = dzer5
  fout->dp   = dper5

  delete([/ter5,zer5,pfix,zsr5,psr5,dzer5,dper5,p2tdz/])

 exit

end 
