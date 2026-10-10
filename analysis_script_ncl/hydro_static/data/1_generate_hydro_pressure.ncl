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
  Rd        = 287.053 

  Groups      = (/"ERA5","EAMv2"/)
  ngrps       = dimsizes(Groups)

 ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  DATA_DIR    = "/global/cfs/cdirs/e3sm/zhan391/local_diag/hydro_static/data/"
  filnam1     = DATA_DIR + "Hydrostatic_analysis_out_"+Groups(0)+".nc"
  filnam2     = DATA_DIR + "Hydrostatic_analysis_out_"+Groups(1)+".nc"
  f1          = addfile(filnam1,"r")
  f2          = addfile(filnam2,"r")

  zso         = f1->ZS
  zsm         = f2->ZS
  zso         = zso / grav
  zsm         = zsm / grav


  zo          = f1->Z
  zm          = f2->Z
  to          = f1->T
  tm          = f2->T
  po          = f1->lev
  pm          = f2->lev

  dzo         = f1->dz
  dzm         = f2->dz

  psm         = f2->PS
  pso         = psm
  pso         = (/f1->PS/)


  zso1        = conform(zo,zso,(/0,2,3/))
  zsm1        = conform(zm,zsm,(/0,2,3/))
  po1         = conform(zo,po,(/1/))
  pm1         = conform(zm,pm,(/1/))

  printVarSummary(zo) 
  printVarSummary(zm)
  printVarSummary(zso)
  printVarSummary(zsm)

  print(min(zo) + " " +max(zo))
  print(min(zm) + " " +max(zm))
  print(min(zso1) + " " +max(zso1))
  print(min(zsm1) + " " +max(zsm1))

  dzo_msk     = where(zo.le.zso1.or.zo.le.zsm1,0.0,dzo)
  dzm_msk     = where(zm.le.zso1.or.zm.le.zsm1,0.0,dzm)
 
  print(num(dzo_msk.eq.0))
  print(num(dzm_msk.eq.0))

  ;; dp = P*g*dz/Rd/T
  dpo  = po1*100.0 * grav * dzo_msk / Rd / to  
  dpm0 = pm1*100.0 * grav * dzm_msk / Rd / tm
  dpm1 = pm1*100.0 * grav * dzm_msk / Rd / to
  dpm2 = pm1*100.0 * grav * dzo_msk / Rd / tm

  pso1 = psm
  psm0 = psm
  psm1 = psm
  psm2 = psm

  pso1 = tofloat(dim_sum_n(dpo,1))
  psm0 = tofloat(dim_sum_n(dpm0,1))
  psm1 = tofloat(dim_sum_n(dpm1,1))
  psm2 = tofloat(dim_sum_n(dpm2,1))

  print(min(pso1)+ " " +max(pso1))
  print(min(psm0)+ " " +max(psm0))
  print(min(psm1)+ " " +max(psm1))
  print(min(psm2)+ " " +max(psm2))

  ;;open netcdf file to save the data for plots;;;;;
  setfileoption("nc", "Format",  "NetCDF4")
  out_file_name = "Hydrostatic_Surface_Pressure_ERA5_EAMv2.nc"
  system("rm " + out_file_name)
  fout = addfile(out_file_name,"c")
  fout->PS_ERA5            = pso
  fout->PS_EAM             = psm
  fout->PS_Hydro_ERA5      = pso1
  fout->PS_Hydro_EAM       = psm0
  fout->PS_Hydro_EAM_TER5  = psm1
  fout->PS_Hydro_EAM_DZER5 = psm2

 exit

end 
