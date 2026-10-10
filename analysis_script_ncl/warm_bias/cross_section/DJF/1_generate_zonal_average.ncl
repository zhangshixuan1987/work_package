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

  ModelOutRoot = "/global/cfs/cdirs/e3sm/zhan391/data/E3SM/climo/"
  ObsOutRoot   = "/global/cfs/cdirs/e3sm/zhan391/data/ERA5/regrid/1981-2010/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Seasons      = (/"DJF"/)
  nseas        = dimsizes(Seasons)
  lats         = -90.0
  latn         = 0.0 
  minz         = 500 
  maxz         = 1000

  ;;;read the experiment names and short names;;;
  exp_list     = "experiment_list.txt"
  data         = asciiread(exp_list,-1,"string")
  delim        = ","
  Casdir       = str_squeeze(str_get_field(data, 1,delim))
  Groups       = str_squeeze(str_get_field(data, 2,delim))
  ngrps        = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList      = (/"T","Q","U","RELHUM","CLOUD"/)
  VobsList     = (/"t","q","u","r","cc"/)
  VnmList      = (/"Temperature","Specific humidity", "Zonal wind","Relative humidity", "Cloud fraction"/)
  UnitList     = (/"K","g kg~S~-1~N~","m s~S~-1~N~","unitless","fraction"/)
  PltNums      = (/"a", "b", "c", "c","d","e","f","g","h","i","j","k"/)
  nvars        = dimsizes(VarList)

 ;;flags to perform significant test;;;;;;;;
 setfileoption("nc", "Format",  "NetCDF4")
 out_file_name = "E3SM_versus_ERA5.nc"
 system("rm " + out_file_name)
 fout = addfile(out_file_name,"c")

do iv = 0,nvars-1, 1

 Varobs = VobsList(iv)
 Varmod = VarList(iv)
 Varunt  = UnitList(iv)
 Vnmstr  = VnmList(iv)

 do itim = 0, nseas-1

 season  = Seasons(itim)
  
 do igp = 0,ngrps-1

  obsfile = ObsOutRoot+ "era5_prelev_analysis_360x180_"+season +".nc"
  modfile = ModelOutRoot+Casdir(igp)+"/"+Casdir(igp)+"_"+season+"_climo.nc"
  File0   = addfile(obsfile,"r")
  File1   = addfile(modfile,"r")
  print(obsfile)
  print(modfile)

  lat             = File0->lat({lats:latn})
  lon             = File0->lon
  lev_p           = (/File0->lev({minz:maxz})/)
  lev_p@long_name = "pressure"      
  lev_p@units     = "hPa"
  print(lev_p)
  x11   = dim_avg_n_Wrap(dim_avg_n_Wrap(File0->$Varobs$(:,{minz:maxz},{lats:latn},:),3),0)

  mlev  = dimsizes(File1->lev) 
  latm  = File1->lat
  lonm  = File1->lon
  P0mb  = 1000.0
  intyp = 2

  yint  = vinth2p(File1->$Varmod$,File1->hyam,File1->hybm,lev_p,File1->PS,intyp,P0mb,1,False)
  yrgd  = linint2_Wrap (lonm, latm, yint, True, lon, lat, 0)
  y11   = x11
  y11   = (/dim_avg_n_Wrap(dim_avg_n_Wrap(yrgd(:,{minz:maxz},{lats:latn},:),3),0)/)

  delete([/yint,mlev,P0mb,latm,lonm,lev_p/])
  ;printVarSummary(x11)
  ;printVarSummary(y11)

  if(Varmod.eq."Q") then
   x11 = x11*1000.0
   y11 = y11*1000.0
  end if
  y11  = y11 - x11
 
  if (igp.eq.0)then
   vout = "ERA5_"+Varmod+"_"+season
   x11@description = Vnmstr+" in " + "ERA5 reanalysis (1981-2010)"
   fout->$vout$    = x11
   fout->lat       = lat
  end if
  vout = Groups(igp)+"-ERA5_"+Varmod+"_"+season
  y11@description  = Vnmstr+" in " + Groups(igp)
  fout->$vout$     = y11

end do 

end do ; end of igp loop

end do ; end of var loop
end
