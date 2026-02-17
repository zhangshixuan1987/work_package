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

  l_create_plot_data = False
  l_create_plot_data = True

  ModelOutRoot = "/compyfs/zhan391/darpa_scratch/post_process/data"
  ObsOutRoot   = "/compyfs/zhan391/darpa_scratch/post_process/data/OBS/" 

 ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir    = (/"UVTQ3_BASE","UVTQPS3_NPBL_NUPP","UVTQ3_CTRL","UVTQ3_GFDL","UVTQ3_NPBL","UVTQ3_NPBL_NUPP"/)
  Groups    = (/"UVTQ3_BASE","UVTQPS3_NPBL_NUPP","UVTQ3_CTRL","UVTQ3_GFDL","UVTQ3_NPBL","UVTQ3_NPBL_NUPP"/)
  ngrps     = dimsizes(Groups)
  plot      = new (ngrps, graphic)

  year      = 2010
  season    = "DJF"
  istr      = 200912
  iend      = 201002
  freq      = "monthly"

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Obsnam    = (/"ERA5",  "ERA5",  "ERA5", "ERA5"/)
  VarList   = (/"PSL",   "U10",   "TS",   "TREFHT"/)
  VnmList   = (/"Sea level pressure", "10m wind speed", \
                "Surface temperature", "2m air temperature"/)
  UnitList  = (/"hPa","m s~S~-1~N~","K","K"/)
  nvars     = dimsizes(VarList)

 do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)

  print("iv : " + iv) 

 do igp = 0, ngrps-1,1

  print(" ")
  print("........................... ")
  print("igp : " + igp) 
  print("........................... ")

  ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  if(Obsnam(iv).eq."ERA5") then 
    fn0    = ObsOutRoot + Obsnam(iv) + "/" + freq + "/" + "2009-2010.nc"
  else 
    fn0    = ObsOutRoot + Obsnam(iv) + "/" + freq + "/" + "2009-2011.nc"
  end if 

  fn1    = ModelOutRoot+ "/" + Casdir(igp)+"/"+freq+"/" + "2010_h0.nc"
  File0  = addfile(fn0,"r")
  File1  = addfile(fn1,"r")

  fn2    = ModelOutRoot+"/landmsk/landmsk.nc"
  fl2    = addfile(fn2, "r")
  lmsk   = fl2->LANDFRAC(0,:,:)
  
  print(" ") 
  print("input : " + fn0) 
  print("input : " + fn1) 
  print(" ") 
  
  if (Obsnam(iv).eq."OAFlux") then 
    tobs   = File0->time
    tobs   = tobs -16.0
    obstim = cd_calendar(tobs,-1) ;File0->time,-1)
  else 
    obstim = cd_calendar(File0->time,-1)
  end if 
  modtim = cd_calendar(File1->time,-1)
  indobs = ind(obstim.ge.istr.and.obstim.le.iend)
  indmod = ind(modtim.ge.istr.and.modtim.le.iend)
  print(obstim(indobs) + "  " +modtim(indmod))

  lat    = File1->lat
  lon    = File1->lon

  if (Varname .eq. "U10" ) then 
    pxcre0300 = wind_speed(File0->U10(indobs,:,:), File0->V10(indobs,:,:))
  else 
    pxcre0300 = File0->$Varname$(indobs,:,:)
  end if 

  wxcre0300 = File1->$Varname$(indmod,:,:)

  if (Varname .eq. "PSL" ) then 
    pxcre0300 = pxcre0300 / 100.0
    wxcre0300 = wxcre0300 / 100.0
  end if

 ;printVarSummary(pxcre0300)
 ;printVarSummary(wxcre0300)

  pncre0300 = dim_avg_n_Wrap(pxcre0300,0)
  wncre0300 = pncre0300
  wncre0300 = dim_avg_n_Wrap(wxcre0300,0)

  wncredf1  = wncre0300
  wncredf1  = wncre0300 - pncre0300
  
  printMinMax(pncre0300,False)
  printMinMax(wncre0300,False)

  pncre0300@long_name = VnmList(iv)   
  pncre0300@units     = UnitList(iv) 
  wncredf1@long_name = VnmList(iv)   
  wncredf1@units     = UnitList(iv) 
  
  wncredf2  = wncre0300 
  wncredf2@long_name = VarList(iv) + "   " + Casdir(igp) + " - " + Casdir(0)
  wncredf2@units     = UnitList(iv)

  mskxx     = pncre0300
  mskxx     = where(pncre0300.ne.0.0,pncre0300,32766)
  mskxx@_FillValue = 32766
  wncredf2  = (wncre0300 - mskxx)*100/mskxx
  wncredf2  = where(ismissing(wncredf2),0.0,wncredf2)
  delete(mskxx)

  if (Varname .eq. "NETCF") then
    wncredf2 = where(abs(wncredf1).le.2.0,0.0,wncredf2)
  end if  
  
  rad     = 4.*atan(1.)/180.    ;
  cost    = cos(lat*rad)       ; cosine weights
  mean    = wgt_areaave(pncre0300,cost,1,0)
  pncre0300@mean = sprintf("%5.2f", mean)

  mean    = wgt_areaave(wncredf1,cost,1,0)
  rmse    = wgt_arearmse(pncre0300,wncre0300,cost,1,0)
  corr    = pattern_cor(pncre0300,wncre0300,cost,0)
  wncredf1@mean = sprintf("%5.2f", mean)
  wncredf1@rmse = sprintf("%5.2f", rmse)
  wncredf1@corr = sprintf("%5.2f", corr)

  wncredf2@mean = sprintf("%5.2f", mean)
  wncredf2@rmse = sprintf("%5.2f", rmse)
  wncredf2@corr = sprintf("%5.2f", corr)

 ;print(mean + " " + rmse + "  "+corr)
  delete([/mean,rmse,corr,rad,cost/])

  outdir = "./"
  system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
  fno    = outdir+ VarList(iv) + "_" + Obsnam(iv) + "_" + Casdir(igp) + ".nc"
  system("rm -f " + fno)
 ;setfileoption("nc", "Format",  "NetCDF4")

  print(" ")
  print("output : " + fno)
  print(" ")

  flo = addfile(fno,"c")

  vna = Varname + "_" + Obsnam(iv)
  vnb = Varname + "_" + Casdir(igp)
  vnc = "d" + Varname

  flo->$vna$=pncre0300
  flo->$vnb$=wncredf1
  flo->$vnc$=wncredf2

  delete([/pncre0300,wncre0300,pxcre0300,wxcre0300,wncredf1,wncredf2/])
  delete([/obstim,modtim,indobs,indmod/])

 end do 

 delete(flo)

end do 

end
