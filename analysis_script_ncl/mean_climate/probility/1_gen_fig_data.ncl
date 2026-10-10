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

  outdir   = "../fig_data/"
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data/Observation"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"

  ystr     = 2008
  yend     = 2017
  nyr      = yend - ystr + 1
  freq     = "6hourly"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  ngrps    = dimsizes(Groups)

  VarList  = (/"TREFHT", "PRECT"/)
  UntList  = (/"K", "mm day~S~-1~N~"/)
  nvars    = dimsizes(VarList)

  ;;;;constant
  pi        = 3.14159265
  rearth    = 6.37122e6
  grav      = 9.806
  sigma     = 5.67e-8 ; Stefan-Boltzmann constant

 do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)

   do igp = 0,ngrps-1,1

    if(igp.eq.0) then
     fl0  = systemfunc("ls " + obsdir + "/"+Casdir(igp)+"/daily/" + Casdir(igp)+"*"+year+"_1x1.nc")
     f0   = addfile(fl0,"r")
     if(isleapyear(year)) then
       tmp    = f0->$Varname$
       utcd   = cd_calendar(f0->time,0)
       indtim = ind(.not.(utcd(:,1).eq.2.and.utcd(:,2).eq.29))
       var    = tmp(indtim,:,:)
       delete([/utcd,indtim,tmp/])
     else
       var  = f0->$Varname$
     end if
     latx = f0->lat
     lonx = f0->lon
    else
     fl0  = systemfunc("ls " + datadir + "/"+Casdir(igp)+"/" + freq + "/" +Varname+"*"+year+"*.nc")
     f0   = addfiles(fl0,"r")
     tmp  = f0[:]->$Varname$
     latx = f0[0]->lat
     lonx = f0[0]->lon
     opt  = True
     opt@nval_crit = 4
     var = calculate_daily_values (tmp, "avg", 0, opt)
     var = var * 24.0 * 3600.0 * 1000.0
     delete([/tmp/])
    end if
    delete([/fl0,f0/])

    lat    = File[0]->lat
    lon    = File[0]->lon

    printVarSummary(pxcre0300)
    pxcre0300@_FillValue = -99999
    pncre1d = ndtooned(pxcre0300)
    pncre1d@_FillValue = -99999
 
    opt     = False
    nbin    = 40
    vpdf1   = pdfx(var1, nbin, opt)

    outdir = "./data/"
    system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
    fno    = outdir+ VarList(iv) + "_" + Casdir(igp) + ".nc"
    system("rm -f " + fno)
   ;setfileoption("nc", "Format",  "NetCDF4")

    print(" ")
    print("output : " + fno)
    print(" ")

    flo = addfile(fno,"c")
 
    vna = Varname + "_PDF_" + Casdir(igp)
    flo->$vna$=vpdf1
    delete([/vpdf1/])

  end do 

 delete(flo)

end do 

end
