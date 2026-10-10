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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Subreg   = (/ "North Pacific",  "CONUS",  "Tropics"/)
  Subrgs   = (/ "NP",             "CONUS",  "Trop"/)
  minlat   = (/  10.0,            30.0,     -10.0/)
  maxlat   = (/  65.0,            48.0,      10.0/)
  minlon   = (/ 110.0,            235.0,     60.0/)
  maxlon   = (/ 260.0,            292.0,    270.0/)
  nsubreg  = dimsizes(Subreg)

  outdir   = "../fig_data/"
  obsdir   = "/global/cfs/cdirs/e3sm/zhan391/data/Observation"
  datadir  = "/pscratch/sd/z/zhan391/DARPA_project/post_process/data/model_output"

  ystr     = 2007
  yend     = 2017
  nyr      = yend - ystr + 1
  freq     = "6hourly"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  Groups   = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ","NDGUVTQ_SRF1","NDGUVTQ_SRF2"/)
  ngrps    = dimsizes(Groups)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"PRECT"/)
  VstList  = (/"PRECT"/)
  UnitList = (/"mm day~S~-1~N~"/)
  nvars    = dimsizes(VarList)

 ;;;loop the regions and plot the PDF

 do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)

  do iy = 0,nyr-1,1
 
   year = ystr + iy 
       
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

    do irg = 0,nsubreg-1,1

     mpMinLatF = minlat(irg)
     mpMaxLatF = maxlat(irg)
     mpMinLonF = minlon(irg)
     mpMaxLonF = maxlon(irg)

     lat = latx({mpMinLatF:mpMaxLatF})
     lon = lonx({mpMinLonF:mpMaxLonF})

     ;;zonal mean;;;;;;;
     rad            = 4.*atan(1.)/180.    ;
     cost           = lat
     cost           = cos(lat*rad)
     pncre          = dim_avg_wgt_n_Wrap(var(:,{mpMinLatF:mpMaxLatF},{mpMinLonF:mpMaxLonF}),cost,1,1)
     pncre@units    = Varunt
     pncre@region   = Subreg(irg)
     pncre@_FillValue = 32766
      
     pncre!0     = "time"
     pncre!1     = "lon"
     pncre&lon   = lon

     outdir = "../fig_data/"
     system("if ! test -d " + outdir +" ; then mkdir -p " + outdir + " ; fi")
     fno    = outdir+ Varname + "_hov_" + Casdir(igp) +"_"+str_sub_str(Subreg(irg)," ","_")+ "_"+year+".nc"
     system("rm -f " + fno)

     print(" ")
     print("output : " + fno)
     print(" ")

     flo = addfile(fno,"c")
     vnb = Varname 
     flo->$vnb$=pncre
     delete([/pncre,lat,lon,rad,cost/])
   end do
   delete([/var,latx,lonx/])
  end do 

 end do ; end of var loop
 
end do ; end of subreg loop

end

