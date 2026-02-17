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

  datadir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/TCS/post_process"
  obsdir  = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/OBS/tc_track"
  outdir  = "./"
  grid    = "./dummy_grid_1x1.nc"

 ;;; flags to control the plot style
  plot_contour      = False ; if True, plot the contour in the figure
  plot_leftstr      = True  ; if True, the lefadd igure captions as left string
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir    = (/"CLIM","NDGUV3_tau06","NDGUVT3_tau06","NDGUVTQ3_tau06"/)
  Groups    = (/"set1",        "set1",         "set1",          "set1"/)
  ngrps     = dimsizes(Groups)
  LabList   = Groups
  plot      = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  timtag    = ystr+"-"+yend 

  ;;construct the 1x1 box 
  fgrid     = addfile(grid,"r")
  latgrd    = fgrid->lat 
  longrd    = fgrid->lon 
  nlat      = dimsizes(latgrd)
  nlon      = dimsizes(longrd)
  lat2d     = conform_dims((/nlat,nlon/),latgrd,0)
  lon2d     = conform_dims((/nlat,nlon/),longrd,1)

  do i = 0,ngrps-1,1

   tdens     = new((/nlat,nlon/),double)
   tdens     = 0.0

   if (Groups(i).eq."OBS") then 
     fna  = obsdir+"/IBTrACS.since1980.v04r00.nc"
     fla  = addfile(fna, "r")
     utc_date = cd_calendar(fla->time(:,0), 0)
     indxx    = ind(utc_date(:,0).ge.ystr.and.utc_date(:,0).le.yend)
     lat  = fla->lat(indxx,:) 
     lon  = fla->lon(indxx,:) 
     stim = fla->time(indxx,:)
     vsmc = fla->wmo_wind(indxx,:)
     vpmc = fla->wmo_pres(indxx,:)
     nst  = dimsizes(lat(:,0))
     do ix = 0,nst-1,1
      indyy = ind(.not.(ismissing(lat(ix,:)).and.ismissing(lon(ix,:))))
      lat1  = lat(ix,indyy) 
      lon1  = lon(ix,indyy)
      lon1  = where(lon1.lt.0.0,lon1+360.0,lon1)
      nm = getind_latlon2d (lat2d,lon2d, lat1, lon1)
      tdens(nm(:,0),nm(:,1)) = tdens(nm(:,0),nm(:,1)) + 1
      delete([/indyy,lat1,lon1,nm/])
     end do 
     delete([/fna,fla,indxx,utc_date/])
   else
     fna  = datadir+"/"+Casdir(i)+"_TCS_"+Groups(i)+"_Track_"+ timtag+ ".txt"
     fla  = asciiread(fna,-1,"string")
     inds = str_match_ind_ic_regex(fla,"start")
     nst  = num(.not.ismissing(inds)) 
     nyy  = 360 
     lat  = new((/nst,nyy/),float)
     lat!0          = "storm"
     lat!1          = "time"
     lat@units      = "degree_north"
     lat@long_name  = "latitude"
     lon            = lat
     lon@units      = "degree_east"
     lon@long_name  = "longitude"
     stim           = lat 
     stim@units     = "hours since 2007-01-01 00:00:00"
     stim@long_name = "date"
     vsmc           = lat
     vsmc@units     = "m/s"
     vsmc@long_name = "maximum sustainable wind speed"
     vpmc           = lat
     vpmc@units     = "Pa"
     vpmc@long_name = "minimum pressure at storm center" 
 
     ntlin = 0
     do ix = 0,nst-1,1
      ntlin = ntlin + 1
      head  = fla(inds(ix))
      ntrk  = stringtointeger(str_get_field(head, 2, "	"))
      ns1   = ntlin
      ne1   = ntlin + ntrk -1
      hsub  = fla(ns1:ne1)
      ntlin = ntlin + ntrk
      lon(ix,0:ntrk-1) = stringtofloat(str_get_field(hsub, 2,"	"))
      lat(ix,0:ntrk-1) = stringtofloat(str_get_field(hsub, 3,"	"))
      vsmc(ix,0:ntrk-1)= stringtofloat(str_get_field(hsub, 5,"	"))
      vpmc(ix,0:ntrk-1)= stringtofloat(str_get_field(hsub, 4,"	"))
 
      yyyy = stringtointeger(str_get_field(hsub, 7,"	"))
      mm   = stringtointeger(str_get_field(hsub, 8,"	"))
      dd   = stringtointeger(str_get_field(hsub, 9,"	"))
      hh   = stringtointeger(str_get_field(hsub, 10,"	"))
      mn   = mm 
      ss   = mm 
      mn   = 0
      ss   = 0 
      date = cd_inv_calendar(yyyy,mm,dd,hh,mn,ss,stim@units, 0)
      stim(ix,0:ntrk-1) = tofloat(date)

      lat1  = lat(ix,0:ntrk-1)
      lon1  = lon(ix,0:ntrk-1)
      lon1  = where(lon1.lt.0.0,lon1+360.0,lon1)
      nm = getind_latlon2d (lat2d,lon2d, lat1, lon1)
      tdens(nm(:,0),nm(:,1)) = tdens(nm(:,0),nm(:,1)) + 1
      delete([/yyyy,mm,dd,hh,mn,ss,date,hsub,ns1,ne1,ntrk,head,lat1,lon1,nm/])
     end do 
     delete([/fna,fla,inds/])

   end if

   print ("total number of TCs in " + Groups(i) + " is " + nst)
   tdens!0 = "lat"
   tdens!1 = "lon"
   tdens&lat = latgrd
   tdens&lon = longrd
   tdens@long_name = "track density (number of sttorms passed the 1x1 grid box)"

   ;;;;open netcdf file and save data for the pfigures;;;;;;;
   setfileoption("nc", "Format",  "NetCDF4")
   out_file_nam1 = outdir+Casdir(i)+"_"+Groups(i)+"_Track_info_"+ timtag+ ".nc"
   out_file_nam2 = outdir+Casdir(i)+"_"+Groups(i)+"_Track_dens_"+ timtag+ ".nc"
   system("rm " + out_file_nam1)
   system("rm " + out_file_nam2)
   fout1 = addfile(out_file_nam1,"cw")
   fout2 = addfile(out_file_nam2,"cw")
   fout1->lat  = lat
   fout1->lon  = lon
   fout1->stim = stim 
   fout1->vsmc = vsmc 
   fout1->vpmc = vpmc

   fout2->lat      = latgrd
   fout2->lon      = longrd
   fout2->density  = tdens  

   delete([/lat,lon,stim,vsmc,vpmc,tdens/])

 end do 

end
