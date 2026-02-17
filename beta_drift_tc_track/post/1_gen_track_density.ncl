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

  datadir = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/out"
  outdir  = "./"

 ;;; flags to control the plot style
  plot_contour      = False ; if True, plot the contour in the figure
  plot_leftstr      = True  ; if True, the lefadd igure captions as left string
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups    = (/"CLIM","NDGUV","NDGUVT","NDGUVTQ"/)
  ngrps     = dimsizes(Groups)
  LabList   = Groups
  plot      = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  nyrs      = yend - ystr +1 
  timtag    = ystr+"-"+yend 

  ;;construct the 5x5 box 
  dlat         = 5.0 
  dlon         = 5.0 
  sgrd         = dlat+"x"+dlon
  nlat         = toint(180.0/dlat)+1 
  nlon         = toint(360.0/dlon)+1 
  latgrd       = fspan(-90,90,nlat)
  longrd       = fspan(0,360,nlon)
  latgrd!0     = "lat"
  longrd!0     = "lon"
  latgrd@units = "degrees_north"
  longrd@units = "degrees_east"
  latgrd&lat   = latgrd
  longrd&lon   = longrd
  lat2d        = conform_dims((/nlat,nlon/),latgrd,0)
  lon2d        = conform_dims((/nlat,nlon/),longrd,1)
  lat2d!0      = "lat"
  lat2d!1      = "lon"
  lon2d!0      = "lat"
  lon2d!1      = "lon"
  lat2d&lat    = latgrd
  lat2d&lon    = longrd
  lon2d&lat    = latgrd
  lon2d&lon    = longrd
  ;out grid info 
  grid         = "./dummy_grid_"+sgrd+".nc"
  system("/bin/rm -f " + grid)
  fgrid        = addfile(grid,"c")
  fgrid->lat   = latgrd
  fgrid->lon   = longrd
  fgrid->lat2d = lat2d
  fgrid->lon2d = lon2d
  
  do i = 0,ngrps-1,1

   do j = 0,nyrs-1,1
    
     year = j + ystr
     fl     =  datadir + "/sythetic_track_forecast_" + Groups(i)+"_"+year+".nc"
     f      = addfile(fl,"r")
      
     nvars  = f->nvars
     ntrks  = f->ntrks
     nsteps = f->nsteps
     vartrk = f->strks
     varnam = str_squeeze(f->sname)
     vartim = f->stime 
  
     ;mark fill values 
     vartrk@_FillValue = -99999 
     vartrk = where(vartrk.eq."nan".or.vartrk.eq."-nan",vartrk@_FillValue,vartrk)

     indlat = ind(varnam.eq."lat")
     indlon = ind(varnam.eq."lon")
     indtim = ind(varnam.eq."time")

     lat            = vartrk(indlat,:,:) 
     lat!0          = "storm"
     lat!1          = "time"
     lat@units      = "degree_north"
     lat@long_name  = "latitude"

     lon            = lat
     lon            = (/vartrk(indlon,:,:)/)
     lon@units      = "degree_east"
     lon@long_name  = "longitude"

     stim           = lat 
     stim           = (/vartrk(indtim,:,:)/)
     stim@units     = "days since "+(year-1)+"-12-31 00:00:00"
     stim@long_name = "date"

     ;ymdh = cd_calendar(stim, 3)
     ;print(vartim(0,:)+ " " +ymdh(0,:))
     ;exit

     ;vsmc           = lat
     ;vsmc@units     = "m/s"
     ;vsmc@long_name = "maximum sustainable wind speed"
     ;vpmc           = lat
     ;vpmc@units     = "Pa"
     ;vpmc@long_name = "minimum pressure at storm center" 
 
     ;printVarSummary(lat)
     ;printVarSummary(lon)

     lat1d = ndtooned(lat)
     lon1d = ndtooned(lon)
     lon1d = where(lon1d.lt.0.0,lon1d+360.0,lon1d)
     ;print (num(.not.ismissing(lat1d)))
     ;print(dimsizes(lat1d))
     ;exit

     if(j.eq.0) then 
       tcounts     = new((/nyrs,nlat,nlon/),double)
       tcounts!0   = "time"
       tcounts!1   = "lat"
       tcounts!2   = "lon"
       tcounts&time = ispan(0,nyrs-1,1) + year 
       tcounts&lat  = latgrd
       tcounts&lon  = longrd
       tcounts@long_name = "TC track counts (number of storms passed the "+sgrd+" grid box)"
       tcounts      = 0.0 
       tdens        = tcounts 
       tdens@long_name = "TC track density (number of storms passed the "+sgrd+" grid box)"
     end if 

     do mm = 0,nlat-1,1
      do nn = 0,nlon-1,1
         xlat1 = lat2d(mm,nn) - 2.5 
         xlat2 = lat2d(mm,nn) + 2.5
         xlon1 = lon2d(mm,nn) - 2.5 
         xlon2 = lon2d(mm,nn) + 2.5 
        ;print(xlat1+ " " +xlat2 + " " + xlon1+ " " +xlon2)
         indnx = ind((lat1d.gt.xlat1).and.(lat1d.le.xlat2).and.(lon1d.gt.xlon1).and.(lon1d.le.xlon2))
         tcounts(j,mm,nn) = tcounts(j,mm,nn) + num(.not.ismissing(indnx))
         delete([/indnx,xlat1,xlat2,xlon1,xlon2/])
      end do 
     end do 
     delete([/lat1d,lon1d,stim,lat,lon/])
   end do 

   ntot  = 5000
   tdens = tcounts*100.0/5000.0

   ;;;;open netcdf file and save data for the pfigures;;;;;;;
   setfileoption("nc", "Format",  "NetCDF4")
   out_file_nam = outdir+Groups(i)+"_Track_dens_"+timtag+ ".nc"
   system("rm " + out_file_nam)
   fout = addfile(out_file_nam,"c")
   fout->lat       = latgrd
   fout->lon       = longrd
   fout->ntrk_tot  = ntot 
   fout->density   = tcounts 
   fout->dens_norm = tdens

 end do 

end
