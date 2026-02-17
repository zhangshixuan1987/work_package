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

 ;;; flags to control the plot style
  l_plot_latlon_box = False
  datdir      = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ARs/post_process/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5", "CLIM", "NDGUV3_tau06", "NDGUVTQ3_tau06"/)
  Groups  = (/"ERA5", "CLIM", "NDGUV3_tau06", "NDGUVTQ3_tau06"/)
  Labels  = (/"ERA5", "E3SM (CLIM)", "E3SM (NDGUV)", "E3SM (NDGUVTQ)"/)
  ngrps   = dimsizes(Groups)
  ystr    = 2007
  yend    = 2017

 ;;construct the 1x1 box
  hres         = 5.0
  nlat         = 36
  nlon         = 72
  latgrd       = fspan(-87.5,87.5,nlat)
  longrd       = fspan(2.5,357.5,nlon)
  latgrd!0     = "lat"
  longrd!0     = "lon"
  latgrd&lat   = latgrd
  longrd&lon   = longrd
  latgrd@units = "degrees_north"
  longrd@units = "degrees_east"
  lat2d        = conform_dims((/nlat,nlon/),latgrd,0)
  lon2d        = conform_dims((/nlat,nlon/),longrd,1)

  do year = ystr, yend 

   do i = 0,0 ;ngrps-1,1

   fna = systemfunc("ls -1 " + datdir+Casdir(i)+"/"+Casdir(i)+"*"+year+"*_nff.nc")
   fla = addfiles(fna, "r")

   fnb = systemfunc("ls -1 " + datdir+Casdir(i)+"/"+Casdir(i)+"*"+year+"*_ivt.nc")
   flb = addfiles(fnb, "r")
   print("working on file: " + fna)

   if (Casdir(i).eq."ERA5") then
     lat   = flb[0]->latitude
     lon   = flb[0]->longitude
     tvq   = flb[:]->VIWVN_PW_AR(:,:,:)
     tuq   = flb[:]->VIWVE_PW_AR(:,:,:)  
    ;tvq   = flb[:]->VIWVN_PW_AR(:,:,:)  + flb[:]->VIWVN_PW_NONAR(:,:,:) 
    ;tuq   = flb[:]->VIWVE_PW_AR(:,:,:)  + flb[:]->VIWVE_PW_NONAR(:,:,:)
     ivt   = tuq
     ivt   = sqrt(tuq*tuq+tvq*tvq)
   else
     lat   = flb[0]->lat
     lon   = flb[0]->lon
     tvq   = flb[:]->VIWVN_PW_AR(:,:)
     tuq   = flb[:]->VIWVE_PW_AR(:,:)
    ;tvq   = tvq + flb[:]->VIWVN_PW_NONAR(:,:)
    ;tuq   = tuq + flb[:]->VIWVE_PW_NONAR(:,:)
     ivt   = tuq
     ivt   = sqrt(tuq*tuq+tvq*tvq)
   end if

   if (Casdir(i).eq."ERA5") then 
     vara   = fla[:]->AR_binary_tag(:,:,:)
   else
     vara   = fla[:]->AR_binary_tag(:,:) 
   end if 
   var0 = byte2flt(vara)
   var0@_FillValue = -9999.0
   tuq@_FillValue  = -9999.0
   tvq@_FillValue  = -9999.0
   ivt@_FillValue  = -9999.0

   ;; do a statistic on a given box 
   vfreq           = new((/nlat,nlon/),double)
   vfreq!0         = "lat"
   vfreq!1         = "lon"
   vfreq&lat       = latgrd
   vfreq&lon       = longrd
   vfreq@long_name = "AR density (number of ARs within the 5x5 grid box)"
   vfreq = 0.0
   vtuq  = vfreq
   vtvq  = vfreq
   vivt  = vfreq
   vtuq@long_name = "AR vapor transport (zonal)"
   vtuq@long_name = "AR vapor transport (meridional)"
   vivt@long_name = "AR vapor transport (total)"

   do ix = 0,nlat-1,1
    do iy = 0,nlon-1,1
     lat1  = latgrd(ix) - hres/2.0
     lat2  = latgrd(ix) + hres/2.0
     lon1  = longrd(iy) - hres/2.0
     lon2  = longrd(iy) + hres/2.0
     if (Casdir(i).eq."ERA5") then 
       indx = ind(lat.ge.lat1.and.lat.le.lat2)
       indy = ind(lon.ge.lon1.and.lon.le.lon2)
       tmp0 = var0(:,indx,indy)
       tmp1 = tuq(:,indx,indy)
       tmp2 = tvq(:,indx,indy)
       tmp3 = ivt(:,indx,indy)
       rad  = 4.0*atan(1.0)/180.0
       clat = cos(lat(indx)*rad)
       vfreq(ix,iy) = dim_sum(ndtooned(tmp0*conform(tmp0,clat,1)))*100.0/dim_sum(ndtooned(conform(tmp0,clat,1)))
       vtuq(ix,iy)  = dim_avg_n_Wrap(wgt_areaave_Wrap(tmp1,clat,1.0,0),0)
       vtvq(ix,iy)  = dim_avg_n_Wrap(wgt_areaave_Wrap(tmp2,clat,1.0,0),0)
       vivt(ix,iy)  = dim_avg_n_Wrap(wgt_areaave_Wrap(tmp3,clat,1.0,0),0)
       delete([/indx,indy,rad,clat,tmp0,tmp1,tmp2,tmp3/])
     else
       indx = ind(lat(0,:).ge.lat1.and.lat(0,:).le.lat2.and.lon(0,:).ge.lon1.and.lon(0,:).le.lon2)
       if(.not.all(ismissing(indx))) then 
        tmp0 = var0(:,indx)
        tmp1 = tuq(:,indx)
        tmp2 = tvq(:,indx)
        tmp3 = ivt(:,indx)
        rad  = 4.0*atan(1.0)/180.0
        clat = cos(lat(:,indx)*rad)
        vfreq(ix,iy) = dim_sum(ndtooned(tmp0*clat))*100.0/dim_sum(ndtooned(clat))
        vtuq(ix,iy)  = dim_sum(ndtooned(tmp1*clat))/dim_sum(ndtooned(clat))
        vtvq(ix,iy)  = dim_sum(ndtooned(tmp2*clat))/dim_sum(ndtooned(clat)) 
        vivt(ix,iy)  = dim_sum(ndtooned(tmp3*clat))/dim_sum(ndtooned(clat))
        delete([/rad,clat,tmp0,tmp1,tmp2,tmp3/])
       end if 
       delete(indx)
     end if 
     delete([/lat1,lat2,lon1,lon2/])
    end do 
   end do 

   out_path    = datdir+"/"+Groups(i)
   system("if ! test -d " + out_path +" ; then mkdir -p " + out_path + " ; fi")
   out_file   = Groups(i)+ "_ARs_annual_freq_" + year+".nc"
   file_ou = out_path+"/"+out_file
   system("/bin/rm -f " + file_ou)
   fou = addfile(file_ou,"cw")
   fou->AR_freq = vfreq
   fou->AR_ivtu = vtuq
   fou->AR_ivtv = vtvq
   fou->AR_ivts = vivt

   delete([/fna,fla,fnb,flb/])
   delete([/vara,var0,tuq,tvq,ivt,lat,lon/])
   delete([/vfreq,vtuq,vtvq,vivt/])

 end do 

end do 

end
