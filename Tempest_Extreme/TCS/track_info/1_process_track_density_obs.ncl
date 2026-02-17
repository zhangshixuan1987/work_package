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


  outdir  = "./"
  datadir = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/OBS/tc_track"
  grid    = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/landmsk_coord/landmsk_coord.nc"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups  = (/"OBS"/)
  Casdir  = (/"IBTrACS"/)
  ngrps   = dimsizes(Groups)
  sset    = "set1"
  ystr    = 1979
  yend    = 2017
  timtag  = ystr+"-"+yend
  nyr     = yend - ystr + 1
  nmon    = 12
  knot2ms = 1/1.94384
  
  ;construct 4x4 global grid to bin storms 
  nlon    = toint(360./4)
  nlat    = toint(180./4)
  latm    = latGlobeFo(nlat,"lat","latitude","degrees_north")
  lonm    = lonGlobeFo(nlon,"lon","longitude","degrees_east")

  metric_vars = (/"trkdens", "tgedens", "tlydens", "tcint", "tcpmin", \
                  "tcltime", "tcspd", "tcace", "tgrate", "tdrate", "tzsfc"/)
  nmetrics    = dimsizes(metric_vars)
  
  do i = 0,ngrps-1,1

   fna     = datadir+"/"+Casdir(i)+"/"+Casdir(i)+".ALL.v04r00.nc"
   fla     = addfile(fna,"r")
   nstorm  = fla->numobs
   xsid    = fla->sid
   xspd    = fla->storm_speed
   xtime   = fla->time
   xbasin  = fla->basin 
   xlat    = fla->lat
   xlon    = fla->lon
   xsfc    = fla->dist2land
   xvmax   = fla->wmo_wind
   xpmin   = fla->wmo_pres
   xsr34   = fla->usa_wind

  ;create temporary array to save diagnostics
   tmpdat   = new((/nmetrics,nyr*nmon,nlat,nlon/),float)
   tmpdat!0 = "metric"
   tmpdat!1 = "time"
   tmpdat!2 = "lat"
   tmpdat!3 = "lon"
   tmpdat&time = yyyymm_time(ystr, yend, "integer")
   tmpdat&lat  = latm
   tmpdat&lon  = lonm
   tmpcount    = new((/nmetrics,nlat,nlon/),float)

   tmpdat = 0.0 
   sdate  = cd_calendar(xtime(:,0),-1)

   do iy = ystr, yend, 1

    do im = 0,nmon-1,1

     itim    = (iy-ystr)*nmon + im
     timestr = sprinti("%04d",iy)+sprinti("%02d",im+1)
     inds    = ind(sdate.eq.timestr)

     if (.not.all(ismissing(inds))) then 

      nst = num(.not.ismissing(inds))

      tmpcount = 0.0      

      do j = 0,nst-1,1

       ih    = inds(j)

       print("working on storm: " + chartostring(xsid(ih,:)))

       indx  = ind((.not.ismissing(xtime(ih,:))).and.(.not.ismissing(xpmin(ih,:))))
       ntrk  = num(.not.ismissing(indx))

       if(.not.all(ismissing(indx))) then 

        slon    = xlon(ih,indx) 
        slat    = xlat(ih,indx) 
        pmin    = short2flt(xpmin(ih,indx))
        vmax    = short2flt(xvmax(ih,indx))
        zsfc    = short2flt(xsfc(ih,indx))
        sspd    = short2flt(xspd(ih,indx))

        vmax    = vmax*knot2ms 
        vmax@units = "m/s"

        units   = "days since "+ystr+"-01-01 00:00:00"
        time    = cd_convert(xtime(ih,indx), "days since 1800-01-01 00:00" )
        printVarSummary(time)
        printVarSummary(vmax)

        ;life time 
        slft    = vmax
        slft    = tofloat(time(ntrk-1)-time(0))
        jj      = maxind(vmax) 
        grat = vmax
        if(.not.ismissing(jj).and.(jj.ne.0)) then
          grat = tofloat((vmax(jj) - vmax(0)) / (time(jj) - time(0)))
        else
          grat = 0 
        end if 

        drat = vmax
        if(.not.ismissing(jj).and.(ntrk-1.ne.jj)) then 
          drat = tofloat((vmax(ntrk-1) - vmax(jj)) / (time(ntrk-1) - time(jj)))
        else
          drat = 0 
        end if 
        ace  = vmax 
        kk   = ind(vmax.ge.(34.0*knot2ms))
        if(.not.all(ismissing(kk)))
          ace = sum(vmax(kk)*vmax(kk))
        else
          ace = 0.0   
        end if 
        delete([/jj,kk/])
 
        ;bin storm tracks        
        slon  = where(slon.lt.0.0,slon+360.0,slon)
        gbin  = new ( (/nlat,nlon/), float ) 
        gknt  = new ( (/nlat,nlon/), integer)

        do jj = 0,dimsizes(metric_vars)-1
          gbin  = new ( (/nlat,nlon/), float )
          gknt  = new ( (/nlat,nlon/), integer)
          gbin  = 0.0                      ; initialization
          gknt  = 0
          if(metric_vars(jj) .eq. "trkdens") then  
            bin_sum(gbin,gknt,lonm,latm,slon,slat,vmax)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + 0 ;gknt 
            printVarSummary(gknt)
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gknt 
          end if 
          if(metric_vars(jj) .eq. "tgedens") then
            bin_sum(gbin,gknt,lonm,latm,slon(0),slat(0),vmax(0))
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + 0 ;gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gknt
          end if
          if(metric_vars(jj) .eq. "tlydens") then
            bin_sum(gbin,gknt,lonm,latm,slon(ntrk-1),slat(ntrk-1),vmax(ntrk-1))
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + 0 ;gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gknt
          end if
          if(metric_vars(jj) .eq. "tcint") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,vmax)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tcpmin") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,pmin)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tcltime") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,slft)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tcspd") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,sspd)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tcace") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,ace)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tgrate") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,grat)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tdrate") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,drat)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          if(metric_vars(jj) .eq. "tzsfc") then
            bin_sum(gbin,gknt,lonm,latm,slon,slat,zsfc)
            tmpcount(jj,:,:)    = tmpcount(jj,:,:)    + gknt
            tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) + gbin
          end if
          delete([/gbin,gknt/])
        end do  
        delete([/vmax,pmin,slft,sspd,ace,grat,drat,slat,slon/])
        delete([/zsfc,units,time/])
       end if 
       delete([/indx/])
      end do ; end of ntrk loop

      ;derive the mean metrics 
      do jj = 0,dimsizes(metric_vars)-1
        if(metric_vars(jj) .eq. "trkdens") then
          tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) ;/ 4.0 / days_in_month (iy,(im+1))
        else if(metric_vars(jj) .eq. "tgedens") then
          tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) ;/ 4.0 / days_in_month (iy,(im+1))
        else if(metric_vars(jj) .eq. "tlydens") then
          tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) ;/ 4.0 / days_in_month (iy,(im+1))
        else 
          gknt  = tmpcount(jj,:,:)
          if(.not.isatt(gknt,"_FillValue"))
            gknt@_FillValue = -9999
          end if
          gknt = where(gknt.eq.0 , gknt@_FillValue, gknt)
          tmpdat(jj,itim,:,:) = tmpdat(jj,itim,:,:) / gknt(:,:)
          delete([/gknt/])
        end if
        end if 
        end if
      end do 
     end if
     delete([/nst,itim,timestr,inds/])       
    end do ; im loop
 
   end do  ; iy loop

   ;;;;open netcdf file and save data for the pfigures;;;;;;;
   ;setfileoption("nc", "Format",  "NetCDF4")
   out_file = outdir+Groups(i)+"_"+sset+"_Track_statistics_"+timtag+".nc"
   system("rm " + out_file)
   fout = addfile(out_file,"cw")
  ;fout->lat  = latm
  ;fout->lon  = lonm
  ;fout->time = yyyymm_time(ystr, yend, "integer")
   printVarSummary(tmpdat)

   do jj = 0,dimsizes(metric_vars)-1
      
     vout = tmpdat(jj,:,:,:) 
     vout!0 = "time"
     vout!1 = "lat"
     vout!2 = "lon"
     vout&time = yyyymm_time(ystr, yend, "integer")
     vout&lat  = latm
     vout&lon  = lonm
     if(metric_vars(jj) .eq. "trkdens") then
       vout@long_name  = "track density (number density of storm tracks passing through a region)"
       vout@short_name = "tdens"
       vout@units     = "transits per unit area equivalent to a 4-degree spherical cap"
     end if
     if(metric_vars(jj) .eq. "tgedens") then
       vout@long_name  = "genesis density (number density of starting points)"
       vout@short_name = "gdens" 
       vout@units     = "transits per unit area equivalent to a 4-degree spherical cap"
     end if
     if(metric_vars(jj) .eq. "tlydens") then
       vout@long_name  = "lysis density (number density of ending points)"
       vout@short_name = "ldens"
       vout@units     = "transits per unit area equivalent to a 4-degree spherical cap"
     end if
     if(metric_vars(jj) .eq. "tcint") then
       vout@long_name  = "mean intensity (vmax)"
       vout@short_name = "int"
       vout@units     = "m/s"
     end if
     if(metric_vars(jj) .eq. "tcpmin") then
       vout@long_name  = "mean minimum seal level pressure"
       vout@short_name = "pmin"
       vout@units     = "m/s"
     end if
     if(metric_vars(jj) .eq. "tcltime") then
       vout@long_name  = "mean lifetime" 
       vout@short_name = "lft"
       vout@units     = "days"
     end if
     if(metric_vars(jj) .eq. "tcspd") then
       vout@long_name  = "mean track speed"
       vout@short_name = "mspd"
       vout@units     = "m/s"
     end if
     if(metric_vars(jj) .eq. "tcace") then
       vout@long_name  = "mean Accumulated Cyclone Energy (ACE)"
       vout@short_name = "ace"
       vout@units     = "m2/s2"
     end if
     if(metric_vars(jj) .eq. "tgrate") then
       vout@long_name  = "mean growth rate"
       vout@short_name = "growth"
       vout@units     = "m/s/day"
     end if
     if(metric_vars(jj) .eq. "tdrate") then
       vout@long_name  = "mean decaying rate"
       vout@short_name = "decay"
       vout@units     = "m/s/day"
     end if
     if(metric_vars(jj) .eq. "tzsfc") then
       vout@long_name  = "surfce terrian height"
       vout@short_name = "zsfc"
       vout@units     = "m"
     end if
     fout->$metric_vars(jj)$ = vout 
     delete(vout)
   end do 

   delete([/tmpdat,fla,fna/])
   delete([/xsid,nstorm,xspd,xtime,xbasin,xlat,xlon,xsfc,xvmax,xvmax,xpmin,xsr34/])

 end do 
end
