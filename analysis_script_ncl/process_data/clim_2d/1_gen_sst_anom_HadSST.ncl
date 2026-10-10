  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/data"
  exp_name  = (/"HadSST"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)

  latS      = -90.  ;negative for southern hemisphere
  latN      =  90.  ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

  ystr      = 1979
  yend      = 2014
  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"
  tobs_tag  = "1870-2020"

  vars      = "sst"
  vout      = "SST"
  vfac      = 1.0
  vunt      = "~S~o~N~C"

 ;declear the data to be used;;;;;

  do k  = 0,nexps-1,1

   fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/"+exp_name(k)+"*"+tobs_tag+"*.nc")
   nens = dimsizes(fl1)
  ;print(fl1)

   do ie = 0,nens-1,1

    print(fl1(ie))

    f1       = addfile(fl1(ie),"r")
    utc_date = cd_calendar(f1->time, 0)
    year     = tointeger(utc_date(:,0))
    indy     = ind(year.ge.ystr.and.year.le.yend)
    obstime  = f1->time(indy)

   ;printVarSummary(obstime)
   ;print(cd_calendar(obstime,-3))
    
    var1 = f1->$vars$(indy,{latS:latN},{lonW:lonE})
    var1@_FillValue = -9999.
    var1 = where(var1.le.-100.0,var1@_FillValue,var1)
    var1       = var1*vfac
    var1@units = vunt
    ttmp       = ispan(0,dimsizes(obstime)-1,1)
    ttmp@units = "months since "+ystr+"-01-01 00:00:00"
    var1&time  = ttmp

   ;printVarSummary(var1)
   ;printMinMax(var1, True)
    delete([/ttmp,utc_date,year,indy/])

    do i = 0, nseas-1,1

      print("working on season "+ seasons(i))
      ;==============================================================
      ; compute climatology and Anomalies
      ;==============================================================

      if (seasons(i).eq."MON") then
        vb1        = var1(time|:,lat|:,lon|:)
        vc1        = var1(time|:,lat|:,lon|:)
        vd1        = var1(time|0:11,lat|:,lon|:)
        var        = vb1(:,:,:)
        Clm        = clmMonTLL(var) ; (12,lat,lon)
        Anom       = calcMonAnomTLL(var, Clm)
        Clm!0      = "month"
        vd1!0      = "month"
        vd1&month  = Clm&month
        vc1(:,:,:) = (/Anom/)
        vd1(:,:,:) = (/Clm/)
        delete([/var,Clm,Anom/])
       ;printVarSummary(vd1)
       ;printMinMax(vd1, True)
      else if (seasons(i).eq."ANN") then
        vb1        = month_to_annual(var1(time|:,lat|:,lon|:), 1)
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(obstime)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp
        vd1        = dim_avg_n_Wrap(vb1,0)
        vc1        = vb1
        vc1        = vb1 - conform(vb1,vd1,(/1,2/))
        delete(ttmp)
      else
        vb1        = month_to_season(var1(time|:,lat|:,lon|:), seasons(i))
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(obstime)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp

        vd1  = dim_avg_n_Wrap(vb1,0)

        vc1  = vb1
        vc1  = vb1 - conform(vb1,vd1,(/1,2/))
        delete(ttmp)
      end if
      end if
 
     ;======================================================================
     ;Save the data for post analysis
     ;==================================================================
      setfileoption("nc", "Format",  "NetCDF4")
      enstr     = "en"+sprinti("%02d",ie)
      out_file  = exp_name(k)+"_anom_"+vout+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
      system("rm " + out_file)
      fout = addfile(out_file,"cw")
      vs1 = vout+"_o"
      vs2 = vout+"_a"
      vs3 = vout+"_c"

      fout->$vs1$  = vb1
      fout->$vs2$  = vc1
      fout->$vs3$  = vd1
      delete([/fout,vb1,vc1,vd1/])

   end do 

  end do ; end exps

 end do 

end
