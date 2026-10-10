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

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/eof_analysis"
  exp_name  = (/"v1.LR.historical"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)

  latS0     = -65.  ;negative for southern hemisphere
  latN0     = -40.  ;negative for southern hemisphere

  ystr      = 1979
  yend      = 2014
  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  vars      = "PSL"
  vfac      = 1.0/100.0 
  vunt      = "hPa"

 ;declear the data to be used;;;;;

  do k  = 0,nexps-1,1

   do i = 0, nseas-1,1

    fl1        = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_en*"+time_tag(k)+".nc")
    f1         = addfiles(fl1,"r")
    ListSetType (f1, "join")
    print(fl1)

    obstime = cd_convert(f1[0]->time,"months since "+ystr+"-01-01 00:00:00") 

    vnm     = vars+"_o"
    slp65   = dim_avg_n(f1[:]->$vnm$(:,:,{latS0},:),2)
    slp40   = dim_avg_n(f1[:]->$vnm$(:,:,{latN0},:),2)

    slp65   = slp65*vfac
    slp40   = slp40*vfac

    slp65!0 = "ens"
    slp40!0 = "ens"

    slp65!1 = "time"
    slp40!1 = "time"

    slp65&time        = obstime
    slp40&time        = obstime

    slp65@_FillValue  = -9999.0
    slp40@_FillValue  = -9999.0

    slp65@units       = vunt
    slp40@units       = vunt

    slp65@description = "mean sea level pressure at 65S"
    slp40@description = "mean sea level pressure at 40S"

    nens  = dimsizes(slp65(:,0))
    nyear = dimsizes(slp65(0,:))/12

    if (seasons(i).eq."MON") then
      var1 = slp40
      var2 = slp65
      var3 = slp65
      var3@units       = "std dev"
      var3@description = "Standardized SAM index"
      a1   = reshape(var1,(/nens,nyear,12/))
      b1   = reshape(var2,(/nens,nyear,12/))
      c1   = dim_standardize_n(a1, 0, 1) - dim_standardize_n(b1, 0, 1)
      var3 = reshape(c1,(/nens,nyear*12/))
      var3 = dim_standardize_n(var3,0,1) 
      delete([/a1,b1,c1/])
    else if (seasons(i).eq."ANN") then
      tmp1 = month_to_annual(slp40(time|:,ens|:),1) 
      tmp2 = month_to_annual(slp65(time|:,ens|:),1)
      time = ispan(0,dimsizes(tmp1(:,0))-1,1)
      time@units = "years since "+ystr+"-01-01 00:00:00"
      tmp1!0     = "time"
      tmp2!0     = "time"
      tmp1&time  = time
      tmp2&time  = time
      var1       = tmp1(ens|:,time|:)
      var2       = tmp2(ens|:,time|:)
      var3       = var1
      var3@units = "std dev"
      var3@description="Standardized SAM index"
      var3 = dim_standardize_n(var1, 0, 1) - dim_standardize_n(var2, 0, 1)
      var3 = dim_standardize_n(var3,0,1)
      delete([/time,tmp1,tmp2/])
    else
      tmp1 = slp40(time|::12,ens|:)
      tmp2 = slp65(time|::12,ens|:)
      do ie = 0,nens-1,1
        tmp1(:,ie) = month_to_season(slp40(ie,:),seasons(i))
        tmp2(:,ie) = month_to_season(slp65(ie,:),seasons(i))
      end do 
      time = ispan(0,dimsizes(tmp1(:,0))-1,1)
      if(seasons(i).ne."DJF") then 
        time@units = "years since "+ystr+"-01-01 00:00:00"
      else
        time@units = "years since "+(ystr-1)+"-01-01 00:00:00"
      end if 
      tmp1!0     = "time"
      tmp2!0     = "time"
      tmp1&time  = time
      tmp2&time  = time
      var1       = tmp1(ens|:,time|:)
      var2       = tmp2(ens|:,time|:)
      var3       = var1
      var3@units = "std dev"
      var3@description="Standardized SAM index"
      var3 = dim_standardize_n(var1, 0, 1) - dim_standardize_n(var2, 0, 1)
      var3 = dim_standardize_n(var3,0,1)
      delete([/time,tmp1,tmp2/])
    end if
    end if

   ;======================================================================
   ;Save the data for post analysis
   ;==================================================================
    out_file  = exp_name(k)+"_SAM_AOI_"+seasons(i)+"_"+ystr+"-"+yend+".nc"
    system("rm " + out_file)
    fout   = addfile(out_file,"c")

    fout->SLP40    = var1
    fout->SLP65    = var2
    fout->SAM_AOI  = var3
    delete([/var1,var2,var3,obstime,slp65,slp40/])

   end do ; end exps

 end do 

end
