  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; ================================================================================================
; compute sam indices (Antarctic Oscillation index)
; Based on the numerical definition of the SAM by Gong and Wang (1999), which is:
; SAM = P*40S -  P*65S
; where P*40S and P*65S are the normalized monthly zonal sea level pressur(SLP) at 40S and 65S
; ================================================================================================
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/polar_diag/pobs/"
  exp_name  = (/"Marshall"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)

  latS0     = -65.  ;negative for southern hemisphere
  latN0     = -40.  ;negative for southern hemisphere

  ;;;read mashall data;;;;
  fx1       = data_dir+"Marshall_mslp_monthly_"+toint(abs(latS0))+"S.txt"
  fx2       = data_dir+"Marshall_mslp_monthly_"+toint(abs(latN0))+"S.txt"
  ncol      = 13
  tslp65    = readAsciiTable(fx1, ncol, "float", (/1/))
  tslp40    = readAsciiTable(fx2, ncol, "float", (/1/))

  yymm      = yyyymm_time(min(toint(tslp65(:,0))), max(toint(tslp65(:,0))), "integer")
  nyear     = dimsizes(yymm)/12
  obstime   = ispan(0,dimsizes(yymm)-1,1)  
  obstime@units = "months since "+ min(toint(tslp65(:,0)))+"-01-01 00:00:00"
  time_tag  = min(toint(tslp65(:,0))) + "-" + max(toint(tslp65(:,0)))

  slp65      = ndtooned(tslp65(:,1:))
  slp40      = ndtooned(tslp40(:,1:))
  slp65!0    = "time"
  slp40!0    = "time"
  slp65&time = obstime
  slp40&time = obstime
  slp65@_FillValue = -9999.0
  slp40@_FillValue = -9999.0
  slp65@units      = "hPa"
  slp40@units      = "hPa"
  slp65@description= "mean sea level pressure at 65S"
  slp40@description= "mean sea level pressure at 40S"

  do k  = 0,nexps-1,1

    do i = 0, nseas-1,1 

      if (seasons(i).eq."MON") then
        var1 = slp40
        var2 = slp65
        var3 = slp65 
        var3@units = "std dev"
        var3@description="Standardized SAM index" 
        a1   = reshape(var1,(/nyear,12/))
        b1   = reshape(var2,(/nyear,12/))
        var3 = ndtooned(dim_standardize_n(a1, 0, 0) - dim_standardize_n(b1, 0, 0))
        var3 = dim_standardize_n(var3,0,0)
        delete([/a1,b1/])
      else if (seasons(i).eq."ANN") then
        var1 = month_to_annual(slp40,1)
        var2 = month_to_annual(slp65,1)
        time = ispan(0,dimsizes(var1)-1,1)
        time@units = "years since "+ min(toint(tslp65(:,0)))+"-01-01 00:00:00"
        var1!0     = "time"
        var2!0     = "time"
        var1&time  = time
        var2&time  = time
        var3       = var1
        var3@units = "std dev"
        var3@description="Standardized SAM index"
        var3 = dim_standardize_n(var1, 0, 0) - dim_standardize_n(var2, 0, 0)
        var3 = dim_standardize_n(var3,0,0)
        delete(time) 
      else
        var1   = month_to_season(slp40,seasons(i))
        var2   = month_to_season(slp65,seasons(i))
        time   = ispan(0,dimsizes(var1)-1,1)
        if(seasons(i).ne."DJF") then
          time@units = "years since "+ min(toint(tslp65(:,0)))+"-01-01 00:00:00"
        else
          time@units = "years since "+ min(toint(tslp65(:,0))-1)+"-01-01 00:00:00"
        end if
        var1!0     = "time"
        var2!0     = "time"
        var1&time  = time
        var2&time  = time
        var3       = var1
        var3@units = "std dev"
        var3@description="Standardized SAM index"
        var3 = dim_standardize_n(var1, 0, 0) - dim_standardize_n(var2, 0, 0)
        var3 = dim_standardize_n(var3,0,0) 
        delete(time)
      end if
      end if 
 
      ;======================================================================
      ;Save the data for post analysis 
      ;==================================================================
     out_file  = exp_name(k)+"_SAM_AOI_"+seasons(i)+"_"+time_tag+".nc"
     system("rm " + out_file)
     fout   = addfile(out_file,"c")
       
     fout->SLP40    = var1
     fout->SLP65    = var2
     fout->SAM_AOI  = var3
     delete([/var1,var2,var3/])

  end do ; end exps

 end do 

 delete([/fx1,fx2,ncol,tslp40,tslp65,slp40,slp65,yymm,obstime/])

end
