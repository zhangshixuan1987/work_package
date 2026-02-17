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

 ;ERA5 -> RD
 ;CLIM -> CR
 ;Nudged -> NU
 ;LSTM(CLIM) -> ML(CR)

 ;;;;;;;strings for the data directory;;;;;;;;
  datadir  = "../fig_data/"
  Groups   = (/"ERA5", "CLIM30", "Nudged10", "CLIM30ML" /)
  Expnams  = (/"RD",   "CR",     "NU",       "ML(CR)"   /)
  ngrps    = dimsizes(Groups)
  syear    = 2007
  eyear    = 2014
 
;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U", "T", "Q"/)
  VstList  = (/"U", "T", "Q"/)
  VnmList  = (/"Zonal wind", "Temperature","Specific humidity"/)
  UnitList = (/"m s~S~-1~N~","K",           "g kg~S~-1~N~"    /)
  nvars    = dimsizes(VarList)

  Seasons  = (/"ANN","DJF","JJA"/) 
  nseas    = dimsizes(Seasons) 

  do is = 0,nseas-1,1

   do iv = 0,nvars-1, 1

    Varname  = VarList(iv)
    Varunt   = UnitList(iv)
    Varstr   = VstList(iv)
    Vnmstr   = VnmList(iv)
 
    do igp = 0,ngrps-1,1

    if(Groups(igp).eq."NDGUVTQ") then 
     fnb = systemfunc("ls " + datadir+ Varname + "_" + Groups(igp) + ".nc")
    else 
      fnb = new(eyear-syear+1,string)
      do year = syear,eyear,1
        fnb(year-syear) = systemfunc("ls " + datadir+ Varname + "_" + Groups(igp) + "_"+year+".nc")
      end do 
    end if 
    print(fna + " "+ fnb)
    exit

    flb   = addfiles(fnb, "r")
    utco  = cd_calendar(fla[:]->time,0)
    utcm  = cd_calendar(flb[:]->time,0)
    indo  = ind(utco(:,0).ge.syear.and.utco(:,0).le.eyear)
    indm  = ind(utcm(:,0).ge.syear.and.utcm(:,0).le.eyear)
    delete([/utco,utcm/])

    Varbias  = "d"+Varname
    if(Seasons(is).eq."ANN") then 
      xtmp  = month_to_annual(fla[:]->$Varname$(indo,:,:),1)
      ytmp  = month_to_annual(flb[:]->$Varname$(indm,:,:),1)
      ztmp  = month_to_annual(flb[:]->$Varbias$(indm,:,:),1)
    else
      xtmp  = month_to_season(fla[:]->$Varname$(indo,:,:),Seasons(is))
      ytmp  = month_to_season(flb[:]->$Varname$(indm,:,:),Seasons(is))
      ztmp  = month_to_season(flb[:]->$Varbias$(indm,:,:),Seasons(is))
    end if

    xtmp!0 = "time"
    ytmp!0 = "time"
    ztmp!0 = "time"
    xtmp@_FillValue = -99999
    ytmp@_FillValue = -99999
    ztmp@_FillValue = -99999
 
    if(Varname.eq."RELHUM") then 
      xtmp = xtmp*100.0
      ytmp = ytmp*100.0
      ztmp = ztmp*100.0
    end if 
 
    printMinMax(xtmp,False)
    printMinMax(ytmp,False)

    ;calculate the global/regional mean values
    latrg = flb[0]->lat
    rad   = 4.*atan(1.)/180.    ;
    wgty  = cos(latrg*rad)      ; cosine weights

    xAve  = dim_avg_Wrap (xtmp(lev|:,lat|:,time|:))     
    yAve  = dim_avg_Wrap (ytmp(lev|:,lat|:,time|:))

    xVar  = dim_variance_Wrap (xtmp(lev|:,lat|:,time|:)) ; calculate variances
    yVar  = dim_variance_Wrap (ytmp(lev|:,lat|:,time|:))

    sigr  = 0.05                                    ; critical sig lvl for r
    xEqv  = equiv_sample_size (xtmp(lev|:,lat|:,time|:), sigr,0)
    yEqv  = equiv_sample_size (ytmp(lev|:,lat|:,time|:), sigr,0)
    xN    = wgt_areaave (xEqv, 1.0, wgty, 0)    ; wgty could be gaussian weights 
    yN    = wgt_areaave (yEqv, 1.0, wgty, 0) 

    iflag = False                        ; population variance similar
    prob  = xAve
    prob  = (/ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)/)
    prob  = where(prob.lt.sigr,prob,prob@_FillValue )

    if(igp.eq.0) then
      zref = ztmp
    else
      ztmp = where(ismissing(zref),ztmp@_FillValue,ztmp) 
    end if 
    vdiff = dim_avg_Wrap (ztmp(lev|:,lat|:,time|:))  ;xAve

    ;normalize figures with obs rms at each level
    ostddev = dim_stddev_n_Wrap(dim_avg_n_Wrap(xtmp,0),1)
   ;printVarSummary(xtmp)
   ;printVarSummary(ostddev)
   ;printVarSummary(vdiff)
   ;vdiff   = vdiff/conform(vdiff,ostddev,(/0/))

    delete([/fnb,xtmp,ytmp,ztmp,xAve,yAve,xVar,yVar,sigr,xN,yN,iflag,indo,indm/])
   end do 
  end do 
 end do ;end of file loop

end
