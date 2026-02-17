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

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/data/ERA5/monthly"
  exp_name  = (/"ERA5"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014
  time_tag  = "1979-2019"
  grid_tag  = "1x1"

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS      = -90.  ;negative for southern hemisphere
  latN      =  90.  ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

  load "./share_var_info.ncl"
  nvars = dimsizes(varList)

  do k  = 0,nexps-1,1

   fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+time_tag(k)+"_"+grid_tag+".nc")
   nens  = dimsizes(fl1)

   do ie = 0,nens-1,1

    print(fl1(ie))
    f1 = addfile(fl1(ie),"r")

    do j = 0,nvars-1,1

    vars   = varList(j)
    varunt = untList(j)
    facmod = facList(j)

    print("working on variable "+ vars)

    lat  = f1->lat
    lon  = f1->lon

    cd_date       = cd_calendar(f1->time, 0)
    year          = tointeger(cd_date(:,0))
    indy          = ind(year.ge.ystr.and.year.le.yend)
    obstime       = f1->time(indy)
    obstime@units = f1->time@units

    if(vars.eq."Z500") then
      var0 = f1->Z3(indy,{500},{latS:latN},{lonW:lonE})
      var0@_FillValue = -9999.
      var0       = var0*facmod
      var0@units = varunt
    else if (vars.eq."Z700") then
      var0 = f1->Z3(indy,{700},{latS:latN},{lonW:lonE})
      var0@_FillValue = -9999.
      var0       = var0*facmod
      var0@units = varunt
    else if (vars.eq."Z200") then
      var0 = f1->Z3(indy,{200},{latS:latN},{lonW:lonE})
      var0@_FillValue = -9999.
      var0       = var0*facmod
      var0@units = varunt
    else if(vars.eq."SST") then
      var0 = f1->SST(indy,{latS:latN},{lonW:lonE})
      var0@_FillValue = -9999.
      var0 = where(var0.le.273.15, var0@_FillValue, var0)
      var0       = var0 - 273.15
      var0@units = "~S~o~N~C"
    else
      var0 = f1->$vars$(indy,{latS:latN},{lonW:lonE})
      var0@_FillValue = -9999.
      var0       = var0*facmod
      var0@units = varunt
    end if
    end if
    end if
    end if

    var0@_FillValue = -9999.
    var0       = var0*facmod
    var0@units = varunt
    delete(var0&time)
    var0&time  = obstime
    printVarSummary(var0)
    printMinMax(var0, True)

    delete([/cd_date,year,indy/])

    ;==============================================================
    ; compute climatology and Anomalies
    ;==============================================================

    do i = 0, nseas-1,1

      print("working on season "+ seasons(i))

      if (seasons(i).eq."MON") then
        vb1   = var0(time|:,lat|:,lon|:)
        vc1   = var0(time|:,lat|:,lon|:)
        vd1   = var0(time|0:11,lat|:,lon|:)

        var  = vb1(:,:,:)
        Clm  = clmMonTLL(var) ; (12,lat,lon)
        Anom = calcMonAnomTLL(var, Clm)
        Clm!0         = "month"
        vd1!0         = "month"
        vd1&month     = Clm&month
        vc1(:,:,:) = (/Anom/)
        vd1(:,:,:) = (/Clm/)
        delete([/var,Clm,Anom/])

        printVarSummary(vd1)
        printMinMax(vd1, True)

      else if (seasons(i).eq."ANN") then

        vb1        = month_to_annual(var0(time|:,lat|:,lon|:), 1)
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(obstime)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp

        vd1  = dim_avg_n_Wrap(vb1,0)

        vc1  = vb1
        vc1  = vb1 - conform(vb1,vd1,(/1,2/))
        delete(ttmp)

      else

        vb1        = month_to_season(var0(time|:,lat|:,lon|:), seasons(i))
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

      setfileoption("nc", "Format",  "NetCDF4")
      enstr     = "en"+sprinti("%02d",ie)
      out_file  = exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
      system("rm " + out_file)
      fout = addfile(out_file,"cw")
      vs1 = vars+"_o"
      vs2 = vars+"_a"
      vs3 = vars+"_c"

      fout->$vs1$  = vb1
      fout->$vs2$  = vc1
      fout->$vs3$  = vd1
      delete([/fout,vb1,vc1,vd1/])

   end do 

  end do 

 end do 

end do 

end
