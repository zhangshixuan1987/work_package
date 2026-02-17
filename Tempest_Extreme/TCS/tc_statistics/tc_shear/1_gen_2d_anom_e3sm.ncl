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

  data_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/nudged_e3smv2_simulation"
  exp_name  = (/"v2.LR.amip"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014
  time_tag  = "1979-2014"
  grid_tag  = "1x1"

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =  90.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  load "./share_var_info.ncl"
  nvars = dimsizes(varList)

  do k  = 0,nexps-1,1

   fl1        = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/"+exp_name(k)+"*"+time_tag(k)+".nc")
   nens  = dimsizes(fl1)

   do ie = 0,nens-1,1

    print(fl1(ie))
    f1 = addfile(fl1(ie),"r")

    do j = 8,10 ;0,nvars-1,1

    vars   = varList(j)
    varunt = untList(j)
    facmod = facList(j)

    print("working on variable "+ vars)

    lat  = f1->lat
    lon  = f1->lon

    time = f1->time
    time = time - 1
    cd_date       = cd_calendar(time, 0)
    year          = tointeger(cd_date(:,0))
    indy          = ind(year.ge.ystr.and.year.le.yend)
    obstime       = time(indy)
    obstime@units = time@units

    if(vars.eq."Z500") then
      xtmp = vinth2p(f1->Z3,f1->hyam,f1->hybm,500.0,f1->PS,1,f1->P0/100.0,1,False)
      var0   = xtmp(indy,0,{latS:latN},{lonW:lonE})
    else if (vars.eq."Z700") then
      xtmp = vinth2p(f1->Z3,f1->hyam,f1->hybm,700.0,f1->PS,1,f1->P0/100.0,1,False)
      var0   = xtmp(indy,0,{latS:latN},{lonW:lonE})
    else if (vars.eq."Z200") then
      xtmp = vinth2p(f1->Z3,f1->hyam,f1->hybm,200.0,f1->PS,1,f1->P0/100.0,1,False)
      var0 = xtmp(indy,0,{latS:latN},{lonW:lonE})
    else if (vars.eq."PRECT") then
      xtmp = f1->PRECC
      xtmp = xtmp + f1->PRECL
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."PRECST") then
      xtmp = f1->PRECSC
      xtmp = xtmp + f1->PRECSL
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."SWCF_SRF") then
      xtmp = f1->FSNS
      xtmp = xtmp - f1->FSNSC
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."LWCF_SRF") then
      xtmp = f1->FLNSC
      xtmp = xtmp - f1->FLNS
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."NETCF_SRF") then
      xtmp = f1->FLNSC
      xtmp = xtmp - f1->FLNS
      xtmp = xtmp + f1->FSNS
      xtmp = xtmp - f1->FSNSC
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."NETCF") then
      xtmp = f1->LWCF
      xtmp = xtmp + f1->SWCF
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."RESSURF") then
      xtmp = f1->FSNS
      xtmp = xtmp - f1->FLNS
      xtmp = xtmp - f1->SHFLX
      xtmp = xtmp - f1->LHFLX
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."RESTOA") then
      xtmp = f1->FSNTOA
      xtmp = f1->FSNTOA - f1->FLUT
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."NET_RAD_FLUX_SRF") then
      xtmp = f1->FSNS
      xtmp = xtmp - f1->FLNS
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."NET_FLUX_SRF") then
      xtmp = f1->FSNS
      xtmp = xtmp - f1->FLNS - f1->SHFLX - f1->LHFLX
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."S10") then
      xtmp = f1->U10
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."U10") then
      nlev = dimsizes(f1->lev)
      xtmp = f1->U(:,nlev-1,:,:)
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else if (vars.eq."V10") then
      nlev = dimsizes(f1->lev)
      xtmp = f1->V(:,nlev-1,:,:)
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    else
      xtmp = f1->$vars$
      var0 = xtmp(:,{latS:latN},{lonW:lonE})
    end if 
    end if 
    end if 
    end if 
    end if 
    end if
    end if
    end if
    end if
    end if
    end if
    end if
    end if
    end if
    end if
    end if 

    var0@_FillValue = -9999.
    var0       = var0*facmod
    var0@units = varunt
    delete(var0&time)
    var0&time  = obstime
    delete([/time,cd_date,year,indy,xtmp/])

    printVarSummary(var0)
    printMinMax(var0, True)

   ;print(cd_calendar(time,-3))


    ;==============================================================
    ; compute climatology and Anomalies
    ;==============================================================

    do i = 0, nseas-1,1

      print("working on season "+ seasons(i))

      if (seasons(i).eq."MON") then
        vb1   = var0(time|:,lat|:,lon|:)
        vc1   = var0(time|:,lat|:,lon|:)
        vd1   = var0(time|0:11,lat|:,lon|:)
        var   = vb1(:,:,:)
        Clm   = clmMonTLL(var) ; (12,lat,lon)
        Anom  = calcMonAnomTLL(var, Clm)
        Clm!0         = "month"
        vd1!0         = "month"
        vd1&month     = Clm&month
        vc1(:,:,:) = (/Anom/)
        vd1(:,:,:) = (/Clm/)
        delete([/var,Clm,Anom/])
      else if (seasons(i).eq."ANN") then
        vb1        = month_to_annual(var0(time|:,lat|:,lon|:), 1)
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(obstime)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp
        vd1        = dim_avg_n_Wrap(vb1,0)
        vc1        = vb1
        vc1        = vb1 - conform(vb1,vd1,(/1,2/))
        delete(ttmp)
      else
        vb1        = month_to_season(var0(time|:,lat|:,lon|:), seasons(i))
        vb1!0      = "time"
        ttmp       = ispan(0,dimsizes(obstime)/12-1,1)
        ttmp@units = "years since "+ystr+"-01-01 00:00:00"
        vb1&time   = ttmp
        vd1        = dim_avg_n_Wrap(vb1,0)
        vc1        = vb1
        vc1        = vb1 - conform(vb1,vd1,(/1,2/))
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
