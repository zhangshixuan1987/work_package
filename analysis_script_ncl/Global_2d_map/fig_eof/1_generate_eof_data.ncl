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

  l_plot_box = True
  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim"
  obs_name   = "ERA5"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP"   , "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  varList    = (/"TS",       "PSL"/)
  varUnit    = (/"~S~o~N~C", "hPa"/)
  nvar       = dimsizes(varList)
  
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  10.
  lonW       =   0.
  lonE       = 360.

  do j = 0,nvar-1,1

  vars   = varList(j)
  vunt   = varUnit(j)
  print("working on " + vars)

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do k = 0,nexps-1

    fl0  = systemfunc("ls -1 "+data_dir+"/"+obs_name+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fm0  = addfile(flmsk,"r")
    f0   = addfiles(fl0,"r")
    f1   = addfiles(fl1,"r")
    ListSetType (f0, "join")
    ListSetType (f1, "join")

    vnm  = vars+"_o"
    vd0  = f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    vd1  = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
    lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})

    dim0 = dimsizes(vd0)
    dim1 = dimsizes(vd1)

    nens = dim0(0)*dim1(0)
    ntim = dim1(1)
    nlat = dim1(2)
    nlon = dim1(3)

    vdat = new((/nens,ntim,nlat,nlon/),typeof(vd0))
    vdat!0 = "ens"
    vdat!1 = "time"
    vdat!2 = "lat"
    vdat!3 = "lon"
    vdat&time = vd1&time
    vdat&lat  = vd1&lat
    vdat&lon  = vd1&lon
    do ii = 0,dim0(0)-1,1
      do jj = 0,dim1(0)-1,1
        vdat(ii*dim1(0)+jj,:,:,:) = vd1(jj,:,:,:) - vd0(ii,:,:,:)
      end do
    end do
    vdat@_FillValue = -9999.0

    lms0 = conform(vdat,lmsk,(/2,3/))
    vdat = where(lms0.le.0,vdat,vdat@_FillValue)
    delete([/dim0,dim1,nens,ntim,nlat,nlon,lmsk,lms0/])

    if(vars.eq."PSL") then
     vdat = vdat/100.0
     vdat@units = "hPa"
     vunt    = "hPa"
    end if 

    print("working on season "+ seasons(i))
    if (seasons(i).eq."ANN") then
      vdx  = vdat(time|:,ens|:,lat|:,lon|:)
      var1 = month_to_annual(vdx,1)
    else
      vdx  = vdat(time|:,ens|:,lat|:,lon|:)
      var1 = month_to_season(vdx,seasons(i))
    end if
    var1!0    = "time"
    var1!1    = "ens"
    var1!2    = "lat"
    var1!3    = "lon"
    var1&time = vdat&time(::12)

    dimx = dimsizes(var1)
    ntim = dimx(0)
    nens = dimx(1)
    nlat = dimx(2)
    nlon = dimx(3)

    ;=================================================================
    ;create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
    ;=================================================================
     rad    = get_d2r("float")
     clat   = var1&lat
     clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid

    ;=================================================================
    ;weight all observations
    ;=================================================================
     xw     = var1*conform(var1, tofloat(clat), 2)
     copy_VarMeta(var1, xw)

    ;=================================================================
    ; Compute EOFs & Standardize time series
    ;=================================================================
     neof    = 20 
     optEOF  = True
     eof     = eofunc_n_Wrap(xw, neof, optEOF, 0)
     eof!0   = "neof"
     eof!1   = "nens"
     eof!2   = "lat"
     eof!3   = "lon"

     eof_ts      = new((/neof,nens,ntim/),typeof(eof))
     do ie = 0,nens-1,1
      optETS  = False
      eof_ts(:,ie,:)  = (/eofunc_ts_n(xw(:,ie,:,:), eof(:,ie,:,:), optETS, 0)/)
     end do 
     eof_ts      = dim_standardize_n(eof_ts, 1, 2)      ; normalize
     eof_ts!0    = "neof"
     eof_ts!1    = "nens"
     eof_ts!2    = "time"
     eof_ts&time = var1&time
    ;printVarSummary(eof)
    ;printVarSummary(eof_ts)
    ;printVarSummary(var1)

    ; =================================================================
    ; Regress
    ; =================================================================
     eof_regres = eof                               ; create an array w meta data
     do ie = 0,neof-1
      do je = 0, nens-1,1
        eof_regres(ie,je,:,:) = (/regCoef_n(eof_ts(ie,je,:),var1(:,je,:,:),0,0)/)
      end do 
     end do 

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file = exp_name(k)+"_"+vars+"_EOF_"+seasons(i)+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
    fout->eof        = eof 
    fout->eof_ts     = eof_ts
    fout->eof_regres = eof_regres

    delete([/rad,clat,xw,eof,eof_ts,eof_regres/])
    delete([/vdx,vdat,var1,vd0,vd1/])
    delete([/fl0,f0,fl1,f1/])

    end do

 end do 

end do 

end

