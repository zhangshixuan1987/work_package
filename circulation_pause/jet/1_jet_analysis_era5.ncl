  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output"
  exp_name  = (/"ERA5"/)
  labList   = (/"ERA5"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014
  time_tag  = ystr+"-"+yend
  out_dir   = "./"
  wgtfile   = "./1deg_to_0.1deg.nc"

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS      = -90. ;negative for southern hemisphere
  latN      = -20. ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

 ;define regions for eof analysis
  mode       = ("JET") 

  varList    = (/"U850"/)
  vnmList    = (/"U850"/)
  facList    = (/1.0/)
  untList    = (/"m/s"/)
  nvars      = dimsizes(varList)

  do i = 0,nvars-1,1

   print("working on season "+ seasons)
   vnam = varList(i)
   vunt = untList(i)
   vfac = facList(i)
   vstr = vnmList(i)

   do k = 0,nexps-1

     if(exp_name(k).eq."ERA5") then
       fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*ens00.nc")
     else
       fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*.nc")
     end if
     f    = addfiles(fl1,"r")

     tim0 = f[:]->time
     yymm = cd_calendar(tim0,-1)
     yms  = ystr+"01"
     yme  = yend+"12"
     indy = ind(yymm.ge.yms.and.yymm.le.yme)
     lat  = f[0]->lat({latS:latN})
     if(vnam.eq."U850") then
       tmp0 = f[:]->U(indy,{850},{latS:latN},{lonW:lonE})
     else
       tmp0 = f[:]->$vnam$(indy,{latS:latN},{lonW:lonE})
     end if
     vp1 = tmp0
     time = tim0(indy)
     vp1&time = time
     delete([/tim0,yymm,yms,yme,indy,lat,tmp0/])

     weightfile = systemfunc("ls "+wgtfile)
     exists = isfilepresent(weightfile)
     if(.not.exists) then
       print("OPeNDAP test unsuccessful.")
       print("Either the file doesn't exist, or NCL does")
       print("not have OPeNDAP cabilities on this system.")
       print("Generate weighting files on the fly")
       Opt                = True
       Opt@SrcFileName    = "src_SCRIP.nc"
       Opt@DstFileName    = "dst_SCRIP.nc"
       Opt@WgtFileName    = "1deg_to_0.1deg.nc"
       Opt@ForceOverwrite = True
       Opt@DstGridType    = "0.1x0.1"
       Opt@DstLLCorner    = (/ latS,   0.25d/)
       Opt@DstURCorner    = (/ latN, 359.75d/)
       vd1                = ESMF_regrid(vp1,Opt)
     else
       vd1 = ESMF_regrid_with_weights(vp1,weightfile,False)
     end if 
     printMinMax (vd1, False) 
     printVarSummary(vd1)

    x1       = dim_avg_n_Wrap(vd1,2)
    x_ins    = dim_max_n(x1,1)
    xmax_ind = dim_maxind(x1(lat|:,time|:), 0)
    x_pos    = x_ins
    ntim     = dimsizes(x_pos(:))

    if(.not.all(ismissing(xmax_ind)))then
      do iy = 0,ntim-1
        x_pos(iy) = tofloat(x1&lat(xmax_ind(iy)))
      end do
    else
      x_pos = -9999.
      x_pos@_FillValue = -9999.
    end if

    ;;calculate the anomalies;;;;
    xa_ins = x_ins
    xa_pos = x_pos

    if(seasons(i).eq."MON") then 
     do ii = 0,11
       xsub = x_ins(ii::12)
       xa_ins(ii::12) = x_ins(ii::12) - avg(xsub)
       delete([/xsub/])
       xsub = x_pos(ii::12)
       xa_pos(ii::12) = x_pos(ii::12) - avg(xsub)
       delete([/xsub/])
     end do
    else
     xsub = x_ins(:)
     xa_ins(:)  = xsub - avg(xsub)
     xsub = x_pos(:)
     xa_pos(:)  = xsub - avg(xsub)
     delete([/xsub/])
    end if 
    delete([/x1,xmax_ind/])

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file  = exp_name(k)+"_Jet_ind_"+vnam+"_"+seasons+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")

    x_ins@description =  "jet strength (original value)"
    x_pos@description =  "jet position (original value)"
    xa_ins@description = "jet strength (anomaly value)"
    xa_pos@description = "jet position (anomaly value)"

    xa_ins!0    = "time" 
    x_ins!0     = "time"
    xa_pos!0    = "time"
    x_pos!0     = "time"

    x_ins&time  = time
    x_pos&time  = time
    xa_ins&time = time
    xa_pos&time = time

    varstr = "jeto_ins"
    fout->$varstr$ = x_ins(time|:)

    varstr = "jeto_pos"
    fout->$varstr$ = x_pos(time|:)

    varstr = "jetd_ins"
    fout->$varstr$ = xa_ins(time|:)

    varstr = "jetd_pos"
    fout->$varstr$ = xa_pos(time|:)

    delete([/time,vp1,vd1,x_ins,x_pos,xa_ins,xa_pos/])

 end do 

end do 

end
