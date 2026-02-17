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
  latN      =  90. ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

 ;define regions for eof analysis
  mode       = ("JET") 

  varList    = (/"SF"/)
  vnmList    = (/"SF"/)
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

     time = f[:]->time
     time = time -1
     yymm = cd_calendar(time,-1)
     yms  = ystr+"01"
     yme  = yend+"12"
     indy = ind(yymm.ge.yms.and.yymm.le.yme)
     lat  = f[0]->lat({latS:latN})
     rad  = 4.0*atan(1.0)/180.0
     wgt  = cos(lat*rad)
    
     lev       = f[0]->lev({50:1000}) 
     lev       = lev * 100.0
     lev!0     = "lev"
     lev&lev   = lev
     lev@units = "Pa"

     umi  = f[:]->U(indy,{50:1000},:,:)
     vmi  = f[:]->V(indy,{50:1000},:,:)
     psm  = f[:]->PS(indy,:,:)

     umi&time = time(indy)
     vmi&time = time(indy)
     psm&time = time(indy)
 
     vmi&lev = lev
     vmi&lev@units = "Pa"
     umi&lev = lev
     umi&lev@units = "Pa"

     if(psm@units.eq."hPa") then 
       psm = psm * 100.0
       psm@units = "Pa"
     end if 

     delete([/time,yymm,yms,yme,indy,lat/])
   
     ;;derive the 3D steam function
     sfvp  = uv2sfvpF (umi,vmi)
     sf    = umi
     vp    = umi
     sf    = (/sfvp(0,:,:,:,:)/)
     vp    = (/sfvp(1,:,:,:,:)/)

     sfzon = dim_avg_n_Wrap(sf,3)
     sfmer = dim_avg_wgt_n_Wrap(sf,wgt,1,2)
     vpzon = dim_avg_n_Wrap(vp,3)
     vpmer = dim_avg_wgt_n_Wrap(vp,wgt,1,2)
     smh   = zonal_mpsi_Wrap(vmi,vmi&lat,lev,psm)

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file  = exp_name(k)+"_stream_function_"+seasons+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")

    ;;;output the variables;;;;
    fout->SF     = smh
    fout->SF_zon = sfzon
    fout->SF_mer = sfmer
    fout->VP_zon = vpzon
    fout->VP_mer = vpmer

    delete([/umi,vmi,smh,psm,wgt,sfvp,sf,vp,sfzon,sfmer,vpzon,vpmer/])

 end do 

end do 

end
