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

  data_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output/"
  exp_name  = (/"ERA5"/)
  labList   = (/"ERA5"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014

  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  ytrend    = 2001
  nseg      = 3  ; 3 trends: 1979-2001, 1979-2014, 2001-2014
  seasons   = (/"SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS      = -90. ;negative for southern hemisphere
  latN      = -60. ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

  gravit    = 9.80616    ; acceleration of gravity ~ m/s^2
  rgrav     = 1.0/gravit ; reciprocal of gravit
  DUfac     = 2.1415e-5  ; 1 DU in kg[O3] per m^2
 ;1 Dobson Unit (DU) is:
 ;2.6867 x 10^20 molecules per meter square
 ;4.4615 x 10^-4 mol/m2
 ;2.1415 x 10^-5 kg[O3]/m2

  varList    = (/"TCO"/)
  vnmList    = (/"Total Column Ozone"/)
  facList    = (/1.0/)
  untList    = (/"DU"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 0,nvars-1,1 

     print("working on season "+ seasons)
     vnam = varList(j)
     vunt = untList(j)
     vfac = facList(j)
     vstr = vnmList(j)

     if(exp_name(k).eq."ERA5") then
       fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*ens00.nc")
     else
       fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*.nc")
     end if
     f1   = addfiles(fl1,"r")

     yymm = cd_calendar(f1[:]->time,-1)
     yms  = ystr+"01"
     yme  = yend+"12"
     indy = ind(yymm.ge.yms.and.yymm.le.yme)
     
     lat  = f1[0]->lat
     lon  = f1[0]->lon
     lev  = f1[0]->lev

     lev0       = lev
     lev0       = lev * 100.0
     lev0@units = "Pa"

     psfc       = 1018*100.
     psfc@units = "Pa"
     ptop       = min(lev0)
     pdel       = dpres_plevel(lev0, psfc, ptop, 0)

     tmp0       = f1[:]->O3(indy,:,:,:) 
     tmp0       = tmp0*conform(tmp0,pdel,1)*rgrav/DUfac
     var0       = dim_sum_n_Wrap(tmp0,1)
     var0@units = vunt
     var0@long_name = vstr 
     delete([/lev0,lev,psfc,ptop,pdel,yymm,yms,yme,indy/])

     rad  = 4.0*atan(1.0)/180.0
     clat = lat
     clat = cos(lat*rad)

     vars = month_to_seasonN(var0,seasons)
     varm = wgt_areaave_Wrap(vars(0,:,{latS:latN},{lonW:lonE}),clat({latS:latN}), 1.0, 1)
     
     out_file  = exp_name(k)+"_"+vnam+"_"+seasons+"_"+time_tag+".nc"
     system("rm " + out_file)
     fout = addfile(out_file,"cw")
     fout->$vnam$  = varm
     delete([/fout,lat,lon,clat,tmp0,vars,varm,var0/])

    end do

  end do 

end
