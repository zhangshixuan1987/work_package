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
  exp_name  = (/"CLIM_30yr","SSP245_30yr","SSP585_30yr"/)
  labList   = (/"AMIP","AMIP+SSP245","AMIP+SSP585"/)
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
  vnmList    = (/"Ozone"/)
  facList    = (/1.0/)
  untList    = (/"DU"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 0,nvars-1,1 

     var = varList(j)
     fl1 = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/monthly/"+"*.nc")
     f1  = addfiles(fl1,"r")
    ;ListSetType (f1, "join")

     lat  = f1[0]->lat
     lon  = f1[0]->lon
     rad  = 4.0*atan(1.0)/180.0
     clat = lat
     clat = cos(lat*rad)

     tmp0 = f1[:]->TCO 
     tmp0 = tmp0+ f1[:]->SCO
     vars = month_to_seasonN(tmp0,seasons)
     varm = wgt_areaave_Wrap(vars(0,:,{latS:latN},{lonW:lonE}),clat({latS:latN}), 1.0, 1)
     
     out_file  = exp_name(k)+"_"+var+"_"+seasons+"_"+time_tag+".nc"
     system("rm " + out_file)
     fout = addfile(out_file,"cw")
     fout->$var$  = varm
     delete([/fout,lat,lon,clat,tmp0,vars,varm/])

    end do

  end do 

end
