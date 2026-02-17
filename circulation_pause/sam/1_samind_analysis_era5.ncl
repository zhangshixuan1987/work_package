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

  data_dir  = "/global/cscratch1/sd/zhan391/DARPA_project/post_process/data/model_output"
  exp_name  = (/"ERA5"/)
  labList   = (/"ERA5"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014
  time_tag  = ystr+"-"+yend
  out_dir   = "./"

  ytrend    = 2001
  nseg      = 3  ; 3 trends: 1979-2001, 1979-2014, 2001-2014
  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS      = -90. ;negative for southern hemisphere
  latN      = -20. ;negative for southern hemisphere
  lonW      =  0.0
  lonE      = 360.0

 ;define regions for eof analysis
  mode      = (/"SAM"/) 
  nmod      = dimsizes(mode)

  varList    = (/"PSL"/)
  vnmList    = (/"SAMIND"/) 
  facList    = (/1.0, 1.0 /) 
  untList    = (/"Pa / std dev"/) 
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


     lat      = f[0]->lat({latS:latN})
     rad      = get_d2r("float")
     clat     = cos(rad*lat)
     clat!0   = "lat"
     clat&lat = lat

     time = f[:]->time
     time = time -1 
     yymm = cd_calendar(time,-1)
     yms  = ystr+"01"
     yme  = yend+"12"
     indy = ind(yymm.ge.yms.and.yymm.le.yme)
     tmp0 = f[:]->$vnam$(indy,{latS:latN},{lonW:lonE})
     tmp0&time = time(indy)
     varm = rmMonAnnCycTLL(tmp0) 
     v65S = wgt_areaave_Wrap(tmp0(:,{-65:-65},:),clat({-65:-65}),1.0,0)
     v40S = wgt_areaave_Wrap(tmp0(:,{-40:-40},:),clat({-40:-40}),1.0,0) 
     delete([/yymm,yms,yme,indy,tmp0,time,clat/])

    ;printMinMax(varm, False) 
    ;printVarSummary(varm)

    ;=================================================================
    ; Compute Standardize time series
    ;=================================================================
    evts  = dim_standardize_n_Wrap(v65S,0,0) 
    evts  = evts - dim_standardize_n_Wrap(v40S,0,0) 
    vreg  = varm(0,:,:)
    vreg  = (/regCoef_n(dim_standardize_n_Wrap(evts,0,0),varm,0,0)/)
    printVarSummary(vreg)

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file  = out_dir+exp_name(k)+"_"+vstr+"_"+seasons+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
  
    vout = mode+"_pcts"
    fout->PSL65    = v65S
    fout->PSL40    = v40S
    fout->$vnam$   = vreg
    fout->$vout$   = evts
    delete([/vreg,evts/])

   delete(fl1)

  end do    

 end do 

end

