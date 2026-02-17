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

  varList    = (/"PSL","Z500"/)
  vnmList    = (/"PSL","Z500"/) 
  facList    = (/1.0, 1.0 /) 
  untList    = (/"Pa / std dev","m/ std dev"/) 
  nvars     = dimsizes(varList)

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
     yymm = cd_calendar(time,-1)
     yms  = ystr+"01"
     yme  = yend+"12"
     indy = ind(yymm.ge.yms.and.yymm.le.yme)

     lat  = f[0]->lat({latS:latN})
     if(vnam.eq."Z500") then 
       tmp0 = f[:]->Z3(indy,{500},{latS:latN},{lonW:lonE})
     else 
       tmp0 = f[:]->$vnam$(indy,{latS:latN},{lonW:lonE})
     end if 
     varm = rmMonAnnCycTLL(tmp0)
     tmp0&time = time(indy)
     delete([/yymm,yms,yme,indy,tmp0,time/])

    ;printMinMax(varm, False) 
    ;printVarSummary(varm)

     dimx  = dimsizes(varm)
     ntime = dimx(0)
     nlat  = dimx(1)
     nlon  = dimx(2)

    ;=================================================================
    ; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
    ;=================================================================
     rad    = get_d2r("float")
     clat   = varm&lat            
     clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid
    ;printVarSummary(clat)
   
    ;=================================================================
    ;weight all observations 
    ;=================================================================
     xw     = varm*conform(varm, tofloat(clat), 1)
     copy_VarMeta(varm, xw)
     xw@long_name = "Wgt: "+vnam

    ;=================================================================
    ; Compute EOFs & Standardize time series
    ;=================================================================
    evcor = eofunc_n_Wrap(xw, 4, 75, 0)
    evcts = eofunc_ts_n_Wrap(xw,evcor,False,0)
    evts  = dim_standardize_n_Wrap(evcts(0,:),0,0)  
    vreg  = varm(0,:,:)
    vreg  = (/regCoef_n(evts,varm,0,0)/)
    if (avg(vreg({-85:-80},{5:10})).ge.0) then
      vreg = vreg * -1.0
      evts = evts * -1.0
    end if
    copy_VarAtts(evcor,evts)
    copy_VarAtts(evcor,vreg)
    printVarSummary(vreg)
    printVarSummary(evcor)
    delete([/evcor,evcts,varm/])

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file  = out_dir+exp_name(k)+"_EOF_PCS_"+vnam+"_"+seasons+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
  
    vout = mode+"_pcts"
    fout->$vnam$   = vreg
    fout->$vout$   = evts
    delete([/vreg,evts/])

   delete(fl1)

  end do    

 end do 

end

