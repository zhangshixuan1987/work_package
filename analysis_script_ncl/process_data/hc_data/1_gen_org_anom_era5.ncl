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

  intyp     = 2

  modvars   = (/"U","V", "OMEGA", "T","PS"/)
  facmod    = (/  1,  1, 100,   1,  1/)
  varunt    = (/"m~S~-1~N~","m~S~-1~N~","hPa s~S~-1~N~","K","Pa","W m~S~-2~N~"/)
  nvars     = dimsizes(modvars)

  do k  = 0,nexps-1,1

   fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+time_tag(k)+"_"+grid_tag+".nc")
   nens  = dimsizes(fl1)

   do ie = 0,nens-1,1

    print(fl1(ie))
    f1 = addfile(fl1(ie),"r")

    cd_date       = cd_calendar(f1->time, 0)
    year          = tointeger(cd_date(:,0))
    indy          = ind(year.ge.ystr.and.year.le.yend)
    obstime       = f1->time(indy)
    obstime@units = f1->time@units

    do i = 0, nseas-1,1

    ;======================================================================
    ;Save the data for analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4Classic")
    enstr      = "en"+sprinti("%02d",ie)
    out_file1  = exp_name(k)+"_HC_analysis_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
    system("rm " + out_file1)
    fout1      = addfile(out_file1,"cw")

    ;;reference pressure levels;;;;
    pinterp           = asciiread("./pres_lev.txt",-1,"float")
    pinterp@longname  = "lev"
    pinterp@units     = "hPa"
    pinterp!0         = "lev"
    pinterp&lev       = pinterp
    pinterp1          = pinterp({10:1000})
    lev               = f1->lev

    ;=================================================================
     print ("Start to process the cross-section of  selected variables ")
     um        = f1->$modvars(0)$(indy,:,:,:)
     vm        = f1->$modvars(1)$(indy,:,:,:)
     psm       = f1->$modvars(4)$(indy,:,:)

     umi       = int2p_n_Wrap(lev,um(:,:,:,:),pinterp1,intyp,1)
     vmi       = int2p_n_Wrap(lev,vm(:,:,:,:),pinterp1,intyp,1)
     delete([/um,vm/])

     ;;;derive the zonal mean for hadley cell;;;;
     umh       = dim_avg_n_Wrap(umi,3)
     vmh       = dim_avg_n_Wrap(vmi,3)

    ;printVarSummary(umh)
    ;printVarSummary(vmh)

     ;;;output the variables;;;;
     fout1->U   = umh
     fout1->V   = vmh
     delete([/umh,vmh/])

     ;;derive the 3D steam function
     lat   = f1->lat
     rad   = 4.0*atan(1.0)/180.0
     wgt   = cos(lat*rad)
     sfvp  = uv2sfvpF (umi,vmi)
     sf    = umi
     vp    = umi
     sf    = (/sfvp(0,:,:,:,:)/)
     vp    = (/sfvp(1,:,:,:,:)/)
     sfzon = dim_avg_n_Wrap(sf,3)
     sfmer = dim_avg_wgt_n(sf,wgt,1,2)
     vpzon = dim_avg_n_Wrap(vp,3)
     vpmer = dim_avg_wgt_n(vp,wgt,1,2)

     ;;;output the variables;;;;
     fout1->SF_zon  = sfzon
     fout1->SF_mer  = sfmer
     fout1->VP_zon  = vpzon
     fout1->VP_mer  = vpmer
     delete([/lat,rad,wgt,sfvp,sf,vp,sfzon,sfmer,vpzon,vpmer/])

     ;;derive the stream function
     levnew     = pinterp1*100.0
     levnew!0   = "lev"
     levnew&lev = levnew
     levnew@units = "Pa"

     vmi&lev    = levnew
     vmi&lev@units = "Pa"

     print(levnew) 
     printVarSummary(vmi)  
     smh        = zonal_mpsi_Wrap(vmi,vmi&lat,levnew,psm)
     smh&lev    = pinterp1
     smh&lev@units = "hPa"

     ;;;output the variables;;;;
     fout1->SF  = smh
     printMinMax(smh,True)
     delete([/umi,vmi,smh,psm/])

     ;;;process the vertical velocity
     tm         = f1->$modvars(3)$(indy,:,:,:)
     omega      = f1->$modvars(2)$(indy,:,:,:)
     pm         = tm
     pm         = tofloat(conform(tm,lev,1))
     pm         = pm * 100.0
     pm@units   = "Pa"
     wm         = omega_to_w(omega, pm, tm)
     wmi        = int2p_n_Wrap(lev,wm(:,:,:,:),pinterp1,intyp,1)
     wmi@units  = "m s~S~-1~N~" ;varunt(2)
     wmh        = dim_avg_n_Wrap(wmi,3)
     ;;;output the variables;;;;
     fout1->W   = wmh
     delete([/wmi,wmh/])
     delete([/tm,wm,pm,omega/])
     delete([/lev,pinterp1/])
  end do ; end exps

  delete([/cd_date,year,indy,obstime/])

 end do 

 end do 

end
