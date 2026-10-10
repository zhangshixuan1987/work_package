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

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/data/E3SMv2/"
  exp_name  = (/"v2.LR.hist-all-xGHG-xaer"/)
  nexps     = dimsizes(exp_name)
  ystr      = 1979
  yend      = 2014
  time_tag  = "1979-2014"
  grid_tag  = "1x1"
  intyp     = 2

  seasons   = (/"MON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

  modvars    = (/"U","V","OMEGA","T","PS","FLUT"/)
  facmod     = (/   1,  1,    100,  1,  1, 1/)
  varunt     = (/"m~S~-1~N~","m~S~-1~N~","hPa s~S~-1~N~","K","Pa","W m~S~-2~N~"/)
  nvars      = dimsizes(modvars)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"/"+exp_name(k)+"*"+time_tag(k)+".nc")
   nens  = dimsizes(fl1)

   do ie = 0,nens-1,1

   print(fl1(ie))
   f1 = addfile(fl1(ie),"r")

   time = f1->time
  ;shift 1-day to be consistent with ERA5;;;
   time = time - 1
  ;printVarSummary(time)
  ;print(cd_calendar(time,-3))

   time_unit = "days since "+sprinti("%04d",ystr)+"-01-01 00:00:00"

   yyyymm       = cd_calendar(time(:),-1)
   yyyymm       = yyyymm - (toint(yyyymm(0))/100)*100  + ystr*100
   yyyy         = toint(yyyymm)/100
   mm           = yyyymm - yyyy*100
   dd           = mm
   hh           = mm
   mn           = mm
   sc           = mm
   dd           = 1
   hh           = 0                                  ; array syntax
   mn           = 0
   sc           = 0

   tnew         = cd_inv_calendar(yyyy,mm,dd,hh,mn,sc,time_unit, 0)
   tnew!0       = "time"

   delete([/time,yyyymm,yyyy,mm,dd,hh,mn,sc,dd,hh,mn,sc/])

   ;;reference pressure levels;;;; 
   pinterp           = asciiread("./pres_lev.txt",-1,"float")
   pinterp@longname  = "lev"
   pinterp@units     = "hPa"
   pinterp!0         = "lev"
   pinterp&lev       = pinterp
   pinterp1          = pinterp({10:1000})
  
   do i = 0, nseas-1,1 

    ;======================================================================
    ;Save the data for eof analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4Classic")
    enstr      = "en"+sprinti("%02d",ie)
    out_file1  = exp_name(k)+"_HC_analysis_"+seasons(i)+"_"+enstr+"_"+ystr+"-"+yend+".nc"
    system("rm " + out_file1)
    fout1      = addfile(out_file1,"cw")

    ;=================================================================
     print ("Start to process the cross-section of  selected variables ")
    ;start to calculate the stream function and velocity potential
     um       = f1->$modvars(0)$(:,:,:,:)
     vm       = f1->$modvars(1)$(:,:,:,:)
     psm      = f1->$modvars(4)$(:,:,:)
     um&time  = tnew
     vm&time  = tnew
     psm&time = tnew

     umi      = vinth2p(um,f1->hyam,f1->hybm,pinterp1,psm,intyp,f1->P0/100.0,1,False)
     vmi      = vinth2p(vm,f1->hyam,f1->hybm,pinterp1,psm,intyp,f1->P0/100.0,1,False)
     umi@units  = varunt(0)
     vmi@units  = varunt(1)
 
     ;;;derive the zonal mean for hadley cell;;;;
     umh        = dim_avg_n_Wrap(umi,3)
     vmh        = dim_avg_n_Wrap(vmi,3)

     printVarSummary(umh)
     printVarSummary(vmh)

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
     printVarSummary(sfvp)
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
     smh        = zonal_mpsi_Wrap(vmi,vmi&lat,levnew,psm)
     smh&lev    = pinterp1
     smh&lev@units = "hPa"

     ;;;output the variables;;;;
     fout1->SF  = smh
     ;printMinMax(smh,True)
     delete([/um,vm,umi,vmi,levnew,smh,psm/])

     ;;;process the vertical velocity
     tm       = f1->$modvars(3)$(:,:,:,:)
     tm&time  = tnew
     pm       = pres_hybrid_ccm (f1->PS,f1->P0,f1->hyam,f1->hybm)  ; Pa [kg/(m s2)]
     pm@units = "Pa"
     copy_VarCoords(tm,pm)

     omega    = f1->$modvars(2)$(:,:,:,:) 
     wm       = omega_to_w(omega, pm, tm)
     copy_VarCoords(tm,wm)
     psm        = f1->$modvars(4)$(:,:,:)
     wmi        = vinth2p(wm,f1->hyam,f1->hybm,pinterp1,psm,intyp,f1->P0/100.0,1,False)
     wmi@units  = "m s~S~-1~N~" ;varunt(2)
     wmh        = dim_avg_n_Wrap(wmi,3)
     ;;;output the variables;;;;
     fout1->W   = wmh
     delete([/wmi,wmh,psm,omega/])
     delete([/tm,wm,pm/])
     delete([/tnew/])
  end do ; end exps

end do 

end do 

end
