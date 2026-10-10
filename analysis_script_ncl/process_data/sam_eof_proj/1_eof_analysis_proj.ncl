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

  yst       = 1979
  yed       = 2014
  time_tag  = yst+"-"+yed

  top_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data"
  data_dir  = top_dir + "/" + "eof_analysis"
  era5_eof  = top_dir + "/sam_eof" 
  out_dir   = "./" 

  exp_name  = (/"ERA5", "v1.LR.amip", "v1.LR.historical", "v2.LR.amip", "v2.LR.historical"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
 
 ;define regions for eof analysis
  mode      = (/"SAM","PSA1","PSA2"/) 
  nmod      = dimsizes(mode)
  vars      = (/"PSL","Z500"/)
  facmod    = (/1.0, 1.0 /)
  varunt    = (/"Pa / std dev","m/ std dev"/)
  nvars     = dimsizes(vars)

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       = -20.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  do iv = 0,nvars-1,1

   do k = 0,nexps-1

    do i = 0, nseas-1,1

     print("working on season "+ seasons(i))

     var   = vars(iv)
     sea   = seasons(i)

     fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_anom_"+var+"_"+sea+"_en"+"*_"+time_tag+".nc")
     nens  = dimsizes(fl1)
     print(fl1)

     do ie = 0,nens-1,1

     f    = addfile(fl1(ie),"r")

     vnm  = var+ "_a"
     vd1  = f->$vnm$(:,{latS:latN},{lonW:lonE})
     time = f->time
     printVarSummary(vd1)
    ;printVarSummary(time)
    ;print(cd_calendar(time,-3))
    ;printMinMax (vd1, False) 
    ;printVarSummary(vd1)

     dimx = dimsizes(vd1)
     ntime = dimx(0)
     nlat  = dimx(1)
     nlon  = dimx(2)

    ; =================================================================
    ; create weights:  sqrt(cos(lat))   [or sqrt(gw) ]
    ; =================================================================
     rad    = get_d2r("float")
     clat   = vd1&lat            
     clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid
    ;printVarSummary(clat)
   
    ; =================================================================
    ; weight all observations 
    ; =================================================================
     xw     = vd1*conform(vd1, tofloat(clat), 1)
     copy_VarMeta(vd1, xw)
     xw@long_name = "Wgt: "+vars(iv)

    ;=================================================================
    ; Compute EOFs & Standardize time series
    ;=================================================================
    eof       = new((/nmod,nlat,nlon/),typeof(vd1))
    eof_ts    = new((/nmod,ntime/),typeof(vd1))
    eof_pcvar = new((/nmod/),typeof(vd1))
    eof_eval  = new((/nmod/),typeof(vd1))

    eof!0     = "mod"
    eof!1     = "lat"
    eof!2     = "lon"
    eof&lat   = vd1&lat
    eof&lon   = vd1&lon

    eof_ts!0    = "mod"
    eof_ts!1    = "time"
    eof_ts&time = time 

    eof_pcvar!0 = "mod"
    eof_eval!0  = "mod"
 
    eof_regres = eof

    optEOF  = True
    eof1    = eofunc_n_Wrap(xw, nmod, optEOF, 0)

    do ieof = 0,nmod-1

       if (ieof.eq.0) then
         avgx = avg(eof1(ieof,{-85:-90},:))
         if(avgx.gt.0.0) then
           eof1(ieof,:,:)  = eof1(ieof,:,:) * -1.0
         end if
         delete([/avgx/])
       end if
    end do 

    ;;read eof pattern obtained by ERA5
    fleof = systemfunc("ls "+era5_eof+"/"+"ERA5_EOF_PCS_"+vars(iv)+"_"+seasons(i)+"_en"+"*_"+time_tag+".nc")
    feof  = addfiles(fleof,"r")
    ListSetType (feof, "join")
    eofpm = dim_avg_n_Wrap(feof[:]->eof(:,:,:,:),0)

    optETS  = False
    eof1_ts = eofunc_ts_n_Wrap (xw, eofpm, optETS, 0)

    do ne=0,nmod-1
      eof(ne,:,:)        = (/eofpm(ne,:,:)/)
      eof_ts(ne,:)       = dim_standardize_n( eof1_ts(ne,:), 0, 0)      ; normalize
      eof_pcvar(ne)      = eof1@pcvar(ne)
      eof_eval(ne)       = eof1@eval(ne)
      eof_regres(ne,:,:) = (/regCoef_n(eof_ts(ne,:), vd1, 0, 0)/)
    end do 

    delete([/eof1,eof1_ts,fleof,feof,eofpm/])

    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    enstr      = "en"+sprinti("%02d",ie)
    out_file  = out_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars(iv)+"_"+seasons(i)+"_"+enstr+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
  
    eof@description    = "EOF pattern"
    eof_ts@description = "EOF PCs (standardized)"

    fout->eof      = eof
    fout->eof_reg  = eof_regres
    fout->eof_ts   = eof_ts
    fout->eof_pcv  = eof_pcvar 
    fout->eof_evl  = eof_eval

    delete([/time,xw,vd1,eof,eof_pcvar,eof_eval,eof_ts,eof_regres/])

   end do 
   end do 
   delete(fl1)

  end do    

 end do 

end

