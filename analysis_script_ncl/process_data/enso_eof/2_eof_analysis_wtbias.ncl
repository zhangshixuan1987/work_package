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
  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/sst_hadsst"
  exp_name  = (/"HadSST","ERA5", "v2.LR.amip", "v2.LR.historical"/)
  nexps     = dimsizes(exp_name)
  
  outdir    = "wtmbias"
  seasons   = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)

 ;define regions for eof analysis
  latS       = -20.  ;negative for southern hemisphere
  latN       =  20.  ;negative for southern hemisphere
  lonW       = 120.0
  lonE       = 300.0
 
 ;define regions for eof analysis
  mode      = (/"ENSO_EP","ENSO_CP"/) 
  nmod      = dimsizes(mode)
  vars      = (/"SST"/)
  facmod    = (/1.0/)
  varunt    = (/"~S~o~N~C / std dev"/)
  nvars     = dimsizes(vars)

  do iv = 0,nvars-1,1

   do k = 0,nexps-1

    do i = 0, nseas-1,1

     var   = vars(iv)
     sea   = seasons(i)

     fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_anom_"+var+"_"+sea+"_en"+"*_"+time_tag+".nc")
     nens  = dimsizes(fl1)
     print(fl1)

     f    = addfiles(fl1,"r")
     ListSetType (f, "join")

     vnm  = var+ "_a"
     vd1  = f[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
     time = f[0]->time
     printVarSummary(vd1)
    ;printVarSummary(time)
    ;print(cd_calendar(time,-3))
    ;printMinMax (vd1, False)
    ;printVarSummary(vd1)

     dimx = dimsizes(vd1)
     ntime = dimx(1)
     nlat  = dimx(2)
     nlon  = dimx(3)

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
     xw     = vd1
     xw     = vd1*conform(vd1, tofloat(clat), 2)
     copy_VarMeta(vd1, xw)
     xw@long_name = "Wgt: "+vars(iv)

    ;=================================================================
    ; Compute EOFs & Standardize time series
    ;=================================================================
    eof       = new((/nens,nmod,nlat,nlon/),typeof(vd1))
    eof_ts    = new((/nens,nmod,ntime/),typeof(vd1))
    eof_pcvar = new((/nens,nmod/),typeof(vd1))
    eof_eval  = new((/nens,nmod/),typeof(vd1))

    eof!0     = "ens"
    eof!1     = "mod"
    eof!2     = "lat"
    eof!3     = "lon"
    eof&lat   = vd1&lat
    eof&lon   = vd1&lon

    eof_ts!0    = "ens"
    eof_ts!1    = "mod"
    eof_ts!2    = "time"
    eof_ts&time = time 

    eof_pcvar!0 = "ens"
    eof_pcvar!1 = "mod"
    eof_eval!0  = "ens"
    eof_eval!1  = "mod"
 
    eof_regres = eof

    do in = 0,nens -1,1

      xw1     = xw(in,:,:,:)
      xa1     = vd1(in,:,:,:)

      optEOF  = True
      eof1    = eofunc_n_Wrap(xw1, nmod, optEOF, 0)

      do ieof = 0,nmod-1

       if (ieof.eq.0) then
         avgx     = avg(eof1(ieof,{-5:5},{180:270}))
         if(avgx.lt.0.0) then
           eof1(ieof,:,:)  = eof1(ieof,:,:) * -1.0
         end if
         delete([/avgx/])
       else
         avgx     = avg(eof1(ieof,{-5:5},{250:270}))
         if(avgx.gt.0.0) then
            eof1(ieof,:,:)  = eof1(ieof,:,:) * -1.0
         end if
         delete([/avgx/])
       end if

      end do 

      optETS  = False
      eof1_ts = eofunc_ts_n_Wrap (xw1, eof1, optETS, 0)

      do ne=0,nmod-1
        eof(in,ne,:,:)        = (/eof1(ne,:,:)/)
        eof_ts(in,ne,:)       = dim_standardize_n( eof1_ts(ne,:), 0, 0)      ; normalize
        eof_pcvar(in,ne)      = eof1@pcvar(ne)
        eof_eval(in,ne)       = eof1@eval(ne)
        eof_regres(in,ne,:,:) = (/regCoef_n(eof_ts(in,ne,:), xa1, 0, 0)/)
      end do 

      delete([/xw1,xa1,eof1,eof1_ts/])

    end do 
  
    ;======================================================================
    ;Save the data for the regression analysis
    ;==================================================================
    setfileoption("nc", "Format",  "NetCDF4")
    out_file  = outdir+"/"+exp_name(k)+"_EOF_PCS_"+vars(iv)+"_"+seasons(i)+"_"+time_tag+".nc"
    system("rm " + out_file)
    fout = addfile(out_file,"cw")
  
    eof@description    = "EOF pattern"
    eof_ts@description = "EOF PCs (standardized)"

    fout->eof      = eof
    fout->eof_reg  = eof_regres
    fout->eof_ts   = eof_ts
    fout->eof_pcv  = eof_pcvar 
    fout->eof_evl  = eof_eval

   ;delete([/tnew,yyyymm,yyyy,mm,dd,hh,mn,sc/])
    delete([/time,xw,vd1,eof,eof_pcvar,eof_eval,eof_ts,eof_regres/])
    delete([/fl1,f/])
   end do 
   
  end do    

 end do 

end

