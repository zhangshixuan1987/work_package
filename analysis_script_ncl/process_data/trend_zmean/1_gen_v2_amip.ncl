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

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/state_var"
  exp_name  = (/"v2.LR.amip"/)
  nexps     = dimsizes(exp_name)

  ystr      = 1979
  yend      = 2014

  time_tag  = ystr+"-"+yend
  grid_tag  = "1x1"

  ytrend    = 2001
  nseg      = 3  ; 3 trends: 1979-2001, 1979-2014, 2001-2014
  seasons   = (/"ANN","DJF","JJA","MAM","SON"/)
  nseas     = dimsizes(seasons)
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  latS       = -90.  ;negative for southern hemisphere
  latN       =   0.  ;negative for southern hemisphere
  lonW       =  0.0
  lonE       = 360.0

  varList    = (/"U","T"/)
  vnmList    = (/"Zonal wind", "Temperature"/)
  facList    = (/  1.0,   1.0/)
  untList    = (/"m s~S~-1~N~","K"/)
  nvars      = dimsizes(varList)
 ;declear the season or annual mean data to be used;;;;;

  do k  = 0,nexps-1,1

   do j = 0,nvars-1,1 

    do i = 0,nseas-1,1

     var = varList(j)
     sea = seasons(i)

     fl1   = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_anom_"+var+"_"+sea+"_en"+"*_"+time_tag(k)+".nc")
     nens  = dimsizes(fl1)
     print(fl1)

     do ie = 0,nens-1,1
      
     print("working on "+ fl1(ie))
     f1    = addfile(fl1(ie),"r")
    ;ListSetType (f1, "join")

     lat    = f1->lat
     lon    = f1->lon
     lev    = f1->lev

     yyyymm = cd_calendar(f1->time, -1)
     yyyy   = yyyymm/100
     ntim   = dimsizes(yyyy)

     rad    = get_d2r("float")
     clat   = lat
     clat   = sqrt( cos(rad*clat) )                 ; gw for gaussian grid
    ;printVarSummary(clat)

     vnam   = var+"_o"
     varx   = dim_avg_n_Wrap(f1->$vnam$,3)
     var0   = dim_avg_wgt_n_Wrap(varx(:,:,{latS:latN}),clat({latS:latN}),0,2)
     delete(varx)

     ;************************************************
     ; Calculate the regression coefficients (slopes)
     ;************************************************

     do kk = 0,nseg-1,1

       if(kk.eq.0) then 
         tstr  = yyyy(0)+"-"+yyyy(ntim-1)
         year  = ispan(yyyy(0), yyyy(ntim-1), 1)
         nyrs  = dimsizes(year)
         var1  = var0(time|:,lev|:)
       else if (kk.eq.1) then 
         indx  = ind(yyyy.eq.ytrend)
         tstr  = yyyy(0)+"-"+yyyy(indx)
         year  = ispan(yyyy(0), yyyy(indx), 1)
         nyrs  = dimsizes(year)
         var1  = var0(time|0:indx,lev|:)
         delete(indx)
       else 
         indx  = ind(yyyy.eq.ytrend)
         tstr  = yyyy(indx)+"-"+yyyy(ntim-1)
         year  = ispan(yyyy(indx), yyyy(ntim-1), 1)
         nyrs  = dimsizes(year)
         var1  = var0(time|indx:,lev|:)
         delete(indx)
       end if 
       end if 
 
      vclimo  = dim_avg_n_Wrap(var1,0)
      vtrend  = vclimo
      vtprob  = vclimo
      vtrend@long_name = "Trend"
      vtrend@units     = var0@units+"/year"
      vtprob@long_name = "probability"

      rc   = regCoef_n(year,var1,0,0)
      tval = onedtond(rc@tval , dimsizes(rc))
      df   = onedtond(rc@nptxy, dimsizes(rc)) - 2
      b    = tval    ; b must be same size as tval (and df)
      b    = 0.5
      prob = betainc(df/(df+tval^2),df/2.0,b)       ; prob(nlat,nlon)

      vtrend = (/tofloat(rc)/)
      vtprob = (/tofloat(prob)*100.0/)
      delete([/year,rc,tval,df,b,prob/])
       
      out_file  = exp_name(k)+"_en"+sprinti("%02d",ie)+"_"+var+"_"+seasons(i)+"_trend_"+tstr+".nc"
      system("rm " + out_file)
      fout = addfile(out_file,"cw")
      vs0  = var+"_climo"
      vs1  = var+"_trend"
      vs2  = var+"_prob"

      fout->$vs0$  = vclimo
      fout->$vs1$  = vtrend
      fout->$vs2$  = vtprob
      delete([/fout,var1,vtprob,vclimo,vs1,vs2,vtrend/])

    end do
    delete(var0)

   end do 

  end do 

 end do 

end do 

end
