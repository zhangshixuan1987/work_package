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

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP","v2 Coupled"/)
  sig_95     = (/0.174, 0.325, 0.250/)
  nexps      = dimsizes(exp_name)

  regions    = (/"Global","SH","NH"/)
  nregs      = dimsizes(regions)

  labels     = (/"Nino3.4 vs GMT", "GMT vs SAM", "Nino3.4 vs SAM"/)
  nlabs      = dimsizes(labels)

 ;define regions for eof analysis
  vars       = (/"TREFHT"/)
  vnam       = (/"temperature anomaly"/)
  facmod     = (/1.0/)
  varunt     = (/"~S~o~N~C"/)

  time_tag   = (/"1981-2010"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1 

  print("working on season "+ seasons(i))

  do k = 0,nexps-1
  
   fl1  = data_dir + "/Global_2d_map/data/"+exp_name(k)+"_GMT_"+seasons(i)+"_"+time_tag+".nc"
   fl2  = data_dir + "/ENSO_EOF/data/"+exp_name(k)+"_Nino34_"+seasons(i)+"_"+time_tag+".nc"
   fl3  = data_dir + "/SAM_EOF/data/"+exp_name(k)+"_EOF_PCS_Z500_"+seasons(i)+"_"+time_tag+".nc"
   f1   = addfile(fl1,"r")
   f2   = addfile(fl2,"r")
   f3   = addfile(fl3,"r")

   gmt  = f1->dGMT
   enso = f2->ONI
   sam  = f3->eof_ts(:,0,:)
   time = yyyymm_to_yyyyfrac(cd_calendar(f1->time,-1),0.0)
 
   if(exp_name(k).eq."v2.LR.amip") then
    tmp = enso
    delete(enso)
    enso = gmt 
    printVarSummary(tmp)
    enso = conform(enso,tmp(0,:),1)
    delete(tmp)
   end if 

   printVarSummary(gmt)
   printVarSummary(sam)
   printVarSummary(enso)

   nens  = dimsizes(gmt(:,0))
   ntime = dimsizes(gmt(0,:))
   rc    = new((/3,ntime/),typeof(enso))

   if(seasons(i).eq."MON") then
     tint = 12 * 5
     do it = tint,ntime-tint-1
       xx1  = ndtooned(gmt(:,it-tint:it+tint))
       xx2  = ndtooned(dim_standardize(enso(:,it-tint:it+tint),1))
       xx3  = ndtooned(sam(:,it-tint:it+tint))
       rc(0,it) = escorc(xx1,xx2)
       rc(1,it) = escorc(xx1,xx3)
       rc(2,it) = escorc(xx2,xx3)
       delete([/xx1,xx2,xx3/])
     end do         
   else 
     tint = 5
     do it = tint,ntime-tint-1
       xx1  = ndtooned(gmt(:,it-tint:it+tint))
       xx2  = ndtooned(dim_standardize(enso(:,it-tint:it+tint),1))
       xx3  = ndtooned(sam(:,it-tint:it+tint))
       rc(0,it) = escorc(xx1,xx2)
       rc(1,it) = escorc(xx1,xx3)
       rc(2,it) = escorc(xx2,xx3)
       delete([/xx1,xx2,xx3/])
     end do
   end if 

   rc@_FillValue = -9999
   delete([/gmt,enso,sam/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = "ind_correlation_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,10,27,19,31/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/0.2,0.2,0.2,0.2,0.2/)
     plot   = new(nexps,graphic)
     dum1   = new(nlabs,graphic)
     dum2   = new(nlabs,graphic)
   end if 

   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.024
    res@tiMainFontHeightF     = FontHeightF
    res@tmYLLabelFontHeightF  = FontHeightF
    res@tmXBLabelFontHeightF  = FontHeightF
    res@gsnStringFontHeightF  = FontHeightF

    if(wtype.eq."png") then
     bthick = 3.0
    else
     bthick = 1.0
    end if

    res@tmBorderThicknessF    = bthick
    res@tmXBMajorThicknessF   = bthick
    res@tmXBMinorThicknessF   = bthick
    res@tmYLMajorThicknessF   = bthick
    res@tmYLMinorThicknessF   = bthick
    res@tmYRMajorThicknessF   = bthick
    res@tmYRMinorThicknessF   = bthick

    res@vpHeightF             = 0.40        ; Changes the aspect ratio
    res@vpWidthF              = 0.80
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@pmLegendDisplayMode    = "Never"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.40                   ; move units right
    res@pmLegendOrthogonalPosF = -0.45                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.12                  ; height of legend.
    res@lgLabelFontHeightF     = .018                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    res@gsnLeftString           = ""
    res@gsnRightString          = ""
    res@tiMainString            = ""
    res@tiYAxisString           = "" 
    res@tiXAxisString           = ""

    res@gsnYRefLine             = (/-sig_95(k),0.,sig_95(k)/)              ; reference line
    res@gsnYRefLineColors       = (/"Grey25","Grey25","Grey25"/)
    res@gsnYRefLineThicknesses  = (/1.0,0.5,1.0/)
    res@gsnYRefLineDashPatterns = (/1,0,1/)
    res@xyMonoLineColor         = False             ; want colored lines

    res@trYMinF                 = -1.0
    res@trYMaxF                 =  1.0
     
    res@trXMinF                 = 1984 ;floor(min(tsm03))
    res@trXMaxF                 = 2006 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************

    res@tiYAxisString         = "Correlation"      ; y-axis label
    res@tiXAxisString         = "Year"
    res@gsnLeftString         = exp_strs(k)
    res@gsnRightString        = seasons(i)
    res@tiMainString          = ""

    do ii = 0,nlabs-1,1

      res1  = res
      res1@xyMonoLineColor         = False             ; want colored lines
      res1@xyLineColors            = colors(ii)
      res1@xyLineThicknesses       = lthick(ii)      ; line thicknesses
      res1@xyDashPatterns          = ldash(ii)      ; make all lines solid
      dum1(ii)     = gsn_csm_xy (wks, time, rc(ii,:), res1)

      res2 = res
      res2@gsnLeftString          = ""
      res2@gsnRightString         = ""
      res2@tiMainString           = ""
      res2@tiYAxisString          = ""
      res2@tiXAxisString          = ""
      res2@gsnXYFillColors        = colors(ii)
      res2@gsnXYFillOpacities     = opcity(ii)
      res2@xyLineColor            = -1          
     ;dum2(ii) = gsn_csm_xy (wks, time, pltdat(1:,ii,:), res2)
     ;overlay(dum1(ii),dum2(ii))
    
      if(ii.eq.0) then
        plot(k) = dum1(ii)
      else
        overlay(plot(k),dum1(ii))
      end if 
    end do 

    delete([/time,rc/])

   if (k.eq.0) then 

   gres = True
   gres@YPosPercent          = 95.     
   gres@XPosPercent          = 15.    
   gres@Position             = "Top"        
   gres@ItemSpacePercent     = 8.          
   gres@LineLabelWhiteSpacePercent = 2 

   lineres = True
   lineres@lgLineColors       = colors
   lineres@lgLineOpacity      = opcity
   lineres@lgLineThicknesses  = lthick
   lineres@lgDashIndexes      = ldash
   lineres@LineLengthPercent  = 5      

   textres                    = True
   textres@lgLabels           = labels
   textres@lgLabelFontHeights = FontHeightF*0.8
   textres@lgLabelColors      = colors
   textres@lgLabelOpacity     = opcity

  ;plot = simple_legend(wks,plot,gres,lineres,textres)

   ;;add text;;
   txid = new(nlabs+1,graphic)                  ; Nine text strings
   txres                  = True             
  ;txres@txFont           = "helvetica-bold" 
   txres@txFontHeightF    = FontHeightF ;*0.8
   txres@txFontThicknessF = 4.0
   txres@txJust           = "BottomLeft"

   do m = 0,nlabs-1,1
     txres@txFontColor   = colors(m)
     txtstr              = labels(m) 
     txid(m) = gsn_add_text (wks,plot(k),txtstr, 1998.0, 0.65-0.2*m ,txres)
   end do 

  end if 

  end do
  
 ;panel plot only resources
  pres1                 = True                ; mods desired
  pres1@gsnFrame        = True               ; save panel until both ready
  gsn_panel(wks,plot,(/nexps,1/),pres1)          ; create first panel

 end do   

end

