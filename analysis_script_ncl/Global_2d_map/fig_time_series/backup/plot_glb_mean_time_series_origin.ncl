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

  data_dir   = "../data"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  ltrend     = new(nexps,string)

  regions    = (/"Global","SH","NH"/)
  nregs      = dimsizes(regions)

  labels     = (/"ERA5", "v2 AMIP", "v2 Coupled"/)
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

  do j = 0,nregs-1,1
  
   if (regions(j).eq."Global") then
     latS = -90.0 
     latN =  90.0
   else if (regions(j).eq."SH") then 
     latS = -90.0
     latN =   0.0
   else 
     latS =   0.0
     latN =  90.0
   end if 
   end if 

  do k = 0,nexps-1
  
   fl   = data_dir+"/"+exp_name(k)+"_anom_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
   f    = addfile(fl,"r")
   time = yyyymm_to_yyyyfrac(cd_calendar(f->time,-1),0.0)
   vstr = "d" + vars
   var  = f->$vstr$(:,:,{latS:latN},:)
   printMinMax(var,False)
   printVarSummary(var)
   exit

   latrg = f->lat({latS:latN})
   rad   = 4.*atan(1.)/180.    ;
   cost  = cos(latrg*rad)       ; cosine weights
   vmn   = wgt_areaave_Wrap(var,cost,1.,0)
   printVarSummary(vmn)
   nens  = dimsizes(vmn(0,:))
   rc    = new(nens,typeof(vmn))
   do ie = 0,nens-1,1
     rc(ie) = regline(time(:),vmn(:,ie)) * 10.0
   end do
   print(avg(rc)+ " " +stddev(rc))
   ltrend(k) = sprintf("%.2f",avg(rc)) + " ~F34~1~N~~F~ " + sprintf("%.3f",stddev(rc))
   delete([/var,rc,cost,rad,latrg,f,fl/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = regions(j)+"_mean_"+vars+"_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,10,27,19,31/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/0.2,0.2,0.2,0.2,0.2/)
     dum1   = new(nexps,graphic)
     dum2   = new(nexps,graphic)
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
    res@vpWidthF              = 0.90
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

    res@gsnYRefLine             = (/0./)              ; reference line
    res@gsnYRefLineColors       = (/"Grey25"/)
    res@gsnYRefLineThicknesses  = (/0.5/)
    res@gsnYRefLineDashPatterns = (/0/)
    res@xyMonoLineColor         = False             ; want colored lines

    res@trYMinF                 = -0.8
    res@trYMaxF                 =  0.8
     
    res@trXMinF                 = 1980 ;floor(min(tsm03))
    res@trXMaxF                 = 2011 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************

    if(k.eq.0) then    
      res@tiYAxisString         = "~F33~D~N~~F~T (" + varunt+")"      ; y-axis label
      res@tiXAxisString         = "Year"
      res@gsnLeftString         = regions(j)+" mean " + vnam 
      res@gsnRightString        = seasons(i)
      res@tiMainString           = ""
    else
      res@gsnLeftString          = ""
      res@gsnRightString         = ""
      res@tiMainString           = ""
      res@tiYAxisString          = ""
      res@tiXAxisString          = ""
    end if 

    pltdat = new((/3,dimsizes(vmn(:,0))/),typeof(vmn))
    pltdat(0,:) = dim_avg_n(vmn,1)
    pltdat(1,:) = dim_avg_n(vmn,1) - dim_stddev_n(vmn,1)
    pltdat(2,:) = dim_avg_n(vmn,1) + dim_stddev_n(vmn,1)

    res1  = res
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = colors(k)
    res1@xyLineThicknesses       = lthick(k)      ; line thicknesses
    res1@xyDashPatterns          = ldash(k)      ; make all lines solid
    dum1(k)     = gsn_csm_xy (wks, time, pltdat(0,:), res1)

    res2 = res
    res2@gsnLeftString          = ""
    res2@gsnRightString         = ""
    res2@tiMainString           = ""
    res2@tiYAxisString          = ""
    res2@tiXAxisString          = ""
    res2@gsnXYFillColors        = colors(k)
    res2@gsnXYFillOpacities     = opcity(k)
    res2@xyLineColor            = -1          
    dum2(k) = gsn_csm_xy (wks, time, pltdat(1:,:), res2)
    overlay(dum1(k),dum2(k))

    if(k.eq.0) then
      plot = dum1(k)
    else
      overlay(plot,dum1(k))
    end if 
 
   delete([/vmn,time,pltdat/])

   end do 
 
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
   txres@txFontHeightF    = FontHeightF*0.8
   txres@txFontThicknessF = 4.0
   txres@txJust           = "BottomLeft"

   do m = 0,nlabs-1,1
    ;if(m.eq.0) then
    ;  txres@txFontColor   = colors(m)
    ;  txtstr = "slope / decade:"
    ;else 
     txres@txFontColor   = colors(m)
     txtstr = labels(m) + " ("+ltrend(m) +")"
    ;end if 
     txid(m) = gsn_add_text (wks,plot,txtstr, 1982.0, 0.6-0.15*m ,txres)
   end do 

  
   ;panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/1,1/),pres1)          ; create first panel

  end do

 end do   

end

