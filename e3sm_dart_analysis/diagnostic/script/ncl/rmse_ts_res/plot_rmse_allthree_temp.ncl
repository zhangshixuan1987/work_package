;*************************************************
; histo_4.ncl
;
; Concepts illustrated:
;   - Drawing a histogram using exact min/max of data
;   - Generating dummy data using "rand"
;   - Drawing X axis labels at a 45 degree angle
;
;************************************************
;
; This file is loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
;************************************************

load "./function_read_dart_diag.ncl"

begin

  l_create_plot_data = True
  Timetag = "2011011500-2011012012"
  nday    = 5
  ilev    = 2
  ilev    = 8
  timestr = "Days since 0000 UTC 15 January 2011"
  datapth = "/compyfs/zhan391/diag_dart/fig_data/"
  ObsStr  = "RADIOSONDE_TEMPERATURE"
  
  hres    = (/"ne11_oQU240","ne16_oQU240","ne30pg2_EC30to60E2r2"/)
  nres    = dimsizes(hres)

  labels  = (/"ne11(2.8deg)","ne16(1.9deg)","ne30(1deg)"/)

  regions = (/"SH","NH","NAmer","Tropics"/)
  regstrs = (/"Southern Hemisphere","Northern Hemisphere","North America","Tropics"/)
  latlons = (/"[0-90S]","[0-90N]","[65~S~o~N~-125~S~o~N~ W, 25~S~o~N~-55~S~o~N~ N]","10S-10N"/)
  nregs   = dimsizes(regions)

  ;variable to be plotted: totalspread,rmse and data rejection rate
  varnam  = (/"RMSE","Spread","Rej. Rate"/)
  varstr  = (/"Root Mean Square Error (RMSE)", \
              "Total Ensemble Spread", \
              "Observation Rejection Rate"/)
  varunt  = (/"K","K","%"/)
  nvars   = dimsizes(varnam)

  ;;;read the variable information;;
  rmse   = new((/nres,nday*4+1/),double)
  rmse@_FillValue = -999.0
  rrate  = rmse 
  spread = rmse 
  
  do k = 0,nregs-1,1 

    region = regions(k)
    regstr = regstrs(k)
    latlon = latlons(k)

    if (region .eq. "NAmer") then 
      ymax    = (/12.0, 8.0, 25.0/)
      ymin    = (/0.0,  0.0, 0.0 /)
    else
     ymax    = (/12.0, 5.0, 15.0/)
     ymin    = (/0.0,  0.0, 0.0 /)
    end if 

    do i = 0,nres-1,1
      fname = datapth + ObsStr+"_"+region+"_"+hres(i)+".nc"
      ff = addfile(fname,"r") 
      rmse(i,:)   = ff->RMSE(:,ilev)
      spread(i,:) = ff->Spread(:,ilev)
      rrate(i,:)  = ff->Rrate(:,ilev)
      pltlev      = ff->plev(ilev)
      time        = ff->time 
    end do 

  ;;start to do the plot ;;;;;
  plot = new(nvars,graphic)
  plot_name = "Fig_dart_diag_rmse_time_series_"+region+"_"+pltlev
  print(" ")
  print("plot figure : " + plot_name+".pdf")
  print(" ")
  
  wtype = "pdf"
  wks = gsn_open_wks(wtype,plot_name)             ; send graphics to PNG file

  gsn_define_colormap(wks,"cb_9step")
  colors  = (/14,28,2/)
  markers = (/16, 4, 7/)
  mkszf   = (/0.8, 1.0, 1.0/)*1.5
  dashind = (/0,0,0/)       ;;line dash pattern
  mkthick = (/1,1,1/)*2.0
  lnthick = (/2,2,2/) ;*1.5   ;;line thichness

  ;---Resources to share between both plots
  res          = True                   ; plot mods desired
  res@vpHeightF = 0.35
  res@vpWidthF =  0.9
  res@gsnDraw  = False
  res@gsnFrame = False

  FontHeight = 0.024
  res@gsnStringFontHeightF = FontHeight*1.2
  res@tmXBLabelFontHeightF = FontHeight*1.2
  res@tmYLLabelFontHeightF = FontHeight*1.2
  res@tmYRLabelFontHeightF = FontHeight*1.2
  res@tiMainFontHeightF    = FontHeight*1.2
  res@tiXAxisFontHeightF   = FontHeight*1.2
  res@tiYAxisFontHeightF   = FontHeight*1.2

  BorderThick = 1.5
  res@tmBorderThicknessF    = BorderThick
  res@tmXBMajorThicknessF   = BorderThick
  res@tmYLMajorThicknessF   = BorderThick
  res@tmXTMajorThicknessF   = BorderThick
  res@tmYRMajorThicknessF   = BorderThick
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0
  res@tmXTMinorThicknessF   = res@tmXTMajorThicknessF/2.0
  res@tmYRMinorThicknessF   = res@tmYRMajorThicknessF/2.0

  res@trYReverse        = False                      ; reverse Y-axis
 ;res@xyYStyle          = "Log"

  res@tmXTOn  = False
  res@tmYROn  = False 
  res@gsnYRefLine             = (/0./)      ; four X reference lines
  res@gsnYRefLineThicknesses  = 1.0
  res@gsnYRefLineDashPatterns = 0
  res@gsnYRefLineColors       = "Grey"

  res@gsnRightString                = region ;pltlev+" (hPa)"   ;VstList(0)
  res@gsnLeftString                 = ""   ;Groups(0)
  res@gsnCenterString               = ""   ;fig_label
  res@gsnRightStringOrthogonalPosF  = 0.02
  res@gsnLeftStringOrthogonalPosF   = 0.02
  res@gsnCenterStringOrthogonalPosF = 0.02

; add a legend
  res@pmLegendDisplayMode    = "Never"              ; turn on legend
  res@pmLegendSide           = "Top"                 ; Change location of
  res@pmLegendWidthF         = 0.12                  ; Change width and
  res@pmLegendHeightF        = 0.10                  ; height of legend.
  res@lgLabelFontHeightF     = FontHeight*0.95       ; change font height
  res@lgPerimOn              = False                 ; no box around

  res@tiMainString         = "" 
  res@tiMainOffsetYF       = 0.02
  res@tiMainFontThicknessF = 3.0
  res@tiXAxisString        = "" 
 ;res@tiXAxisOffsetYF      = -0.015

  res@trXMaxF         = nday+0.5
  res@trXMinF         = 0
  res@tmXBMode        = "Explicit"	
  res@tmXBValues      = ispan(0,nday,1) 
  res@tmXBLabels      = "" + res@tmXBValues 
  res@tmXBMinorValues = fspan(0,nday,2*nday+1)

  ;left axis
  resL  = res
  resL@tiYAxisString          := ""
  resL@tiXAxisString          := timestr 
 ;resL@tmXBLabelsOn           := False
  
  resL@xyMarkLineMode         := "MarkLines"
  resL@xyLineColors           := colors(0:2)
  resL@xyMarkers              := markers(0:2)
  resL@xyMarkerColors         := colors(0:2)
  resL@xyMarkerSizes          := 0.015*mkszf(0:2)
  resL@xyMarkerThicknesses    := mkthick
  resL@xyDashPatterns         := dashind(0:2)
  resL@xyLineThicknesses      := lnthick(0:2)

; add a legend
  resL@xyExplicitLegendLabels := "  " +labels
  resL@pmLegendDisplayMode    := "Always" ;"Never"
 ;resL@pmLegendSide           := "Top"                 ; Change location of
  resL@pmLegendParallelPosF   :=  0.75                   ; move units right
  resL@pmLegendOrthogonalPosF := -0.95                ; more neg = down
  resL@pmLegendSide           = "Top"                 ; Change location of
  resL@pmLegendWidthF         = 0.16                  ; Change width and
  resL@pmLegendHeightF        = 0.18                  ; height of legend.
  resL@lgLabelFontHeightF     = FontHeight*1.05       ; change font height
  resL@lgPerimOn              = False                 ; no box around
  
  res0 = resL
  res0@pmLegendDisplayMode    := "Always" 
  res0@gsnLeftString  = varstr(0)
  res0@trYMaxF       := ymax(0) ;0.21
  res0@trYMinF       := ymin(0) ;0.0    
  res0@tiYAxisString := varnam(0) + " ("  +varunt(0)+ ")" 
  plot(0) = gsn_csm_xy(wks,tofloat(time),tofloat(rmse),res0)

  res1 = resL
  res1@pmLegendDisplayMode    := "Never"
  res1@gsnLeftString  = varstr(1)
  res1@trYMaxF       := ymax(1) ;0.21
  res1@trYMinF       := ymin(1) ;0.0  
  res1@tiYAxisString := varnam(1) + " ("  +varunt(1)+ ")"
  plot(1) = gsn_csm_xy(wks,tofloat(time),tofloat(spread),res1)

  res2 = resL
  res2@pmLegendDisplayMode    := "Never"
  res2@gsnLeftString  = varstr(2) 
  res2@trYMaxF       := ymax(2) ;0.21
  res2@trYMinF       := ymin(2) ;0.0 
  res2@tiYAxisString := varnam(2) + " ("  +varunt(2)+ ")"
  ;printVarSummary(time)
  ;printVarSummary(rmse)
  ;printVarSummary(rrate)
  ;printVarSummary(spread)   
  plot(2) = gsn_csm_xy(wks,tofloat(time),tofloat(rrate),res2)

;--Draw both plots in a panel
  pnres                         = True
  pnres@gsnMaximize             = False
  pnres@gsnPanelMainString      = " "
  pnres@gsnPanelMainFontColor   = "Black"
  pnres@gsnPanelMainFontHeightF = 0.020
  pnres@gsnPanelLabelBar        = False
  pnres@gsnPanelBottom          = 0.05
  pnres@gsnPanelYWhiteSpacePercent = 2
  pnres@gsnPanelYF = (/0.80,0.53,0.26/)    ; Adjust rightmost plots
  gsn_panel(wks,(/plot/),(/3,1/),pnres)
  delete([/pnres/])

 end do

end
