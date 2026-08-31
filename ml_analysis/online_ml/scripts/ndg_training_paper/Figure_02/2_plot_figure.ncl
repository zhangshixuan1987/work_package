;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Figure 5 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  ModelOutRoot  = "./"
  file_name     = "PCC_and_RMSE_1year_mean_with_respect_to_obs.nc"
  plotFileLabel = "Figure_02"  ;; "global_mean_rmse_bias_clim"

  plot     = new(2,graphic)
  dum1     = new((/10/),graphic)
  dum2     = new((/10/),graphic)
  dum3     = new((/10/),graphic)
  dum4     = new((/10/),graphic)

  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+plotFileLabel)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  leftstring    = (/"(a) Global Pattern Correlation Coefficient (PCC)","(b) Relative difference in global pattern"/)
  right_string = ""
  Groups       = (/"CLIM", "DNDG-UV", "DNDG-UVT", "DNDG-UVQ", "DNDG-UVTQ","RNDG-UVTQ"/)
  ngrps        = dimsizes(Groups)
  labels       = (/"CTRL            ", "DNDG-UV      ", "DNDG-UVT   ", "DNDG-UVQ   ", "DNDG-UVTQ","RNDG-UVTQ"/)
  load  "./fig_setup.ncl"

  ; Assign the new colormap
  gsn_define_colormap(wks, "NCV_bright")
  ; Set your plot color indices (gray is now at index 0)
  ;lpcols = (/0, 254, 171, 120, 70, 95, 2/)

  gsn_define_colormap(wks, "t2m_29lev")
  ; Set your plot color indices (gray is now at index 0)
  lpcols = (/15, 24, 22, 19, 2, 27/)

  colors = lpcols
  BarChartFillScale = 0.3
 
  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  fin = addfile(ModelOutRoot+"/"+file_name,"r")

  xlabnam             = fin->Evaluated_variables
  nvars               = dimsizes(xlabnam)
 
 ;;;processe the zonal mean data;;;;
  plotbias  = new((/ngrps,nvars/),float)
  plotrmse  = new((/ngrps,nvars/),float)
  plotxlab  = new((/ngrps,nvars/),float)

;;......................................................
;; loop over cases (D, 6h, 3h, 1h) 
;;......................................................
  do igp  = 0,ngrps-1
     vnmbias             = "PCC_"+Groups(igp)
     vnmrmse             = "RMSE_"+Groups(igp)
     plotxlab(igp,:)     = ispan(1,nvars,1)+0.5*(igp-ngrps/2.0)
     plotbias(igp,:)     = fin->$vnmbias$
     plotrmse(igp,:)     = fin->$vnmrmse$
  end do 

;;......................................................
;; plot 
;;......................................................
  res               = True           ; plot mods desired
  res@vpXF          = 0.17
  res@vpYF          = 0.8
  res@vpWidthF      = 0.9                      ; set width of plot
  res@vpHeightF     = 0.2                      ; set height of plot
  res@gsnDraw  = False
  res@gsnFrame = False

  res@tmYROn = False            ; Turn off right tickmarks.
  res@tmXTOn = False            ; Turn off top tickmarks
  FontHeight = 0.016
  res@gsnStringFontHeightF  = FontHeight*1.15
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight*0.3
  res@tiYAxisFontHeightF    = FontHeight


  ;res@tiMainFontThicknessF  = 8.0
  ;res@tiYAxisFontThicknessF = 5.0
  ;res@tiXAxisFontThicknessF = 5.0

  bthick = 1.0
  res@tmBorderThicknessF    = bthick
  res@tmXBMajorThicknessF   = bthick
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMajorThicknessF   = bthick
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0

  res@gsnYRefLine           = 0.0             ; create a reference line   
  ;res@gsnYRefLineDashPattern= 2             
  ;res@gsnYRefLineColor      = "Black"
  res@gsnYRefLineThicknessF = 0.25

  res@tiMainString            = ""
  res@gsnLeftString           = "" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor  = ""
  res@gsnRightString          = ""
  res@gsnRightStringFontColor = ""

  
  res@xyMarkLineModes         = "Lines"

  res@trXMinF               = 0.1 ;min(latrgd)
  res@trXMaxF               = 2*nvars+1.25 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = fspan(1.075,2*nvars+.075,nvars) ;ispan(1,nvars,1)
  res@tmXBMinorValues       = fspan(1.075,2*nvars+.075,nvars) ;ispan(1,nvars,1)
  res@tmXBLabels            = xlabnam
  res@tmXBLabelDirection    = "Across"
  res@tmXBLabelAngleF       = 90.0
  res@tmXBLabelJust         = "CenterRight"
  res@tmYLMinorPerMajor     = 2

   res@tmYMajorGrid            = False        ; Add horizontal grid lines
  ;res@tmYMinorGrid            = False
  ;res@tmYMajorGridLineColor   = "LightGrey"
  ;res@tmYMinorGridLineColor   = "White"
   res@tmYMajorGridLineDashPattern = 0
  ;res@tmYMinorGridLineDashPattern = 5
   res@tmYMajorGridThicknessF  = 0.3
   res@tmGridDrawOrder         = "PreDraw"    ; new resource added in NCL V6.5.0

  res@gsnXYBarChart                   = True ; create bar chart
  res@gsnXYBarChartOutlineThicknessF  = 0.5;1.008
  res@gsnXYBarChartFillLineThicknessF = 1.0

  ;res@xyLineThicknessF     = 10.0            ; make thicker
  ;res@xyDashPatterns       = dashind(0:nvars)
  ;res@xyMarkerSizes        = 0.012
  ;res@xyLineColors         = colors(0:nvars)
  ;res@xyMarkers            = markers(0:nvars)
  ;res@xyMarkerColors       = colors(0:nvars)

  ;res@tmXTAutoPrecision    = False
  ;res@tmYLPrecision        =  2
   res@tmXBOn               = True            ; Turn off top tickmarks
   res@tmXBLabelsOn         = True            ; Turn off top tickmarks

   res@tmYLMajorLengthF = 0.01   ; Shorter major ticks (default is usually 0.015)
   res@tmYLMinorLengthF = 0.005  ; Shorter minor ticks (default is usually 0.0075)
   res@tmXBMajorLengthF = 0.01  ; Shorter major ticks (default is usually 0.015)
   res@tmXBMinorLengthF = 0.005  ; Shorter minor ticks (default is usually 0.0075)

  do igp = 0,ngrps-1,1

   res@xyLineColors              = (/-1,-1/);(/colors(0),colors(0)/)
   res@gsnXYBarChartBarWidth     = barWidth(igp)          ; change bar widths
   res@gsnXYBarChartColors       = colors(igp)
   res@gsnXYBarChartFillOpacityF = bopc(igp)
   res@gsnXYBarChartPatterns     = 0 ;pattns(igp)
   res@gsnXYBarChartFillScaleF   = BarChartFillScale 

   xx  = fspan(1.25,2*nvars+0.25,nvars)+0.25*(igp-ngrps/2.0)

   ;;;;plot global mean biases
   res0 = res
   res0@tmXBLabelsOn         = False            ; Turn off top tickmarks
   res0@trYMinF              =  0.75  ;floor(min(pltobs)) ;min(ynflx)-0.1 ;You could also just use
   res0@trYMaxF              =  1.02  ;ceil(max(pltobs)) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
   res0@tmYLMode             = "Explicit"
   res0@tmYLValues          := fspan(res0@trYMinF,1.0,6)
   res0@tmYLMinorValues     := fspan(res0@trYMinF,1.0,11)
   res0@tmYLLabels          := sprintf("%5.2f",res0@tmYLValues)
   res0@tmYLLabels          := (/"0.75","0.80","0.85","0.90", "0.95", "1.0"/)

   if(igp.eq.0)then
     res0@tiMainString                = ""
     res0@tiYAxisString               = "PCC" ;(NDG - OBS) / OBS"
     res0@tiXAxisString               = "";Varname(iv)+ " ("+var@units+")"
     res0@gsnLeftString               = leftstring(0)
     res0@gsnLeftStringOrthogonalPosF = 0.04
     res0@gsnLeftStringFontColor      = "Black"
     res0@gsnRightString              = right_string
     res0@gsnRightStringFontColor     = "Black"
   else
     res0@tiMainString                = ""
     res0@tiYAxisString               = ""
     res0@tiXAxisString               = ""
     res0@gsnLeftString               = ""
     res0@gsnRightString              = ""
   end if 

   dum1(igp) = gsn_csm_xy(wks,xx,plotbias(igp,:),res0)

   ;; for bars with fill pattern 
   res0@gsnXYBarChartPatterns         = pattns(igp)  
   res0@gsnXYBarChartColors           = lpcols(igp) ;1
   res0@gsnXYBarChartFillScaleF       = BarChartFillScale 
   res0@gsnXYBarChartFillOpacityF     = 0.5 
   res0@gsnXYBarChartFillLineThicknessF = 0.01
   dum3(igp) = gsn_csm_xy(wks,xx,plotbias(igp,:),res0)

   ;;;;plot global mean centered RMSD with respect to RMS of CLIM
   res1 = res
   res1@tmXBOn               = True            ; Turn off top tickmarks
   res1@tmXBLabelsOn         = True            ; Turn off top tickmarks

   res1@tiMainString         = ""
   res1@trYMinF              = 0.0
   res1@trYMaxF              = 2.0
   res1@tmYLMode             = "Explicit"
   res1@tmYLValues          := fspan(res1@trYMinF,res1@trYMaxF,6)
   res1@tmYLMinorValues     := fspan(res1@trYMinF,res1@trYMaxF,11)
   res1@tmYLLabels          := (/"0.0","0.4","0.8","1.2","1.6", "2.0"/)

   if(igp.eq.0)then
     res1@tiMainString                = ""
     res1@tiYAxisString               = "RMSE / RMS~B~OBS~N~"
     res1@tiXAxisString               = "";Varname(iv)+ " ("+var@units+")"
     res1@gsnLeftString               = leftstring(1)
     res1@gsnLeftStringOrthogonalPosF = 0.04
     res1@gsnLeftStringFontColor      = "Black"
     res1@gsnRightString              = right_string
     res1@gsnRightStringFontColor     = "Black"
   else
     res1@tiMainString                = ""
     res1@tiYAxisString               = ""
     res1@tiXAxisString               = ""
     res1@gsnLeftString               = ""
     res1@gsnRightString              = ""
   end if 

   print("plot pattern : " + igp + "   " + pattns(igp)) 

   dum2(igp) = gsn_csm_xy(wks,xx,plotrmse(igp,:),res1)

   ;; for bars with fill pattern 

   res1@gsnXYBarChartPatterns           = pattns(igp)
   res1@gsnXYBarChartFillScaleF         = BarChartFillScale 
   res1@gsnXYBarChartFillLineThicknessF = 0.01 
   res1@gsnXYBarChartColors             = lpcols(igp) ; 1
   res1@gsnXYBarChartFillOpacityF       = 0.5

   dum4(igp) = gsn_csm_xy(wks,xx,plotrmse(igp,:),res1)

;;..........................................
;; overlay bars from the 2nd one  
;;..........................................
   if(igp.eq.0)then
    plot(0) = dum1(igp)
    plot(1) = dum2(igp)
   else
    overlay (plot(0),dum1(igp))
    overlay (plot(1),dum2(igp))
   end if

  ;;..........................................
  ;; for UVT, use different fill pattern  
  ;;..........................................
  ;overlay (plot(0),dum3(igp))
  ;overlay (plot(1),dum4(igp))

  end do  ; end of the "group" loop, corresponding to different nudging configs.
  ;----------------------

  ;----------------------
  ;Attach a legend
  ;----------------------
  lgres                    = True
  lgres@vpWidthF           = 0.60          ; width of legend (NDC)
  lgres@vpHeightF          = 0.06          ; height of legend (NDC)
  lgres@lbBoxMajorExtentF  = 0.12         ; puts space between color boxes
  lgres@lbBoxMinorExtentF  = 0.20
  lgres@lbMonoFillPattern  = False          ; Solid fill pattern
  lgres@lgLabelFontHeightF = FontHeight*0.4     ; legend label font thickness
  lgres@lgMonoDashIndex    = False
  lgres@lgMonoLineThickness= False
  lgres@lgPerimOn          = True
  lgres@lgPerimThicknessF  = 1.0
  lgres@lbLabelJust        = "CenterLeft"  ; left justify labels

  do i = 0,ngrps-1,1
    ;lgres@lbBoxLineDashPattern = dashind(i) 
     lgres@lbBoxLineColor       = "black"
     lgres@lbBoxLineThicknessF  = 0.5
     lgres@lbFillColors         = colors(i)
     lgres@lbFillOpacityF       = bopc(1)*3.0
     lgres@lbFillPatterns       = pattns(i)
     lgres@gsnXYBarChartFillScaleF = 0.4
     lgres@lbFillPatterns         = pattns(i)

     print("ndguv pattern : " + i + "   " + pattns(i)) 

     lgres@vpWidthF           = 0.12          ; width of legend (NDC)
     ypos = 0.85 + (ngrps-i-1)*0.018
     gsn_labelbar_ndc(wks,1,labels(i),0.85,ypos,lgres)
     ;ypos = 0.618 + (ngrps-i-1)*0.018
     ;gsn_labelbar_ndc(wks,1,labels(i),0.75,ypos,lgres)

  end do  ; do i = 0,ngrps-1,1

 panres = True                                         ; panel resource list
 panres@gsnMaximize                = True
 panres@gsnPaperMargin             = 0.3
 panres@gsnPanelBottom             = 0.38     ; move bottom up so size is 0.3
 panres@gsnPanelYWhiteSpacePercent = 5.0       ;set spacing vertically between 2 panels
;panres@gsnPanelYF          = (/ -1.0,.64/)      ; Adjust bottommost plots.
 gsn_panel(wks,(/plot(0),plot(1)/),(/3,1/),panres)
;frame(wks)

end



