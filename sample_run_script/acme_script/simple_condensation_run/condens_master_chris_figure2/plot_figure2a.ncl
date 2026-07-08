;****************************************************
; xy_25.ncl
;
; Concepts illustrated:
;   - Drawing an XY plot with two different Y axes
;   - Adding curves to an existing XY plot that contains two different Y axes
;   - Removing trailing zeros from tickmark labels
;   - Maximizing plots after they've been created
;
;****************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
;****************************************************
begin
;***************************
; read in data
;***************************
  Datadir = "./"
 
  dtime   = 4  ; time step size is 4s 

  filnam1 = "SGR_P000_RKZ_fac_1069852.txt"   ; cloud fraction
  filnam2 = "SGR_P000_RKZ_term_A_1069852.txt"
  filnam3 = "SGR_P000_RKZ_term_B_1069852.txt"
  filnam4 = "SGR_P000_RKZ_term_C_1069852.txt"

  strfac  = asciiread(Datadir+filnam1, -1, "float")   ; -1 means readd all rows
  strtma  = asciiread(Datadir+filnam2, -1, "float") 
  strtmb  = asciiread(Datadir+filnam3, -1, "float") 
  strtmc  = asciiread(Datadir+filnam4, -1, "float") 
  time    = strfac(0:1800)
  fac     = strfac(1801:)
  tma     = strtma(1801:)*dtime 
  tmb     = strtmb(1801:)*dtime
  tmc     = strtmc(1801:)*dtime
  tmabc   = (tma + tmb + tmc)*1e8
  tmc     = tmc * 1e8
  labels  = (/"  (~F14~A~N~~F~+~F14~B~N~~F~+~F14~C~N~~F~)~F33~D~F10~t~F~","  ~F14~C~N~~F~~F33~D~F10~t~F~"," ~F14~f"/)

  tmin = -3.0 ;min((/tmabc,tmc/))          
  tmax =  3.0 ;max((/tmabc,tmc/))             ; of axes in advance.
  pmin = 0.07 ;min((/fac/))
  pmax = 0.13 ;max((/fac/))

;---Start the graphics

  wks = gsn_open_wks("pdf","SGR-P[0,0,0]_problem_cell_a")    ; send graphics to PNG file

  gsn_define_colormap(wks,"StepSeq25")
  color   = (/4,9,18/); (/"palegreen4","orange4","navy"/)
  markers = (/16,11,8/)
  mkszf   = (/0.8,0.8,0.8/)
  dashind = (/0,2,7/)
  lnthick = (/3.0,3.0,3.0/)

;---Set common resources for all plots
  res                  = True
  res@gsnMaximize      = False
 ;res@vpWidthF         = 0.8            ; Change the aspect ratio, but 
 ;res@vpHeightF        = 0.4            ; make plot as large as possible.

; Uncomment these if you don't want to draw the first two plots
  res@gsnDraw          = False
  res@gsnFrame         = False

  ;res@tmYROn = False            ; Turn off right tickmarks.
  ;res@tmXTOn = False            ; Turn off top tickmarks
  FontHeight = 0.025
  res@gsnStringFontHeightF  = FontHeight
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight
  res@tiYAxisFontHeightF    = FontHeight
  ;res@tiMainFontThicknessF  = 8.0
  ;res@tiYAxisFontThicknessF = 5.0
  ;res@tiXAxisFontThicknessF = 5.0

  res@tiMainString      = ""
  res@trXMinF          = 0           ; Set min/max of left Y axis
  res@trXMaxF          = 7200
  res@tmXBMode         = "Explicit"
  res@tmXBValues       = ispan(res@trXMinF,res@trXMaxF,1800)
  res@tmXBLabels       = ""+res@tmXBValues 
  res@tmXBMinorValues  = ispan(res@trXMinF,res@trXMaxF,300)
  res@tiXAxisString    = "Time (s)"

  res@gsnYRefLine            = 0.0             ; create a reference line
  res@gsnYRefLineDashPattern = 5
  res@gsnYRefLineColor       = "Grey"
  res@gsnYRefLineThicknessF  = 1.0


; resources for "left" variable
  resL                  = res
  resL@trYMinF          = tmin           ; Set min/max of left Y axis
  resL@trYMaxF          = tmax
  resL@xyLineColor      = color(0)
  resL@xyLineThicknessF = lnthick(0)
  resL@xyLineThicknessF = lnthick(0)
  resL@tiYAxisString    = "Concentration (10~S~-8~N~ kg kg~S~-1~N~)"
  ;resL@tiYAxisFontColor = resL@xyLineColor
  
; resources for "right" variable
  resR                    = res                      
  resR@trYMinF            = pmin           ; Set min/max of right Y axis
  resR@trYMaxF            = pmax
  resR@xyLineColor        = color(2)
  resR@xyDashPattern      = dashind(2)              ; Dashed line for 2nd curve
  resR@xyLineThicknessF   = lnthick(2)
  resR@tiYAxisString      = "Cloud fraction"
  resR@tiYAxisFontColor   = resR@xyLineColor
  resR@tmYRLabelFontColor = resR@xyLineColor

 ;resR@tmYRFormat       = "f"            ; Remove trailing zeros from labels

  plot = gsn_csm_xy2(wks,time,tmabc,fac,resL,resR)

;---Add the third line to the plots by using overlay
  res@xyLineColor       = color(1)
  res@xyDashPattern     = dashind(1)                    ; different dash line but same color
  res@xyLineThicknessF  = lnthick(1)
  overlay_plot = gsn_csm_xy(wks,time,tmc,res)
  
  overlay(plot,overlay_plot)

;Attach a legend
;;manual legend
   lgres                    = True
   lgres@lgItemType         = "Lines"
   lgres@lgLineLabelFontThicknessF= 1.0
  ;lgres@lgLineDashSegLenF  = 0.4
  ;lgres@lgBoxMinorExtentF  = 0.8
   lgres@vpWidthF           = 0.22           ; width of legend (NDC)
   lgres@vpHeightF          = 0.16           ; height of legend (NDC)
   lgres@lgLabelFontHeightF = 0.16            ; legend label font thickness
   lgres@lgMonoDashIndex    = False
   lgres@lgMonoLineThickness= False
   lgres@lgPerimOn          = True
  ;lgres@lgRightMarginF     = 0.5
  ;lgres@lgLeftMarginF      = 0.5
  ;lgres@lgTopMarginF       = 0.5
  ;lgres@lgBottomtMarginF   = 0.5
  
   lgres@lgPerimThicknessF   = 1.0
   lgres1 = lgres
   lgres1@lgLineColors       = color
   lgres1@lgDashIndexes      = dashind
   lgres1@lgLineThicknesses  = lnthick
   legend1 = gsn_create_legend (wks,dimsizes(labels), labels,lgres1)

;;;;;;;;;;;;
  amres                  = True
  amres@amJust           = "TopRight"    ; Use bottom right corner of box
  amres@amParallelPosF   =  0.2             ; Move legend to right
  amres@amOrthogonalPosF =  -0.45              ; Move legend down.
  annoid0 = gsn_add_annotation(plot,legend1,amres)  ; add legend to the top figure




  draw(plot)    ; This will draw all three curves
  frame(wks)

end
