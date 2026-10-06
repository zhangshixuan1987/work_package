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
begin

  ldatalocal = True
  ldatalocal = True

  if (ldatalocal) then
     ModelOutRoot = "/pic/projects/E3SM_DA/zhan391/DART_Paper_Materials/data_scripts_for_report/data_for_plot/e3smv1_vs_e3smv0_obs_diag/"
  else
     ModelOutRoot = "/pic/projects/E3SM_DA/zhan391/e3sm_dart_NE30/diagnostics/EAM_NE30_Diags_2010-12_fnl/"
  end if

  Pltfil   = (/"EAMv0_vs_EAMv1_20ens_rmse_bias_with_respect_to_radiosonde.nc"/)
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Casdir       = (/"FV20",  "NE30"/)
  Groups       = (/"EAMv0_En20", "EAMv1_En20"/)
  ngrps        = dimsizes(Groups)
  labels       = (/"  EAMv0","  EAMv1"/)

  VarNameList = (/"U","V","T","Q"/)
  VarStrList  = (/"Zonal wind","Meridional wind", "Temperature","Specific humidity"/)
  UnitList    = (/"m s~S~-1~N~","m s~S~-1~N~","K","g kg~S~-1~N~"/)
  nvars       = dimsizes(VarNameList)

  regions  = (/"Tropics"/)
  nregs    = dimsizes(regions)

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 ;;There are four regions, but we only plot Northern Hemisphere in this script;;;
 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
 ;;;;x-aix range for each of the variable in the plot;;; 
  ;biases
  xmin1  = new(nvars, double)
  xmax1  = new(nvars, double)
  xmin1  = (/-1.2, -0.5, -2.0, -0.2/)
  xmax1  = (/ 1.8,  2.0,  1.0,  0.6/)

  ;rmse
  xmin2  = new(nvars, double)
  xmax2  = new(nvars, double)
  xmin2  = (/ 1.5,  1.5, 1.0,  0.0/)
  xmax2  = (/ 4.5,  4.5, 3.0,  1.0/)

  ;data rejection rate 
  xmin3  = new(nvars, double)
  xmax3  = new(nvars, double)
  xmin3  = 0.0
  xmax3  = 100.0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  plot        = new(3,graphic)

  do irg   = 0,nregs-1

  do iv = 0,nvars-1,1

  plotlabel =   VarNameList(iv)+"_"+str_sub_str(str_strip(regions(irg))," ","_")+"_profile"

  wks = gsn_open_wks("pdf",plotlabel)             ; send graphics to PNG file
  gsn_define_colormap(wks,"temp1")              ; choose colormap
  gsn_define_colormap(wks,"cb_9step")
  ;colors  = (/1,2,6,28,59/)
  ;gsn_merge_colormaps(wks,"StepSeq25","cb_9step")
  ;colors  = (/1,2,17,85,55/) ;(/1,55,2,17,85/)       ;;(/18,4,9,18,4,9/); (/"palegreen4","orange4","navy"/)
  colors  = (/2,6,28,1,36/)

  markers = (/16, 4, 3, 6, 16/)
  mkszf   = (/0.7, 0.9, 0.8, 0.8/)*2.0
  dashind = (/0,0,0,0,0/)       ;;line dash pattern
  lnthick = (/3,3,3,3,3/)/1.5   ;;line thichness
  mkthk   = mkszf
  mkthk   = 1.0

  fxx   = addfile(ModelOutRoot+Pltfil,"r")

  plevel  = fxx->plev
  nlevs = dimsizes(plevel)
  pltdat= new((/ngrps,3,nlevs/),double)

  do igp = 0, ngrps-1, 1

   varbias = Groups(igp)+"_"+VarNameList(iv)+"_bias"+"_"+str_sub_str(str_strip(regions(irg))," ","_")
   varrmse = Groups(igp)+"_"+VarNameList(iv)+"_rmse"+"_"+str_sub_str(str_strip(regions(irg))," ","_")
   varrate = Groups(igp)+"_"+VarNameList(iv)+"_rejection_rate"+"_"+str_sub_str(str_strip(regions(irg))," ","_")
   pltdat(igp,0,:)= (/fxx->$varbias$/)
   pltdat(igp,1,:)= (/fxx->$varrmse$/)
   pltdat(igp,2,:)= (/fxx->$varrate$/)
   timrange = fxx->$varrate$@timerange
  end do 

  mdat   = dim_avg_n(pltdat,2)
  pltdat@_FillValue = 32766
 
;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
  res          = True                   ; plot mods desired
  res@vpHeightF = 0.8
  res@vpWidthF =  0.5
  res@gsnDraw  = False
  res@gsnFrame = False

  FontHeight = 0.028
  res@gsnStringFontHeightF = FontHeight
  res@tmXBLabelFontHeightF = FontHeight
  res@tmYLLabelFontHeightF = FontHeight
  res@tmYRLabelFontHeightF = FontHeight
  res@tiMainFontHeightF    = FontHeight*1.1
  res@tiXAxisFontHeightF   = FontHeight
  res@tiYAxisFontHeightF   = FontHeight

  res@tmBorderThicknessF    = 1.5 
  res@tmXBMajorThicknessF   = 1.5 
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0

; add a legend
  res@pmLegendDisplayMode    = "Never"              ; turn on legend
  res@pmLegendSide           = "Top"                 ; Change location of
  res@pmLegendParallelPosF   = .75                   ; move units right
  res@pmLegendOrthogonalPosF = -0.35                  ; more neg = down
  res@pmLegendWidthF         = 0.14                  ; Change width and
  res@pmLegendHeightF        = 0.08                  ; height of legend.
  res@lgLabelFontHeightF     = .022                   ; change font height
  res@lgPerimOn              = False                 ; no box around

  res@trYReverse        = True                       ; reverse Y-axis
 ;res@xyYStyle          = "Log"
  res@tmYLMode          = "Explicit"                ; explicit labels
  res@tmYLValues        = (/ 1000,850,700,600,500,400,300,200,100,10/)
  res@tmYLLabels        = ""+res@tmYLValues         ; make strings
  res@tmXTOn  = False
  res@trYMaxF = 1000.0
  res@trYMinF = 10 ;-200.0

  res@gsnXRefLine             = (/0./)      ; four X reference lines
  res@gsnXRefLineThicknesses  = 1.0
  res@gsnXRefLineDashPatterns = 2
  res@gsnXRefLineColors       = "Black"

  res@xyMarkLineMode         = "MarkLines"
  res@xyLineColors           = colors(0:1)
  res@xyMarkers              = markers(0:1)
  res@xyMarkerColors         = colors(0:1)
  res@xyMarkerSizes          = 0.015*mkszf(0:1)
  res@xyMarkerThicknesses    = lnthick(0:1)
  res@xyDashPatterns         = dashind(0:1)
  res@xyLineThicknesses      = lnthick(0:1)

  res@tiMainString           = ""; Groups(igp) 
  res@tiMainOffsetYF         = 0.02
  res@tiMainFontThicknessF   = 3.0

  res@gsnRightString         = "" ;VstList(0)
  res@gsnLeftString          = "" ;VarStrList(iv); "Temperature (K)" ;Groups(0)

  res1 = res
  res1@tmYLOn                 = True
  res1@tiYAxisString          = "Pressure (hPa)";"log~B~10~N~(~F8~D~F10~t~F~)"
  res1@tiXAxisString          = "        "+VarStrList(iv) + "~C~" \
                                + ""+timrange
  res1@tiXAxisOffsetYF        = -0.015
  res1@gsnLeftString          = "Mean bias"; "Temperature (K)" ;Groups(0)
  res1@trXMaxF               := xmax1(iv)
  res1@trXMinF               := xmin1(iv)
  plot(0)  = gsn_csm_xy (wks,pltdat(:,0,:),plevel,res1)   ; create plot

  res2 = res
  res2@tmYLOn          = False
  res2@tiYAxisString   = ""
  res2@tiXAxisString   = "        "+VarStrList(iv) + "~C~" \
                           + ""+timrange
  res2@tiXAxisOffsetYF  = -0.015
  res2@gsnLeftString          = "Root mean square error"; "Temperature (K)" ;Groups(0)
  res2@trXMaxF                := xmax2(iv)
  res2@trXMinF                := xmin2(iv)
  plot(1)  = gsn_csm_xy (wks,pltdat(:,1,:),plevel,res2)   ; create plot

  res3 = res
  res3@xyLineColors           = colors(2:3)
  res3@xyMarkerColors         = colors(2:3)
  res3@tmYLOn                 = False
  res3@tiYAxisString          = ""
  res3@tiXAxisString          = "    Data rejection rate (%)"+ "~C~" \
                                  + ""+timrange
  res3@tiXAxisOffsetYF        = -0.015
  res3@gsnLeftString          = "Data rejection" ; "Temperature (K)" ;Groups(0)
  res3@trXMaxF                := xmax3(iv)
  res3@trXMinF                := xmin3(iv)
  plot(2)  = gsn_csm_xy (wks,pltdat(:,2,:),plevel,res3)   ; create plot
 

;Attach a legend
;;manual legend
   lgres                    = True
   lgres@lgItemType         = "MarkLines"
   lgres@lgLineLabelFontThicknessF= 1.0
 ; lgres@lgLineDashSegLenF  = 0.01
   lgres@lgBoxMinorExtentF  = 0.25
   lgres@lgBoxMajorExtentF  = 0.1
   lgres@lgMonoDashIndex    = False
   lgres@lgMonoLineThickness= False
   lgres@lgMonoMarkerThickness = False
   lgres@lgPerimOn          = False ;True
  ;lgres@lgRightMarginF     = 0.5
  ;lgres@lgLeftMarginF      = 0.5
  ;lgres@lgTopMarginF       = 0.5
  ;lgres@lgBottomtMarginF   = 0.5

   lgres@lgPerimThicknessF  = 1.0

  lgres1 = lgres
  lgres1@vpWidthF           = 0.30           ; width of legend (NDC)
  lgres1@vpHeightF          = 0.12           ; height of legend (NDC)
  lgres1@lgLineColors       = colors(0:1) ;(/"Black","Black","Black"/)
  lgres1@lgDashIndexes      = dashind(0:1)
  lgres1@lgLineThicknesses  = lnthick(0:1)
  lgres1@lgMarkerColors     = colors(0:1)
  lgres1@lgMarkerIndexes    = markers(0:1)
  lgres1@lgMarkerThicknesses= mkthk(0:1)

  lgres2 = lgres
  lgres2@vpWidthF           = 0.30           ; width of legend (NDC)
  lgres2@vpHeightF          = 0.12           ; height of legend (NDC)
  lgres2@lgLineColors       = colors(2:3) ;(/"Black","Black","Black"/)
  lgres2@lgDashIndexes      = dashind(0:1)
  lgres2@lgLineThicknesses  = lnthick(0:1)
  lgres2@lgMarkerColors     = colors(2:3)
  lgres2@lgMarkerIndexes    = markers(0:1)
  lgres2@lgMarkerThicknesses= mkthk(0:1)


  ;Labels for the legend
  labls = labels
  do ii= 0,1
   labls(ii) = labels(ii) + " ("+sprintf("%5.2f",mdat(ii,0))+")"
  end do
  legend1 = gsn_create_legend (wks,dimsizes(labls), labls,lgres1)

  labls = labels
  do ii= 0,1
   labls(ii) = labels(ii) + " ("+sprintf("%5.2f",mdat(ii,1))+")"
  end do
  legend2 = gsn_create_legend (wks,dimsizes(labls), labls,lgres1)

  labls = labels
  do ii= 0,1
   labls(ii) = labels(ii) + " ("+sprintf("%5.2f",mdat(ii,2))+")"
  end do
  legend3 = gsn_create_legend (wks,dimsizes(labls), labls,lgres2)

  amres1                  = True
  if(VarNameList(iv).eq."T") then 
    amres1@amJust           = "TopLeft"
    amres1@amParallelPosF   =  -0.48
    amres1@amOrthogonalPosF =  -0.30
  end if 
  if(VarNameList(iv).eq."Q") then
    amres1@amJust           = "TopRight"
    amres1@amParallelPosF   =  0.45
    amres1@amOrthogonalPosF =  -0.45
  end if

  if(VarNameList(iv).eq."U") then
    amres1@amJust           = "TopRight"
    amres1@amParallelPosF   =   0.48
    amres1@amOrthogonalPosF =  -0.35
  end if

  if(VarNameList(iv).eq."V") then
    amres1@amJust           = "TopRight"
    amres1@amParallelPosF   =   0.48
    amres1@amOrthogonalPosF =  -0.35
  end if

  annoid1 = gsn_add_annotation(plot(0),legend1,amres1)  ; add legend to the top figure

  amres2                  = True
  if(VarNameList(iv).eq."T") then
    amres2@amJust           = "BottomRight"
    amres2@amParallelPosF   =  0.40
    amres2@amOrthogonalPosF =  0.1
  end if
  if(VarNameList(iv).eq."Q") then
    amres2@amJust           = "BottomLeft"
    amres2@amParallelPosF   =  -0.40
    amres2@amOrthogonalPosF =   0.45
  end if

  if(VarNameList(iv).eq."U") then
    amres2@amJust           = "TopLeft"
    amres2@amParallelPosF   = -0.45
    amres2@amOrthogonalPosF = -0.30
  end if

  if(VarNameList(iv).eq."V") then
    amres2@amJust           = "TopLeft"
    amres2@amParallelPosF   = -0.45
    amres2@amOrthogonalPosF = -0.30
  end if

  annoid2 = gsn_add_annotation(plot(1),legend2,amres2)  ; add legend to the top figure

  amres3                  = True
  amres3@amJust           = "BottomRight"
  amres3@amParallelPosF   =  0.40
  amres3@amOrthogonalPosF =  0.1
  annoid3 = gsn_add_annotation(plot(2),legend3,amres3)  ; add legend to the top figure

;---Draw both plots in a panel
  pnres                         = True
  pnres@gsnMaximize             = False
  pnres@gsnPanelMainString      = "    "+VarNameList(iv) + " ("+ regions(irg)+")" 
  pnres@gsnPanelMainFontColor   = "Black"
  pnres@gsnPanelMainFontHeightF = 0.020
  pnres@gsnPanelLabelBar        = False
  pnres@gsnPanelBottom          = 0.05
  pnres@gsnPanelYWhiteSpacePercent = 2
  pnres@gsnPanelXF = (/0.15,0.45,0.75/)    ; Adjust rightmost plots
  gsn_panel(wks,(/plot/),(/1,3/),pnres)
  delete([/pnres/])

end do

end do 
end
