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

  fig_grp    = "nombias"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/state_var"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  LabList    = (/"SAM+","SAM-"/)
 
  vars       = (/"U","T"/)
  facmod     = (/1.0,1.0/)
  varunt     = (/"m s~S~-1~N~","~S~o~N~C"/)
  nvars      = dimsizes(vars)

  statistical_siglev = 0.05
  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  d2r        = get_d2r("float")

 ;declear the season or annual mean data to be used;;;;;
 do i = 0,nseas-1,1

  do k = 0,nexps-1,1

    print("working on season "+ seasons(i))

    fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+vars(0)+"*en*"+time_tag+".nc")
    fl2  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+vars(1)+"*en*"+time_tag+".nc")
    f1   = addfiles(fl1,"r")
    f2   = addfiles(fl2,"r")
    ListSetType (f1, "join")
    ListSetType (f2, "join")

    up  = dim_avg_n_Wrap(f1[:]->U_p,0)
    tp  = dim_avg_n_Wrap(f2[:]->T_p,0)
    pp  = dim_avg_n_Wrap(f1[:]->prob_p,0)

    pp  = where(pp.le.statistical_siglev,pp,pp@_FillValue)
    pp  = (1.0 - pp)*100.0

    un  = dim_avg_n_Wrap(f1[:]->U_n,0)
    tn  = dim_avg_n_Wrap(f2[:]->T_n,0)
    pn  = dim_avg_n_Wrap(f1[:]->prob_n,0)
    pn  = where(pn.le.statistical_siglev,pn,pn@_FillValue)
    pn  = (1.0 - pn)*100.0

    printVarSummary(up)
    
  ;----------------------------------------------------------------------
  ; Plotting section
  
  if (k.eq.0 ) then
     wtype           = "png"
     wtype@wkWidth   = 2500
     wtype@wkHeight  = 2500
     fignam          = "SAM_composite_"+seasons(i)+"_"+time_tag + "_"+fig_grp
     wks             = gsn_open_wks(wtype,fignam)
     plot      = new((/nexps*2/),graphic)
     dum       = new((/nexps*2/),graphic)
  end if 

  load "./var_share_colorbar_preshgt.ncl"

  ;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False

  res@vpWidthF             = 0.9
  res@vpHeightF            = 0.6
  res@gsnMajorLatSpacing   =  10             ;-- change major lat tickmark spacing
  res@gsnMinorLatSpacing   =  2.5            ;-- change major lat tickmark spacing

  FontHeightF = 0.032
  res@tiMainFontThicknessF    = 1.0
  res@tiMainFontHeightF       = FontHeightF*1.4
  res@tmYLLabelFontHeightF    = FontHeightF
  res@tmXBLabelFontHeightF    = FontHeightF
  res@tiXAxisFontHeightF      = FontHeightF*1.1
  res@tiYAxisFontHeightF      = FontHeightF*1.1
  res@gsnStringFontHeightF    = FontHeightF*1.1
  res@cnLineLabelFontHeightF  = FontHeightF
  res@lbLabelFontHeightF      = FontHeightF
  res@lbTitleFontHeightF      = FontHeightF   

   BorderThick = 3.0
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"

  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
 ;res@lbLabelStride           = 1                 ; skip every other label
  res@lbLabelFontHeightF      = FontHeightF*0.9
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.6
  res@pmLabelBarWidthF        = 0.12
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15

  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYLOn                  = True
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"
  res@tmYLOn                  = True

  res@trYMinF                 = 10.0                   ; min value on y-axis
  res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/1000, 700, 500, 300, 200, 100, 50, 30, 20, 10/) ;(/ 1000, 850, 700,500,400,300,200,100,10/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings

  res@trYLog                  = True

  res@gsnLeftString           = "" 
  res@gsnRightString          = "" 
  res@tiMainString            = ""

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  res1=res
  res1@cnLevelSelectionMode           = "ExplicitLevels"  
  res1@cnFillPalette                  = "cmp_b2r" ;"cmocean_balance" ;"BlueDarkRed18"
  res1@cnLevels                       = fspan(-1.0,1.0,21)
  res1@cnLinesOn                      = True
  res1@cnLineLabelsOn                 = False            ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@gsnContourZeroLineThicknessF   = 0.2
  res1@gsnContourNegLineDashPattern   = 0
  res1@cnLineThicknessF               = 0.2
  res1@lbTitleString                  = "" ;"("+ Varunt +")"

  res1@gsnLeftString  = exp_strs(k)
  res1@gsnRightString = "positive SAM"
  plot(k*2) = gsn_csm_pres_hgt(wks,tp,res1)

  res1@gsnLeftString  = exp_strs(k)
  res1@gsnRightString = "negative SAM"
  plot(k*2+1) = gsn_csm_pres_hgt(wks,tn,res1)

  res11 = res
  res11@gsnLeftString                  = "" 
  res11@gsnRightString                 = ""
  res11@gsnCenterString                = ""
  res11@cnFillOn                       = False
  res11@cnLinesOn                      = True
  res11@cnLineLabelsOn                 = True            ; turn off contour labels
  res11@cnLevelSelectionMode           = "ExplicitLevels"
  res11@cnLevels                       = fspan(-5,5,21)              ; set a different set of contour levels
  res11@cnInfoLabelOn                  = False
  res11@cnLineLabelDensityF            = 2.0
  res11@cnLineLabelFontHeightF         = FontHeightF 
  res11@cnLineLabelBackgroundColor     = -1
  res11@gsnContourLineThicknessesScale = 1.0
  res11@cnLineThicknessF               = 6.0
  res11@cnLineColor                    = "Black" ;Magenta"
  res11@cnLineLabelFontColor           = "Black"
  res11@gsnContourZeroLineThicknessF   = 12.0
  res11@gsnContourNegLineDashPattern   = 1
  plot_cont1 = gsn_csm_pres_hgt(wks,up,res11)
  plot_cont2 = gsn_csm_pres_hgt(wks,un,res11)

  overlay(plot(k*2),plot_cont1)
  overlay(plot(k*2+1),plot_cont2)

  rescn2 = True
  rescn2@lbLabelBarOn    = False
  rescn2@gsnTickMarksOn  = False     ; no tickmarks
  rescn2@gsnDraw         = False     ; don't draw
  rescn2@gsnFrame        = False     ; don't advance frame
  rescn2@tiMainString    = ""
  rescn2@gsnLeftString   = ""
  rescn2@gsnRightString  = ""
  rescn2@cnMinLevelValF  = (1.0 - statistical_siglev)*100.0
  rescn2@cnMaxLevelValF  = 100.0
 ;rescn2@cnLevelSpacingF = 1.
  rescn2@cnMonoFillPattern = True
  rescn2@cnMonoFillColor = True
  rescn2@cnFillOn        = True      ; color fill
  rescn2@cnLinesOn       = False
  rescn2@cnLineLabelsOn  = False        ; True is default
  rescn2@cnInfoLabelOn   = False                ; turn off info label
  rescn2@cnFillMode      = "AreaFill" ; raster fill significance patterns
  rescn2@cnFillPattern   = 17
  rescn2@cnFillDotSizeF  = 0.002
  rescn2@cnFillColor     = "black"
  rescn2@cnFillOpacityF  = 0.5

  plotA = gsn_csm_contour(wks,pp,rescn2) ; contours are at 95,100
  plotB = gsn_csm_contour(wks,pn,rescn2) ; contours are at 95,100

  overlay(plot(k*2),plotA)
  overlay(plot(k*2+1),plotB)

  delete([/res11,res/])
  delete([/fl1,fl2,f1,f2,up,tp,pp,un,tn,pn/])
 end do ;end of exp loop

;---Draw both plots in a panel
 pnres                         = True
 pnres@gsnFrame                = True
 pnres@gsnMaximize             = False
 pnres@gsnPanelLabelBar        = True
 pnres@gsnPanelMainString      = " "; "Surface Temperature ("+ts@units+")"
 pnres@gsnPanelLabelBar        = False                ; add common colorbar
 pnres@gsnPanelMainFontColor   = "Black"
 pnres@gsnPanelMainFontHeightF = 0.030
;pnres@lbLabelFontHeightF      = 0.018               ; make labels smaller
;pnres@pmLabelBarParallelPosF  = 0.026
;pnres@pmLabelBarOrthogonalPosF= 0.015
;pnres@pmLabelBarHeightF       = 0.1
;pnres@pmLabelBarWidthF        = 0.57
;pnres@gsnPanelBottom          = 0.05
 pnres@gsnPanelYWhiteSpacePercent = 2
 pnres@gsnPanelXWhiteSpacePercent = 1
 gsn_panel(wks,(/plot(:)/),(/nexps,2/),pnres)

end do 

end
