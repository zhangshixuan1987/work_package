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

  data_dir   = "../data"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
 
  vars       = (/"U","T"/)
  facmod     = (/1.0,1.0/)
  varunt     = (/"m s~S~-1~N~","~S~o~N~C"/)
  nvars      = dimsizes(vars)

  statistical_siglev = 0.05
  time_tag   = (/"1981-2010"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  d2r        = get_d2r("float")

 ;declear the season or annual mean data to be used;;;;;
 do i = 0,nseas-1,1

  do k = 1,nexps-1,1

    print("working on season "+ seasons(i))

    fl1  = systemfunc("ls "+data_dir+"/"+exp_name(0)+"*"+vars(0)+"*en*"+time_tag+".nc")
    fl2  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+vars(0)+"*en*"+time_tag+".nc")
    fl3  = systemfunc("ls "+data_dir+"/"+exp_name(0)+"*"+vars(1)+"*en*"+time_tag+".nc")
    fl4  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"*"+vars(1)+"*en*"+time_tag+".nc")

    f1   = addfiles(fl1,"r")
    f2   = addfiles(fl2,"r")
    f3   = addfiles(fl3,"r")
    f4   = addfiles(fl4,"r")

    ListSetType (f1, "join")
    ListSetType (f2, "join")
    ListSetType (f3, "join")
    ListSetType (f4, "join")

    ur  = dim_avg_n_Wrap(f1[:]->$vars(0)$,0)
    up  = dim_avg_n_Wrap(f2[:]->$vars(0)$,0)
    tr  = dim_avg_n_Wrap(f3[:]->$vars(1)$,0)
    tp  = dim_avg_n_Wrap(f4[:]->$vars(1)$,0)
    tp  = tp - 273.15
    tr  = tr - 273.15

    ud  = ur
    ud  = up - ur
    td  = tr
    td  = tp - tr

    printVarSummary(up)
    delete([/fl1,fl2,fl3,fl4,f1,f2,f3,f4/])

  ;----------------------------------------------------------------------
  ; Plotting section
  
  if (k.eq.1 ) then
     wtype           = "png"
     wtype@wkWidth   = 2500
     wtype@wkHeight  = 2500
     fignam          = "zonal_wind_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     plot      = new((/nexps/),graphic)
     dum       = new((/nexps/),graphic)
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
  res1@cnLinesOn                      = True
  res1@cnLineLabelsOn                 = False            ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@gsnContourZeroLineThicknessF   = 0.2
  res1@gsnContourNegLineDashPattern   = 0
  res1@cnLineThicknessF               = 0.2
  res1@lbTitleString                  = "" ;"("+ Varunt +")"

  if(k.eq.1) then 
    res1@cnFillPalette                = "MPL_rainbow" ;"cmocean_balance" ;"BlueDarkRed18"
    res1@cnLevels                     = (/-65,-60,-55,-50,-45,-40,-30,-20,-10,0,10,20/);fspan(-60.0,20.0,21)
    res1@gsnLeftString                = exp_strs(0)
    res1@gsnRightString               = "Zonal wind"
    plot(0) = gsn_csm_pres_hgt(wks,tr,res1)
  end if 

  res1@cnLevelSelectionMode           = "ExplicitLevels"
  res1@cnFillPalette                  = "cmp_b2r" ;"cmocean_balance" ;"BlueDarkRed18"
  res1@cnLevels                      := fspan(-5.0,5.0,11)
  res1@gsnLeftString                  = exp_strs(k)
  res1@gsnRightString                 = "Model biases"
  plot(k) = gsn_csm_pres_hgt(wks,td,res1)


  res2 = res
  res2@gsnLeftString                  = ""
  res2@gsnRightString                 = ""
  res2@gsnCenterString                = ""
  res2@cnFillOn                       = False
  res2@cnLinesOn                      = True
  res2@cnLineLabelsOn                 = True            ; turn off contour labels
  res2@cnLevelSelectionMode           = "ExplicitLevels"
  res2@cnInfoLabelOn                  = False
  res2@cnLineLabelDensityF            = 2.0
  res2@cnLineLabelFontHeightF         = FontHeightF
  res2@cnLineLabelBackgroundColor     = -1
  res2@gsnContourLineThicknessesScale = 1.0
  res2@cnLineThicknessF               = 6.0
  res2@cnLineColor                    = "Black" ;Magenta"
  res2@cnLineLabelFontColor           = "Black"
  res2@gsnContourZeroLineThicknessF   = 12.0
  res2@gsnContourNegLineDashPattern   = 1

  if(k.eq.1) then
    res2@cnLevels                     = fspan(-30.0,30.0,21)
    res2@gsnLeftString                = ""
    res2@gsnRightString               = ""
    plot1 = gsn_csm_pres_hgt(wks,ur,res2)
    overlay(plot(0),plot1)
  end if

  res2@cnLevels                      := fspan(-5.0,5.0,11)
  res2@gsnLeftString                  = ""
  res2@gsnRightString                 = ""
  plot2 = gsn_csm_pres_hgt(wks,ud,res2)
  overlay(plot(k),plot2)

  delete([/res,res1,res2/])
  delete([/ur,up,ud,tr,tp,td/])
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
 gsn_panel(wks,(/plot(:)/),(/nexps,1/),pnres)

end do 

end
