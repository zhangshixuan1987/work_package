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

  statistical_siglev = 0.05

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/trend_zonal"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
 
  varList    = (/"T","U"/)
  varFact    = (/10.0,10.0/)
  varUnit    = (/"~S~o~N~C per decade","m s~S~-1~N~ per decade"/)
  nvars      = dimsizes(varList)
  
  seasons    = (/"ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  d2r        = get_d2r("float")

  time_tag   = (/"1979-2014","1979-2001","2001-2014"/)
  panel_str  = (/"Depletion period", "Recovery period", "Change"/)
  np         = dimsizes(panel_str)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0,nseas-1,1

 do j = 0,nvars-1,1 

  varnam = varList(j)
  varunt = varUnit(j)
  varfac = varFact(j)

  do k = 0,nexps-1,1

    print("working on season "+ seasons(i))

    fl0  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_"+"*en*"+"_"+varnam+"_"+seasons(i)+"_trend_"+time_tag(0)+".nc")
    fl1  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_"+"*en*"+"_"+varnam+"_"+seasons(i)+"_trend_"+time_tag(1)+".nc")
    fl2  = systemfunc("ls "+data_dir+"/"+exp_name(k)+"_"+"*en*"+"_"+varnam+"_"+seasons(i)+"_trend_"+time_tag(2)+".nc")

    f0   = addfiles(fl0,"r")
    f1   = addfiles(fl1,"r")
    f2   = addfiles(fl2,"r")
    enstr= "enmn"

    ListSetType (f0, "join")
    ListSetType (f1, "join")
    ListSetType (f2, "join")

    vs0  = varnam+"_climo"
    vs1  = varnam+"_trend"
    vs2  = varnam+"_prob"
    
    uc1  = dim_avg_n_Wrap(f0[:]->$vs0$,0)
    ut1  = dim_avg_n_Wrap(f1[:]->$vs1$,0)
    ut2  = dim_avg_n_Wrap(f2[:]->$vs1$,0)

    ut1  = ut1 * varfac
    ut2  = ut2 * varfac

    ut1@units = varunt
    ut2@units = varunt
    lbstr     = varnam + "(" + varunt+")"

    utd  = ut1
    utd  = ut2 - ut1
    printMinMax(ut1,False)
    printMinMax(ut2,False)
    printMinMax(uc1,False)
 
  ;----------------------------------------------------------------------
  if (k.eq.0 ) then
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fignam          = "fig_trend_"+enstr+"_"+seasons(i)+"_"+varnam+"_"+time_tag(0)
    wks             = gsn_open_wks(wtype,fignam)
    plot            = new((/nexps*np/),graphic)
    dum             = new((/nexps*np/),graphic)
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
  res@tmYLValues              = (/1000, 700, 500, 300, 200, 100, 50, 30, 20, 10/) 
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
  res1@lbTitleString                  = lbstr
  res1@cnFillPalette                  = "cmp_b2r" ;"cmocean_balance" ;"BlueDarkRed18"
  res1@cnLevels                      := cnlev2

  res1@gsnLeftString                  = panel_str(0)
  res1@gsnRightString                 = exp_strs(k)
  plot(k*np) = gsn_csm_pres_hgt(wks,ut1,res1)

  res1@gsnLeftString                  = panel_str(1)
  res1@gsnRightString                 = exp_strs(k)
  plot(k*np+1) = gsn_csm_pres_hgt(wks,ut2,res1)

  res1@gsnLeftString                  = panel_str(2)
  res1@gsnRightString                 = exp_strs(k)
  plot(k*np+2) = gsn_csm_pres_hgt(wks,utd,res1)

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

  res2@cnLevels                      := cnlev1
  res2@gsnLeftString                  = ""
  res2@gsnRightString                 = ""
  dum(k*np)   = gsn_csm_pres_hgt(wks,uc1,res2)
  dum(k*np+1) = gsn_csm_pres_hgt(wks,uc1,res2)
  dum(k*np+2) = gsn_csm_pres_hgt(wks,uc1,res2)

  overlay(plot(k*np),dum(k*np))
  overlay(plot(k*np+1),dum(k*np+1))
  overlay(plot(k*np+2),dum(k*np+2))

  delete([/res,res1,res2/])
  delete([/uc1,ut1,ut2,utd/])
  delete([/fl0,fl1,fl2/])

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
  gsn_panel(wks,(/plot(:)/),(/nexps,np/),pnres)

end do 
end do 

end
