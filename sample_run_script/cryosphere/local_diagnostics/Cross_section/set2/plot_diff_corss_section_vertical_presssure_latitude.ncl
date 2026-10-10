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

;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

  ;;;the time for the plot;;;
  Plot_time = (/"ANN","JJA","DJF"/)
  nplttim   = dimsizes(Plot_time)
  
  TimeString = " (2010-2012)"

  ;;;flags to determine which ERA5 datasets will be compared;;;;
  l_era5_single_level   = False
  l_era5_pressure_level = True
  l_era5_surface_level  = False
  load "../../share_script/NCL_experiment_information.ncl"

 ;;;;;;;strings for the data directory;;;;;;;;
  Expnames    = Exp_short_names
  Groups      = Exp_long_names
  ngrps       = dimsizes(Groups)
  LabList     = Exp_short_names
  LabList(1:) = LabList(1:) + " ~F19~:~N~~F~ " + LabList(0)
  LabList(0)  = LabList(0)+TimeString

  Expnames    = Exp_short_names
  Groups      = Exp_long_names
  ngrps       = dimsizes(Groups)
  LabList     = Exp_short_names
  LabList(0)  = Exp_short_names(0)
  LabList(1:) = Exp_short_names(1:) + " ~F19~:~N~~F~ " + Exp_short_names(0)
  LabList(0)  = LabList(0)+TimeString

  plot        = new(ngrps,graphic)

;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList      = (/"T","Q","U"/)
  VstList      = (/"T","Q","U"/)
  VnmList      = (/"Temperature", "Specific humidity", "Zonal Wind"/)
  UnitList     = (/"K","g kg~S~-1~N~","m s~S~-1~N~"/)
  nvars        = dimsizes(VarList)

 ;;; flags to control the plot style
  plotcntl          = True   ; if True, plot the CNTL origianl field on the leftmost panel
  plot_contour      = True   ; if True, plot the contour in the figure  
  plot_leftstr      = True   ; if True, the lefadd igure captions as left string 
  l_studentttest    = False  ; if student-t test is applied
  l_equiv_size      = False  
  plot_ttest_hatch  = False  ; if True, plot the significant region as hatch region
  plot_ttest_hatrv  = False  ; if True, plot the non-significant region as hatch region

do ison = 0,nplttim-1

 SEASON_STRING = Plot_time(ison)

do iv = 0,nvars-1, 1

 Varname  = VarList(iv)
 Varunt   = UnitList(iv)
 Varstr   = VstList(iv)
 Vnmstr   = VnmList(iv)

 plotFileLabel = "pres_height_"+SEASON_STRING+"_"+Varname

 ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  OUTDir  = "./"
  wtype = "pdf"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks1 = gsn_open_wks(wtype,OUTDir+"Fig_2d_"+plotFileLabel)

do igp = 1,ngrps-1,1

 ;;;set up the plot information;;;;
 load "../../share_script/var_share_colorbar_preshgt.ncl"

;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  File0 = addfile(ERA5_HRES,"r")
  File1 = addfile(Exp_data_files(igp),"r")

  lat   = File0->lat
  lon   = File0->lon
  lev   = File0->level
  P0mb  = 1000.0
  x11   = File0->$str_lower(Varname)$(:,:,:,:) 
  y11   = vinth2p(File1->$Varname$,File1->hyam,File1->hybm,lev,File1->PS,1,P0mb,1,False)

  ;x11!0 = "time"
  ;x11!1 = "lev"
  ;x11!2 = "lat"
  ;x11!3 = "lon"
  
  ;x11&lat = lat
  ;x11&lon = lon
  printVarSummary(x11)
  printVarSummary(y11)

  if(Varname.eq."RELHUM") then
   x11 = x11/100.0
   y11 = y11/100.0
  end if

  if(Varname.eq."Q") then
   x11 = x11*1000.0
   y11 = y11*1000.0
  end if

  if (SEASON_STRING .eq. "ANN")then
   pncre0300 = dim_avg_n_Wrap(dim_avg_n_Wrap(x11,0),2)
   wncre0300 = dim_avg_n_Wrap(dim_avg_n_Wrap(y11,0),2)
  else
   xx1 = month_to_seasonN(x11,SEASON_STRING)
   xx2 = month_to_seasonN(y11,SEASON_STRING)
   pncre0300 = dim_avg_n_Wrap(dim_avg_n_Wrap(xx1(0,:,:,:,:),0),2)
   wncre0300 = dim_avg_n_Wrap(dim_avg_n_Wrap(xx2(0,:,:,:,:),0),2)
   delete([/xx1,xx2/])
  end if
   
  wncredf1  = wncre0300
  wncredf1  = wncre0300 - pncre0300
  wncredf2  = wncre0300

  sigr = 0.05                      ; critical sig lvl for r
  if (l_studentttest) then
    prob = wncredf1
    prob = 0.0
    rad  = 4.*atan(1.)/180.    ;
    wgty = cos(lat*rad)
    x12 = month_to_annual(x11,1)
    y12 = month_to_annual(y11,1)
    x1   = dim_avg_Wrap(x12)
    y1   = dim_avg_Wrap(y12)
    xtmp = x1(lev|:,lat|:,year|:)
    ytmp = y1(lev|:,lat|:,year|:)
    xAve = dim_avg (xtmp)              ; calculate means at each grid point
    yAve = dim_avg (ytmp)
    xVar = dim_variance (xtmp)         ; calculate variances
    yVar = dim_variance (ytmp)
    if(l_equiv_size)
     xEqv    = equiv_sample_size (xtmp, sigr,0)
     yEqv    = equiv_sample_size (ytmp, sigr,0)
     xN      = dim_avg_wgt_n(xEqv, wgty, 1, 0)    ; wgty could be gaussian weights
     yN      = dim_avg_wgt_n(yEqv, wgty, 1, 0)
     delete([/xEqv,yEqv/])
    else
     xN      = dimsizes(xtmp(0,0,:))
     yN      = dimsizes(ytmp(0,0,:))
    end if
    iflag   = False                        ; population variance similar
    prob    = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)
    delete([/xtmp,ytmp,x12,y12,x11,y11,x1,y1,xAve,yAve,xVar,yVar,xN,yN/])

    v1mnts1 = wncredf1 ;smth9_Wrap(wncredf2, 0.50,  0.25, True)
    v1mnts1 = where(prob.le.sigr,wncredf1, 32766)
    if(plot_ttest_hatrv) then
       v1mnts1 = wncredf1
       v1mnts1 = where(prob.gt.sigr,wncredf1, 32766)
    end if
    v1mnts1@_FillValue = 32766
    delete([/prob/])
 
  end if

  strlabels  = (/Varstr ,"~F33~D~N~~F~"+ Varname,"~F33~D~N~~F~"+ Varname/)
  nlabls  = dimsizes(strlabels)

;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False
  res@vpHeightF            = 0.4               ; change aspect ratio of plot
  res@vpWidthF             = 0.8

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"


  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
  res@lbLabelStride           = 1                 ; skip every other label
  res@lbLabelFontHeightF      = 0.024
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.4
  res@pmLabelBarWidthF        = 0.10
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15

  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"

 ;res@trYMinF                 = 100.0                   ; min value on y-axis
 ;res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/1000, 900, 800, 700, 600, 500, 400, 300, 200, 100/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings
  res@trYLog                  = False

  res@tiMainFontThicknessF    = 3.0
  res@tiMainFontHeightF       = 0.025
  res@tmYLLabelFontHeightF    = 0.023
  res@tmXBLabelFontHeightF    = 0.023
  res@tiXAxisFontHeightF      = 0.025
  res@tiYAxisFontHeightF      = 0.025
  res@gsnStringFontHeightF    = 0.025
  res@cnLineLabelFontHeightF  = 0.023
  res@lbLabelFontHeightF      = 0.023
  res@lbTitleFontHeightF      = 0.023                             ; font height

  leftstr                     = "" ;"Surface Temperature ("+ts@units+")"
  res@gsnLeftString           = Varstr+" ("+ Varunt +")" ; leftstr ;var_name0    ; long_name is too long!
  res@gsnRightString          = "" ;" ("+Varunt+")";""
 ;res@tiMaingString           = "Annual Mean Response to Present-day climate"

;---Resources for plotting original data
  if(plotcntl.and.igp.eq.1)then
    res1=res
    res1@cnFillPalette                  = ocname
    res1@cnLevelSelectionMode           = "ExplicitLevels"  
    res1@cnLevels                       = orglev
   ;res1@cnFillColors                   = orgcolors
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False  
   ;res1@lbLabelStride                  = 15
    res1@gsnLeftString                  = LabList(0)
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.2
    res1@gsnRightString                 =  strlabels(0)+" ("+ Varunt +")"
    res1@gsnRightStringOrthogonalPosF   = 0.02
    res1@tiMainString                   = "";"CNTL (PD1800)"
    res1@cnLinesOn                      = plot_contour       
    res1@cnLineLabelsOn                 = False     
    res1@cnInfoLabelOn                  = False
    res1@cnLineLabelDensityF            = 2.0
    res1@cnLineLabelFontHeightF         = 0.018
    res1@cnLineLabelBackgroundColor     = -1
    res1@gsnContourLineThicknessesScale = 1.0
    res1@cnLineThicknessF               = 1.0
    res1@cnLineColor                    = "Black"
    res1@gsnContourZeroLineThicknessF   = 2.
   ;res1@gsnContourNegLineDashPattern   = 1
    res1@lbTitleString                  = "" ;"("+ Varunt+")" 
   ;plot(0)  = gsn_csm_contour(wks1,pncre0300,res1)
    plot(0)  = gsn_csm_pres_hgt(wks1,pncre0300,res1)
    delete([/res1/])
  end if 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@cnLevelSelectionMode           = "ExplicitLevels"  
    res2@cnFillPalette                  = dcname
   ;res2@cnFillColors                   = colorgroup1
    res2@cnLevels                      := diffgrp1
    res2@cnLinesOn                      = False 
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@gsnContourZeroLineThicknessF   = 2.
    res2@gsnContourNegLineDashPattern   = 1
    res2@tiMainString                   = "";"CNTL (PD1800)"
    res2@gsnLeftString                  = LabList(igp)
    res2@gsnLeftStringOrthogonalPosF    = 0.02
   ;res2@gsnLeftStringParallelPosF      = 0.2 
    res2@gsnRightString                 = strlabels(1)+" ("+ Varunt +")"
    res2@gsnRightStringOrthogonalPosF   = 0.02
    res2@lbTitleString                  = "" ;"("+ Varunt +")"
    ;plot(igp)                  = gsn_csm_contour(wks1,wncredf1,res2)
     plot(igp)                  = gsn_csm_pres_hgt(wks1,wncredf1,res2)

    if(plot_contour)
     res21 = res2
     res21@gsnLeftString                  = "" 
     res21@gsnRightString                 = ""
     res21@gsnCenterString                = ""
     res21@cnFillOn                       = False
     res21@cnLinesOn                      = True
     res21@cnLineLabelsOn                 = False            ; turn off contour labels
     res21@cnInfoLabelOn                  = False
     res21@cnLineLabelDensityF            = 2.0
     res21@cnLineLabelFontHeightF         = 0.018
     res21@cnLineLabelBackgroundColor     = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF               = 1.0
     res21@cnLineColor                    = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
    ;plot_cont1                  = gsn_csm_contour(wks1,wncredf1,res21)
     plot_cont1                  = gsn_csm_pres_hgt(wks1,wncredf1,res21)
     overlay(plot(igp),plot_cont1)
     delete([/res21/])
    end if 
   ;overlay(mpid2,plot_rgd1)

   ;---Add test shading to plot
   if(plot_ttest_hatch.or.plot_ttest_hatrv)then
    res3=res
    res3@cnFillOn              = True     ; color plot desired
    res3@lbLabelBarOn          = False    ; turn off contour labels
    res3@cnLinesOn             = False    ; turn off contour lines
    res3@cnLineLabelsOn        = False    ; turn off contour labels
    res3@cnMissingValFillColor = -1
    res3@cnMonoFillColor       = True     ; default color is fground [black]
    res3@cnMonoFillPattern     = True     ; want multiple patterns
    res3@cnMonoFillScale       = True     ; want patterns w/ diff densities
    res3@cnFillColor           = FillColor
    res3@cnFillPattern         = FillPattern
    res3@cnFillScaleF          = FillScaleF
    res3@cnFillDotSizeF        = FillDotSizeF
    res3@cnFillOpacityF        = FillOpacityF
    res3@tiMainString          = "" ;label
    res3@cnInfoLabelOn         = False
    str1=unique_string("tmp1")
    str2=unique_string("tmp2")
   ;wrf_smooth_2d(v1mnts1,3)
    plot@str1                = gsn_csm_pres_hgt(wks1,wncredf1,res3) 
    plot@str2                = gsn_csm_pres_hgt(wks1,wncredf1,res3) ;gsn_csm_contour(wks1,v1mnts1,res3)
    overlay(plot(igp),plot@str1)
    overlay(plot(igp),plot@str2)
    delete([/res3/])
   end if
   delete([/res/])

 end do

;---Draw both plots in a panel
 pnres                         = True
 pnres@gsnFrame                = False
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
 pnres@gsnPanelXWhiteSpacePercent = 5
;pnres@gsnPanelXF = (/0.10,0.38,0.66/)    ; Adjust rightmost plots
 gsn_panel(wks1,(/plot/),(/5,1/),pnres)

;;************************************************
;; Draw text on plot using plot coordinates.
;;************************************************
; tres                      = True                ; text mods desired
; tres@txFontHeightF        = 0.011               ; make smaller
; gsn_text_ndc(wks1,LabList(0),0.19,0.95,tres)
; gsn_text_ndc(wks1,LabList(1),0.48,0.95,tres)
; gsn_text_ndc(wks1,LabList(2),0.78,0.95,tres)

 frame(wks1)                                     ; now advance frame

 end do ;end of variable loop

end do ; end of season loop

end
