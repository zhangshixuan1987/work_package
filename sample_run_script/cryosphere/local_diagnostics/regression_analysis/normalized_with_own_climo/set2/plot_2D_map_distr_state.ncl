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
  Regions   = (/"Polar_SH"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    = -20.
  MinLon    =   0.
  MaxLon    = 360.

  ;;;flags to determine which ERA5 datasets will be compared;;;; 
  l_era5_single_level   = True
  l_era5_pressure_level = False
  l_era5_surface_level  = False
  l_plot_era5_spread    = True ;; determine if we plot era5 ensemble spread
  l_plot_reference      = True ;; determine if we plot a reference climo difference between observations

  load "./NCL_experiment_information.ncl"

 ;;;;;;;date strings for the data directory;;;;;;;;
  Expnames  = Exp_short_names
  Groups    = Exp_long_names
  ngrps     = dimsizes(Groups)

  LabList     = Exp_short_names
  plot        = new(ngrps,graphic)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  load "./NCL_variable_list_state.ncl"
  VarList      = VarList3
  VstList      = VstList3
  VnmList      = VnmList3
  UnitList     = UnitList3
  nvars        = dimsizes(VarList)

 ;;; flags to control the plot style
  plotcntl          = False   ; if True, plot the CNTL origianl field on the leftmost panel

  plot_ttest_mskout = False  ; if True, mask the non-significant region as missing values

  plot_ttest_hatch  = False ; if True, plot the significant region as hatch region

  plot_ttest_hatrv  = False  ; if True, plot the non-significant region as hatch region

  plot_contour      = False   ; if True, plot the contour in the figure

  plot_leftstr      = True   ; if True, the lefadd igure captions as left string

  plot_msktrop      = False   ; if True, mask out the tropics region [-20,20]

 ;;color for the map boundary line
  lndcolor = "Grey28"

do iv = 0,nvars-1, 1

 Varname = VarList(iv)
 Varunt  = UnitList(iv)
 Varstr  = VstList(iv)

 plotFileLabel = Varname+""

 ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
 OUTDir  = "./"
 wtype = "png"
 wtype@wkWidth = 2500
 wtype@wkHeight = 2500 
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks1 = gsn_open_wks(wtype,OUTDir+"fig_regression_2d_map_"+plotFileLabel)

 ;;;set up the plot information;;;;
 load "./var_share_colorbar.ncl"

 do igp = 0, ngrps-1
 
  print(Exp_data_files(igp))
  File    = addfile(Exp_data_files(igp),"r")
  time    = tofloat(File->time)
  lat     = File->lat
  lon     = File->lon

  ntime   = dimsizes(time) ;; 3 years = 36 months
  nlat    = dimsizes(lat)  ;; latitude point
  nlon    = dimsizes(lon)  ;; longitutde point
 
  nyear   = ntime/12
  ncol    = 13 
  data    = asciiread(SAM_indx_files(igp),(/nyear,ncol/),"float")     ; 2 x 70 x 70
  samindx =  ndtooned(data(:,1:))
  if(igp.eq.0) then
   Missing_Value = File->mlspr@_FillValue
  end if 

 Var1   = new((/ntime,nlat,nlon/),double,Missing_Value)
 Var1!0 = "time"
 Var1!1 = "lat"
 Var1!2 = "lon"
 Var1&lat = lat
 Var1&lon = lon

 ;;directly read monthly mean data ;;;;
 if(igp.eq.0)then

  if (Varname.eq."U10")then
   Var1(:,:,:)    = (/File->si10/)
  end if 

  if (Varname.eq."TREFHT")then
   Var1(:,:,:)    = (/File->t2m/)
  end if

  if (Varname.eq."PSL")then
   Var1(:,:,:)    = (/File->msl/)
  end if 

  if (Varname.eq."LHFLX")then
   Var1(:,:,:)    = (/-File->slhf/)
   Var1           = Var1/86400.0
  end if

  if (Varname.eq."SHFLX")then
   Var1(:,:,:)    = (/-File->sshf/)
   Var1           = Var1/86400.0
  end if

  if (Varname.eq."PS")then
   Var1(:,:,:)    = (/File->sp/)
  end if

  if(Varname.eq."TS")then
   Var1(:,:,:)    = (/File->skt/)
  end if 

  if(Varname.eq."TMQ")then
   Var1(:,:,:)    = (/File->tcwv/)
  end if

  if(Varname.eq."TAUX")then
   Var1(:,:,:)    = (/-File->metss/)
  end if

  if(Varname.eq."TAUY")then
   Var1(:,:,:)    = (/-File->mntss/)
  end if

  if (Varname.eq."TAUXY")then
   Var1   = sqrt((/File->metss*File->metss  + File->mntss*File->mntss/))
  end if 

  if (Varname.eq."PRECT")then
   Var1   = (/File->mtpr/) 
   Var1   = Var1   *24.0*3600.0
  end if 

  if (Varname.eq."CLDTOT")then
   Var1   = (/File->tcc/)
  end if

 else

  ;;;;read E3SM data;;;;;;;;;
  if (Varname.eq."PRECT")then
   Var1   = (/File->PRECC  + File->PRECL/)
   Var1   = Var1   *24.0*3600.0*1000.0
  else if (Varname.eq."TAUXY")then
   Var1   = sqrt((/File->TAUX*File->TAUX  + File->TAUY*File->TAUY/))
  else
   Var1   = (/File->$Varname$/)
  end if 
  end if 

 end if  

 Var1@_FillValue = Missing_Value


 xClm   = clmMonTLL(Var1)                             ; (12,lat,lon)
 printVarSummary(xClm)

 xAnom  = calcMonAnomTLL ( Var1, xClm)                ; (time, lat,lon)
 printVarSummary(xAnom)
 printMinMax(xAnom, True)

; =================================================================
; Regress
; =================================================================

 rc   = regCoef_n(samindx, xAnom, 0, 0)
  printVarSummary(rc)
 df   = rc@nptxy-2   ; degrees of freedom
 tval = rc@tval      ; t-statistic
 b    = tofloat(df)          ; b must be same size as tval (and df)
 b    = 0.5
 prob  = betainc(df/(df+tval^2),df/2.0,b)
 printVarSummary(prob)

 mpdmpi = xClm(0,:,:)
 mpdmpi = (/rc/)
 printVarSummary(mpdmpi)
 printMinMax(mpdmpi, True)

 printVarSummary(mpdmpi)
 print(min(prob)+ " "+max(prob))

 if (plot_ttest_mskout.or.plot_ttest_hatch.or.plot_ttest_hatrv)then

   sigr = 0.05

   ;;;mask out non-significant point with FillValue;;;;;
   v1mnts1 = mpdmpi

   v1mnts1 = where(prob.le.sigr,mpdmpi, mpdmpi@_FillValue)

   if(plot_ttest_hatrv) then
     v1mnts1 = mpdmpi
     v1mnts1 = where(prob.gt.sigr,mpdmpi, mpdmpi@_FillValue)
   end if
   v1mnts1@_FillValue = mpdmpi@_FillValue

   delete([/prob/])

 end if 


;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False

    FontHeightF = 0.022
    res@tiMainFontHeightF        = FontHeightF*1.2
    res@tmYLLabelFontHeightF     = FontHeightF
    res@tmXBLabelFontHeightF     = FontHeightF

    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF/1.2
    res@gsnStringFontHeightF = FontHeightF*1.2

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61 
    res@pmLabelBarOrthogonalPosF = 0.15
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.55
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = Varstr+ " ("+Varunt+")";
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height


    leftstr                  = "" ;"Surface Temperature ("+ts@units+")"
    res@gsnLeftString        = "" ; leftstr ;var_name0    ; long_name is too long!
    res@gsnRightString       = "" ;" ("+Varunt+")";""
 ;  res@tiMaingString        = "Annual Mean Response to Present-day climate"

    ;---Control appearance of map.
   mpres@gsnPolarLabelFontHeightF = FontHeightF
   mpres@gsnPolar               = "SH"
   mpres@mpFillOn               = False        ; turn off map fill
   mpres@mpMaxLatF              = MaxLat
   mpres@mpCenterLonF            = 0
   mpres@mpLabelsOn              = False
   mpres@mpPerimOn               = True
   mpres@mpFillOn                = False
   mpres@mpGeophysicalLineColor  = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
   mpres@mpOutlineOn              = True
   mpres@mpOutlineDrawOrder       = "PostDraw"
   mpres@mpFillDrawOrder          = "Predraw"
   mpres@mpGridAndLimbOn          = True
   mpres@mpGridLineThicknessF     = 3.0
;  mpres@mpOceanFillColor         = "lightskyblue1"
;  mpres@mpLandFillColor          = "gray33"

    ;if (igp.eq.0) then
    ; mpres@tmYLLabelsOn           = True
    ; mpres@tmYLOn                 = True
    ; plot1(igp) = gsn_csm_map(wks1,mpres)
    ;else if (igp.lt.ngrps-1) then
    ; mpres@tmYLLabelsOn           = False
    ; mpres@tmYLOn                 = False
    ; plot1(igp) = gsn_csm_map(wks1,mpres)
    ;else
    ; mpres@tmYLLabelsOn           = False
    ; mpres@tmYLOn                 = False
    ; mpres@tmYROn                 = True
    ; plot(igp) = gsn_csm_map(wks1,mpres)
    ;end if
    ;end if

    mpres@tmYLLabelsOn           = True
    mpres@tmYLOn                 = True
    plot(igp) = gsn_csm_map(wks1,mpres)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = dcname
   ;res2@cnLevelSelectionMode          = "ManualLevels"
   ;res2@cnMinLevelValF                = mnmxint(0)
   ;res2@cnMaxLevelValF                = mnmxint(1)
   ;res2@cnLevelSpacingF               = mnmxint(2)
    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                       := diffgrp1
    res2@lbLabelStrings                 := sprintf("%5.2f",diffgrp1)
  ; res2@cnFillColors                   := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@tiMainString                   = LabList(igp) ;Figstr + LabList1(igp) ;"PD300 - CNTL"
    res2@tiMainOffsetYF                 = 0.01
    res2@gsnLeftString                  = " " ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
    res2@gsnLeftStringOrthogonalPosF     = 0.04
    ;res2@gsnLeftStringParallelPosF      = 0.05 
    res2@gsnRightString                  = " " ;bias+ "~C~" + rmse + "~C~"+pcrr
    res2@gsnRightStringOrthogonalPosF    = 0.04
    res2@gsnRightStringParallelPosF     = 0.95
    res2@gsnCenterString                 = "" ;"Global mean: "+ mdfstr 
    res2@gsnCenterStringOrthogonalPosF     = 0.04
    if(plot_ttest_mskout) then
     plot_rgd1                  = gsn_csm_contour(wks1,v1mnts1,res2)
    else
     plot_rgd1                  = gsn_csm_contour(wks1,mpdmpi,res2)  
    end if 

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
     plot_cont1                  = gsn_csm_contour(wks1,mpdmpi,res21)
     overlay(plot_rgd1,plot_cont1)
     delete([/res21/])
    end if
    overlay(plot(igp),plot_rgd1)
    delete([/res2/])

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
     str=unique_string("tmp1")
    ;wrf_smooth_2d(v1mnts1,3)
     plot@str                = gsn_csm_contour(wks1,v1mnts1,res3)
     overlay(plot(igp),plot@str)
     delete([/res3/])
    end if


   ; ;figstring = ""+rmse + "~C~"+ "    "+bias+ "~C~" + "    "+pcrr
   ; tres                      = True                ; text mods desired
   ; tres@txFontHeightF        = FontHeightF         ; make smaller
   ; txid1 = gsn_create_text(wks1, figstring, tres)

   ; amres = True
   ; amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
   ; amres@amOrthogonalPosF =  0.6    ; This is the bottom edge of the plot.
   ; annoid1 = gsn_add_annotation(plot(igp+2), txid1, amres)
   ; delete([/res,time/])
   delete([/mpres,res/])
 end do ;end of experiment group loop 

;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = True
    pnres@gsnPanelLabelBar        = False
    pnres@gsnPanelMainFontHeightF = 0.016
    pnres@gsnPanelMainString      = VnmList(iv) ;+ " ("+SEASON_STRING+")"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@lbOrientation           = "vertical"        ; vertical label bar
    pnres@pmLabelBarHeightF       = 0.12
    pnres@pmLabelBarWidthF        = 0.03
    pnres@lbLabelFontHeightF      = 0.006               ; make labels smaller
    pnres@pmLabelBarParallelPosF  = -0.005
    pnres@pmLabelBarOrthogonalPosF= 0.005
    pnres@lbTitleString           = Varstr +"" + "  ("+ Varunt + ")" 
    pnres@lbTitlePosition         = "Right"                           ; title location
    pnres@lbTitleDirection        = "Across"                          ; letter angle
    pnres@lbTitleAngleF           = 90.                               ; title angle
    pnres@lbTitleFontHeightF      = 0.006                              ; font height
    pnres@gsnPanelBottom          = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5

    pnres@gsnPanelBottom      = 0.4        ; move bottom up from 0.0 to 0.4
    gsn_panel(wks1,plot,(/1,3/),pnres)

   delete([/pnres/])

end do ; end of var loop

end

