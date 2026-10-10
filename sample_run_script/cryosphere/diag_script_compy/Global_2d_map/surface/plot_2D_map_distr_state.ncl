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
load "../../share_script/NCL_experiment_information.ncl"

begin
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Polar_SH"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    = -20.
  MinLon    =   0.
  MaxLon    = 360.

 ;;;;;;;date strings for the data directory;;;;;;;;
  Expnames  = Exp_short_Names
  Groups    = Exp_long_names
  ngrps     = dimsizes(Groups)

  LabList     = Exp_short_Names
  LabList(0)  = Exp_short_Names(0)
  LabList(1:) = Exp_short_Names(1:) + " ~F19~:~N~~F~ " + Exp_short_Names(0)
  plot        = new(ngrps+1,graphic)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  load "./NCL_variable_list_state.ncl"
  VarList      = VarList3
  VstList      = VstList3
  VnmList      = VnmList3
  UnitList     = UnitList3
  nvars        = dimsizes(VarList)

 ;;; flags to control the plot style
  plotcntl          = True   ; if True, plot the CNTL origianl field on the leftmost panel

  plot_ttest_mskout = False  ; if True, mask the non-significant region as missing values

  plot_ttest_hatch  = True  ; if True, plot the significant region as hatch region

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
  wks1 = gsn_open_wks(wtype,OUTDir+"fig_ANN_2d_map_"+plotFileLabel)

 ;;;set up the plot information;;;;
 load "./var_share_colorbar.ncl"

 FERA5_SFC_ENM  = addfile(ERA5_SFC_LEV,"r")
 FERA5_SIG_ENM  = addfile(ERA5_SIG_LEV,"r")
 FERA5_PRE_ENM  = addfile(ERA5_PRE_LEV,"r")
 FERA5_SIG_E10  = addfile(ERA5_SIG_LEV_ENS,"r")
 FERA5_PRE_E10  = addfile(ERA5_PRE_LEV_ENS,"r")
 FE3SM_NE30     = addfile(E3SM_NE30,"r")
 FE3SM_NE120    = addfile(E3SM_NE120,"r")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 time   = FE3SM_NE30->time
 lat    = FE3SM_NE30->lat
 lon    = FE3SM_NE30->lon

;;;;construct  the array to save the data;;;;;;;;;
 nens   = 11   ;; 1 hirs analysis + 10 ensembles
 ntime  = dimsizes(time) ;; 3 years = 36 months
 nlat   = dimsizes(lat)  ;; latitude point
 nlon   = dimsizes(lon)  ;; longitutde point 
 Var1   = new((/ntime,nens,nlat,nlon/),double)
 Var2   = new((/ntime,nlat,nlon/),double)
 Var3   = new((/ntime,nlat,nlon/),double)
 Var1!0 = "time"
 Var1!1 = "nens"
 Var1!2 = "lat"
 Var1!3 = "lon"
 
 Var1&lat = lat
 Var1&lon = lon

 ;;directly read monthly mean data ;;;;
 if (Varname.eq."U10")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->si10/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->si10/)
 end if 

 if (Varname.eq."TREFHT")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->t2m/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->t2m/)
 end if

 if (Varname.eq."PSL")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->msl/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->msl/)
 end if 

 if (Varname.eq."LHFLX")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->mslhf/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->mslhf/)
 end if

 if (Varname.eq."SHFLX")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->msshf/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->msshf/)
 end if

 if (Varname.eq."PS")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->sp/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->sp/)
 end if

 if(Varname.eq."TS")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->skt/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->skt/)
 end if 

 if(Varname.eq."TMQ")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->tcw/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->tcw/)
 end if


  if (Varname.eq."PRECT")then
   Var2   = (/FE3SM_NE30->PRECC  + FE3SM_NE30->PRECL/)
   Var3   = (/FE3SM_NE120->PRECC + FE3SM_NE120->PRECL/)
  else
   Var2   = (/FE3SM_NE30->$Varname$/)
   Var3   = (/FE3SM_NE120->$Varname$/)
  end if 
 
  if (Varname.eq."PSL".or.Varname.eq."PS")then
   Var1 = Var1/100.0
   Var2 = Var2/100.0
   Var3 = Var3/100.0 
  end if 

 ;;;;;; calculate the  3-year mean climo;;;;
 Clim1  = dim_avg_n_Wrap(Var1,0) 
 Clim2  = dim_avg_n_Wrap(Var2,0)
 Clim3  = dim_avg_n_Wrap(Var3,0)

 difclim1 = Clim1
 difclim2 = Clim1
 difclim1 = (/conform(Clim1,Clim2,(/1,2/)) - Clim1/)
 difclim2 = (/conform(Clim1,Clim3,(/1,2/)) - Clim1/)
 ;printVarSummary(difclim1)
 ;printVarSummary(difclim2)

 ;;calculate the ensemble mean difference;;;;;
 mpivar   = dim_avg_n_Wrap(Clim1,0)
 spivar   = dim_stddev_n_Wrap(Clim1,0)
 difmn1   = dim_avg_n_Wrap(difclim1,0)
 difsd1   = dim_stddev_n_Wrap(difclim1,0)
 difmn2   = dim_avg_n_Wrap(difclim2,0)
 difsd2   = dim_stddev_n_Wrap(difclim2,0)

 do igp = 1, ngrps-1,1

  Figstr  = "" ;PltNums(iv*(ngrps-1)+ igp-1)

  if (igp.eq.1)then
   pivar   = Clim1
   pdvar   = Clim1
   pdvar   = (/conform(Clim1,Clim2,(/1,2/))/)
   pdmpi   = difclim1
   mpdmpi  = difmn1
   spdmpi  = difsd1
  else
   pivar   = Clim1
   pdvar   = Clim1
   pdvar   = (/conform(Clim1,Clim3,(/1,2/))/)
   pdmpi   = difclim2
   mpdmpi  = difmn2
   spdmpi  = difsd2
  end if 

  ;;do a statistical test on the difference;;;;;;
  sigr  = 0.05                        ; critical sig lvl for r
  prob  = (/difmn1/)
  prob  = sigr
  rad   = 4.*atan(1.)/180.    ;
  wgty  = cos(lat*rad)
  xtmp  = pivar(lat|:,lon|:,nens|:)
  ytmp  = pdvar(lat|:,lon|:,nens|:)
  xAve  = dim_avg (xtmp)              ; calculate means at each grid point
  yAve  = dim_avg (ytmp)
  xVar  = dim_variance (xtmp)         ; calculate variances
  yVar  = dim_variance (ytmp)
  xEqv  = equiv_sample_size (xtmp, sigr,0)
  yEqv  = equiv_sample_size (ytmp, sigr,0)
  iflag = False                        ; population variance similar
  prob  = ttest(xAve,xVar,xEqv, yAve,yVar,yEqv, iflag, False)
  delete([/xtmp,ytmp,xAve,yAve,xVar,yVar,xEqv,yEqv/])

  ;;;;calculate the mean and spread for ERA5;;;;;;;
  latrg     = lat({MinLat:MaxLat})
  rad       = 4.*atan(1.)/180.    ;
  cost      = cos(latrg*rad)       ; cosine weights
  mpimean   = wgt_areaave (mpivar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  spimean   = wgt_areaave (spivar({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mpistr    = "Mean: " + sprintf("%5.2f", mpimean)
  spistr    = "Mean: " + sprintf("%5.2f", spimean)
  
  ;;;;;;;;;calculate mean bias, RMSE and pattern correlation;;;;;;;
  pdmpim    = wgt_areaave (mpdmpi({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpis    = wgt_areaave (spdmpi({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpirmse = wgt_arearmse(pivar(:,{MinLat:MaxLat},{MinLon:MaxLon}),pdvar(:,{MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpicorr = pattern_cor(pivar(:,{MinLat:MaxLat},{MinLon:MaxLon}),pdvar(:,{MinLat:MaxLat},{MinLon:MaxLon}),cost,0)
  bias      = "Bias: "+sprintf("%5.2f", pdmpim) + " ~F34~1~N~~F~ " + sprintf("%5.2f", pdmpis)   
  rmse      = "RMSE: "+sprintf("%5.2f", dim_avg(pdmpirmse)) + " ~F34~1~N~~F~ " + sprintf("%5.2f", dim_stddev(pdmpirmse)) 
  pcrr      = "Corr: "+sprintf("%5.2f", dim_avg(pdmpicorr)) + " ~F34~1~N~~F~ " + sprintf("%5.2f", dim_stddev(pdmpicorr)) 

  ;;;mask out non-significant point with FillValue;;;;;
  v1mnts1 = mpdmpi 
  v1mnts1 = where(prob.le.sigr,mpdmpi, 32766)

  if(plot_ttest_hatrv) then
    v1mnts1 = mpdmpi
    v1mnts1 = where(prob.gt.sigr,mpdmpi, 32766)
  end if
  v1mnts1@_FillValue = 32766

  if(isvar("prob"))then
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
   mpres@mpCenterLonF            = 180
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
     plot(igp+1) = gsn_csm_map(wks1,mpres)
    ;end if
    ;end if

    mpres@tmYLLabelsOn           = True
    mpres@tmYLOn                 = True
    plot(igp+1) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
  if(plotcntl)then
    plot(0)                   = gsn_csm_map(wks1,mpres)
    plot(1)                   = gsn_csm_map(wks1,mpres)
    res1=res
    res1@gsnAddCyclic         = True
    res1@cnFillPalette                   = ocname
   ;res1@cnLevelSelectionMode            = "ManualLevels"
   ;res1@cnMinLevelValF                  = mnmxint(0)
   ;res1@cnMaxLevelValF                  = mnmxint(1)
   ;res1@cnLevelSpacingF                 = mnmxint(2)
    res1@cnLevelSelectionMode           := "ExplicitLevels"
    res1@cnLevels                       := orglev
    res1@lbLabelStrings                 := sprintf("%5.2f",orglev)
   ;res1@cnFillColors                   := orgcolors
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = LabList(0) + " (ensemble mean)" ;Figstr + LabList0(igp)
    res1@tiMainOffsetYF                 = 0.01
    res1@gsnLeftString                  = ""
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.05
    res1@gsnRightString                 = ""; mpistr ;
    res1@gsnRightStringOrthogonalPosF   = 0.02
    res1@gsnRightStringParallelPosF     = 0.98
    res1@gsnCenterString                = "" ; mpistr
    res1@gsnCenterStringOrthogonalPosF  = 0.04
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
    res1@gsnContourNegLineDashPattern   = 1
    plot_org1  = gsn_csm_contour(wks1,mpivar,res1)
    overlay(plot(0),plot_org1)
    delete([/res1/])
  
    res1x=res
    res1x@gsnAddCyclic                   = True
    res1x@cnFillPalette                  = xcname
    res1x@cnLevelSelectionMode           = "ExplicitLevels"
    res1x@cnLevels                      := diffgrp0
   ;res1x@lbBoxLinesOn                   = True
   ;res1x@lbBoxSeparatorLinesOn          = False
   ;res1x@lbLabelStride                  = 15
    res1x@cnLinesOn                      = False             ; turn off contour lines
    res1x@cnLineLabelsOn                 = False            ; turn off contour labels
    res1x@cnInfoLabelOn                  = False
    res1x@tiMainString                   = LabList(0)+" (ensemble spread)" ;Figstr + LabList1(igp) ;"PD300 - CNTL"
    res1x@tiMainOffsetYF                 = 0.01
    res1x@gsnLeftString                  = "" ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
    res1x@gsnLeftStringOrthogonalPosF     = 0.04
   ;res1x@gsnLeftStringParallelPosF      = 0.05
    res1x@gsnRightString                  = "" ;spistr
    res1x@gsnRightStringOrthogonalPosF    = 0.04
   ;res1x@gsnRightStringParallelPosF     = 0.9
    res1x@gsnCenterString                 = "" ;"Global mean: "+ mdfstr
    res1x@gsnCenterStringOrthogonalPosF     = 0.04
   ;res1x@lbLabelStrings                 := sprintf("%5.2f",diffgrp0)
    plot_org2  = gsn_csm_contour(wks1,spivar,res1x)
    overlay(plot(1),plot_org2)
    delete([/res1x/])

    tres                      = True                ; text mods desired
    tres@txFontHeightF        = FontHeightF               ; make smaller
    txidx1 = gsn_create_text(wks1, mpistr, tres)
    txidx2 = gsn_create_text(wks1, spistr, tres)
    amres = True
    amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
    amres@amOrthogonalPosF =  0.5    ; This is the bottom edge of the plot.
    annoidx1 = gsn_add_annotation(plot(0), txidx1, amres)
    annoidx2 = gsn_add_annotation(plot(1), txidx2, amres)

   end if 

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
    overlay(plot(igp+1),plot_rgd1)
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
     overlay(plot(igp+1),plot@str)
     delete([/res3/])
    end if


    figstring = ""+rmse + "~C~"+ "    "+bias+ "~C~" + "    "+pcrr
    tres                      = True                ; text mods desired
    tres@txFontHeightF        = FontHeightF         ; make smaller
    txid1 = gsn_create_text(wks1, figstring, tres)

    amres = True
    amres@amParallelPosF   =  0.5    ; This is the right edge of the plot.
    amres@amOrthogonalPosF =  0.6    ; This is the bottom edge of the plot.
    annoid1 = gsn_add_annotation(plot(igp+1), txid1, amres)
    delete([/res/])
 end do ;end of experiment group loop 

;---Draw both plots in a panel
    pnres                         = True
    pnres@gsnMaximize             = True
    pnres@gsnPanelLabelBar        = False
    pnres@gsnPanelMainString      = "";Varstr; "Surface Temperature ("+ts@units+")"
    pnres@gsnPanelMainFontColor   = "Black"
    pnres@lbOrientation           = "vertical"        ; vertical label bar
    pnres@pmLabelBarHeightF       = 0.12
    pnres@pmLabelBarWidthF        = 0.03
    pnres@gsnPanelMainFontHeightF = 0.022
    pnres@lbLabelFontHeightF      = 0.006               ; make labels smaller
    pnres@pmLabelBarParallelPosF  = -0.005
    pnres@pmLabelBarOrthogonalPosF= 0.005
    pnres@lbTitleString           = Varstr +"" + "  ("+ Varunt + ")" 
    pnres@lbTitlePosition         = "Right"                           ; title location
    pnres@lbTitleDirection        = "Across"                          ; letter angle
    pnres@lbTitleAngleF           = 90.                               ; title angle
    pnres@lbTitleFontHeightF      = 0.006                              ; font height
    pnres@gsnPanelBottom          = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 5
    pnres@gsnPanelXWhiteSpacePercent = 5

    pnres@gsnPanelBottom      = 0.4        ; move bottom up from 0.0 to 0.4
    gsn_panel(wks1,plot,(/2,2/),pnres)

   delete([/pnres/])

end do ; end of var loop
end

