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
load "../share_script/NCL_experiment_information.ncl"

begin
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Polar_SH"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    = -55.
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
  load "./NCL_variable_list_cloud.ncl"
  VarList      = VarList1
  VstList      = VstList1
  VnmList      = VnmList1
  UnitList     = UnitList1
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
  wks1 = gsn_open_wks(wtype,OUTDir+"fig_pd_minus_pi_"+plotFileLabel)

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
 if (Varname.eq."CLDTOT")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->tcc/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->tcc/)
 end if 

 if (Varname.eq."CLDLOW")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->lcc/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->lcc/)
 end if

 if (Varname.eq."CLDHGH")then
  Var1(:,0,:,:)    = (/FERA5_SIG_ENM->hcc/)
  Var1(:,1:,:,:)   = (/FERA5_SIG_E10->hcc/)
 end if 

 Var2   = (/FE3SM_NE30->$Varname$/)
 Var3   = (/FE3SM_NE120->$Varname$/)
 
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
 difmn2   = dim_avg_n_Wrap(difclim2,0)


 do igp = 1, ngrps-1,1

  Figstr  = "" ;PltNums(iv*(ngrps-1)+ igp-1)

  if (igp.eq.1)then
   pivar   = Clim1
   pdvar   = Clim1
   pdvar   = (/conform(Clim1,Clim2,(/1,2/))/)
   pdmpi   = difclim1
   mpdmpi  = difmn1
  else
   pivar   = Clim1
   pdvar   = Clim1
   pdvar   = (/conform(Clim1,Clim3,(/1,2/))/)
   pdmpi   = difclim2
   mpdmpi  = difmn2
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

  ;;;;calculate the mean and standard deviation;;;;;;;
  latrg     = lat({MinLat:MaxLat})
  rad       = 4.*atan(1.)/180.    ;
  cost      = cos(latrg*rad)       ; cosine weights
  pimean    = wgt_areaave (pivar(:,{MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  pdmpim    = wgt_areaave (pdmpi(:,{MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mpistr    = sprintf("%5.2f", dim_avg_n(pimean,0)) 
  mdfstr    = sprintf("%5.2f", dim_avg_n(pdmpim,0)) 

  print(mpistr)
  print(mdfstr)

 
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
    res@tiMainFontHeightF    = 0.022
    res@tmYLLabelFontHeightF = 0.022
    res@tmXBLabelFontHeightF = 0.022

    mpres                    = res

    res@lbLabelFontHeightF   = 0.018
    res@gsnStringFontHeightF = 0.022

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
    ;res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
    ;res@lbTitlePosition         = "Right"                           ; title location
    ;res@lbTitleDirection        = "Across"                          ; letter angle
    ;res@lbTitleAngleF           = 90.                               ; title angle
    ;res@lbTitleFontHeightF      = 0.022                              ; font height


    leftstr                  = "" ;"Surface Temperature ("+ts@units+")"
    res@gsnLeftString        = "" ; leftstr ;var_name0    ; long_name is too long!
    res@gsnRightString       = "" ;" ("+Varunt+")";""
 ;  res@tiMaingString        = "Annual Mean Response to Present-day climate"

    ;---Control appearance of map.
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
     plot(igp) = gsn_csm_map(wks1,mpres)
    ;end if
    ;end if

    mpres@tmYLLabelsOn           = True
    mpres@tmYLOn                 = True
    plot(igp) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
  if(plotcntl)then
    plot(0)                   = gsn_csm_map(wks1,mpres)
    res1=res
    res1@gsnAddCyclic         = True
    res1@cnFillPalette                   = ocname
   ;res1@cnLevelSelectionMode            = "ManualLevels"
   ;res1@cnMinLevelValF                  = mnmxint(0)
   ;res1@cnMaxLevelValF                  = mnmxint(1)
   ;res1@cnLevelSpacingF                 = mnmxint(2)
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                       := orglev
   ;res1@cnFillColors                   := orgcolors
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = LabList(0) ;Figstr + LabList0(igp)
    res1@gsnLeftString                  = ""
    res1@gsnLeftStringOrthogonalPosF    = 0.02
   ;res1@gsnLeftStringParallelPosF      = 0.05
    res1@gsnRightString                 = "";
    res1@gsnRightStringOrthogonalPosF   = 0.02
   ;res1@gsnRightStringParallelPosF     = 0.9
    res1@gsnCenterString                 = "Global mean: "+ mpistr
    res1@gsnCenterStringOrthogonalPosF     = 0.04
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

    txres                  =  True
    txres@txFontHeightF    = 0.014
    txres@txFontColor      = resb@gsLineColor
    txres@txFontThicknessF = 10 ;resb@gsLineThicknessF
    txres@txJust           = "CenterCenter"
    overlay(plot(0),plot_org1)
    delete([/res1/])
  

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
  ; res2@cnFillColors                   := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@tiMainString                   = LabList(igp) ;Figstr + LabList1(igp) ;"PD300 - CNTL"
    res2@gsnLeftString                  = "" ;Varstr +"" + "  ("+ Varunt + ")" ;Groups(igp) + " ~F19~:~N~~F~ " + Groups(0)
    res2@gsnLeftStringOrthogonalPosF     = 0.04
    ;res2@gsnLeftStringParallelPosF      = 0.05 
    res2@gsnRightString                  = "Mean: "+mdfstr 
    res2@gsnRightStringOrthogonalPosF    = 0.04
    ;res2@gsnRightStringParallelPosF     = 0.9
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

    tres                      = True                ; text mods desired
    tres@txFontHeightF        = 0.020               ; make smaller
    txid1 = gsn_create_text(wks1, Varunt, tres)
    txid2 = gsn_create_text(wks1, Varunt, tres)
    txid3 = gsn_create_text(wks1, "%", tres)
    amres = True
    amres@amParallelPosF   =  0.43    ; This is the right edge of the plot.
    amres@amOrthogonalPosF = -0.5    ; This is the bottom edge of the plot.
    ;annoid1 = gsn_add_annotation(mpid1, txid1, amres)
    ;annoid2 = gsn_add_annotation(mpid2, txid2, amres)
    ;annoid3 = gsn_add_annotation(mpid3, txid3, amres)
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
    pnres@gsnPanelYWhiteSpacePercent = 2

    ;if (plot_zonal_mean) then
    ;  pnres@gsnPanelXF = (/0.05,0.24,0.43,0.495/)    ; Adjust rightmost plots
    ;else
    ;  pnres@gsnPanelXF = (/0.05,0.34,0.63/)    ; Adjust rightmost plots
    ;end if

    pres1  = pnres
    pres1@gsnFrame  = False      ; don't advance frame yet
    pres1@gsnPanelBottom      = 0.4        ; move bottom up from 0.0 to 0.4
    gsn_panel(wks1,plot,(/1,ngrps/),pres1)

  ;  pres2  = pnres
  ;  pres2@gsnFrame  = False 
  ; ;pres2@gsnPanelTop         = 0.4        ; draw up to the bdry of upper plot
  ;  pres2@gsnPanelBottom      = 0.1        ; move bottom up so size is 0.3
  ;  gsn_panel(wks1,plot,(/1,ngrps/),pres2)

  ; now advance frame for all plots
    frame(wks1)
   
   delete([/pnres,pres1/])

end do ; end of var loop
end

