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

  l_plot_cntl       = True  ; if True, plot the CNTL origianl field on the leftmost panel
  plot_contour      = False   ; if True, plot the contour in the figure
  plot_leftstr      = True    ; if True, the lefadd igure captions as left string
  plot_sig_test     = True

  datadir     = "../fig_data/"
  mskfile     = "../fig_data/CLIM30_QBOT_1979-2014.nc"
  fmsk        = addfile(mskfile,"r")
  lmsk        = where(ismissing(fmsk->QBOT),0.0,1.0)

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
 ;ERA5 -> RD
 ;CLIM -> CR
 ;Nudged -> NU
 ;LSTM(CLIM) -> ML(CR)

  Groups      = (/"ERA5","CLIM30","CLIM30ML"/)
  ngrps       = dimsizes(Groups)
  LabList     = (/"RD","CR","ML(CR)"/)
  LabList(1:) = LabList(1:) + " ~F19~:~N~~F~ " + LabList(0)
  yst         = 1979
  yed         = 2014

  if((yed-yst+1).lt.5) then 
    plot_sig_test  = False 
  end if 

  PltNums     = (/"(a)", "(b)", "(c)", "(d)"/)
  VarList     = (/"UBOT"/)
  VarName     = "UBOT"
  VarUnit     = "m s~S~-1~N~"
  VarStr      = "Zonal wind at bottom model level (UBOT, m s~S~-1~N~)"
  VarSt1      = "Bias in UBOT (m s~S~-1~N~)"
  nvars       = dimsizes(VarList)
  rstr        = "Annual mean" ;VarName+ "(" + VarUnit+")"
  plotname    = "fig_2d_map_clim_bias_"+VarName

  do igp = 1, ngrps-1,1

   model = Groups(igp)
   obs   = Groups(0)
   fl0   = datadir+obs+"_"+VarList+"_"+yst+"-"+yed+".nc"
   f0    = addfile(fl0, "r")
   fl1   = datadir+model+"_"+VarList+"_"+yst+"-"+yed+".nc"
   f1    = addfile(fl1, "r")
 
   uboto = f0->$VarName$
   ubotm = f1->$VarName$
  ;printMinMax(uboto,False)
  ;printMinMax(ubotm,False)

   uboto = where(ismissing(ubotm),uboto@_FillValue,uboto)
   ubotm = where(ismissing(uboto),ubotm@_FillValue,ubotm)

   ubotd = ubotm
   ubotd = ubotm - uboto

  ;print(num(ismissing(uboto)))
  ;print(num(ismissing(ubotm)))
  
   xAve  = f0->$VarName$
   yAve  = f1->$VarName$
   xVar  = f0->variance
   yVar  = f1->variance

   sigr  = 0.05                                   ; critical sig lvl for r
   xN    = (yed - yst +1)*12.0
   yN    = (yed - yst +1)*12.0

   iflag = False                        ; population variance similar
   prob  = xAve
   prob  = (/ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)/)
   prob  = where(prob.lt.sigr,prob,prob@_FillValue )
   delete([/xAve,yAve,xVar,yVar,sigr,xN,yN,iflag/])
 
   ubotm = where(lmsk.eq.0.0,ubotm@_FillValue,ubotm)
   uboto = where(lmsk.eq.0.0,uboto@_FillValue,uboto)
   prob  = where(lmsk.eq.0.0,prob@_FillValue,prob)
   ubotd = where(lmsk.eq.0.0,ubotd@_FillValue,ubotd)
  
   ;calculate the global/regional mean values
   latrg = f1->lat
   rad   = 4.*atan(1.)/180.    ;
   wgty  = cos(latrg*rad)      ; cosine weights
   grmn  = wgt_areaave_Wrap(uboto,   wgty, 1.0, 0)
   gbmn  = wgt_areaave_Wrap(ubotd,   wgty, 1.0, 0)
   grms  = wgt_arearmse(uboto,ubotm, wgty, 1.0, 0)


   rstr1  = "Mean: " +sprintf("%.2f",grmn) + " " + VarUnit;+ VarName+ " (" +VarUnit +")" 
   rstr2  = "Bias/RMSE: "+ sprintf("%.2f",gbmn) + "/" + sprintf("%.2f",grms)
   rstr2  = "RMSE: " + sprintf("%.2f",grms) + " " + VarUnit
   delete([/latrg,rad,wgty,gbmn,grms/])

   guess     = 1                ; use zonal means
   is_cyclic = True             ; cyclic [global]
   nscan     = 1500             ; usually much less than this
   eps       = 1.e-2            ; variable dependent
   relc      = 0.6              ; relaxation coefficient
   opt       = 0                ; not used
   poisson_grid_fill(uboto, is_cyclic, guess, nscan, eps, relc, opt)

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(igp.eq.1)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
   ;wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
    plot     = new(ngrps,"graphic")
    ;;colors for the original plot
    ;omap0      = read_colormap_file("prcp_3")
    ;orgcolors := (/omap0(3,:),omap0(4,:),omap0(5,:),omap0(6,:), omap0(7,:), \
    ;               omap0(10,:),omap0(11,:),omap0(12,:),omap0(13,:), \
    ;               omap0(14,:),omap0(15,:), omap0(16,:), omap0(17,:)/)

    ;colors for the difference plot
    omap       = read_colormap_file("BlAqGrYeOrReVi200")
    orglev    := fspan(-9,12,21) ;(/-4,-2,0,2,4,6,8,10,12,15,20/)
    orgcolors := span_color_rgba(omap,dimsizes(orglev)+1)

    cmap       = read_colormap_file("NCV_blu_red") ;"NCV_blu_red") ;"cmp_b2r")    ; returns 254 x 4 array
    cmap(dimsizes(cmap(:,0))/2-7:dimsizes(cmap(:,0))/2+6,:) = 0
    cnlevels  := (/-4,-3.0,-2.5,-2.0,-1.5,-1.0,-0.3,0.3,1.0,1.5,2.0,2.5,3.0,4.0/)
                
    cncolors   = span_color_rgba(cmap,dimsizes(cnlevels)+1) ;(/72,62,55,5,7,10,0,120,14,16,91,85,77/)

   end if

   ;;;set up the plot information;;;;
   ;load "./var_share_colorbar.ncl"

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 1.0
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.024
   res@tiMainFontHeightF    = FontHeightF*0.95
   res@tmYLLabelFontHeightF = FontHeightF*0.95
   res@tmXBLabelFontHeightF = FontHeightF*0.95

   mpres     = res
   res@gsnStringFontHeightF  = FontHeightF*1.4
   res@lbLabelFontHeightF    = FontHeightF*1.1
   res@lbTitleFontHeightF    = FontHeightF*1.15

   res@cnFillOn              = True     ; color plot desired
  ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = False   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.20
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.12 ;0.60
   res@pmLabelBarWidthF         = 0.65 ;0.06
   res@lbTitleString            = "~F33~D~N~~F~U (850hPa, m/s)"; 
  ;res@lbTitlePosition          = "Right"         ; title location
  ;res@lbTitleDirection         = "Across"        ; letter angle
  ;res@lbTitleAngleF            = 90.          ; title angle

   res@gsnCenterString    = ""
   res@gsnLeftString      = ""
   res@gsnRightString     = "" 
   res@tiMainString       = "" 

;---Control appearance of map.
    mpres@gsnMajorLatSpacing     = 30
    mpres@gsnMajorLonSpacing     = 45
    mpres@gsnMinorLatSpacing     = 15
    mpres@gsnMinorLonSpacing     = 15
   ;mpres@mpProjection           = "Robinson"       ; choose projection
    mpres@mpLabelsOn             = False
    mpres@mpPerimOn              = False
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = "Grey23"  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 1.0          ; thickness of outlines
   ;mpres@mpOutlineOn            = True
   ;mpres@mpOutlineDrawOrder     = "PostDraw"
   ;mpres@mpFillDrawOrder        = "Predraw"
   ;mpres@mpGridAndLimbOn        = True
   ;mpres@mpGridLatSpacingF      = 90.0
   ;mpres@mpGridLonSpacingF      = 180.0
   ;mpres@mpGridLineThicknessF   = 3.0
   ;mpres@mpGridAndLimbDrawOrder="PreDraw"
   ;mpres@mpMinLonF              = 0  ;-180.
   ;mpres@mpMaxLonF              = 360;180.
    mpres@mpCenterLonF           = 180.0
;   mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(0)   = gsn_csm_map(wks1,mpres)
   plot(igp) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                = True
   res1@cnLevelSelectionMode        = "ExplicitLevels"  
   res1@cnLevels                    = orglev
   res1@cnFillColors                = orgcolors
   res1@cnFillPalette               = "NCV_jaisnd" ;"WhBlGrYeRe" ;"BlAqGrYeOrReVi200"
   res1@lbBoxLinesOn                = True
   res1@lbBoxSeparatorLinesOn       = False  
   res1@lbLabelStride               = 5
   res1@gsnLeftString               = PltNums(0) + "  " +LabList(0) + " (" + VarName + ")"
   res1@gsnLeftStringOrthogonalPosF = 0.02
  ;res1@gsnLeftStringParallelPosF   = 0.2
   res1@gsnRightString              = rstr1
   res1@gsnRightStringOrthogonalPosF= 0.02
   res1@cnLinesOn                   = False ;plot_contour  
   res1@cnLineLabelsOn              = False  
   res1@cnInfoLabelOn               = False
   res1@cnLineLabelDensityF         = 2.0
   res1@cnLineLabelFontHeightF      = 0.018
   res1@cnLineLabelBackgroundColor  = -1
   res1@gsnContourLineThicknessesScale = 1.0
   res1@cnLineThicknessF            = 0.001
   res1@cnLineColor                 = "Black"
   res1@gsnContourZeroLineThicknessF= 0.2
   res1@gsnContourNegLineDashPattern= 0
   plot_org1  = gsn_csm_contour(wks1,uboto,res1)
   overlay(plot(0) ,plot_org1)

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@gsnAddCyclic                 = True
    res2@cnLevelSelectionMode         = "ExplicitLevels"  
    res2@lbBoxLinesOn                 = True
    res2@lbBoxSeparatorLinesOn        = False
    res2@lbLabelStride                = 2
    res2@cnFillColors                 = cncolors
    res2@cnLevels                     = cnlevels
    res2@cnLinesOn                    = False    ; turn off contour lines
    res2@cnLineLabelsOn               = False   ; turn off contour labels
    res2@cnInfoLabelOn                = False
    res2@gsnLeftString                = PltNums(igp)+"  " + LabList(igp) + " (" + VarName + ")"
    res2@gsnLeftStringOrthogonalPosF  = 0.02
   ;res2@gsnLeftStringParallelPosF    = 0.2 
    res2@gsnRightString               = rstr2 
    res2@gsnRightStringOrthogonalPosF = 0.02
   ;wrf_smooth_2d(ubotd,3)
    plot_rgd1 = gsn_csm_contour(wks1,ubotd,res2)  
    overlay(plot(igp),plot_rgd1)

    if(plot_contour)
     res21 =    res2
     res21@gsnLeftString      = ""
     res21@gsnRightString     = ""
     res21@gsnCenterString    = ""
     res21@cnFillOn           = False
     res21@cnLinesOn          = True
     res21@cnLineLabelsOn     = False   ; turn off contour labels
     res21@cnInfoLabelOn      = False
     res21@cnLineLabelDensityF      = 2.0
     res21@cnLineLabelFontHeightF   = 0.018
     res21@cnLineLabelBackgroundColor     = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF      = 1.0
     res21@cnLineColor     = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
     plot_cont1      = gsn_csm_contour(wks1,ubotd,res21)
     overlay(plot(igp),plot_cont1)
     delete([/res21/])
    end if

    if(plot_sig_test) then
      res21 = res2
      res21@gsnLeftString                = ""
      res21@gsnRightString               = ""
      res21@gsnCenterString              = ""
      res21@cnMissingValFillColor        = "Transparent"
      res21@cnFillOn                     = True
      res21@cnLinesOn                    = False
      res21@cnLineLabelsOn               = False       ; turn off contour labels
      res21@cnInfoLabelOn                = False
      res21@cnMonoFillColor              = True        ; default color is fground [black]
      res21@cnMonoFillPattern            = True
      res21@cnMonoFillScale              = True
      res21@cnFillPattern                = 17
      res21@cnFillScaleF                 = 0.5
      res21@cnFillDotSizeF               = 0.002
      res21@lbLabelBarOn                 = False            
      plot_cont2 = gsn_csm_contour(wks1,prob,res21)
      overlay(plot(igp),plot_cont2)
      delete([/res21/])
    end if
    delete([/ubotm,ubotd/])
    delete([/res,res2,mpres/])
  end do ; end of experiment group loop

;---Specify y positions of each row of plots.
  nrows       = 1
  plot_height = 1./tofloat(nrows)

  ypos1  = 0.96     ; top row
  ypos2  = 0.64     ; middle row
  ypos   = (/0.97,0.69,0.46,0.23/) + 0.01
  xpos   = (/0.05,0.315,0.57/)     ; top row

  dist_from_title_to_plot = 0.013
  ytitl1 = ypos1 + dist_from_title_to_plot
  ytitl2 = ypos2 + dist_from_title_to_plot

  pres                          = True
; pres@gsnPanelDebug            = True      ; prints debug info
  pres@gsnFrame                 = False
  pres@gsnMaximize              = False
  pres@gsnPanelMainFontHeightF  = 0.015
  pres@gsnPanelFigureStrings    = ""
  pres@gsnPanelFigureStringsFontHeightF  = 0.01

;---Top row
  pres1 = pres
  pres1@gsnPanelLabelBar         = True
  pres1@pmLabelBarParallelPosF   = 0.0
  pres1@pmLabelBarOrthogonalPosF = -0.01
  pres1@pmLabelBarHeightF        = 0.04 ;0.60
  pres1@pmLabelBarWidthF         = 0.32 ;0.06
  pres1@lbBoxLineThicknessF      = 0.1
  pres1@lbBoxLinesOn             = True
  pres1@lbBoxSeparatorLinesOn    = False
  pres1@lbLabelStride            = 2
  pres1@lbLabelFontHeightF       = FontHeightF*0.35
  pres1@lbTitleFontHeightF       = FontHeightF*0.35
  pres1@lbTitleString            = VarStr
  pres1@gsnPanelTop              = 1.0
  pres1@gsnPanelBottom           = 0.715
  pres1@gsnPanelRight            = 0.5
  pres1@gsnPanelYF              := ypos(0)
 ;pres1@gsnPanelXF              := xpos(0)
  pres1@gsnPanelMainString       = ""
  pres1@gsnPanelMainPosYF        = ytitl1
  gsn_panel(wks1,plot(0),(/1,1/),pres1)

;---Bottom row
  pres2 = pres
  pres2@gsnPanelLabelBar         = True
  pres2@pmLabelBarParallelPosF   = 0.00
  pres2@pmLabelBarOrthogonalPosF = -0.016
  pres2@pmLabelBarHeightF        = 0.04 ;0.60
  pres2@pmLabelBarWidthF         = 0.33 ;0.06
  pres2@lbBoxLineThicknessF      = 0.1
  pres2@lbBoxLinesOn             = True
  pres2@lbBoxSeparatorLinesOn    = False
  pres2@lbLabelStride            = 1
  pres2@lbLabelFontHeightF       = FontHeightF*0.35
  pres2@lbTitleFontHeightF       = FontHeightF*0.35
  pres2@lbTitleString            = VarSt1
  pres2@gsnPanelTop              = 0.75
  pres2@gsnPanelBottom           = 0.0
  pres2@gsnPanelYF              := ypos(1:)
 ;pres2@gsnPanelXF              := xpos(1:)
  pres2@gsnPanelRight            = 0.5
 ;pres2@gsnPanelLeft             = 0.5                
  pres2@gsnPanelMainString       = ""
  pres2@gsnPanelMainPosYF        = ytitl2    ; Y position of title
  gsn_panel(wks1,plot(1:),(/3,1/),pres2)

  frame(wks1)


end
