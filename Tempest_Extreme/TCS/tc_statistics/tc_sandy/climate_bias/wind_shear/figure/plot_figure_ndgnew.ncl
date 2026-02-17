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

  l_create_plot_data = False
  l_create_plot_data = True

  l_plot_cntl       = True  ; if True, plot the CNTL origianl field on the leftmost panel
  l_plot_absdif     = True  ; if True, plot biases/differences
  l_plot_reldif     = False ; if True, plot relative differences
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region

  outdir    = "../../data/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Obsnam    = (/"CERES_EBAF","ERA5","CERES_EBAF","CERES_EBAF"/)
  Casdir    = (/"Obs.","CLIM","UVTQ3_BASE", "UVTQ3_GFDL","UVTQ3_NPBL_NUPP","UVTQPS3_NPBL_NUPP"/)
  Groups    = Casdir
  ngrps     = dimsizes(Groups)

  LabList     = Groups
  LabList(1:) = Groups(1:) + " ~F19~:~N~~F~ " + Groups(0)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList   = (/"CLDTOT", "PRECT", "SWCF", "FLUT"/) ;, "LWCF", "NETCF"/)
  VstList   = (/"CLDTOT", "PRECT", "SWCF", "FLUT"/) ;, "LWCF", "NETCF"/)
  UnitList  = (/"fraction","mm day~S~-1~N~",\
                "W m~S~-2~N~", "W m~S~-2~N~", "W m~S~-2~N~"/)
  nvars     = dimsizes(VarList)
  PltNums   = (/"a", "b", "c", "d", "e", "f", "g", "h"/)

 ;;; flags to control the plot style
  l_plot_cntl    = True   ; if True, plot the CNTL origianl field on the leftmost panel
  plot_contour   = False   ; if True, plot the contour in the figure
  plot_leftstr   = True    ; if True, the lefadd igure captions as left string

  if(l_plot_cntl) then 
    plot        = new (ngrps*nvars, graphic)
  else
    plot        = new ((ngrps-1)*nvars, graphic)
  end if 

 do igp = 1, ngrps-1,1

  do iv = 0,nvars-1, 1

   Varname = VarList(iv)
   Varunt  = UnitList(iv)
   Varstr  = VstList(iv)

   print("iv : " + iv)

   fna = outdir+ Varname + "_" + Casdir(igp) + ".nc"
   fla = addfile(fna, "r")

   print(" ")
   print("........................... ")
   print("igp : " + igp) 
   print("........................... ")

   vna = Varname + "_" + Obsnam(iv)
   vnb = Varname + "_" + Casdir(igp)
   vnc = "d" + Varname

   pncre0300=fla->$vna$
   wncredf1=fla->$vnb$
   wncredf2=fla->$vnc$

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(igp.eq.1.and.iv.eq.0)then

    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "fig_2d_map_ndgnew"
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")

   end if 

   ;;;set up the plot information;;;;
   load "./var_share_colorbar.ncl"

   Figstr1 = "" ;PltNums(0)
   Figstr2 = "" ;PltNums(igp)
   Figstr3 = "" ;PltNums(igp)


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
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF

   mpres     = res
   resx      = res
   resx1     = res

   res@gsnStringFontHeightF  = FontHeightF*1.5
   res@lbLabelFontHeightF    = FontHeightF*0.6
   res@lbTitleFontHeightF    = FontHeightF*0.8

   res@cnFillOn              = True     ; color plot desired
  ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = True   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.20
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.12 ;0.60
   res@pmLabelBarWidthF         = 0.65 ;0.06
  ;res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
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
    mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
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

   plot(iv) = gsn_csm_map(wks1,mpres)
   plot(igp*nvars+iv)     = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                = True
   res1@cnLevelSelectionMode        = "ExplicitLevels"  
   res1@cnLevels                   := orglev
   res1@cnFillColors               := orgcolors
   res1@lbBoxLinesOn                = True
  ;res1@lbBoxSeparatorLinesOn       = False  
  ;res1@lbLabelStride               = 15
   res1@gsnLeftString               = Figstr1 + "  " +LabList(0)+ " ("+Obsnam(iv)+ ")"
   res1@gsnLeftStringOrthogonalPosF = 0.04
  ;res1@gsnLeftStringParallelPosF  = 0.2
   res1@gsnRightString              = Varstr+ " ("+Varunt+")";
   res1@gsnRightStringOrthogonalPosF= 0.055
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
  ;res1@sfXArray     = lon
  ;res1@sfYArray     = lat
   plot_org1  = gsn_csm_contour(wks1,pncre0300,res1)
   overlay(plot(iv) ,plot_org1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    res2=res
    res2@gsnAddCyclic            = True
    res2@cnLevelSelectionMode     = "ExplicitLevels"  
    res2@cnFillColors             = cncolors
    res2@cnLevels                 = cnlevels
    res2@cnLinesOn                = False    ; turn off contour lines
    res2@cnLineLabelsOn           = False   ; turn off contour labels
    res2@cnInfoLabelOn            = False
    res2@gsnLeftString            = Figstr2+"  " + LabList(igp)
    res2@gsnLeftStringOrthogonalPosF = 0.04
   ;res2@gsnLeftStringParallelPosF   = 0.2 
    res2@gsnRightString     = "~F33~D~N~~F~"+ Varstr + "  ("+Varunt+")"
    res2@gsnRightStringOrthogonalPosF   = 0.055
    plot_rgd1      = gsn_csm_contour(wks1,wncredf1,res2)  

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
     plot_cont1      = gsn_csm_contour(wks1,wncredf1,res21)
     overlay(plot_rgd1,plot_cont1)
     delete([/res21/])
    end if

    if(l_plot_absdif) then
      overlay(plot(igp*nvars+iv),plot_rgd1)
    end if

    ;;;;;add a box for selected region;;;;
    
    if(l_plot_latlon_box.and.Groups(igp).eq."RNDG_UVT6") then
       regx = new((/4/),graphic)
       resb                  = True ; polyline mods desired
       resb@gsLineColor      = "Magenta" ; color of lines
       resb@gsLineThicknessF = 2.5 ; thickness of lines
       do i = 0,3
        regx(i)=gsn_add_polyline(wks1,plot_rgd1,xptsa(i:i+1),yptsa(i:i+1),resb)
       end do
    end if

   ;;;;;;;;;;;;;;;;;;;;;;;;
   ;relative difference 
   ;;;;;;;;;;;;;;;;;;;;;;;;
    res2@cnFillColors       = cncolors
    res2@cnLevels           = cnlevels
    res2@tiMainString       = "" ;
    res2@gsnLeftString      = Figstr3+"  " +LabList(igp) 
    res2@gsnLeftStringOrthogonalPosF = 0.04
   ;res2@gsnLeftStringParallelPosF   = 0.2
    res2@gsnRightString     = "~F33~D~N~~F~"+ Varstr + ", Relative  (%)"
    res2@gsnRightStringOrthogonalPosF = 0.055
    plot_rgd2      = gsn_csm_contour(wks1,wncredf2,res2)

    if(plot_contour)
     res21 =    res2
     res21@gsnLeftString      = ""
     res21@gsnRightString     = ""
     res21@gsnCenterString    = ""
     res21@cnFillOn        = False
     res21@cnLinesOn       = True
     res21@cnLineLabelsOn     = False   ; turn off contour labels
     res21@cnInfoLabelOn      = False
     res21@cnLineLabelDensityF   = 2.0
     res21@cnLineLabelFontHeightF   = 0.018
     res21@cnLineLabelBackgroundColor  = -1
     res21@gsnContourLineThicknessesScale = 1.0
     res21@cnLineThicknessF      = 1.0
     res21@cnLineColor     = "Black"
     res21@gsnContourZeroLineThicknessF   = 2.
     res21@gsnContourNegLineDashPattern   = 1
    ;res21@sfXArray     = lon
    ;res21@sfYArray     = lat
     plot_cont1      = gsn_csm_contour(wks1,wncredf2,res21)
     overlay(plot_rgd2,plot_cont1)
     delete([/res21/])
    end if
 
    if(l_plot_reldif) then 
      overlay(plot(igp*nvars+iv),plot_rgd2)
    end if 

 delete([/res,res2,resx,mpres/])

end do 

end do ; end of experiment group loop

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
  pres1                          = True
  pres1@gsnFrame                 = True      ; don't advance frame yet
  pres1@gsnMaximize              = False
  pres1@gsnPanelLabelBar         = False
  pres1@pmLabelBarParallelPosF   = -0.010
  pres1@pmLabelBarOrthogonalPosF = 0.0
  pres1@lbLabelFontHeightF       = FontHeightF * 0.45
  pres1@lbOrientation            = "vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.21
  pres1@pmLabelBarWidthF         = 0.05
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.028

  if (.not.l_plot_cntl) then 
    gsn_panel(wks1,(/plot(nvars:)/),(/ngrps-1,nvars/),pres1)
  else
    gsn_panel(wks1,(/plot(:)/),(/ngrps,nvars/),pres1)
  end if 

end
