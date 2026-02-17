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

 ;;; flags to control the plot style
  l_plot_latlon_box = False
  datdir      = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ARs/post_process/"
  outdir      = "./"
  wgtFileDir  = "./"
  wgtFileName = "map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5"/)
  Groups  = (/"ERA5"/)
  Labels  = (/"ERA5"/)
  ngrps   = dimsizes(Groups)
  year    = 2017
  timplt  = 2017

  varList = (/"VIWVN_PW_AR","VIWVN_PW_NONAR","VIWVE_PW_AR","VIWVE_PW_NONAR"/)
  nvars   = dimsizes(varList)

  do i = 0,ngrps-1,1

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "png"
    wtype@wkWidth  = 2500
    wtype@wkHeight = 2500 
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "global_map" 
    wks1 = gsn_open_wks(wtype,plotname)
    plot = new (ngrps, graphic)
    dum1 = new (ngrps, graphic)
    dum2 = new (ngrps, graphic)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
   end if 

   ;;;set up the plot information;;;;
   load "./var_share_colorbar.ncl"

   Figstr = "" 

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 2.0
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
   res@cnFillPalette         = cmap ;"precip3_16lev" ;"WhBlGrYeRe" ;"precip_11lev" ;"precip3_16lev"

   res@cnLinesOn             = False    ; turn off contour lines
   res@cnLineLabelsOn        = False    ; turn off contour labels
  ;res@cnFillMode            = "RasterFill"   ; turn raster on
   res@cnMissingValFillColor = "White"

   res@lbLabelBarOn             = True   ; Will turn on in panel later
   res@pmLabelBarParallelPosF   = 0.50 ;61
   res@pmLabelBarOrthogonalPosF = 0.20
  ;res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   res@pmLabelBarHeightF        = 0.12 ;0.60
   res@pmLabelBarWidthF         = 0.62 ;0.06
  ;res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
  ;res@lbTitlePosition          = "Right"         ; title location
  ;res@lbTitleDirection         = "Across"        ; letter angle
  ;res@lbTitleAngleF            = 90.          ; title angle

   res@gsnCenterString    = ""
   res@gsnLeftString      = ""
   res@gsnRightString     = "" 
   res@tiMainString       = "" 

;---Control appearance of map.
    ;mpres@gsnMajorLatSpacing     = 30
    ;mpres@gsnMajorLonSpacing     = 45
    ;mpres@gsnMinorLatSpacing     = 15
    ;mpres@gsnMinorLonSpacing     = 15
    mpres@mpProjection            = "Robinson"       ; choose projection
    mpres@mpLabelsOn              = False
    mpres@mpPerimOn               = False
    mpres@mpFillOn                = True   ;-- turn map fill on
    mpres@mpFillDrawOrder         = "PreDraw"   ;-- draw filled map first
    mpres@mpOceanFillColor        = (/ 0.824, 0.961, 1.0 /)   ;-- ocean color
    mpres@mpInlandWaterFillColor  = (/ 0.824, 0.961, 1.0 /)   ;-- inland water color
    mpres@mpLandFillColor         = (/ 0.7, 0.7, 0.7 /)   ;-- land color
    mpres@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 2.0 ;0.1          ; thickness of outlines
    mpres@mpOutlineOn            = True
    mpres@mpOutlineDrawOrder     = "PostDraw"
    mpres@mpFillDrawOrder        = "Predraw"
    mpres@mpGridAndLimbOn        = True
   ;mpres@mpGridLatSpacingF      = 180.0
   ;mpres@mpGridLonSpacingF      = 720.0
    mpres@mpGridLineThicknessF   = 1.0
    mpres@mpGridAndLimbDrawOrder = "PostDraw"
    mpres@mpLimbLineThicknessF   = 2.0
    mpres@mpGridLineColor        = "Transparent"

   ;mpres@mpMinLonF              = 0  ;-180.
   ;mpres@mpMaxLonF              = 360;180.
    mpres@mpCenterLonF           = 240.0
;   mpres@mpCenterLatF           = -90.     ; This is necessary to get the correct map
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(i) = gsn_csm_map(wks1,mpres)

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
  gsn_panel(wks1,(/plot/),(/1,1/),pres1)

end
