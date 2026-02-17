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
  l_plot_latlon_box = False ; if True, add lat-lon box for SEP region
  datdir = "/global/cscratch1/sd/zhan391/DARPA_project/Track_model/post"

  ;********************************
  ; select sub-regions
  ;********************************
   minlat =  0.0 ;0.0    ; deg N
   maxlat =  70.0 ;57.0    ; deg N
   minlon = 200.0 ;-100.0    ; - deg E = deg W
   maxlon = 350.0 ;-20.0    ; - deg E = deg W
   cenlon = 200.0 ;

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups    = (/"CLIM","NDGUV"/) ;,"NDGUVT","NDGUVTQ"/)
  LabList   = (/"E3SM (CLIM)","E3SM (NDGUV)"/) ;, "E3SM (NDGUVTQ)"/)
  ngrps     = dimsizes(Groups)
  plot      = new (ngrps, graphic)
  ystr      = 2007
  yend      = 2017
  timtag    = ystr+"-"+yend 

  do i = 0,ngrps-1,1

   fna = datdir+"/"+Groups(i)+"_Track_dens_"+ timtag+ ".nc"
   fla = addfile(fna, "r")

   print(" ")
   print("........................... ")
   print("igp : " + i) 
   print("........................... ")

   lat   = fla->lat
   lon   = fla->lon
   var   = dim_avg_n_Wrap(fla->density,0)
   nyear = yend - ystr + 1
   printVarSummary(var)
   printMinMax(var,False)

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "track_density_sythentic" 
    wks1 = gsn_open_wks(wtype,plotname)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")


   end if 

   ;;;set up the plot information;;;;
   ;sub region for sythentic analysis
   lats  =  0.         ;21. -2  ;-30
   late  = 50.         ;36. -2  ;-10
   lons  = -100 + 360. ;-137. -8.0 + 360
   lone  = -40 + 360.  ;-122. -8.0 + 360

   ;;draw box for the region
   yptsa = (/late, lats, lats, late, late/) ;define lat points for each line
   xptsa = (/lons, lons, lone, lone, lons/) ; define lon points for each line

   ;;color for the map boundary line
   lndcolor = "Black" ;"Grey50"

  ;colorbar setup for difference plots
   cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
   cnlevels  := ispan(0,2000,100) ;(/0.1,1,2,3,4,5,6,7,8,9,10/)
   cnlevels(0) = 1

   Figstr = "" 

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 0.5
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
   res@cnFillPalette         = "precip3_16lev"

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
    mpres@mpLimitMode  = "LatLon"
    mpres@mpMinLatF    = minlat
    mpres@mpMaxLatF    = maxlat
    mpres@mpMinLonF    = minlon
    mpres@mpMaxLonF    = maxlon
    mpres@mpCenterLonF = cenlon
 
    mpres@gsnMajorLatSpacing     = 30
    mpres@gsnMajorLonSpacing     = 45
    mpres@gsnMinorLatSpacing     = 15
    mpres@gsnMinorLonSpacing     = 15
   ;mpres@mpProjection           = "Robinson"       ; choose projection
    mpres@mpLabelsOn             = False
    mpres@mpPerimOn              = False
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 0.1          ; thickness of outlines
   ;mpres@mpOutlineOn            = True
   ;mpres@mpOutlineDrawOrder     = "PostDraw"
   ;mpres@mpFillDrawOrder        = "Predraw"
   ;mpres@mpGridAndLimbOn        = True
   ;mpres@mpGridLatSpacingF      = 90.0
   ;mpres@mpGridLonSpacingF      = 180.0
   ;mpres@mpGridLineThicknessF   = 1.0
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

   plot(i) = gsn_csm_map(wks1,mpres)

;---Resources for plotting original data
   res1=res
   res1@gsnAddCyclic                = False
   res1@cnLevelSelectionMode        = "ExplicitLevels"  
   res1@cnLevels                   := cnlevels
  ;res1@cnFillColors               := cncolors
   res1@lbBoxLinesOn                = True
  ;res1@lbBoxSeparatorLinesOn       = False  
  ;res1@lbLabelStride               = 15
   res1@gsnLeftString               = Figstr + "  " +LabList(i)
   res1@gsnLeftStringOrthogonalPosF = 0.04
  ;res1@gsnLeftStringParallelPosF   = 0.2
   res1@gsnRightString              = "TC track counts" 
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
   plot_org1  = gsn_csm_contour(wks1,var,res1)
   overlay(plot(i) ,plot_org1)

   ;;;;;add a box for selected region;;;;
   
   if(l_plot_latlon_box) then 
      regx = new((/4/),graphic)
      resb                  = True ; polyline mods desired
      resb@gsLineColor      = "Magenta" ; color of lines
      resb@gsLineThicknessF = 2.5 ; thickness of lines
      do i = 0,3
       regx(i)=gsn_add_polyline(wks1,plot_rgd1,xptsa(i:i+1),yptsa(i:i+1),resb)
      end do
   end if

 delete([/res,res1,mpres/])

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
  gsn_panel(wks1,(/plot/),(/2,2/),pres1)

end
