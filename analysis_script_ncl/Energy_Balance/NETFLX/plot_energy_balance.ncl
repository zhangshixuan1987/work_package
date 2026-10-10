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

  DATA_DIR    = "./"
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

  ;;;;;;;date strings for the data directory;;;;;;;;
  Groups      = (/"ERA5","EAMv2"/)
  ngrps       = dimsizes(Groups)

  npanel      = 3 
  LabList     = new(npanel,string)
  LabList(0)  = "ERA5 (1981-2010)"
  LabList(1)  = "EAMv2 (1981-2010)"
  LabList(2)  = "EAMv2 - ERA5" 
  Rstring     = "Surface enerngy fluxes (W m~S~-2~N~)"

  ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  fnam1  = DATA_DIR + "Surface_energy_fluxes_" + Groups(0)+".nc" 
  fnam2  = DATA_DIR + "Surface_energy_fluxes_" + Groups(1)+".nc" 
  f1     = addfile(fnam1,"r")
  f2     = addfile(fnam2,"r")
  
  lat    = f1->lat({MinLat:MaxLat})
  lon    = f1->lon({MinLon:MaxLon})
  var1   = dim_avg_n_Wrap(f1->NETFLX1(:,{MinLat:MaxLat},{MinLon:MaxLon}),0)
  var2   = dim_avg_n_Wrap(f2->NETFLX2(:,{MinLat:MaxLat},{MinLon:MaxLon}),0)
  vard   = var1
  vard   = var2 - var1

  ;;;;calculate the global mean and correlation
  rad     = 4.*atan(1.)/180.    ;
  cost    = cos(lat*rad)       ; cosine weights
  mvar1   = wgt_areaave_Wrap (var1,cost,1.,0)
  mvar2   = wgt_areaave_Wrap (var2,cost,1.,0)
  mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
  mstr2   = "Mean: " + sprintf("%5.2f", mvar2)

  ;----------------------------------------------------------------------
  ; Plotting section
  ;----------------------------------------------------------------------
   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_e3sm_vs_era5_PSL")
   plot        = new(npanel,graphic)

   cnlevs = fspan(-200,200,41)
   cndiff = fspan(-75,75,31)
    
  ;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@tmXTOn               = False
    res@tmYROn               = False

    FontHeightF              = 0.024
    res@tiMainFontHeightF    = FontHeightF*1.2
    res@tmYLLabelFontHeightF = FontHeightF
    res@tmXBLabelFontHeightF = FontHeightF
    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF*0.6
    res@gsnStringFontHeightF = FontHeightF*1.5

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = True      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61
    res@pmLabelBarOrthogonalPosF = 0.05
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.65
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = ""
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height

    ;---Control appearance of map.
    mpres@gsnMajorLatSpacing      = 30
    mpres@gsnMajorLonSpacing      = 45
    mpres@gsnMinorLatSpacing      = 15
    mpres@gsnMinorLonSpacing      = 15
    mpres@mpFillOn                = False        ; turn off map fill
    mpres@mpCenterLonF            = 180
    mpres@mpLabelsOn              = False
    mpres@mpPerimOn               = True
    mpres@mpFillOn                = False
    mpres@mpGeophysicalLineColor  = "Black";lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
    mpres@mpOutlineOn              = True
    mpres@mpOutlineDrawOrder       = "PostDraw"
    mpres@mpFillDrawOrder          = "Predraw"
    mpres@mpGridAndLimbOn          = False
    mpres@mpGridLineThicknessF     = 1.0
;   mpres@mpOceanFillColor         = "lightskyblue1"
;   mpres@mpLandFillColor          = "gray33"
    mpres@tmYLLabelsOn             = True
    mpres@tmYLOn                   = True
    do igp = 0,npanel-1
      plot(igp) = gsn_csm_map(wks1,mpres)
    end do

    res@gsnLeftString                  = ""
    res@gsnLeftStringOrthogonalPosF    = 0.02
   ;res@gsnLeftStringParallelPosF      = 0.05
    res@gsnRightString                 = ""
    res@gsnRightStringOrthogonalPosF   = 0.02
    res@gsnRightStringParallelPosF     = 0.98
    res@gsnCenterString                = ""
    res@gsnCenterStringOrthogonalPosF  = 0.04

;---Resources for plotting original data
    res1=res
    res1@gsnAddCyclic                  = True
    res1@cnFillPalette                 = "NCV_jaisnd"
    res1@cnLevelSelectionMode          = "ExplicitLevels"
    res1@cnLevels                     := cnlevs
    res1@lbBoxLinesOn                   = True
    res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = ""
    res1@tiMainOffsetYF                 = 0.01
    res1@gsnRightString                 = Rstring
    res1@cnLinesOn                      = False
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

    res1@gsnLeftString                  = LabList(0)
    plot_org1  = gsn_csm_contour(wks1,var1,res1)
    overlay(plot(0),plot_org1)

    res1@gsnLeftString                  = LabList(1)
    plot_org2  = gsn_csm_contour(wks1,var2,res1)
    overlay(plot(1),plot_org2)
    delete([/res1/])

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = "MPL_PRGn"
    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                      := cndiff
  ; res2@cnFillColors                  := colorgroup1
    res2@lbBoxLinesOn                   = True
    res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@gsnRightString                 = "~F33~D~N~~F~"+ Rstring
    res2@gsnCenterString                = "" ;"Global mean: "+ mdfstr
    res2@gsnCenterStringOrthogonalPosF  = 0.04

    res2@gsnLeftString                  = LabList(2)
    plot_rgd1                           = gsn_csm_contour(wks1,vard,res2)
    overlay(plot(2),plot_rgd1)

;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = True
    pnres@gsnPanelMainFontHeightF    = 0.016
    pnres@gsnPanelMainString         = ""
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4
    gsn_panel(wks1,plot,(/3,1/),pnres)

   delete([/pnres/])


end 
