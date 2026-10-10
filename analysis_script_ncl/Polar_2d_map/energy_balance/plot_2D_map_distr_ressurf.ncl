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

  DATA_DIR    = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  ERA5_File   = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  E3SM_File   = DATA_DIR + "20210727.v2.LR.F20TR_eam.h0_monthly_1981-2010.nc"
  DiagTime  = "1981-2010"

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

  
 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups      = (/"ERA5","v2.LR.F20TR"/)
  ngrps       = dimsizes(Groups)
  LabList     = new(4,string)
  LabList(0)  = Groups(0) + " reanalysis"
  LabList(1)  = "v2.LR.F20TR"
  LabList(2)  = "v2.LR.F20TR - ERA5"

 ;;; flags to control the plot style
  lndcolor = "Grey28"

  f1     = addfile(ERA5_File,"r")
  f2     = addfile(E3SM_File,"r")
  
  time   = f2->time
  lat    = f2->lat
  lon    = f2->lon

  ;;;energy balance at surface;;;;
  ;;;RESSURF  = FSNS - FLNS - SHFLX - LHFLX;;;
 
  ;;; ERA5 ;; 
  er5_ressurf = f2->ssr
  er5_ressurf = f2->ssr + f2->str - f2->slhf - f2->sshf

  ;;;E3SM ;;;
  eam_ressurf = f2->FSNS
  eam_ressurf = FSNS
 
  var1   = f1->z(0,:,:) 
  var2   = f2->PHIS(0,:,:)
  var2   = f2->PHIS(0,:,:) - f1->z(0,:,:)
  var1   = var1 / 9.8065
  var2   = var2 / 9.8065

  Rstring1     = "PHIS (m)"
  Rstring2     = "~F33~D~N~~F~PHIS (m)"

  ;;;;calculate the global mean and correlation 
  latrg   = lat({MinLat:MaxLat})
  rad     = 4.*atan(1.)/180.    ;
  cost    = cos(latrg*rad)       ; cosine weights
  mvar1   = wgt_areaave_Wrap (var1({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar2   = wgt_areaave_Wrap (var2({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)

  mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
  mstr2   = "Mean: " + sprintf("%5.2f", mvar2)
  
  ;----------------------------------------------------------------------
  ; Plotting section

   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_e3sm_vs_era5_PHIS")
  plot        = new(4,graphic)

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
    mpres@mpCenterLonF            = 0
    mpres@mpLabelsOn              = False
    mpres@mpPerimOn               = True
    mpres@mpFillOn                = False
    mpres@mpGeophysicalLineColor  = lndcolor  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
    mpres@mpOutlineOn             = True 
    mpres@mpOutlineDrawOrder      = "PostDraw"
    mpres@mpFillDrawOrder         = "Predraw"
    mpres@mpGridAndLimbOn         = False 
    mpres@mpGridLineThicknessF    = 1.0
;   mpres@mpOceanFillColor        = "lightskyblue1"
;   mpres@mpLandFillColor         = "gray33"
    mpres@tmYLLabelsOn            = True
    mpres@tmYLOn                  = True
    do igp = 0,3 
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
    res1@gsnAddCyclic                   = True
    res1@cnFillPalette                  = "WhBlGrYeRe"
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                      := (/10,100,200,400,800,1200,1600,2000,2400,2800,3200,3600,4000,5000/)
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = ""
    res1@tiMainOffsetYF                 = 0.01
    res1@cnLinesOn                      = True 
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
    res1@gsnRightString                 = Rstring1
    plot_org1  = gsn_csm_contour(wks1,var1,res1)
    overlay(plot(0),plot_org1)

    delete([/res1/])

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                = True
    res2@cnFillPalette               = "CBR_coldhot"
    res2@cnLevelSelectionMode        = "ExplicitLevels"
    res2@cnLevels                   := (/-800,-500,-300,-100,-10,10,100,300,500,800/)
  ; res2@cnFillColors               := colorgroup1
    res2@lbBoxLinesOn                = True
   ;res2@lbBoxSeparatorLinesOn       = False
   ;res2@lbLabelStride               = 15
    res2@cnLinesOn                   = False             ; turn off contour lines
    res2@cnLineLabelsOn              = False            ; turn off contour labels
    res2@cnInfoLabelOn               = False
    res2@gsnLeftString               = LabList(1)
    res2@gsnRightString              = Rstring2
    plot_rgd2                        = gsn_csm_contour(wks1,var2,res2)
    overlay(plot(1),plot_rgd2)
      
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
    gsn_panel(wks1,plot,(/1,2/),pnres)

   delete([/pnres/])

end

