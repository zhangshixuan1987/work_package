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

  E3SM_dir  = "/compyfs/zhan391/cryosphere_scratch/post_process/data/180x360/" 
  ERA5_dir  = "/compyfs/zhan391/cryosphere_scratch/post_process/data/era5_anl/"
  CaseName  = "e3sm_v2_compy_F20TR_ne30pg2_r05_oECv3_e3smv2_99a485645"
  DiagTime  = "2010-01"

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

  
 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups    = (/"ERA5","EAMv2"/)
  ngrps     = dimsizes(Groups)

  LabList     = new(4,string)
  LabList(0)  = Groups(0) + " reanalysis"
  LabList(1)  = "E3SM - ERA5"
  LabList(2)  = Groups(0) + " reanalysis"
  LabList(3)  = "E3SM - ERA5"

 ;;; flags to control the plot style
  lndcolor = "Grey28"

  file_era5  = ERA5_dir + "era5_siglev_"+DiagTime+".nc"
  file_e3sm  = E3SM_dir + CaseName+".eam.h0."+DiagTime+".nc"

  f1         = addfile(file_era5,"r")
  f2         = addfile(file_e3sm,"r")

  
  time   = f2->time
  lat    = f2->lat
  lon    = f2->lon

  var1   = f1->skt(0,:,:) 
  var2   = f2->TS(0,:,:)
  var3   = f1->skt(0,:,:)
  var4   = f2->TBOT(0,:,:)
  
  var2   = f2->TS(0,:,:) - f1->skt(0,:,:)
  var4   = f2->TREFHT(0,:,:) - f1->t2m(0,:,:)
  Rstring1     = "TS  (K)"
  Rstring2     = "T2m (K)"

   ;;;;calculate the global mean and correlation 
   latrg   = lat({MinLat:MaxLat})
   rad     = 4.*atan(1.)/180.    ;
   cost    = cos(latrg*rad)       ; cosine weights
   mvar1   = wgt_areaave_Wrap (var1({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
   mvar2   = wgt_areaave_Wrap (var2({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
   mvar3   = wgt_areaave_Wrap (var3({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
   mvar4   = wgt_areaave_Wrap (var4({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)

   mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
   mstr2   = "Mean: " + sprintf("%5.2f", mvar2)
   mstr3   = "Mean: " + sprintf("%5.2f", mvar3)
   mstr4   = "Mean: " + sprintf("%5.2f", mvar4)
  
  ;----------------------------------------------------------------------
  ; Plotting section

   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_e3sm_vs_era5_TS")
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
    mpres@mpOutlineOn              = True 
    mpres@mpOutlineDrawOrder       = "PostDraw"
    mpres@mpFillDrawOrder          = "Predraw"
    mpres@mpGridAndLimbOn          = False 
    mpres@mpGridLineThicknessF     = 1.0
;   mpres@mpOceanFillColor         = "lightskyblue1"
;   mpres@mpLandFillColor          = "gray33"
    mpres@tmYLLabelsOn             = True
    mpres@tmYLOn                   = True
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
    res1@gsnAddCyclic                  = True
    res1@cnFillPalette                 = "WhBlGrYeRe"
    res1@cnLevelSelectionMode          = "ExplicitLevels"
    res1@cnLevels                     := (/210,220,230,240,250,260,270,275,280,285,290,295,300,305,310/)
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

    res1@gsnLeftString                  = LabList(2)
    res1@gsnRightString                 = Rstring2
    plot_org2  = gsn_csm_contour(wks1,var3,res1)
    overlay(plot(2),plot_org2)

    delete([/res1/])

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = "CBR_coldhot"
    res2@cnLevelSelectionMode           = "ExplicitLevels"
    res2@cnLevels                      := (/-20,-10,-5,-2,-1,1,2,5,10,20/)
  ; res2@cnFillColors                  := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False

    res2@gsnLeftString                  = LabList(1)
    res2@gsnRightString                 = "~F33~D~N~~F~"+ Rstring1
    plot_rgd1                           = gsn_csm_contour(wks1,var2,res2)  
    overlay(plot(1),plot_rgd1)

    res2@gsnLeftString                  = LabList(3)
    res2@gsnRightString                 = "~F33~D~N~~F~"+ Rstring2
    plot_rgd4                           = gsn_csm_contour(wks1,var4,res2)
    overlay(plot(3),plot_rgd4)
      
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
    gsn_panel(wks1,plot,(/2,2/),pnres)

   delete([/pnres/])

end

