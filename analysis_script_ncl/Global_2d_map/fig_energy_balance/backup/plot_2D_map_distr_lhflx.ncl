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

  DATA_DIR  = "/global/cfs/cdirs/e3sm/zhan391/data/polar_diag/"
  ERA5_File = DATA_DIR + "ERA5_analysis_360x180_monthly_1981-2010.nc"
  E3SM_File = DATA_DIR + "20210727.v2.LR.F20TR_eam.h0_monthly_1981-2010.nc"
  DiagTime  = "1981-2010"

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.

  Season    = (/"ANN","DJF","JJA"/)
  nseason   = dimsizes(Season)
  
 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups    = (/"ERA5_analysis_360x180","20210727.v2.LR.F20TR"/)
  ngrps     = dimsizes(Groups)
  LabList   = (/"ERA5 reanalysis", "v2.LR.F20TR(land model) - ERA5", "New formula - Land model"/)
  nlabs     = dimsizes(LabList)

  Rstring1  = "LHFLX (W m~S~-2~N~)"
  Rstring2  = "~F33~D~N~~F~LHFLX (W m~S~-2~N~)"

  ;;; flags to control the plot style
  lndcolor = "Grey28"

 do igp = 1,ngrps-1,1
 
  file1  = DATA_DIR + Groups(0)+"_monthly_"+DiagTime+".nc"
  file2  = DATA_DIR + Groups(igp)+"_eam.h0_monthly_"+DiagTime+".nc"
  
  f1     = addfile(file1,"r")
  f2     = addfile(file2,"r")
  
  time   = f2->time
  lat    = f2->lat
  lon    = f2->lon

  ;;;energy balance at surface;;;;
  ;;;RESSURF  = FSNS - FLNS - SHFLX - LHFLX;;;
 
  ;;; ERA5 ;; 
  lhflx_er5  = f1->slhf
  lhflx_er5  = -lhflx_er5 / 86400.0

  lhflx_eam1 = f2->LHFLX
  
  ;;;e3sm lhflx recalucated with the QFLX 
  lvap = 2.501e6
  lice = 3.337e5
  lsub = lvap + lice

  lhflx_eam2 = f2->QFLX
  tmp1 = f2->QFLX
  tmp2 = f2->PRECC + f2->PRECL - f2->PRECSC- f2->PRECSL
  lhflx_eam2 = lsub*tmp1 - lice*1e3*tmp2
 
  ;;;prepare the data for plot 

  do isea = 0,nseason-1,1

  if (Season(isea).eq."ANN") then 
    vtp1  = month_to_annual(lhflx_er5,1)
    vtp2  = month_to_annual(lhflx_eam1,1)
    vtp3  = month_to_annual(lhflx_eam2,1)
  else
    vtp1  = month_to_season(lhflx_er5,Season(isea))
    vtp2  = month_to_season(lhflx_eam1,Season(isea))
    vtp3  = month_to_season(lhflx_eam2,Season(isea))
  end if 

  var1  = dim_avg_n_Wrap(vtp1,0)
  var2  = dim_avg_n_Wrap(vtp2,0)
  var3  = dim_avg_n_Wrap(vtp3,0)

  ;;differences between E3SM and ERA5;;;
  ;;and the two methods in E3SM
  var3  = var3 - var2
  var2  = var2 - var1
  delete([/vtp1,vtp2,vtp3/])
  
  ;;;;calculate the global mean and correlation 
  latrg   = lat({MinLat:MaxLat})
  rad     = 4.*atan(1.)/180.    ;
  cost    = cos(latrg*rad)       ; cosine weights

  mvar1   = wgt_areaave_Wrap (var1({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar2   = wgt_areaave_Wrap (var2({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
  mvar3   = wgt_areaave_Wrap (var3({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)

  mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
  mstr2   = "Mean: " + sprintf("%5.2f", mvar2)
  mstr3   = "Mean: " + sprintf("%5.2f", mvar3)

  ;----------------------------------------------------------------------
  ; Plotting section

   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   fignam = OUTDir+"fig_2d_map_"+Groups(igp)+"_vs_era5_"+Season(isea)

   print("")
   print("plotting " + fignam)
   print("")

   wks1 = gsn_open_wks(wtype,fignam)
   plot = new(nlabs,graphic)

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
    do m = 0,nlabs-1 
      plot(m) = gsn_csm_map(wks1,mpres)
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
    res1@cnLevels                      := (/-5,-2,0,2,5,10,20,30,60,90,120,150,180/)
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

    res1@lbTitleString                  = Rstring1
    res1@gsnLeftString                  = LabList(0)
    res1@gsnRightString                 = mstr1
    plot_org1  = gsn_csm_contour(wks1,var1,res1)
    overlay(plot(0),plot_org1)

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                = True
    res2@cnFillPalette               = "CBR_coldhot"
    res2@cnLevelSelectionMode        = "ExplicitLevels"
    res2@cnLevels                   := (/-30,-20,-10,-5,-1,1,5,10,20,30/)
  ; res2@cnFillColors               := colorgroup1
    res2@lbBoxLinesOn                = True
   ;res2@lbBoxSeparatorLinesOn       = False
   ;res2@lbLabelStride               = 15
    res2@cnLinesOn                   = False             ; turn off contour lines
    res2@cnLineLabelsOn              = False            ; turn off contour labels
    res2@cnInfoLabelOn               = False

    res2@lbTitleString               = Rstring2
    res2@gsnLeftString               = LabList(1)
    res2@gsnRightString              = mstr2
    plot_rgd2                        = gsn_csm_contour(wks1,var2,res2)
    overlay(plot(1),plot_rgd2)

    res2@lbTitleString               = Rstring2
    res2@gsnLeftString               = LabList(2)
    res2@gsnRightString              = mstr3
    plot_rgd3                        = gsn_csm_contour(wks1,var3,res2)
    overlay(plot(2),plot_rgd3)
      
;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = True
    pnres@gsnPanelMainFontHeightF    = 0.012
    pnres@gsnPanelMainString         = Season(isea)+" (1981 - 2010)"
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4        
    gsn_panel(wks1,plot,(/3,1/),pnres)

    delete([/res1,res,mpres,pnres/])
 end do 

end do 

end

