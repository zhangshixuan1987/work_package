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

  DATA_DIR  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/energy_balance/data/"
  DiagTime  = "1981-2010"

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Polar"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    = -55.
  MinLon    =   0.
  MaxLon    = 360.

  Season    = (/"ANN","DJF","JJA"/)
  nseason   = dimsizes(Season)
  
 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups    = (/"ERA5","EAMv2"/)
  ngrps     = dimsizes(Groups)
  LabList   = (/"ERA5", "v2.LR.F20TR"/)
  nlabs     = dimsizes(LabList)

  Rstring1     = "LHFLX (W m~S~-2~N~)"
  Rstring2     = "~F33~D~N~~F~LHFLX (W m~S~-2~N~)"

  ;;; flags to control the plot style
  lndcolor = "Grey28"

 do igp = 0,ngrps-1,1
 
  fil1  = DATA_DIR + "Surface_energy_fluxes_"+Groups(igp)+".nc"
  f     = addfile(fil1,"r")
  
  lat    = f->lat
  lon    = f->lon

  lhflx0 = f->LHF_LAND
  lhflx1 = f->LHF

  ;;;prepare the data for plot 

  do isea = 0,nseason-1,1

    if (Season(isea).eq."ANN") then 
      vtp1  = month_to_annual(lhflx0,1)
      vtp2  = month_to_annual(lhflx1,1)
    else
      vtp1  = month_to_season(lhflx0,Season(isea))
      vtp2  = month_to_season(lhflx1,Season(isea))
    end if 

    var1  = dim_avg_n_Wrap(vtp1,0)
    var2  = dim_avg_n_Wrap(vtp2,0)
    ;;differences between two methods
    var2  = var2 - var1
    delete([/vtp1,vtp2/])
  
  ;;;;calculate the global mean and correlation 
    latrg   = lat({MinLat:MaxLat})
    rad     = 4.*atan(1.)/180.    ;
    cost    = cos(latrg*rad)       ; cosine weights
    mvar1   = wgt_areaave_Wrap (var1({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
    mvar2   = wgt_areaave_Wrap (var2({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
    mstr1   = "Mean: " + sprintf("%5.2f", mvar1)
    mstr2   = "Mean: " + sprintf("%5.2f", mvar2)

    delete([/mvar1,mvar2,latrg,rad,cost/])

  ;----------------------------------------------------------------------
  ; Plotting section

   if (isea.eq.0) then 
     OUTDir  = "./"
     wtype = "png"
     wtype@wkWidth = 2500
     wtype@wkHeight = 2500
    ;wtype@wkPaperSize     = "A4"
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
    ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
     fignam = OUTDir+"fig_polar_map_"+Groups(igp) ;+Season(isea)

     print("")
     print("plotting " + fignam)
     print("")
  
     wks1 = gsn_open_wks(wtype,fignam)
     plot = new(nseason*2,graphic)
     dum1 = new(nseason*2,graphic)
     dum2 = new(nseason*2,graphic)

    end if 

  ;---Resources to share between both plots
    res                      = True              ; Plot modes desired.
    res@gsnDraw              = False
    res@gsnFrame             = False
    res@gsnMaximize          = False              ; Maximize plot

    res@gsnPolarLabelFontHeightF      = 0.02
    res@gsnPolar                      = "SH"
    res@mpFillOn                      = False        ; turn off map fill
    res@mpMaxLatF                     = MaxLat
    res@mpGeophysicalLineThicknessF   = 2.0

    res@tmXTOn                        = False
    res@tmYROn                        = False
    res@tmBorderThicknessF            = 4.0
    res@tmXBMajorThicknessF           = 4.0
    res@tmXBMinorThicknessF           = 4.0
    res@tmYLMajorThicknessF           = 4.0
    res@tmYLMinorThicknessF           = 4.0
    res@tmYRMajorThicknessF           = 4.0
    res@tmYRMinorThicknessF           = 4.0

    FontHeightF = 0.020
    res@tiMainFontThicknessF          = 3.0
    res@tiMainFontHeightF             = FontHeightF*1.2
    res@tmYLLabelFontHeightF          = FontHeightF
    res@tmXBLabelFontHeightF          = FontHeightF
    res@tiXAxisFontHeightF            = FontHeightF
    res@tiYAxisFontHeightF            = FontHeightF
    res@gsnStringFontHeightF          = FontHeightF
    res@cnLineLabelFontHeightF        = FontHeightF
    res@lbLabelFontHeightF            = FontHeightF*0.8
    res@lbTitleFontHeightF            = FontHeightF*0.8

    res@cnFillOn                 = True              ; color plot desired
   ;res@cnFillPalette            = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"

    res@lbLabelBarOn             = True
    res@lbOrientation            = "Vertical"
    res@lbPerimOn                = False             ; default has box
   ;res@pmLabelBarParallelPosF   = 0.61
    res@pmLabelBarOrthogonalPosF = 0.15
    res@lbLabelFontHeightF       = 0.02               ; make labels smaller
   ;res@pmLabelBarHeightF        = 0.38
   ;res@pmLabelBarWidthF         = 0.1
    res@lbTitleString            = ""
    res@lbTitlePosition          = "Right"                           ; title location
    res@lbTitleDirection         = "Across"                          ; letter angle
    res@lbTitleAngleF            = 90.                               ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2                              ; font height

    res@gsnLeftString                  = ""
    res@gsnRightString                 = ""
    res@gsnCenterString                = ""

;---Resources for plotting original data
    res1=res
    res1@gsnAddCyclic                   = True
    res1@cnFillPalette                  = "WhBlGrYeRe"
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                      := (/-5,-2,0,2,5,8,10,12,16,20,30,40,50,60/)
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

    res1@tiMainString                   = Season(isea)+" ("+DiagTime+")"
    res1@lbTitleString                  = Rstring1
    res1@gsnLeftString                  = LabList(igp) + " (Land)" 
    res1@gsnRightString                 = mstr1
    plot(isea*2) = gsn_csm_contour_map(wks1,var1,res1)

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                = True
    res2@cnFillPalette               = "CBR_coldhot"
    res2@cnLevelSelectionMode        = "ExplicitLevels"
    res2@cnLevels                   := (/-8,-4,-2,-1,-0.1,0.1,1,2,4,8/)
  ; res2@cnFillColors               := colorgroup1
    res2@lbBoxLinesOn                = True
   ;res2@lbBoxSeparatorLinesOn       = False
   ;res2@lbLabelStride               = 15
    res2@cnLinesOn                   = False             ; turn off contour lines
    res2@cnLineLabelsOn              = False            ; turn off contour labels
    res2@cnInfoLabelOn               = False

    res2@tiMainString                = Season(isea)+" ("+DiagTime+")"
    res2@lbTitleString               = Rstring2
    res2@gsnLeftString               = LabList(igp) + " (New - Land)"
    res2@gsnRightString              = mstr2
    plot(isea*2+1) = gsn_csm_contour_map(wks1,var2,res2)
    delete([/res,res1,res2/])
    delete([/var1,var2/])
   end do ; end of season

;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = True
    pnres@gsnPanelMainFontHeightF    = 0.012
    pnres@gsnPanelMainString         = "" ;Season(isea)+" (1981 - 2010)"
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4        
    gsn_panel(wks1,plot,(/3,2/),pnres)
    delete([/dum1,dum2,plot/])

    delete([/lat,lon,lhflx0,lhflx1/])

end do 

end

