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

  Diag_dir  = "/global/cfs/cdirs/e3sm/zhan391/diag_out/cvdp_e3sm_amip_nudging_2008-2017"
  Groups    = (/"ERA5-PSL"/)
  LabList   = (/"ERA5"/)
  ngrps     = dimsizes(Groups)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"NA"/)
  nregs     = dimsizes(Regions)
  MinLat    =   0.
  MaxLat    =  90.
  MinLon    = 210
  MaxLon    = 360.

  Season    = (/"DJF","MAM","JJA","SON"/)
  VarList   = "nao_pattern_"+str_lower(Season)
  nvars     = dimsizes(VarList)

 ;;; flags to control the plot style
  lndcolor = "Grey28"

 ;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  do iv = 0,nvars-1,1
   
   do igp = 0,ngrps-1,1

    if(igp.eq.0) then 
      fl     = Diag_dir + "/"+Groups(igp)+"_en00.cvdp_data.psl.nam_nao.2008-2017.nc"
    else
      fl     = Diag_dir + "/"+Groups(igp)+"_en00.cvdp_data.2008-2017.nc"
    end if 

    f      = addfile(fl,"r")
    lat    = f->lat
    lon    = f->lon
    var    = f->$VarList(iv)$
    mstr   = var@pcvar
  
  ;----------------------------------------------------------------------
  ; Plotting section

    if (igp.eq.0) then 
      OUTDir  = "./"
      wtype = "png"
      wtype@wkWidth   = 2500
      wtype@wkHeight  = 2500
     ;wtype@wkPaperSize     = "A4"
     ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
     ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
      wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_" + VarList(iv))
      plot = new(ngrps,graphic)
      dum = new(ngrps,graphic)
    end if 

  ;---Resources to share between both plots
    res                       = True              ; Plot modes desired.
    res@gsnDraw               = False
    res@gsnFrame              = False
    res@gsnMaximize           = False              ; Maximize plot
    res@tmXTOn                = False
    res@tmYROn                = False

    BorderThick = 8.0
    res@tmBorderThicknessF    = BorderThick
    res@tmXBMajorThicknessF   = BorderThick
    res@tmXBMinorThicknessF   = BorderThick*0.5
    res@tmYLMajorThicknessF   = BorderThick
    res@tmYLMinorThicknessF   = BorderThick*0.5
    res@tmXTMajorThicknessF   = BorderThick
    res@tmXTMinorThicknessF   = BorderThick*0.5
    res@tmYRMajorThicknessF   = BorderThick
    res@tmYRMinorThicknessF   = BorderThick*0.5

    FontHeightF              = 0.024
    res@tiMainFontHeightF    = FontHeightF*1.2
    res@tmYLLabelFontHeightF = FontHeightF
    res@tmXBLabelFontHeightF = FontHeightF
    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF*0.6
    res@gsnStringFontHeightF = FontHeightF*1.5

    res@cnFillOn                 = True              ; color plot desired
    res@cnLinesOn                = False             ; turn off contour lines
    res@cnLineLabelsOn           = False             ; turn off contour labels
   ;res@cnFillMode               = "RasterFill"      ; turn raster on
    res@cnMissingValFillColor    = "White"
    res@lbLabelBarOn             = False      ; Will turn on in panel later
   ;res@pmLabelBarParallelPosF   = 0.61 
    res@pmLabelBarOrthogonalPosF = 0.05
    res@lbOrientation            = "vertical"        ; vertical label bar
    res@pmLabelBarHeightF        = 0.65
    res@pmLabelBarWidthF         = 0.12
    res@lbTitleString            = "" 
    res@lbTitlePosition          = "Right"         ; title location
    res@lbTitleDirection         = "Across"        ; letter angle
    res@lbTitleAngleF            = 90.             ; title angle
    res@lbTitleFontHeightF       = FontHeightF/1.2 ; font height

   ;---Control appearance of map.
    mpres@mpLimitMode            = "LatLon"
    mpres@mpMinLonF              = MinLon
    mpres@mpMaxLonF              = MaxLon
    mpres@mpMinLatF              = MinLat
    mpres@mpMaxLatF              = MaxLat
    mpres@mpCenterLatF           = 0.5*(MinLat +MaxLat)
    mpres@mpCenterLonF           = 0.5*(MinLon +MaxLon)

    mpres@gsnMajorLatSpacing     = 30
    mpres@gsnMajorLonSpacing     = 30
    mpres@gsnMinorLatSpacing     = 10
    mpres@gsnMinorLonSpacing     = 10

    mpres@mpDataBaseVersion      = "MediumRes"         ; better resolution
    mpres@mpDataSetName          = "Earth..4"
    mpres@mpOutlineBoundarySets  = "GeophysicalAndUSStates"; turn on states
   ;mpres@mpOutlineBoundarySets  = "AllBoundaries"     ; more outlines

    mpres@mpFillOn               = False        ; turn off map fill
    mpres@mpPerimOn              = True
    mpres@mpPerimLineThicknessF  = BorderThick 
    mpres@mpFillOn               = False
    mpres@mpGeophysicalLineColor = "Black"  ; color of cont. outlines
    mpres@mpGeophysicalLineThicknessF = 3.0          ; thickness of outlines
    mpres@mpOutlineOn            = True
;   mpres@mpOceanFillColor       = "lightskyblue1"
;   mpres@mpLandFillColor        = "gray33"

    mpres@tmYLLabelsOn           = True
    mpres@tmXBLabelsOn           = True
    mpres@tmYLOn                 = True
    mpres@tmYROn                 = False
    mpres@tmXBOn                 = True

    plot(igp) = gsn_csm_map(wks1,mpres)

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
    res1@cnFillPalette                  = "BlueDarkRed18"
    res1@cnLevelSelectionMode           = "ExplicitLevels"
    res1@cnLevels                      := (/-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4/) 
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
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

    res1@gsnRightString                 = mstr
    res1@gsnLeftString                  = "(a) NAO regress on PSL (" + Season(iv)+")"
    dum(igp)  = gsn_csm_contour(wks1,var,res1)
    overlay(plot(igp),dum(igp))

    delete([/res1,res/])

  end do 

;---Draw both plots in a panel
  pnres                            = True
  pnres@gsnMaximize                = True
  pnres@gsnPanelLabelBar           = True 
  pnres@pmLabelBarParallelPosF     = 0.02
  pnres@pmLabelBarOrthogonalPosF   = -0.02
  pnres@lbLabelFontHeightF         = FontHeightF * 0.5
  pnres@lbOrientation              = "Horizontal" 
  pnres@pmLabelBarHeightF          = 0.06
  pnres@pmLabelBarWidthF           = 0.44
  pnres@lbTitleString              = "PSL Anomalies (hPa)";
  pnres@lbTitlePosition            = "Top" ;"Right"         ; title location
  pnres@lbTitleFontHeightF         = FontHeightF * 0.6
 ;pnres@lbTitleDirection           = "Across"        ; letter angle
 ;pnres@lbTitleAngleF              = 90.          ; title angle
  pnres@gsnPanelMainString         = "" ;"NAO pattern (" + Season(iv)+")"
  pnres@gsnPanelMainFontColor      = "Black"
  pnres@gsnPanelMainFontHeightF    = FontHeightF*0.95
  pnres@gsnPanelYWhiteSpacePercent = 2
  pnres@gsnPanelXWhiteSpacePercent = 2
  gsn_panel(wks1,plot,(/2,2/),pnres)

 delete([/pnres/])

end do 

end

