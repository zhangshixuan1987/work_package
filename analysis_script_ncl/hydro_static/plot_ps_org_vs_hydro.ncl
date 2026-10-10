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

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = -90.
  MaxLat    =  90.
  MinLon    =   0.
  MaxLon    = 360.
  grav      = 9.80665
  Rd        = 287.053 

 ;;;;;;;date strings for the data directory;;;;;;;;
  Groups      = (/"ERA5","EAMv2"/)
  ngrps       = dimsizes(Groups)
  DATA_DIR    = "/global/cfs/cdirs/e3sm/zhan391/local_diag/hydro_static/data/"
  filname     = DATA_DIR + "Hydrostatic_Surface_Pressure_"+Groups(0)+"_"+Groups(1)+".nc"

  npanel      = 6 
  LabList     = new(npanel,string)
  LabList(0)  = "ERA5 reanalysis"
  LabList(1)  = "E3SM"
  LabList(2)  = "Hydrostatic - Original  (ERA5)"
  LabList(3)  = "Hydrostatic - Original  (E3SM)"
  LabList(4)  = "E3SM - ERA5 (Original)"
  LabList(5)  = "E3SM - ERA5 (Hydrostatic)"
  Rstring     = "PS (hPa)"

 ;;;;;;;begin to work with the simulations;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  filname     = DATA_DIR + "Hydrostatic_Surface_Pressure_"+Groups(0)+"_"+Groups(1)+".nc"
  f           = addfile(filname,"r")

  pso         = dim_avg_n_Wrap(f->PS_ERA5(:,{MinLat:MaxLat},:),0)
  psm         = dim_avg_n_Wrap(f->PS_EAM(:,{MinLat:MaxLat},:),0)  

  pso0        = dim_avg_n_Wrap(f->PS_Hydro_ERA5(:,{MinLat:MaxLat},:),0)
  psm0        = dim_avg_n_Wrap(f->PS_Hydro_EAM(:,{MinLat:MaxLat},:),0)  
  
  pso         = pso / 100.0
  psm         = psm / 100.0

  pso0        = pso0 / 100.0
  psm0        = psm0 / 100.0

  var1        = pso
  var2        = psm

  ;;hyro - org 
  var3        = pso0 
  var3        = pso0 - pso
  var4        = psm0
  var4        = psm0 - psm

  ;;org - org 
  var5        = pso
  var5        = psm - pso 
  var6        = pso
  var6        = psm0 - pso0
  
  ;----------------------------------------------------------------------
  ; Plotting section

   OUTDir  = "./"
   wtype = "png"
   wtype@wkWidth = 2500
   wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   wks1 = gsn_open_wks(wtype,OUTDir+"fig_2d_map_e3sm_vs_era5_PS_vs_hydro")
   plot        = new(npanel,graphic)

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

    res@gsnRightString                 = Rstring

;---Resources for plotting original data
    res1=res
    res1@gsnAddCyclic                  = True
    res1@cnFillPalette                 = "WhBlGrYeRe"
    res1@cnLevelSelectionMode          = "ExplicitLevels"
    res1@cnLevels                     := (/850,900,950,960,970,980,990,1000,1010,1015,1020/);
   ;res1@cnLevels                     := (/965,970,975,985,990,995,1000,1005,1010,1015,1020/);
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
    plot_org1  = gsn_csm_contour(wks1,var1,res1)
    overlay(plot(0),plot_org1)

    res1@gsnLeftString                  = LabList(1)
    plot_org2  = gsn_csm_contour(wks1,var2,res1)
    overlay(plot(1),plot_org2)

    delete([/res1/])

   ;;;;difference between E3SM and ERA5
    res2=res
    res2@gsnAddCyclic                   = True
    res2@cnFillPalette                  = "BlueDarkRed18" ;"NCV_blu_red"
    res2@cnLevelSelectionMode           = "ExplicitLevels"
  ; res2@cnFillColors                  := colorgroup1
    res2@lbBoxLinesOn                   = True
   ;res2@lbBoxSeparatorLinesOn          = False
   ;res2@lbLabelStride                  = 15
    res2@cnLinesOn                      = False             ; turn off contour lines
    res2@cnLineLabelsOn                 = False            ; turn off contour labels
    res2@cnInfoLabelOn                  = False
    res2@gsnRightString                 = "~F33~D~N~~F~"+ Rstring
    res2@gsnCenterString                = "" ;"Global mean: "+ mdfstr
    res2@gsnCenterStringOrthogonalPosF  = 0.04

    res2@cnLevels                      := (/-12,-8,-4,-2,0,2,4,8,12/)
    res2@gsnLeftString                  = LabList(2)
    plot_rgd1                           = gsn_csm_contour(wks1,var3,res2)
    overlay(plot(2),plot_rgd1)
    res2@gsnLeftString                  = LabList(3)
    plot_rgd2                           = gsn_csm_contour(wks1,var4,res2)
    overlay(plot(3),plot_rgd2)

    res2@cnLevels                      := (/-12,-8,-4,-2,0,2,4,8,12/)
    res2@gsnLeftString                  = LabList(4)
    plot_rgd3                           = gsn_csm_contour(wks1,var5,res2)
    overlay(plot(4),plot_rgd3)
    res2@gsnLeftString                  = LabList(5)
    plot_rgd4                           = gsn_csm_contour(wks1,var6,res2)
    overlay(plot(5),plot_rgd4)

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
    gsn_panel(wks1,plot,(/3,2/),pnres)

   delete([/pnres/])


end 
