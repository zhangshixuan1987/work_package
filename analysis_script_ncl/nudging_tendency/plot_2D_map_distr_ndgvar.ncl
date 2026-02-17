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

  outdir    = "../data/"

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions   = (/"Global"/)
  nregs     = dimsizes(Regions)
  MinLat    = (/-90./)
  MaxLat    = (/90./)
  MinLon    = (/0./)
  MaxLon    = (/360./)
  tim_str   = "2010-01"
  plev      = 850

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir    = (/"NDGUVT_tau6"/)
  Groups    = Casdir
  ngrps     = dimsizes(Groups)
  LabList   = (/"NDGUVT (tau = 6h)"/)

  ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList   = (/"U","V","T","Q"/)

  VarName   = VarList
  VarUnit   = (/"m s~S~-1~N~", "m s~S~-1~N~", "K", "g kg~S~-1~N~"/)
  VarFrac   = (/1, 1,1,1e3/)
  VtdName   = VarList+ " tend"
  VtdUnit   = (/"m s~S~-2~N~", "m s~S~-2~N~", "K s~S~-1~N~", "g kg~S~-1~N~ s~S~-1~N~"/)

  VarFrac   = (/1, 1,1,1e3/)

  VtdFrac   = (/1e5, 1e5,1e5,1e8/)
  VtdUnit   = "10~S~-5~N~ " + VtdUnit

  VarMin    = (/-5,-5,-5,-1/)
  VarMax    = (/ 5, 5, 5, 1/)

  VtdMin    = (/-2,-2,-5,-1/)
  VtdMax    = (/ 2, 2, 5, 1/)

  nvars     = dimsizes(VarList)
  
  do j = 0,nvars-1,1
  
   do i = 1,ngrps,1
 
    Varnam  = VarList(j)
    Varunt  = VarUnit(j)
    Varstr  = VarName(j)
    Vtdnam  = VtdName(j)
    Vtdunt  = VtdUnit(j)
    Vtdstr  = VtdName(j)

    print("j : " + j)

    fref = outdir + "ERA5" + "_" + tim_str+".nc" 
    fclm = outdir + "CLIM" + "_" + tim_str+".nc"
    ften = outdir + "E3SMv2_" + Casdir(i-1) + "_"+tim_str+".nc"
    f0   = addfile(fref, "r")
    f1   = addfile(fclm, "r")
    f2   = addfile(ften, "r")

    print(" ")
    print("........................... ")
    print("i : " + i)
    print("........................... ")

    lev   = sprintf("%.2f",f0->lev({plev}))
    vr    = f0->$Varnam$(0,{plev},:,:)
    vs    = f1->$Varnam$(0,{plev},:,:)
    vb    = vr
    vb    = vs - vr 
    vb    = vb * VarFrac (j)

    varnm = "Nudge_" + Varnam
    var   = f2->$varnm$(0,{plev},:,:)
    var   = var * VtdFrac (j)

    print(min(var) + " " +max(var))
    printVarSummary(var)

    ;;;;calculate the global mean and correlation 
    latrg = var&lat
    latrg!0 = "lat"
    latrg&lat = latrg
    rad   = 4.*atan(1.)/180.    ;
    cost  = cos(latrg({MinLat:MaxLat})*rad)       ; cosine weights
    mvar  = wgt_areaave_Wrap (vb({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
    mstr  = "Mean: " + sprintf("%5.2f", mvar)
    mvar  = wgt_areaave_Wrap (var({MinLat:MaxLat},{MinLon:MaxLon}),cost,1.,0)
    mstd  = "Mean: " + sprintf("%5.2f", mvar)
    print(mstr)
    delete([/latrg,rad,cost/])
  ;----------------------------------------------------------------------
  ; Plotting section
   OUTDir  = "./"
   wtype = "pdf"
  ;wtype@wkWidth = 2500
  ;wtype@wkHeight = 2500
  ;wtype@wkPaperSize     = "A4"
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
   if (i.eq.1)then 
     fignam = OUTDir+"fig_2d_ndg_tend_"+Varnam+"_"+ tim_str
     wks1 = gsn_open_wks(wtype,fignam)
     plot = new(ngrps+1,graphic)
     dum  = new(ngrps+1,graphic)
     lndcolor = "Grey28"
   end if 

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

    lthick = 1.0
    res@tmBorderThicknessF    = lthick/2.0
    res@tmXBMajorThicknessF   = lthick/2.0
    res@tmYLMajorThicknessF   = lthick/2.0
    res@tmXBMinorThicknessF   = lthick/2.0
    res@tmYLMinorThicknessF   = lthick/2.0

    mpres                    = res

    res@lbLabelFontHeightF   = FontHeightF*0.6
    res@gsnStringFontHeightF = FontHeightF*1.5

    res@cnFillOn                 = True              ; color plot desired
    res@cnFillPalette            = "temp_19lev" ;"ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
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
    mpres@mpGeophysicalLineThicknessF = 1.0          ; thickness of outlines
    mpres@mpOutlineOn             = True 
    mpres@mpOutlineDrawOrder      = "PostDraw"
    mpres@mpFillDrawOrder         = "Predraw"
    mpres@mpGridAndLimbOn         = False 
    mpres@mpGridLineThicknessF    = 1.0
;   mpres@mpOceanFillColor        = "lightskyblue1"
;   mpres@mpLandFillColor         = "gray33"
    mpres@tmYLLabelsOn            = True
    mpres@tmYLOn                  = True
    plot(i) = gsn_csm_map(wks1,mpres)

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
    res1@cnFillMode                     = "RasterFill"
    res1@cnFillPalette                  = "cmp_b2r" ;"MPL_rainbow"
    res1@lbBoxLinesOn                   = True
   ;res1@lbBoxSeparatorLinesOn          = False
   ;res1@lbLabelStride                  = 15
    res1@gsnAddCyclic                   = True
    res1@tiMainString                   = ""
    res1@tiMainOffsetYF                 = 0.01
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

    res11 = res1
    res11@cnFillPalette                 := "NCV_jaisnd" ;"MPL_rainbow"
    res11@cnLevelSelectionMode          := "ExplicitLevels"
    res11@cnLevels                      := fspan(VtdMin(j),VtdMax(j),21) 
    res11@gsnLeftString                  = LabList(i-1)
    res11@gsnRightString                 = Vtdnam + " (" + Vtdunt + ")"
    res11@gsnCenterString                = "" ;mstd
    dum(i)  = gsn_csm_contour(wks1,var,res11)
    overlay(plot(i),dum(i))

    if(i.eq.1)then
     res12 = res1
     res12@cnFillPalette                 := "temp_diff_18lev" ;"MPL_rainbow"
     res12@cnLevelSelectionMode          := "ExplicitLevels"
     res12@cnLevels                      := fspan(VarMin(j),VarMax(j),11) 
     res12@gsnLeftString                  = "Mean model bias (CLIM - ERA5)" 
     res12@gsnRightString                 = Varnam + " (" + Varunt + ")"
     res12@gsnCenterString                = "" ;mstr
     plot(i-1) = gsn_csm_map(wks1,mpres)
     dum(i-1)  = gsn_csm_contour(wks1,vb,res12)
     overlay(plot(i-1),dum(i-1))
    end if 

    delete([/res1,res/])

   end do 

;---Draw both plots in a panel
    pnres                            = True
    pnres@gsnMaximize                = False
    pnres@gsnPanelMainFontHeightF    = 0.016
    pnres@gsnPanelMainString         = ""
    pnres@gsnPanelMainFontColor      = "Black"
    pnres@gsnPanelBottom             = 0.05
    pnres@gsnPanelYWhiteSpacePercent = 10
    pnres@gsnPanelXWhiteSpacePercent = 5
    pnres@gsnPanelBottom             = 0.4        
    gsn_panel(wks1,plot,(/1,ngrps+1/),pnres)

   delete([/pnres/])

  end do

end

