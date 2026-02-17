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

  ldatalocal = True
  ldatalocal = True

  outdir    = "../data/"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Groups    = (/"NDGUVT_tau6","NDGUVT_tau6", "NDGUVT_tau3"/)
  Expnams   = Groups
  ngrps     = dimsizes(Groups)
  LabList   = (/"NDGUV (tau = 6h)","NDGUVT (tau = 6h)", "NDGUVT (tau = 3h)"/)
  tim_str   = "2010-07"

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

  VarMin    = (/-5,-5,-4,-1/)
  VarMax    = (/ 5, 5, 4, 1/)

  VtdMin    = (/-1,-1,-2,-1/)
  VtdMax    = (/ 1, 1, 2, 1/)

  nvars     = dimsizes(VarList)

  do j = 0,nvars-1,1

    Varnam  = VarList(j)
    Varunt  = VarUnit(j)
    Varstr  = VarName(j)
    Vtdnam  = VtdName(j)
    Vtdunt  = VtdUnit(j)
    Vtdstr  = VtdName(j)

   do i = 1,ngrps,1

    fref = outdir + "ERA5" + "_" + tim_str+".nc"
    fclm = outdir + "CLIM" + "_" + tim_str+".nc"
    fndg = outdir + "NDGUV" + "_" + tim_str+".nc"
    ften = outdir + "E3SMv2_" + Groups(i-1) + "_"+tim_str+".nc"
    f0   = addfile(fref, "r")
    f1   = addfile(fclm, "r")
    f2   = addfile(ften, "r")

    vr    = dim_avg_n_Wrap(f0->$Varnam$(0,:,:,:),2)
    vs    = dim_avg_n_Wrap(f1->$Varnam$(0,:,:,:),2)
    vb    = vr
    vb    = vs - vr
    vb    = vb * VarFrac (j)

    varnm = "Nudge_" + Varnam
    var   = dim_avg_n_Wrap(f2->$varnm$(0,:,:,:),2)
    var   = var * VtdFrac (j)

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
     fignam = OUTDir+"fig_lat_hgt_ndg_tend_"+Varnam+"_"+ tim_str
     wks1 = gsn_open_wks(wtype,fignam)
     plot = new(ngrps+1,graphic)
     dum  = new(ngrps+1,graphic)
  end if

;----------------------------------------------------------------------
; Plotting section
;---Resources to share between both plots
  res                       = True              ; Plot modes desired.
  res@gsnDraw               = False
  res@gsnFrame              = False
  res@gsnMaximize           = False              ; Maximize plot
  res@gsnAddCyclic          = False
  res@vpHeightF             = 0.6
  res@vpWidthF              = 0.8

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"


  FontHeightF = 0.028
  res@tiMainFontThicknessF    = 1.0
  res@tiMainFontHeightF       = FontHeightF*1.4
  res@tmYLLabelFontHeightF    = FontHeightF
  res@tmXBLabelFontHeightF    = FontHeightF
  res@tiXAxisFontHeightF      = FontHeightF*1.1
  res@tiYAxisFontHeightF      = FontHeightF*1.1
  res@gsnStringFontHeightF    = FontHeightF*1.1
  res@cnLineLabelFontHeightF  = FontHeightF
  res@lbLabelFontHeightF      = FontHeightF
  res@lbTitleFontHeightF      = FontHeightF   

  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
  res@lbLabelStride           = 1                 ; skip every other label
  res@lbPerimOn               = False             ; default has box
  res@lbOrientation           = "vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.6
  res@pmLabelBarWidthF        = 0.12
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15
  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"
  res@tmYLOn                  = True
  res@tiYAxisString           = "Pressure (hPa)"

  res@trYMinF                 = 100.0                   ; min value on y-axis
 ;res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/975, 850, 700,600, 500,400,300,200,100/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings
  res@trYLog                  = False

   BorderThick = 1.0
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   res2=res
   res2@cnLineLabelsOn                 = False            ; turn off contour labels
   res2@cnInfoLabelOn                  = False
   res2@cnLineLabelDensityF            = 2.0
   res2@cnLineLabelFontHeightF         = 0.018
   res2@cnLineLabelBackgroundColor     = -1
   res2@gsnContourLineThicknessesScale = 1.0
   res2@cnLineThicknessF               = 1.0
   res2@cnLineColor                    = "Black"
   res2@gsnContourZeroLineThicknessF   = 2.
   res2@gsnContourNegLineDashPattern   = 1
   res2@gsnContourZeroLineThicknessF   = 2.
   res2@gsnContourNegLineDashPattern   = 1
   res2@tiMainString                   = "";"CNTL (PD1800)"
   res2@gsnLeftString                  = "~F33~D~N~~F~"+Varstr+" ("+ Varunt +")"  
   res2@gsnLeftStringOrthogonalPosF    = 0.02
  ;res2@gsnLeftStringParallelPosF      = 0.2 
   res2@gsnRightString                 = "" ;Varstr ;"~F33~D~N~~F~"+ Varstr ;Varunt
   res2@gsnRightStringOrthogonalPosF   = 0.02

   res11 = res2
   res11@cnFillPalette                 := "NCV_jaisnd" ;"MPL_rainbow"
   res11@cnLevelSelectionMode          := "ExplicitLevels"
   res11@cnLevels                      := fspan(VtdMin(j),VtdMax(j),11)
   res11@gsnLeftString                  = LabList(i-1)
   res11@gsnRightString                 = Vtdnam + " (" + Vtdunt + ")"
   res11@gsnCenterString                = "" ;mstd
   res11@lbTitleString                  = "" ;"("+ Varunt +")"
   plot(i)   = gsn_csm_pres_hgt(wks1,var,res11)

   if(i.eq.1)then
    res12 = res2
    res12@cnFillPalette                 := "temp_diff_18lev" ;"MPL_rainbow"
    res12@cnLevelSelectionMode          := "ExplicitLevels"
    res12@cnLevels                      := (/-4,-3,-2,-1,-0.5,0,0.5,1,2,3,4/) ;fspan(VarMin(j),VarMax(j),11)
    res12@gsnLeftString                  = "Mean model bias (CLIM - ERA5)"
    res12@gsnRightString                 = "~F33~D~N~~F~"+Varnam + " (" + Varunt + ")"
    res12@gsnCenterString                = "" ;mstr
    res12@lbTitleString                  = "" ;"("+ Varunt +")"
    plot(i-1) = gsn_csm_pres_hgt(wks1,vb,res12)
   end if

   delete([/res,res2/])

 end do

;---Draw both plots in a panel
 pnres                         = True
 pnres@gsnFrame                = True
 pnres@gsnMaximize             = False
 pnres@gsnPanelLabelBar        = True
 pnres@gsnPanelMainString      = " "; "Surface Temperature ("+ts@units+")"
 pnres@gsnPanelLabelBar        = False                ; add common colorbar
 pnres@gsnPanelMainFontColor   = "Black"
 pnres@gsnPanelMainFontHeightF = 0.030
;pnres@lbLabelFontHeightF      = 0.018               ; make labels smaller
;pnres@pmLabelBarParallelPosF  = 0.026
;pnres@pmLabelBarOrthogonalPosF= 0.015
;pnres@pmLabelBarHeightF       = 0.1
;pnres@pmLabelBarWidthF        = 0.57
;pnres@gsnPanelBottom          = 0.05
 pnres@gsnPanelYWhiteSpacePercent = 4
 pnres@gsnPanelXWhiteSpacePercent = 4
 gsn_panel(wks1,(/plot(:)/),(/1,ngrps+1/),pnres)

end do

end
