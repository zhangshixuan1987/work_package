;read ascii data (OMI & GOME2) and write to netcdf file
; 5th Nov. 2014
; by Jianfeng Li

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

  Explist   = (/"true", "pred"/)
  Labels    = (/"Truth (Nudging)", "ML prediction (Unet)"/)
  ngrps     = dimsizes(Explist)
  timstr    = "2015_N20"
  ilev      = 48 

  fnum      = (/"(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"/)
  rstr      = "" 
  lev_p     = 500 
  varList   = (/"Nudge_U","Nudge_V","Nudge_T","Nudge_Q"/)
  varLsub   = (/"Nudge_U","Nudge_V","Nudge_T","Nudge_Q"/)
  varUnit   = (/"m s~S~-1~N~ day~S~-1~N~","m s~S~-1~N~ day~S~-1~N~","K day~S~-1~N~","g kg~S~-1~N~ day~S~-1~N~"/)
  varDecs   = (/"U tend", "V tend","T tend","Q tend"/)
  varFacs   = (/1.0,1.0,1.0,1.0/)
  nvars     = dimsizes(varList)

  do i = 0,ngrps-1,1
  do j = 0,nvars-1,1
   varnam    = varList(j)
   varunt    = varUnit(j)
   varfac    = varFacs(j)
   vardec    = varDecs(j)

   flnm0     = systemfunc("ls " + "./rgd_model_"+Explist(0)+"*.nc")
   f0        = addfile (flnm0, "r")
   flnm1     = systemfunc("ls " + "./rgd_model_"+Explist(1)+"*.nc")
   f1        = addfile (flnm1, "r")

   vtmp0     = dim_avg_n_Wrap(f0->$varnam$(:,ilev,:,:),0)
   vtmp1     = dim_avg_n_Wrap(f1->$varnam$(:,ilev,:,:),0)
   printVarSummary(vtmp0)
   printVarSummary(vtmp1)
   printMinMax(vtmp0,False)
   printMinMax(vtmp1,False)
   clat = cos(0.01745329*f0->lat)

   lfstr  = fnum(i*nvars+j) +" " + Labels(i) ;+ "(" + vardec+")"
   if (i.eq.0) then 
     vtend = vtmp0 
     mean = wgt_areaave(vtend, clat, 1.0, 1)
     colmap = "BlAqGrYeOrReVi200"
     rstr   = "Mean: " +sprintf("%.2f",mean)
   else
     vtend = vtmp1 - vtmp0
     rmse = wgt_arearmse(vtmp0,vtmp1,clat,1.0,0)
     pcor = pattern_cor(vtmp0,vtmp1,clat,1)
     colmap = "BlueWhiteOrangeRed" ;"cmocean_haline"
     rstr = "RMSE/CORR: " + sprintf("%.2f",rmse) + "/" + sprintf("%.2f",pcor)
   end if 

   delete([/f0,flnm0,f1,flnm1/])

   plot_sig_test = False 
   lbstr  = vardec + " " + " (" + varunt+")"

   if (varnam.eq."Nudge_U") then 
    if (i.eq.0) then
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/) 
    else
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/)
    end if 
   end if 
   if (varnam.eq."Nudge_V") then
    if (i.eq.0) then
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/) 
    else
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/)
    end if
   end if
   if (varnam.eq."Nudge_T") then
    if (i.eq.0) then
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/) 
    else
      cnlev    = (/-3.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,3.0/)
    end if
   end if
   if (varnam.eq."Nudge_Q") then
    if (i.eq.0) then
      cnlev    = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.8,1.2/)
    else
      cnlev    = (/-1.2,-0.8,-0.4,-0.2,-0.1,0.1,0.2,0.4,0.8,1.2/)
    end if
   end if

   if(i.eq.0.and.j.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
   ;wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
    plotname = "./diag_figure/fig_2d_map_ndgtend_" + varnam 
    wks = gsn_open_wks(wtype,plotname)
    ;;colors for the original plot
    orgcolors := read_colormap_file("precip3_16lev")
    ;colors for the difference plot
    gsn_merge_colormaps(wks,"temp_diff_18lev","WhViBlGrYeOrReWh")
    cncolors  := (/72,62,55,5,7,10,0,120,14,16,91,85,77/)
    plot = new(nvars*ngrps,graphic)
   end if

   res                  = True
   res@gsnDraw          = False
   res@gsnFrame         = False
   res@gsnMaximize      = False
   res@tmXTOn           = False
   res@tmYROn           = False

   res@gsnLeftString    = ""
   res@gsnRightString   = ""
   res@tiMainString     = ""

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

   FontHeightF = 0.022
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF
   res@gsnStringFontHeightF = FontHeightF*5.0
   res@tiXAxisFontHeightF   = FontHeightF
   res@tiYAxisFontHeightF   = FontHeightF

   mresP = res
  ;mresP@pmTickMarkDisplayMode       = "Always"; use NCL default lat/lon labels
   mresP@gsnMajorLatSpacing          = 30
   mresP@gsnMajorLonSpacing          = 60
   mresP@gsnMinorLatSpacing          = 10
   mresP@gsnMinorLonSpacing          = 10
   mresP@mpCenterLonF                = 240.0       ; choose center lon
  ;mresP@mpCenterLatF                = 30.         ; choose center lat
  ;mresP@mpGridSpacingF              = 10.
   mresP@mpLabelsOn                  = False
   mresP@mpPerimOn                   = False             ; turn off box around plot(i)
   mresP@mpFillOn                    = False ;True
   mresP@mpLimitMode                 = "LatLon"
   mresP@mpMinLonF                   = -180.0
   mresP@mpMaxLonF                   =  180.0
   mresP@mpMinLatF                   = -90.0
   mresP@mpMaxLatF                   = 90.0
   mresP@mpDataBaseVersion           = "MediumRes"    ; better map outlines
   mresP@mpGridAndLimbOn             = True           ; turn on limb and grid lines
   mresP@mpGridLineColor             = -1             ; ...but don't draw grid lines...
   mresP@mpGridAndLimbDrawOrder      = "PostDraw"   ; draw lines first, so
   mresP@mpOutlineOn                 = True
   mresP@mpOutlineDrawOrder          = "PostDraw"      ; draw continental outline last
   mresP@mpOutlineBoundarySets       = "GeophysicalAndUSStates" ; state boundaries
   mresP@mpGeophysicalLineThicknessF = 1.0 
   mresP@mpUSStateLineThicknessF     = 1.0 
   mresP@mpNationalLineThicknessF    = 1.0 
   mresP@mpLimbLineThicknessF        = 1.0
   mresP@mpPerimLineColor            = 0 ;-1
   mresP@mpPerimOn                   = True ;False          ; turn off map perimeter
   mresP@mpOceanFillColor            = "skyblue1"
   mresP@mpLandFillColor             = "transparent"
   mresP@mpInlandWaterFillColor      = "skyblue1" 
   mresP@mpFillDrawOrder             = "PostDraw"
 
   plot(i*nvars+j)                   = gsn_csm_map(wks,mresP)

   ref_res                          = res
   ref_res@cnFillOn                 = True
   ref_res@cnLevelSelectionMode     = "ExplicitLevels"
  ;ref_res@cnFillMode               = "RasterFill"
   ref_res@cnRasterSmoothingOn      = True
   ref_res@cnOutOfRangeFillColor    = "white"
   ref_res@cnInfoLabelOn            = False
   ref_res@cnLineLabelsOn           = False
   ref_res@cnLinesOn                = False
   ref_res@cnLineLabelBackgroundColor = -1    ; transparent
   ref_res@cnFillDrawOrder            = "Draw";"PostDraw"

   ref_res@lbLabelBarOn             = True   ; Will turn on in panel later
   ref_res@pmLabelBarParallelPosF   = 0.50 ; 0.50
   ref_res@pmLabelBarOrthogonalPosF = 0.15 ; 0.04
  ;ref_res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
   ref_res@pmLabelBarHeightF        = 0.12 ;0.60
   ref_res@pmLabelBarWidthF         = 0.60 ;0.06

   ref_res@lbLabelFontHeightF       = FontHeightF*0.75
   ref_res@lbTitleFontHeightF       = FontHeightF*0.8
   ref_res@lbTitleString            = lbstr 
  ;ref_res@lbTitlePosition          = "Right"         ; title location
  ;ref_res@lbTitleDirection         = "Across"        ; letter angle
  ;ref_res@lbTitleAngleF            = 90.          ; title angle

   ref_res@gsnRightString               := rstr 
   ref_res@gsnRightStringOrthogonalPosF  = 0.03
   ref_res@gsnLeftString                := lfstr
   ref_res@gsnLeftStringOrthogonalPosF   = 0.02
   ref_res@gsnCenterString               = "" ;timstr
   ref_res@gsnCenterStringOrthogonalPosF = -0.05
   ref_res@gsnCenterStringFontColor      = "Black"
   ref_res@gsnLeftStringFontHeightF      = FontHeightF*1.5
   ref_res@gsnRightStringFontHeightF     = FontHeightF*1.5
   ref_res@gsnCenterStringFontHeightF    = FontHeightF*1.5
  ;ref_res@gsnCenterStringParallelPosF   = 0.90
  ;ref_res@gsnRightStringParallelPosF    = 0.85

    if(.not.isvar("contou2")) then
      contou2 = new(100,graphic)
    end if
    res2                             = ref_res
    res2@gsnAddCyclic                = True
    res2@cnLevelSelectionMode       := "ExplicitLevels" ;"ManualLevels"
    res2@cnLineThicknessF           := 1.0
    res2@cnLevels                    = cnlev
   ;res2@cnFillColors                = cncolors
    res2@cnFillPalette               = colmap
   ;res2@lbLabelStride               = 2
   ;res2@lbBoxLinesOn                = True
   ;res2@lbBoxSeparatorLinesOn       = False
    contou2(i*nvars+j) = gsn_csm_contour(wks,vtend,res2)
    overlay(plot(i*nvars+j),contou2(i*nvars+j))
    delete([/res,res2,ref_res,mresP/])
  end do 
  end do 

;---Specify y positions of each row of plots.
  pres                          = True
; pres@gsnPanelDebug            = True      ; prints debug info
  pres@gsnFrame                 = False
  pres@gsnMaximize              = False
  pres@gsnPanelMainFontHeightF  = 0.015
  pres@gsnPanelFigureStrings    = ""
  pres@gsnPanelFigureStringsFontHeightF  = 0.01

  pres2 = pres
; pres2@vpXF                     = 0.15
; pres2@vpYF                     = 0.3
; pres2@vpWidthF                 = 0.7
; pres2@vpHeightF                = 0.18
  pres2@gsnPanelLabelBar         = False
  pres2@pmLabelBarParallelPosF   = 0.02
  pres2@pmLabelBarOrthogonalPosF = -0.02
  pres2@pmLabelBarHeightF        = 0.05 ;0.60
  pres2@pmLabelBarWidthF         = 0.85 ;0.06
  pres2@lbLabelFontHeightF       = FontHeightF*0.40
  pres2@lbTitleFontHeightF       = FontHeightF*0.45
  pres2@lbTitleString            = lbstr
  pres2@gsnPanelTop              = 0.5
  pres2@gsnPanelBottom           = 0.0
 ;pres2@gsnPanelYF              := (/ypos2,ypos2,ypos2/)
 ;pres2@gsnPanelMainString       = ""
 ;pres2@gsnPanelMainPosYF        = ytitl2    ; Y position of title
  gsn_panel(wks,plot,(/2,nvars/),pres2)

  frame(wks)


end
