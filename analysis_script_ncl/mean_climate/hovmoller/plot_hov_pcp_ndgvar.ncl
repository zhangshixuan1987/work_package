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

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  Subreg   = (/ "North Pacific",  "CONUS",  "Tropics"/)
  Subrgs   = (/ "NP",             "CONUS",  "Trop"/)
  minlat   = (/  10.0,            30.0,     -10.0/)
  maxlat   = (/  65.0,            48.0,      10.0/)
  minlon   = (/ 110.0,            235.0,     60.0/)
  maxlon   = (/ 260.0,            292.0,    270.0/)
  nsubreg  = dimsizes(Subreg)
  outdir   = "../fig_data/"

 ;;;;;;strings for the data directory;;;;;;;;
  Groups    = (/"TRMM","CLIM","NDGUV","NDGUVT","NDGUVTQ"/)
  Expnams   = Groups
  ngrps     = dimsizes(Groups)

  season    = (/"DJF"/)
  days      = (/31,31,28/)
  plot      = new(ngrps,graphic)

 ;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList   = (/"PRECT"/)
  UnitList  = (/"mm day~S~-1~N~"/)
  PltNums   = (/"a", "b", "c", "d", "e", "f", "g", "h"/)
  nvars     = dimsizes(VarList)

 do irg = 0,nsubreg-1,1

  do iv = 0,nvars-1, 1

  Varname = VarList(iv)
  Varunt  = UnitList(iv)

  do igp = 0,ngrps-1,1

   fna   = outdir+ Varname + "_" + Groups(igp) +"_"+str_sub_str(Subreg(irg)," ","_")+ ".nc"
   fla   = addfile(fna,"r")

   vnb    = Varname + "_" + Groups(igp)
   precm   = fla->$vnb$
   Ptimes  = fla->time
   datestr = cd_calendar(Ptimes, -2)
  ;print(datestr)

   tprecm = precm
   do im = 0,dimsizes(days)-1
     if(im.eq.0)then 
       ids0 = 0 
       ide0 = days(im) -1
     else
       ids0 = sum(days(0:im-1))
       ide0 = sum(days(0:im))-1 
     end if 
     print(ids0 + "  "+ide0)
     indmon = ind(Ptimes.ge.ids0.and.Ptimes.le.ide0)
     avgm   = dim_avg_n(precm(indmon,:),0)
     tprecm(indmon,:) = precm(indmon,:) - conform(precm(indmon,:),avgm,1)
     delete([/indmon,avgm,ids0,ide0/])
   end do

   if(igp.eq.0)then
     precm@pc   = -9999
     precm@rmse = -9999
     anmobs     = tprecm
     refobs     = precm
   else
     precm@pc   = sprintf("%3.2f",pattern_cor(tprecm, anmobs, 1.0, 1))
     precm@rmse = sprintf("%3.2f",wgt_arearmse(precm,refobs, 1.0, 1, 0))
   end if
   delete([/tprecm/])
   
  ;----------------------------------------------------------------------
  ; Plotting section
  ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
  wtype = "png"
  wtype@wkWidth = 2500
  wtype@wkHeight = 2500
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length

  if(igp.eq.0)then
    fignam  = "Fig_hov_" +Varname+"_"+ str_sub_str(Subreg(irg)," ","_")+"_ndgvar"
    print("plotname : " + fignam + ".pdf" )
    wks1 = gsn_open_wks(wtype,fignam)
  end if

  Figstr  = "" ; "("+PltNums (igp)+") "

  ;;;set up the plot information;;;;
  load "./var_share_colorbar.ncl"

;---Resources to share between both plots
  res                           = True                   ; plot mods desired
  res@vpHeightF                 = vpHeightF
  res@vpWidthF                  =  vpWidthF
 ;res@vpXF                      = 0.02
  res@gsnDraw                   = False
  res@gsnFrame                  = False
  res@gsnStringFontHeightF      = FontHeight
  res@gsnRightStringFontHeightF = FontHeight
  res@gsnLeftStringFontHeightF  = FontHeight
  res@tmXBLabelFontHeightF      = FontHeight
  res@tmYLLabelFontHeightF      = FontHeight
  res@tmYRLabelFontHeightF      = FontHeight
  res@tiMainFontHeightF         = FontHeight*1.1
  res@tiXAxisFontHeightF        = FontHeight
  res@tiYAxisFontHeightF        = FontHeight

  res@lbLabelBarOn            = False      ; Will turn on in panel later
 ;res@pmLabelBarParallelPosF  = 0.55
  res@pmLabelBarOrthogonalPosF= 0.15
  res@lbOrientation           = "Horizontal";"vertical"        ; vertical label bar
  res@pmLabelBarHeightF       = 0.25 ;0.65
  res@pmLabelBarWidthF        = 0.7 ;0.1
 ;res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
 ;res@lbTitlePosition         = "Right"                           ; title location
 ;res@lbTitleDirection        = "Across"                          ; letter angle
 ;res@lbTitleAngleF           = 90.                               ; title angle
 ;res@lbTitleFontHeightF      = 0.022       

  BorderThick = 2.0
  res@tmBorderThicknessF    = BorderThick
  res@tmXBMajorThicknessF   = BorderThick
  res@tmXBMinorThicknessF   = BorderThick*0.5
  res@tmYLMajorThicknessF   = BorderThick
  res@tmYLMinorThicknessF   = BorderThick*0.5
  res@tmXTMajorThicknessF   = BorderThick
  res@tmXTMinorThicknessF   = BorderThick*0.5
  res@tmYRMajorThicknessF   = BorderThick
  res@tmYRMinorThicknessF   = BorderThick*0.5

  ;res@gsnLeftStringParallelPosF    = -0.15
  ;res@gsnLeftStringOrthogonalPosF  =  0.06
  ;res@gsnRightStringOrthogonalPosF =  0.06
  ;res@gsnXYBarChart             = True              ; Create bar plot
  ;res@gsnXYBarChartOutlineOnly  = True
  ;res@trXLog                    = False
  ;res@trYLog                    = True

  res@tmXTOn                    = False
  res@tmYROn                    = True

  res@tmYLMode                  = "Explicit"
  res@tmYLMinorOn               = True
  res@tmYLMajorLengthF          = 0.015
  res@tmYLMajorOutwardLengthF   = 0.015
  res@tmYLMinorLengthF          = 0.01
  res@tmYLMinorOutwardLengthF   = 0.01
  res@tmYLMinorPerMajor         = 4
 ;res@tmYLMode                  = "Manual"
 ;res@tmYLTickSpacingF          = 2
 ;res@tmYLFormat                = "f"
 ;res@tmYLTickSpacingF          = yint
 ;res@trYMaxF                   = ymax
 ;res@trYMinF                   = ymin
  res@trYMinF                   = min(Ptimes)
  res@trYMaxF                   = max(Ptimes)
  res@tmYLValues                = Ptimes(::6)
  res@tmYLLabels                = datestr(::6)
  res@tmYLMinorValues           = Ptimes(::3)

  res@tmXBMode                  = "Explicit"
  res@tmXBMinorOn               = True
  res@tmXBMajorLengthF          = 0.015
  res@tmXBMajorOutwardLengthF   = 0.015
  res@tmXBMinorLengthF          = 0.01
  res@tmXBMinorOutwardLengthF   = 0.01
  res@tmXBMinorPerMajor         = 4
  res@tmXBLabelFontHeightF      = FontHeight
  res@tmXBLabelDeltaF           = -0.6
  ;res@tmXBMode          = "Manual"
  ;res@tmXBFormat        = "f"
  ;res@tmXBTickSpacingF  = xint
  ;res@trXMaxF           = xmax
  ;res@trXMinF           = xmin
  res@gsnMajorLonSpacing = 60.
  res@gsnMinorLonSpacing = 30. ;res@gsnMajorLonSpacing
  res@trXMinF            = precm&lon(0)
  res@trXMaxF            = precm&lon(dimsizes(precm&lon)-1)
  res@trYReverse         = True

  res@tmYLLabelJust     = "CenterRight"
  res@tmYRLabelJust     = "CenterRight"
  res@tmYUseLeft        = True
  res@tmYRLabelsOn      = False
  res@tiXAxisString     = "Longitude"
  res@tiYAxisString     = ""
  res@tiXAxisOffsetYF   = 0.
  res@tiYAxisOffsetXF   = 0.005
 ;res@xyComputeYMin     = True

; plotting
  ;res@gsnYRefLine            = (/-100.0/)            ; create a reference line
  ;res@gsnYRefLineDashPattern = 0
  ;res@gsnYRefLineColor       = "Grey";"Black"
  ;res@gsnYRefLineThicknessF  = 1.0

  res@tmYLLabelsOn           = False
  res@tmYLOn                 = False

  res=res
  res@cnFillOn             = True               ; turn on color fill
 ;res@cnFillPalette        = "WhiteBlueGreenYellowRed" ;"WhBlGrYeRe"
 ;res@cnLevelSelectionMode = "ManualLevels"     ; manual contour levels
 ;res@cnMinLevelValF       = ymin1
 ;res@cnMaxLevelValF       = ymax1
 ;res@cnLevelSpacingF      = yint1                 ; contour level spacing
  res@cnLevelSelectionMode  = "ExplicitLevels"
  res@cnLevels             := ylevs
  res@cnFillColors         := colorgroup1
  res@cnLinesOn                      = False ;plot_contour
  res@cnLineLabelsOn                 = False
  res@cnInfoLabelOn                  = False
  res@cnLineLabelDensityF            = 1.0
  res@cnLineLabelFontHeightF         = 0.018
  res@cnLineLabelBackgroundColor     = -1
  res@gsnContourLineThicknessesScale = 1.0
  res@cnLineThicknessF               = 1.0
  res@cnLineColor                    = "Black"
  res@gsnContourZeroLineThicknessF   = 1.
  res@gsnContourNegLineDashPattern   = 2

  res1 = res
  if(igp.eq.0) then 
    res1@tiMainString         = Figstr+ Expnams(igp) + " (" + Subreg(irg)+")"
  else
    res1@tiMainString         = Figstr+ Expnams(igp) 
  end if 

  if(precm@pc.eq.-9999)then
   res1@gsnLeftString        = "  "
   res1@gsnRightString       = "  "
   res1@tiMainOffsetYF       = 0.005
   ;res1@gsnRightStringOrthogonalPosF = 0.04
   ;res1@gsnLeftStringOrthogonalPosF = 0.04
  else
   res1@tiMainOffsetYF       = 0.0
   res1@gsnLeftString        = "Corr: "+precm@pc
   res1@gsnRightString       = "RMSE: "+precm@rmse
   ;res1@gsnRightStringOrthogonalPosF = 0.02
   ;res1@gsnLeftStringOrthogonalPosF = 0.02
  end if 

  if(igp.eq.0)then 
   res1@tmYLLabelsOn         = True
   res1@tmYLOn               = True
   res1@gsnRightStringOrthogonalPosF = 0.01
   res1@gsnLeftStringOrthogonalPosF  = 0.01
  else
   res1@tmYLLabelsOn         = False
   res1@tmYLOn               = False
   res1@gsnRightStringOrthogonalPosF = 0.02
   res1@gsnLeftStringOrthogonalPosF = 0.02
  end if 

  plot(igp) = gsn_csm_hov(wks1,precm,res1)

  delete([/precm/])

 end do 

;---Draw both plots in a panel
  pnres                         = True
  pnres@gsnMaximize             = False
  pnres@gsnPanelMainString      = "" ;"Hovmoller diagrams of "+Varname+" over " + Figstr 
  pnres@gsnPanelMainFontColor   = "Black"
  pnres@gsnPanelMainFontHeightF = 0.012
  pnres@gsnPanelLabelBar        = True
  pnres@lbOrientation           = "vertical"        ; vertical label bar
  pnres@pmLabelBarHeightF       = 0.24
  pnres@pmLabelBarWidthF        = 0.03
  pnres@lbLabelFontHeightF      = 0.006               ; make labels smaller
  pnres@pmLabelBarParallelPosF  = 0.005
 ;pnres@pmLabelBarOrthogonalPosF= 0.01
  pnres@lbTitleString           = Varname+ " ("+Varunt+")"
  pnres@lbTitlePosition         = "Right"                           ; title location
  pnres@lbTitleDirection        = "Across"                          ; letter angle
  pnres@lbTitleAngleF           = 90.                               ; title angle
  pnres@lbTitleFontHeightF      = 0.006                              ; font height
 ;pnres@gsnPanelBottom          = 0.05
  pnres@gsnPanelYWhiteSpacePercent = 2
  xpos = fspan(0.1,0.68,ngrps)
  pnres@gsnPanelXF = xpos ; Adjust rightmost plots
  gsn_panel(wks1,(/plot(:)/),(/1,ngrps/),pnres)

  system ("convert "+fignam+".png" + " -trim " +fignam+".pdf")
  system ("rm " + fignam+".png")

 end do ; end of var loop

 delete([/anmobs,refobs/])

end do ; end of subreg loop

end

