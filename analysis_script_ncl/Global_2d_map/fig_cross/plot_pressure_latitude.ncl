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

  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/state_var"
  obs_name   = "ERA5"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP"   , "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"
  ystr       = 2001
  yend       = 2014
  stat_slev  = 0.1
  l_plot_ttest = True

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11
  d2r        = get_d2r("float")

  varList    = (/"U","T"/)
  varFact    = (/1.0,1.0/)
  varUnit    = (/"m s~S~-1~N~","~S~o~N~C"/)
  nvars      = dimsizes(varList)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =   0.
  lonW       =   0.
  lonE       = 360.

  do j = 0,nvars-1,1

  vars   = varList(j)
  vunt   = varUnit(j)
  print("working on " + vars)

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do k = 0,nexps-1

    fl0  = systemfunc("ls -1 "+data_dir+"/"+obs_name+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
    fm0  = addfile(flmsk,"r")
    f0   = addfiles(fl0,"r")
    f1   = addfiles(fl1,"r")
    ListSetType (f0, "join")
    ListSetType (f1, "join")
 
    tim0 = f0[0]->time
    cd_date       = cd_calendar(tim0, 0)
    year          = tointeger(cd_date(:,0))
    indy          = ind(year.ge.ystr.and.year.le.yend)
   ;obstime       = f0[0]->time(indy)
   ;obstime@units = f0[0]->time@units
   ;print(cd_calendar(obstime,0))

    vnm  = vars+"_o"
    vd0  = dim_avg_n_Wrap(f0[:]->$vnm$(:,indy,:,{latS:latN},{lonW:lonE}),0)
    vd1  = dim_avg_n_Wrap(f1[:]->$vnm$(:,indy,:,{latS:latN},{lonW:lonE}),0)
    lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})

    vd0@_FillValue = -9999.0
    vd1@_FillValue = -9999.0
    printVarSummary(lmsk)
    printVarSummary(vd0)
    lms0 = conform(vd0,lmsk,(/2,3/))
    vd0  = where(lms0.le.0,vd0,vd0@_FillValue)
    vd1  = where(lms0.le.0,vd1,vd1@_FillValue)
    
    vdm0 = dim_avg_n_Wrap(vd0,0)
    vdm1 = dim_avg_n_Wrap(vd1,0)
    rad  = 4.0*atan(1.0)/180.0
    wgty = cos(f1[0]->lat({latS:latN})*rad)
    glm1 = wgt_areaave(vdm0,wgty,1,0)
    glm2 = wgt_areaave(vdm1,wgty,1,0)

    vpt1 = dim_avg_n_Wrap(vd0,3)
    vpt2 = dim_avg_n_Wrap(vd1,3)

    if(l_plot_ttest) then 
      xtmp  = vpt1(lev|:,lat|:,time|:)       ; reorder but do it only once [temporary]
      ytmp  = vpt2(lev|:,lat|:,time|:)
      xAve  = dim_avg (xtmp)              ; calculate means at each grid point
      yAve  = dim_avg (ytmp)
      xVar  = dim_variance (xtmp)         ; calculate variances
      yVar  = dim_variance (ytmp)
      xEqv  = equiv_sample_size (xtmp, stat_slev,0)
      yEqv  = equiv_sample_size (ytmp, stat_slev,0)
      xN    = dim_avg_n(dim_avg_wgt_n(xEqv, wgty, 0, 1),0)
      yN    = dim_avg_n(dim_avg_wgt_n(yEqv, wgty, 0, 1),0)
      iflag = True                       ; population variance similar
      vprob = ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)
      delete([/xtmp,ytmp,xAve,yAve,xVar,yVar,xEqv,yEqv,xN,yN/])
    end if 

    vor1 = dim_avg_n_Wrap(vpt1,0)
    vor2 = dim_avg_n_Wrap(vpt2,0)
    vdat = vor1
    vdat = vor2 - vor1
    delete([/vpt1,vpt2,vd0,vd1,lmsk,lms0/])
    delete([/fl0,fl1,fm0,f0,f1/])
  ;----------------------------------------------------------------------
  ; Plotting section
  
   if (k.eq.0 ) then
     wtype           = "png"
     wtype@wkWidth   = 2500
     wtype@wkHeight  = 2500
     fignam          = "fig_clim_bias_U_T_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     plot      = new((/nexps/),graphic)
     dum       = new((/nexps/),graphic)
   end if 

  load "./var_share_colorbar_preshgt.ncl"

  ;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False

  res@vpWidthF             = 0.9
  res@vpHeightF            = 0.6
  res@gsnMajorLatSpacing   =  10             ;-- change major lat tickmark spacing
  res@gsnMinorLatSpacing   =  2.5            ;-- change major lat tickmark spacing

  FontHeightF = 0.032
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

   BorderThick = 3.0
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "White"

  res@lbLabelBarOn            = True      ; Will turn on in panel later
  res@lbAutoManage            = False             ; we control label bar
 ;res@lbLabelStride           = 1                 ; skip every other label
  res@lbLabelFontHeightF      = FontHeightF*0.9
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
  res@tmYLOn                  = True
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"
  res@tmYLOn                  = True

  res@trYMinF                 = 10.0                   ; min value on y-axis
  res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/1000, 700, 500, 300, 200, 100, 50, 30, 20, 10/) ;(/ 1000, 850, 700,500,400,300,200,100,10/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings

  res@trYLog                  = True

  res@gsnLeftString           = "" 
  res@gsnRightString          = "" 
  res@tiMainString            = ""

 ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  res1=res
  res1@cnLevelSelectionMode           = "ExplicitLevels"  
  res1@cnLinesOn                      = True
  res1@cnLineLabelsOn                 = False            ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@gsnContourZeroLineThicknessF   = 0.2
  res1@gsnContourNegLineDashPattern   = 0
  res1@cnLineThicknessF               = 0.2
  res1@lbTitleString                  = "" ;"("+ Varunt +")"

  if(k.eq.0) then 
    res1@cnFillPalette                = "MPL_rainbow" ;"cmocean_balance" ;"BlueDarkRed18"
    res1@cnLevels                     = (/-65,-60,-55,-50,-45,-40,-30,-20,-10,0,10,20/);fspan(-60.0,20.0,21)
    res1@gsnLeftString                = exp_strs(0)
    res1@gsnRightString               = "Zonal wind"
    plot(0) = gsn_csm_pres_hgt(wks,vor1,res1)
  end if 

  res1@cnLevelSelectionMode           = "ExplicitLevels"
  res1@cnFillPalette                  = "cmp_b2r" ;"cmocean_balance" ;"BlueDarkRed18"
  res1@cnLevels                      := fspan(-5.0,5.0,11)
  res1@gsnLeftString                  = exp_strs(k)
  res1@gsnRightString                 = "Model biases"
  plot(k) = gsn_csm_pres_hgt(wks,vdat,res1)


  res2 = res
  res2@gsnLeftString                  = ""
  res2@gsnRightString                 = ""
  res2@gsnCenterString                = ""
  res2@cnFillOn                       = False
  res2@cnLinesOn                      = True
  res2@cnLineLabelsOn                 = True            ; turn off contour labels
  res2@cnLevelSelectionMode           = "ExplicitLevels"
  res2@cnInfoLabelOn                  = False
  res2@cnLineLabelDensityF            = 2.0
  res2@cnLineLabelFontHeightF         = FontHeightF
  res2@cnLineLabelBackgroundColor     = -1
  res2@gsnContourLineThicknessesScale = 1.0
  res2@cnLineThicknessF               = 6.0
  res2@cnLineColor                    = "Black" ;Magenta"
  res2@cnLineLabelFontColor           = "Black"
  res2@gsnContourZeroLineThicknessF   = 12.0
  res2@gsnContourNegLineDashPattern   = 1

  if(k.eq.0) then
    res2@cnLevels                     = fspan(-30.0,30.0,21)
    res2@gsnLeftString                = ""
    res2@gsnRightString               = ""
    plot1 = gsn_csm_pres_hgt(wks,vor1,res2)
    overlay(plot(0),plot1)
  end if

  res2@cnLevels                      := fspan(-5.0,5.0,11)
  res2@gsnLeftString                  = ""
  res2@gsnRightString                 = ""
  plot2 = gsn_csm_pres_hgt(wks,vdat,res2)
  overlay(plot(k),plot2)

  delete([/res,res1,res2/])
  delete([/vdat,vor1,vor2/])
 end do ;end of exp loop

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
 pnres@gsnPanelYWhiteSpacePercent = 2
 pnres@gsnPanelXWhiteSpacePercent = 1
 gsn_panel(wks,(/plot(:)/),(/nexps,1/),pnres)

end do 
end do 

end
