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

 ldatalocal    = True
 ldatalocal    = True
 plot_sig_test = True
 plot_contour  = True   ; if True, plot the contour in the figure
 datadir       = "../fig_data/"

;ERA5 -> RD
;CLIM -> CR
;Nudged -> NU
;LSTM(CLIM) -> ML(CR)

;;;;;;;strings for the data directory;;;;;;;;
 Obsnam   = "ERA5"
 Obsstr   = "RD"
 Groups   = (/"CLIM30", "CLIM30ML" /)
 Expnams  = (/    "CR",   "ML(CR)"   /)
 ngrps    = dimsizes(Groups)
 LabList  = Expnams + " ~F19~:~N~~F~ " + Obsstr 
 yst      = 1979
 yed      = 2014

;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
 VarList  = (/"U", "T", "Q"/)
 VstList  = (/"U", "T", "Q"/)
 VnmList  = (/"Zonal wind", "Temperature","Specific humidity"/)
 UnitList = (/"m s~S~-1~N~","K",           "g kg~S~-1~N~"    /)
 PltNums  = (/"a", "b",  "d", "e", "g", "h"/)
 PltNums  = "(" + PltNums+")"

 PltNumx  = (/"c","f","i"/)
 PltNumx  = "(" + PltNumx+")"

 nvars    = dimsizes(VarList)

 Seasons  = (/"ANN"/) 
 nseas    = dimsizes(Seasons) 

 do is = 0,nseas-1,1

  do iv = 0,nvars-1, 1

   Varname  = VarList(iv)
   Varunt   = UnitList(iv)
   Varstr   = VstList(iv)
   Vnmstr   = VnmList(iv)

   do igp = 0,ngrps-1,1
    model = Groups(igp)
    fna   = datadir+Obsnam+"_"+Varname+"_"+yst+"-"+yed+".nc"
    fnb   = datadir+model+"_"+Varname+"_"+yst+"-"+yed+".nc"
    fla   = addfile(fna, "r")
    flb   = addfile(fnb, "r")
 
    xtmp  = fla->$Varname$
    ytmp  = flb->$Varname$
    prob  = flb->prob 

    xtmp@_FillValue = -99999
    ytmp@_FillValue = -99999
    prob@_FillValue = -99999

    if(Varname.eq."RELHUM") then 
      xtmp = xtmp*100.0
      ytmp = ytmp*100.0
    end if 

    if(Varname.eq."Q") then
      xtmp = where(xtmp.le.0.005.or.ytmp.le.0.005,xtmp@_FillValue,xtmp)
      ytmp = where(xtmp.le.0.005.or.ytmp.le.0.005,ytmp@_FillValue,ytmp)
    end if

    vdiff     = xtmp
    vdiff     = ytmp - xtmp 
    vdiff&lev = ispan(1,dimsizes(vdiff&lev),1)
    prob&lev  = ispan(1,dimsizes(vdiff&lev),1)

   ;normalize figures with obs rms at each level
    ostddev = dim_stddev_n_Wrap(xtmp,1)
    ostddev&lev = ispan(1,dimsizes(vdiff&lev),1)
    vdiff   = vdiff/conform(vdiff,ostddev,(/0/))
     
    delete([/fna,fnb,xtmp,ytmp/])
   
;----------------------------------------------------------------------
 ; Plotting section
  if(igp.eq.0.and.iv.eq.0)then 
    wtype = "pdf"
    OUTDir = "./"
    outnam = OUTDir + "Fig_pres_lat_ndg_"+Seasons(is)
    wks1   = gsn_open_wks(wtype,outnam)
    plot   = new((/nvars,ngrps/),graphic)
    pltx   = new((/nvars,ngrps/),graphic)
    plox   = new((/nvars/),graphic)
  end if 

  load "./var_share_colorbar_nom.ncl"
  labels = (/Varstr ,"~F33~D~N~~F~"+ Varstr/)
  nlabls = dimsizes(labels)

;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@vpWidthF             = 0.8
  res@vpHeightF            = 0.6

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"

  res@trYMinF                 = 10.0                   ; min value on y-axis
  res@trYMaxF                 = 72.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = ispan(0,80,10) ;(/1000, 850, 700,500,400,300,200,100,1/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings
  res@trYLog                  = False
  res@tmYLOn                  = True
  res@trYReverse              = True                       ; reverse Y-axis

  FontHeightF = 0.032
  res@tiMainFontThicknessF    = 1.0
  res@tiMainFontHeightF       = FontHeightF*1.4
  res@tmYLLabelFontHeightF    = FontHeightF*1.1
  res@tmXBLabelFontHeightF    = FontHeightF*1.1
  res@tiXAxisFontHeightF      = FontHeightF*1.1
  res@tiYAxisFontHeightF      = FontHeightF*1.1
  res@gsnStringFontHeightF    = FontHeightF*1.1

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

  res0 = res

  res@cnFillOn                = True              ; color plot desired
 ;res@cnFillPalette           = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn               = False             ; turn off contour lines
  res@cnLineLabelsOn          = False             ; turn off contour labels
 ;res@cnFillMode              = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn           = False
  res@cnMissingValFillColor   = "Transparent"
  res@cnLineLabelFontHeightF  = FontHeightF
  res@lbLabelFontHeightF      = FontHeightF
  res@lbTitleFontHeightF      = FontHeightF                             ; font height

  res@lbLabelBarOn            = True            ; Will turn on in panel later
  res@lbAutoManage            = False           ; we control label bar
  res@lbLabelStride           = 1               ; skip every other label
  res@lbPerimOn               = False           ; default has box
  res@lbOrientation           = "vertical"      ; vertical label bar
  res@pmLabelBarHeightF       = 0.6
  res@pmLabelBarWidthF        = 0.12
; res@pmLabelBarParallelPosF  = 0.026
; res@pmLabelBarOrthogonalPosF= 0.15

  res@lbTitleString           = "" ;Varunt;Varstr+ " ("+Varunt+")";
  res@lbTitlePosition         = "Right"                           ; title location
  res@lbTitleDirection        = "Across"                          ; letter angle
  res@lbTitleAngleF           = 90.                               ; title angle
  res@gsnLeftString           = PltNums(iv*ngrps+igp)+ " " +LabList(igp) 
  res@gsnRightString          = "~F33~D~N~~F~"+Varstr+" ("+ Varunt +")" 
  res@gsnLeftStringOrthogonalPosF    = 0.02
 ;res@gsnLeftStringParallelPosF      = 0.2
  res@gsnRightStringOrthogonalPosF   = 0.02
  res@tiMainString                   = ""

  if(igp.eq.0) then
    res@tiYAxisString          = "Sigma level"
    res@tmYLLabelsOn           = True
  else
    res@tiYAxisString          = ""
    res@tmYLLabelsOn           = False
  end if

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  res2=res
  res2@cnLevelSelectionMode           = "ExplicitLevels"  
  res2@cnFillColors                  := colorgroup1
  res2@cnLevels                      := diffgrp1
  res2@cnLineThicknessF               = 0.5
  res2@cnLinesOn                      = plot_contour
  res2@cnLineLabelsOn                 = False            ; turn off contour labels
  res2@cnInfoLabelOn                  = False
  res2@gsnContourZeroLineThicknessF   = 0.5
  res2@gsnContourNegLineDashPattern   = 1
  if(igp.eq.1) then
     res2@lbLabelBarOn          = True            ; Will turn on in panel later
  else
     res2@lbLabelBarOn          = False           ; Will turn on in panel later
  end if
  plot(iv,igp) = gsn_csm_contour(wks1,vdiff,res2)

  if(plot_sig_test) then 
    res21 = res2
    res21@gsnLeftString                = "" 
    res21@gsnRightString               = ""
    res21@gsnCenterString              = ""
    res21@cnFillOn                     = True
    res21@cnLinesOn                    = False
    res21@cnLineLabelsOn               = False    
    res21@cnInfoLabelOn                = False
    res21@cnMonoFillColor              = True       
    res21@cnMonoFillPattern            = True 
    res21@cnMonoFillScale              = True 
    res21@cnFillPattern                = 17
    res21@cnFillScaleF                 = 0.8
    res21@cnFillDotSizeF               = 0.002 
    res21@lbLabelBarOn                 = False    
    pltx(iv,igp) = gsn_csm_contour(wks1,prob,res21)
    overlay(plot(iv,igp),pltx(iv,igp))
    delete([/res21/])
  end if 

  delete([/res/])

 end do

 ;plot rms for observations
 res0@vpWidthF           = 0.42
 res0@vpHeightF          = 0.6
 res0@trYMinF            = 10.0                   ; min value on y-axis
 res0@trYMaxF            = 72.0                   ; max value on y-axis
 res0@tmYLMode           = "Explicit"                ; explicit labels
 res0@tmYLValues         = ispan(0,80,10) ;(/1000, 850, 700,500,400,300,200,100,1/)
 res0@tmYLLabels         = ""+res0@tmYLValues         ; make strings
 res0@trYLog             = False
 res0@tmYLOn             = True 
 res0@tiYAxisString      = "" ;"Pressure (hPa)" ;""
 res0@tiXAxisString      = ""
 res0@tmYLLabelsOn       = False
 res0@tmXBMinorPerMajor  = 1
 res0@trXMinF            = floor(min(ostddev))                   ; min value on y-axis
 res0@trXMaxF            = ceil(max(ostddev))                   ; max value on y-axis
 res0@tmXBMode           = "Explicit"                ; explicit labels
 xmin = toint(res0@trXMinF)
 xmax = toint(res0@trXMaxF)
 xint = (xmax - xmin) / 4
 res0@tmXBValues        := ispan(xmin,xmax,xint) 
 res0@tmXBLabels        := ""+res0@tmXBValues         ; make strings
 res0@gsnLeftString                  = PltNumx(iv)+ " "+Obsstr+ " (~F33~s~N~~F~~B~"+Varstr+"~N~, "+ Varunt +")"
 res0@gsnRightString                 = "" ;Obsnam
 res0@gsnLeftStringOrthogonalPosF    = 0.02
;res0@gsnLeftStringParallelPosF      = 0.2
 res0@gsnRightStringOrthogonalPosF   = 0.02
 res0@tiMainString                   = ""

 plox(iv) = gsn_csm_xy (wks1,ostddev,vdiff&lev,res0)           ; create plot

 end do ;end of variable loop

 res_R                            = True
 res_R@gsnFrame                   = False
 res_R@gsnPanelLabelBar           = False
 res_R@lbAutoManage               = False             ; we control label bar
 res_R@lbLabelStride              = 1                 ; skip every other label
 res_R@lbLabelFontHeightF         = FontHeightF*0.30
 res_R@lbPerimOn                  = False             ; default has box
;res_R@lbOrientation              = "vertical"        ; vertical label bar
 res_R@pmLabelBarHeightF          = 0.04
 res_R@pmLabelBarWidthF           = 0.50
 res_R@gsnPanelMainString         = ""
 res_R@gsnPanelMainFont           = ""
 res_R@pmLabelBarParallelPosF     = 0.025

; plot_height = 1./tofloat(ngrps)
; do igp = 0,ngrps-1,1
;   res_R@gsnPanelTop               = 1.0 - igp * 0.25
;   res_R@gsnPanelBottom            = res_R@gsnPanelTop - plot_height
;   gsn_panel(wks1,plot(:,igp),(/1,nvars/),res_R)
; end do

 plot_height = 1./tofloat(nvars)
 do iv = 0,nvars-1,1
   res_R@gsnPanelTop       = 1.0 - iv * 0.27
   res_R@gsnPanelBottom    = res_R@gsnPanelTop - plot_height
   res_R@gsnPanelXF        = (/0.2,0.50,0.85/)-0.12 ; Adjust rightmost plots
   dumx = new(ngrps+1,graphic)
   do igp = 0,ngrps-1
     dumx(igp) = plot(iv,igp)
   end do 
   dumx(ngrps) = plox(iv)
   gsn_panel(wks1,(/dumx/),(/1,ngrps+1/),res_R)
 end do

 frame(wks1)

 end do ;end of file loop

end
