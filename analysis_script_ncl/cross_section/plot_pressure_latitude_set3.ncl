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

 ;;;;;;;strings for the data directory;;;;;;;;
  Obsnam   = "ERA5"
  Groups   = (/"CLIM30","NDGUVTQ","CLIM30ML"/)
  Expnams  = Groups
  ngrps    = dimsizes(Groups)
  LabList  = Expnams + " ~F19~:~N~~F~ " + Obsnam

  syear    = 2010
  eyear    = 2014
 
;;;;;;variables to be plotted;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  VarList  = (/"U", "T", "Q", "RELHUM"/)
  VstList  = (/"U", "T", "Q", "RH"/)
  VnmList  = (/"Zonal wind", "Temperature", "Cloud fraction", "Relative humidity"/)
  UnitList = (/"m s~S~-1~N~","K",           "g/kg",           "%"                /)
  PltNums  = (/"a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l"/)
  PltNums  = "(" + PltNums+")"
  nvars    = dimsizes(VarList)

  Seasons  = (/"ANN","DJF","JJA"/) 
  nseas    = dimsizes(Seasons) 

  do is = 0,nseas-1,1

   do iv = 0,nvars-1, 1

    Varname  = VarList(iv)
    Varunt   = UnitList(iv)
    Varstr   = VstList(iv)
    Vnmstr   = VnmList(iv)
 
    do igp = 0,ngrps-1,1

    fna = new(5,string)
    do year = syear,eyear,1
      fna(year-syear) = systemfunc("ls " + datadir+ Varname+ "_" + Obsnam + "_"+year+".nc")
    end do 

    if(Groups(igp).eq."NDGUVTQ") then 
     fnb = systemfunc("ls " + datadir+ Varname + "_" + Groups(igp) + ".nc")
    else 
      fnb = new(5,string)
      do year = syear,eyear,1
        fnb(year-syear) = systemfunc("ls " + datadir+ Varname + "_" + Groups(igp) + "_"+year+".nc")
      end do 
    end if 

    fla   = addfiles(fna, "r")
    flb   = addfiles(fnb, "r")
    utco  = cd_calendar(fla[:]->time,0)
    utcm  = cd_calendar(flb[:]->time,0)
    indo  = ind(utco(:,0).ge.syear.and.utco(:,0).le.eyear)
    indm  = ind(utcm(:,0).ge.syear.and.utcm(:,0).le.eyear)
    delete([/utco,utcm/])

    Varbias  = "d"+Varname
    if(Seasons(is).eq."ANN") then 
      xtmp  = month_to_annual(fla[:]->$Varname$(indo,:,:),1)
      ytmp  = month_to_annual(flb[:]->$Varname$(indm,:,:),1)
      ztmp  = month_to_annual(flb[:]->$Varbias$(indm,:,:),1)
    else
      xtmp  = month_to_season(fla[:]->$Varname$(indo,:,:),Seasons(is))
      ytmp  = month_to_season(flb[:]->$Varname$(indm,:,:),Seasons(is))
      ztmp  = month_to_season(flb[:]->$Varbias$(indm,:,:),Seasons(is))
    end if

    xtmp!0 = "time"
    ytmp!0 = "time"
    ztmp!0 = "time"
    xtmp@_FillValue = -99999
    ytmp@_FillValue = -99999
    ztmp@_FillValue = -99999
 
    if(Varname.eq."RELHUM") then 
      xtmp = xtmp*100.0
      ytmp = ytmp*100.0
      ztmp = ztmp*100.0
    end if 
 
    printMinMax(xtmp,False)
    printMinMax(ytmp,False)

    ;calculate the global/regional mean values
    latrg = flb[0]->lat
    rad   = 4.*atan(1.)/180.    ;
    wgty  = cos(latrg*rad)      ; cosine weights

    xAve  = dim_avg_Wrap (xtmp(lev|:,lat|:,time|:))     
    yAve  = dim_avg_Wrap (ytmp(lev|:,lat|:,time|:))

    xVar  = dim_variance_Wrap (xtmp(lev|:,lat|:,time|:)) ; calculate variances
    yVar  = dim_variance_Wrap (ytmp(lev|:,lat|:,time|:))

    sigr  = 0.05                                    ; critical sig lvl for r
    xEqv  = equiv_sample_size (xtmp(lev|:,lat|:,time|:), sigr,0)
    yEqv  = equiv_sample_size (ytmp(lev|:,lat|:,time|:), sigr,0)
    xN    = wgt_areaave (xEqv, 1.0, wgty, 0)    ; wgty could be gaussian weights 
    yN    = wgt_areaave (yEqv, 1.0, wgty, 0) 

    iflag = False                        ; population variance similar
    prob  = xAve
    prob  = (/ttest(xAve,xVar,xN, yAve,yVar,yN, iflag, False)/)
    prob  = where(prob.lt.sigr,prob,prob@_FillValue )
    vdiff = dim_avg_Wrap (ztmp(lev|:,lat|:,time|:))  ;xAve
    delete([/fna,fnb,xtmp,ytmp,ztmp,xAve,yAve,xVar,yVar,sigr,xN,yN,iflag,indo,indm/])
    
;----------------------------------------------------------------------
; Plotting section
  if(igp.eq.0.and.iv.eq.0)then 
    wtype = "pdf"
    OUTDir = "./"
    outnam = OUTDir + "Fig_pres_lat_"+Seasons(is)
    wks1   = gsn_open_wks(wtype,outnam)
    plot   = new((/nvars,ngrps/),graphic)
  end if 

  load "./var_share_colorbar2.ncl"
  labels = (/Varstr ,"~F33~D~N~~F~"+ Varstr/)
  nlabls = dimsizes(labels)

;---Resources to share between both plots
  res                      = True              ; Plot modes desired.
  res@gsnDraw              = False
  res@gsnFrame             = False
  res@gsnMaximize          = False              ; Maximize plot
  res@gsnAddCyclic         = False
  res@vpWidthF             = 0.8
  res@vpHeightF            = 0.6

  res@cnFillOn              = True              ; color plot desired
 ;res@cnFillPalette         = "ViBlGrWhYeOrRe";"testcmap" ;"NCV_jaisnd";"temp_19lev"
  res@cnLinesOn             = False             ; turn off contour lines
  res@cnLineLabelsOn        = False             ; turn off contour labels
 ;res@cnFillMode            = "RasterFill"      ; turn raster on
  res@cnInfoLabelOn         = False
  res@cnMissingValFillColor = "Transparent"

  res@lbLabelBarOn            = True            ; Will turn on in panel later
  res@lbAutoManage            = False           ; we control label bar
  res@lbLabelStride           = 1               ; skip every other label
  res@lbLabelFontHeightF      = 0.028
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

  res@tmXTOn                  = False
  res@tmYROn                  = False
  res@tmYRMode                = "Automatic"
  res@tiYAxisString           = "Pressure (hPa)"

  res@trYMinF                 = 100.0                   ; min value on y-axis
  res@trYMaxF                 = 1000.0                   ; max value on y-axis
  res@tmYLMode                = "Explicit"                ; explicit labels
  res@tmYLValues              = (/1000, 850, 700,500,400,300,200,100/)
  res@tmYLLabels              = ""+res@tmYLValues         ; make strings
  res@trYLog                  = False

  FontHeightF = 0.032
  res@tiMainFontThicknessF    = 1.0
  res@tiMainFontHeightF       = FontHeightF*1.4
  res@tmYLLabelFontHeightF    = FontHeightF*1.1
  res@tmXBLabelFontHeightF    = FontHeightF*1.1
  res@tiXAxisFontHeightF      = FontHeightF*1.1
  res@tiYAxisFontHeightF      = FontHeightF*1.1
  res@gsnStringFontHeightF    = FontHeightF*1.1
  res@cnLineLabelFontHeightF  = FontHeightF
  res@lbLabelFontHeightF      = FontHeightF
  res@lbTitleFontHeightF      = FontHeightF                             ; font height

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

  res@gsnLeftString                  = PltNums(iv*ngrps+igp)+ " " +LabList(igp) 
  res@gsnRightString                 = "~F33~D~N~~F~"+Varstr+" ("+ Varunt +")" ;" ;" ("+Varunt+")";""
  res@gsnLeftStringOrthogonalPosF    = 0.02
  ;res@gsnLeftStringParallelPosF      = 0.2
  res@gsnRightStringOrthogonalPosF   = 0.02
  res@tiMainString                   = ""
  res@tmYLLabelsOn                   = False 
  res@tiYAxisString                  = ""

  if(igp.gt.0) then
     res@tmYLOn                 = True
     res@tiYAxisString          = ""
     res@tmYLLabelsOn           = True
  else
     res@tmYLOn                 = True
     res@tmYLLabelsOn           = True
     res@tiYAxisString          = "Pressure (hPa)"
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

  plot(iv,igp) = gsn_csm_pres_hgt(wks1,vdiff,res2)

  if(plot_sig_test) then 
    res21 = res2
    res21@gsnLeftString                = "" 
    res21@gsnRightString               = ""
    res21@gsnCenterString              = ""
    res21@cnFillOn                     = True
    res21@cnLinesOn                    = False
    res21@cnLineLabelsOn               = False            ; turn off contour labels
    res21@cnInfoLabelOn                = False
    res21@cnMonoFillColor              = True                 ; default color is fground [black]
    res21@cnMonoFillPattern            = True 
    res21@cnMonoFillScale              = True 
    res21@cnFillPattern                = 17
    res21@cnFillScaleF                 = 0.8
    res21@cnFillDotSizeF               = 0.002 
    plot_cont1 = gsn_csm_pres_hgt(wks1,prob,res21)
    overlay(plot(iv,igp),plot_cont1)
    delete([/res21/])
  end if 

  delete([/res/])

 end do

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

 plot_height = 1./tofloat(ngrps)
 do igp = 0,ngrps-1,1
   res_R@gsnPanelTop               = 1.0 - igp * 0.25
   res_R@gsnPanelBottom            = res_R@gsnPanelTop - plot_height
   gsn_panel(wks1,plot(:,igp),(/1,nvars/),res_R)
 end do

 frame(wks1)

 end do ;end of file loop

end
