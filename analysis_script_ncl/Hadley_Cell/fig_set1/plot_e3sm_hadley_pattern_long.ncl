  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
  load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; Purpose: A tempolate to make a eof plot for AAO patterns 
;          Data is vertically interpolated to pressure levels
;          1st panel: EOFs for the first three modes
;          2nd panel: time series plot
; Author : Shixuan Zhang
; Date   : Dec 06, 2020

begin

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/EOF_analysis/v2_4CO2/data"
  exp_name  = (/"v2.LR.piControl", \
                "v2.LR.piControl.CRYOnoismf", \
                "v2.LR.piControl.CRYO", \
                "v2.LR.abrupt-4xCO2", \
                "v2.LR.CRYOnoismf.abrupt-4xCO2", \
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(exp_name)

  main_str  = (/"v2.LR.piControl            ", "v2.LR.abrupt-4xCO2           ", "4xCO2 - piControl (WC)             ", \
                "v2.LR.piControl.CRYOnoismf ", "v2.LR.CRYOnoismf.abrupt-4xCO2", "4xCO2 - piControl (CRYOnoismf)     ", \
                "v2.LR.piControl.CRYO       ", "v2.LR.abrupt-4xCO2.CRYO      ", "4xCO2 - piControl (CRYO)           ", \
                "CRYOnoismf - WC (piControl)", "CRYOnoismf - WC (4xCO2)      ", "CRYOnoismf - WC (4xCO2 - piControl)", \
                "CRYO - WC (piControl)      ", "CRYO - WC (4xCO2)            ", "CRYO - WC (4xCO2 - piControl)      " /)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)

  facsfm    = 1e-9
  facwm     = 100.0 

  latS      = -90.0
  latN      =  90.0
 
  time_tag  = (/"0051-0150"/)
  time_str  = (/"Year 51-150"/)
  time_prd  = (/"0051-0100", "0101-0150"/)
  nprd      = dimsizes(time_prd)

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    do k = 0,nexps-1,1

     fls = new(nprd,string)
     do ii = 0,nprd-1,1
       fls(ii)  = systemfunc ("ls "+data_dir+"/"+exp_name(k)+"_HC_analysis_Monthly_"+time_prd(ii)+".nc")
     end do
     f    = addfiles(fls,"r")
   ;print(fls)
    delete(fls)

    rad      = 4.0*atan(1.0)/180.0
    wgty     = cos(f[0]->lat*rad)

    if (k.eq.0) then 
      ud1    = dim_avg_n_Wrap(f[:]->U,0)
      vd1    = dim_avg_n_Wrap(f[:]->V,0)
      wd1    = dim_avg_n_Wrap(f[:]->W,0)
      sfd1   = dim_avg_n_Wrap(f[:]->SF,0)
      sfd1   = sfd1 * facsfm
      wd1    = wd1 * facwm
    end if 

    if (k.eq.1) then
      ud2    = dim_avg_n_Wrap(f[:]->U,0)
      vd2    = dim_avg_n_Wrap(f[:]->V,0)
      wd2    = dim_avg_n_Wrap(f[:]->W,0)
      sfd2   = dim_avg_n_Wrap(f[:]->SF,0)
      sfd2   = sfd2 * facsfm
      wd2    = wd2 * facwm
    end if

    if (k.eq.2) then
      ud3    = dim_avg_n_Wrap(f[:]->U,0)
      vd3    = dim_avg_n_Wrap(f[:]->V,0)
      wd3    = dim_avg_n_Wrap(f[:]->W,0)
      sfd3   = dim_avg_n_Wrap(f[:]->SF,0)
      sfd3   = sfd3 * facsfm
      wd3    = wd3 * facwm
    end if

    if (k.eq.3) then
      up1    = dim_avg_n_Wrap(f[:]->U,0)
      vp1    = dim_avg_n_Wrap(f[:]->V,0)
      wp1    = dim_avg_n_Wrap(f[:]->W,0)
      sfp1   = dim_avg_n_Wrap(f[:]->SF,0)
      sfp1   = sfp1 * facsfm
      wp1    = wp1 * facwm
    end if

    if (k.eq.4) then
      up2    = dim_avg_n_Wrap(f[:]->U,0)
      vp2    = dim_avg_n_Wrap(f[:]->V,0)
      wp2    = dim_avg_n_Wrap(f[:]->W,0)
      sfp2   = dim_avg_n_Wrap(f[:]->SF,0)
      sfp2   = sfp2 * facsfm
      wp2    = wp2 * facwm
    end if

    if (k.eq.5) then
      up3    = dim_avg_n_Wrap(f[:]->U,0)
      vp3    = dim_avg_n_Wrap(f[:]->V,0)
      wp3    = dim_avg_n_Wrap(f[:]->W,0)
      sfp3   = dim_avg_n_Wrap(f[:]->SF,0)
      sfp3   = sfp3 * facsfm
      wp3    = wp3 * facwm
    end if

    end do  

    ;;difference due to 4CO2
    udp1     = ud1
    udp1     = up1 - ud1
    vdp1     = vd1
    vdp1     = vp1 - vd1
    wdp1     = wd1
    wdp1     = wp1 - wd1
    sfdp1    = sfd1
    sfdp1    = sfp1 - sfd1
    glCordp1 = pattern_cor(sfp1,sfd1,1.0,0)
    glAvedp1 = wgt_arearmse(sfp1,sfd1, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp1)
    mstrdp1  = mstrx1 + "~C~" + mstrx2

    udp2     = ud2
    udp2     = up2 - ud2
    vdp2     = vd2
    vdp2     = vp2 - vd2
    wdp2     = wd2
    wdp2     = wp2 - wd2
    sfdp2    = sfd2
    sfdp2    = sfp2 - sfd2
    glCordp2 = pattern_cor(sfp2,sfd2,1.0,0)
    glAvedp2 = wgt_arearmse(sfp2,sfd2, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp2)
    mstrdp2  = mstrx1 + "~C~" + mstrx2

    udp3     = ud3
    udp3     = up3 - ud3
    vdp3     = vd3
    vdp3     = vp3 - vd3
    wdp3     = wd2
    wdp3     = wp3 - wd3
    sfdp3    = sfd3
    sfdp3    = sfp3 - sfd3
    glCordp3 = pattern_cor(sfp3,sfd3,1.0,0)
    glAvedp3 = wgt_arearmse(sfp3,sfd3, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp3)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp3)
    mstrdp3  = mstrx1 + "~C~" + mstrx2

    ;;difference due to configuration
    udd1     = ud1
    udd1     = ud2 - ud1
    vdd1     = vd1
    vdd1     = vd2 - vd1
    wdd1     = wd1
    wdd1     = wd2 - wd1
    sfdd1    = sfd1
    sfdd1    = sfd2 - sfd1
    glCordd1 = pattern_cor(sfd1,sfd2,1.0,0)
    glAvedd1 = wgt_arearmse(sfd1,sfd2, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd1)
    mstrdd1  = mstrx1 + "~C~" + mstrx2

    udd2     = ud1
    udd2     = ud3 - ud1
    vdd2     = vd1
    vdd2     = vd3 - vd1
    wdd2     = wd1
    wdd2     = wd3 - wd1
    sfdd2    = sfd1
    sfdd2    = sfd3 - sfd1
    glCordd2 = pattern_cor(sfd1,sfd3,1.0,0)
    glAvedd2 = wgt_arearmse(sfd1,sfd3, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd2)
    mstrdd2  = mstrx1 + "~C~" + mstrx2

    upp1     = up1
    upp1     = up2 - up1
    vpp1     = vp1
    vpp1     = vp2 - vp1
    wpp1     = wp1
    wpp1     = wp2 - wp1
    sfpp1    = sfp1
    sfpp1    = sfp2 - sfp1
    glCorpp1 = pattern_cor(sfp1,sfp2,1.0,0)
    glAvepp1 = wgt_arearmse(sfp1,sfp2, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorpp1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvepp1)
    mstrpp1  = mstrx1 + "~C~" + mstrx2
   
    upp2     = up1
    upp2     = up3 - up1
    vpp2     = vp1
    vpp2     = vp3 - vp1
    wpp2     = wp1
    wpp2     = wp3 - wp1
    sfpp2    = sfp1
    sfpp2    = sfp3 - sfp1
    glCorpp2 = pattern_cor(sfp1,sfp3,1.0,0)
    glAvepp2 = wgt_arearmse(sfp1,sfp3, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorpp2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvepp2)
    mstrpp2  = mstrx1 + "~C~" + mstrx2

    uppdd1    = upp1
    uppdd1    = upp1 - udd1
    vppdd1    = vpp1
    vppdd1    = vpp1 - vdd1
    wppdd1    = wpp1
    wppdd1    = wpp1 - wdd1
    sfppdd1   = sfpp1
    sfppdd1   = sfpp1 - sfdd1
    glCorppdd1= pattern_cor(sfpp1,sfdd1,1.0,0)
    glAveppdd1= wgt_arearmse(sfpp1,sfdd1, 1.0, 1.0,0)
    mstrx1    = "CORR. w/r " + sprintf("%.2f", glCorppdd1)
    mstrx2    = "RMSD  w/r " + sprintf("%.2f", glAveppdd1)
    mstrppdd1 = mstrx1 + "~C~" + mstrx2

    uppdd2    = upp2
    uppdd2    = upp2 - udd2
    vppdd2    = vpp2
    vppdd2    = vpp2 - vdd2
    wppdd2    = wpp2
    wppdd2    = wpp2 - wdd2
    sfppdd2   = sfpp2
    sfppdd2   = sfpp2 - sfdd2
    glCorppdd2= pattern_cor(sfpp2,sfdd2,1.0,0)
    glAveppdd2= wgt_arearmse(sfpp2,sfdd2, 1.0, 1.0,0)
    mstrx1    = "CORR. w/r " + sprintf("%.2f", glCorppdd2)
    mstrx2    = "RMSD  w/r " + sprintf("%.2f", glAveppdd2)
    mstrppdd2 = mstrx1 + "~C~" + mstrx2

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype             = "png"
    wtype@wkWidth     = 2500
    wtype@wkHeight    = 2500
    fignam            = "Hadley_circulation_"+seasons(i)+"_"+time_tag
    wks               = gsn_open_wks(wtype,fignam)

    plot              = new((/15/),graphic)
    tmpA              = new((/15/),graphic)
    tmpB              = new((/15/),graphic)

    leftstring        = "Zonal mean stream fuction (10~S~9~N~ kg s~S~-1~N~)"
    XBValues          = ispan(toint(latS),toint(latN),30)
    XBLabels          = (/"90S","60S","30S","0","30N","60N","90N"/)
    XBMinorValues     = ispan(toint(latS),toint(latN),15)
    YLValues          = (/ 1000, 850, 700, 500, 300, 200, 100, 10, 1/)

    Contour_color_bar = "BlueDarkRed18"      ;"MPL_PRGn"
    Diff_corlor_bar   = "amwg_blueyellowred" ;"BlueDarkOrange18" ;"GreenMagenta16"
    
    cnsflev           = (/-200,-100,-50,-30,-10,-1,0,1,10,30,50,100,200/)
    cnwdlev           = (/-50,-40,-30,-20,-10,-5,0,5,10,20,30,40,50/)
    cnsflevd          = (/-60,-40,-20,-10,-5,-1,0,1,5,10,20,40,60/)
    cnwdlevd          = (/-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10/)
    cnsfleva          = (/-10,-8,-6,-4,-2,-1,0,1,2,4,6,8,10/)

    if(seasons(i).eq."ANN") then 
      vcRefMag   = 1.0
      vcRefMagd1 = 0.5
      vcRefMagd2 = 0.2       
    else if (seasons(i).eq."DJF") then 
      vcRefMag   = 2.0   
      vcRefMagd1 = 1.0
      vcRefMagd2 = 0.2  
    else 
      vcRefMag   = 2.0   
      vcRefMagd1 = 1.0
      vcRefMagd2 = 0.2 
    end if 
    end if 

    ;============================================================
    ; PLOTS
    ;============================================================
     res                          = True         
     res@gsnDraw                  = False        ; don't draw yet
     res@gsnFrame                 = False        ; don't advance frame yet
     res@vpWidthF                 = 1.0
     res@vpHeightF                = 0.5
 
     FontHeightF = 0.028
     res@tiMainFontThicknessF     = 3.0
     res@tiMainFontHeightF        = FontHeightF*1.2
     res@tmYLLabelFontHeightF     = FontHeightF
     res@tmXBLabelFontHeightF     = FontHeightF
     res@tiXAxisFontHeightF       = FontHeightF
     res@tiYAxisFontHeightF       = FontHeightF
     res@gsnStringFontHeightF     = FontHeightF
     res@cnLineLabelFontHeightF   = FontHeightF
     res@lbLabelFontHeightF       = FontHeightF*0.8
     res@lbTitleFontHeightF       = FontHeightF*0.8

     res@cnFillOn                 = True         ; turn on color fill
     res@cnLinesOn                = False        ; True is default
     res@cnLineLabelsOn           = False        ; True is default
     res@cnMissingValFillColor    = "White" ;"gray76"               ; black = 1

     if(wtype.eq."png") then
      lthick = 3.0
     else
      lthick = 1.0
     end if
  
     res@tmBorderThicknessF       = lthick
     res@tmXBMajorThicknessF      = lthick
     res@tmXBMinorThicknessF      = lthick
     res@tmYLMajorThicknessF      = lthick
     res@tmYLMinorThicknessF      = lthick
     res@tmYRMajorThicknessF      = lthick
     res@tmYRMinorThicknessF      = lthick

     res@tmXTOn                   = False
     res@tmYROn                   = False
     res@tmXBMode                 = "Explicit"
     res@tmXBValues               = XBValues
     res@tmXBLabels               = XBLabels
     res@tmXBMinorValues          = XBMinorValues
     res@tmYLMode                 = "Explicit"                
     res@tmYLValues               = YLValues
     res@tmYLLabels               = ""+res@tmYLValues         

     res@cnInfoLabelOn                = False
     res@trYReverse                   = True
     res@tiYAxisString                = "Pressure (hPa)"
     res@gsnLeftStringOrthogonalPosF  = 0.02
     res@gsnRightStringOrthogonalPosF = 0.02

    ;*******************************************
    ; contour plot for zonal wind 
    ;*******************************************
     res0                              = res
     res0@cnFillOn                     = False                  ; turn on color fill
     res0@cnLinesOn                    = True                 ; no contour lines
     res0@cnLineLabelsOn               = True
     res0@cnMissingValFillColor        = "White" ;"gray76"               ; black = 1
     res0@cnLineColor                  = "maroon1";"orangered";"Magenta";"Blue3"
     res0@cnLineThicknessF             = lthick
     res0@cnInfoLabelOn                = False
     res0@cnLevelSelectionMode         = "ExplicitLevels"       ; manual levels
     res0@cnLineLabelInterval          = 1             ; default = 2
     res0@cnLineLabelFontColor         = res0@cnLineColor
     res0@cnLineLabelPlacementMode     = "constant"    ; choose constant label method
     res0@tiMainString                 = ""
     res0@gsnLeftString                = ""
     res0@gsnRightString               = ""

    ;*******************************************
    ; vector plot for (u,w)
    ;*******************************************
     vecres                           = True            ; vector only resources
     vecres@gsnDraw                   = False           ; don't draw
     vecres@gsnFrame                  = False           ; don't advance frame
     vecres@vcGlyphStyle              = "LineArrow"     ;"CurlyVector"   ; curly vectors
     vecres@vcRefLengthF              = 0.06            ; define length of vec ref
     vecres@vcMinDistanceF            = 0.015
     vecres@vcLineArrowThicknessF     = lthick
     vecres@gsnRightString            = " "             ; turn off right string
     vecres@gsnLeftString             = " "             ; turn off left string
     vecres@tiXAxisString             = " "             ; turn off axis label
     vecres@vcMapDirection            = False
     vecres@vcRefAnnoOrthogonalPosF   = -0.40    ; move ref vector into plot
     vecres@vcRefAnnoPerimOn          = False
     vecres@vcRefAnnoParallelPosF     = 0.99    ; move ref vector into plot
     vecres@vcRefAnnoFontHeightF      = FontHeightF*0.5
     vecres@vcRefAnnoSide             = "Top"
     vecres@vcRefMagnitudeF           = vcRefMag         ; define vector ref mag
     vecres@vcRefAnnoString1          = vecres@vcRefMagnitudeF
     vecres@vcRefAnnoString2          = "(v,w*100)"

    ;*******************************************
    ; shading plot for stream function 
    ;*******************************************
     res1 = res
     res1@cnLinesOn                    = True                 ; no contour lines
     res1@cnLineLabelsOn               = False
     res1@cnLineColor                  = "black"
     res1@cnLineThicknessF             = 0.2
     res1@lbLabelBarOn                 = True
     res1@lbOrientation                = "Vertical"
     res1@lbTitleString                = ""
     res1@lbTitlePosition              = "Right"     ; title location
     res1@lbTitleDirection             = "Across"    ; letter angle
     res1@lbTitleAngleF                = 90.         ; title angle
     res1@lbBoxLineColor               = "black"
    ;res1@lbBoxLinesOn                 = False
    ;res1@pmLabelBarParallelPosF       = 0.61
     res1@pmLabelBarOrthogonalPosF     = 0.02
     res1@pmLabelBarWidthF             = 0.15
     res1@pmLabelBarHeightF            = 0.52
     res1@cnInfoLabelOn                = False
     res1@cnLevelSelectionMode         = "ExplicitLevels"       ; manual levels
     res1@gsnLeftString                = leftstring
     res1@tiMainOffsetYF               = 0.0

     res0@cnLevels                    := cnwdlev
     tmpA(0)                           = gsn_csm_contour(wks,ud1,res0)
     vecres@vcRefMagnitudeF            = vcRefMag    ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(0)                           = gsn_csm_vector(wks,vd1,wd1,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfd1(:,:)))
     res1@tiMainString                 = main_str(0)
     res1@cnLevels                    := cnsflev
     plot(0)                           = gsn_csm_contour(wks,sfd1,res1)
     overlay(plot(0),tmpA(0))           
     overlay(plot(0),tmpB(0))      

     res0@cnLevels                    := cnwdlev
     tmpA(1)                           = gsn_csm_contour(wks,up1,res0)
     vecres@vcRefMagnitudeF            = vcRefMag              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(1)                           = gsn_csm_vector(wks,vp1,wp1,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfp1(:,:)))
     res1@tiMainString                 = main_str(1)
     res1@cnLevels                    := cnsflev
     plot(1)                           = gsn_csm_contour(wks,sfp1,res1)
     overlay(plot(1),tmpA(1))
     overlay(plot(1),tmpB(1))

     res0@cnLevels                    := cnwdlevd
     tmpA(2)                           = gsn_csm_contour(wks,udp1,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd1              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(2)                           = gsn_csm_vector(wks,vdp1,wdp1,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdp1(:,:)))
     res1@tiMainString                 = main_str(2)
     res1@cnLevels                    := cnsflevd
     plot(2)                           = gsn_csm_contour(wks,sfdp1,res1)
     overlay(plot(2),tmpA(2))
     overlay(plot(2),tmpB(2))

     res0@cnLevels                    := cnwdlev
     tmpA(3)                           = gsn_csm_contour(wks,ud2,res0)
     vecres@vcRefMagnitudeF            = vcRefMag              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(3)                           = gsn_csm_vector(wks,vd2,wd2,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfd2(:,:)))
     res1@tiMainString                 = main_str(3)
     res1@cnLevels                    := cnsflev
     plot(3)                           = gsn_csm_contour(wks,sfd2,res1)
     overlay(plot(3),tmpA(3))
     overlay(plot(3),tmpB(3))

     res0@cnLevels                    := cnwdlev
     tmpA(4)                           = gsn_csm_contour(wks,up2,res0)
     vecres@vcRefMagnitudeF            = vcRefMag            ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(4)                           = gsn_csm_vector(wks,vp2,wp2,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfp2(:,:)))
     res1@tiMainString                 = main_str(4)
     res1@cnLevels                    := cnsflev
     plot(4)                           = gsn_csm_contour(wks,sfp2,res1)
     overlay(plot(4),tmpA(4))
     overlay(plot(4),tmpB(4))

     res0@cnLevels                    := cnwdlevd
     tmpA(5)                           = gsn_csm_contour(wks,udp2,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd1             ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(5)                           = gsn_csm_vector(wks,vdp2,wdp2,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdp2(:,:)))
     res1@tiMainString                 = main_str(5)
     res1@cnLevels                    := cnsflevd
     plot(5)                           = gsn_csm_contour(wks,sfdp2,res1)
     overlay(plot(5),tmpA(5))
     overlay(plot(5),tmpB(5))

     res0@cnLevels                    := cnwdlev
     tmpA(6)                           = gsn_csm_contour(wks,ud3,res0)
     vecres@vcRefMagnitudeF            = vcRefMag              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(6)                           = gsn_csm_vector(wks,vd3,wd3,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfd3(:,:)))
     res1@tiMainString                 = main_str(6)
     res1@cnLevels                    := cnsflev
     plot(6)                           = gsn_csm_contour(wks,sfd3,res1)
     overlay(plot(6),tmpA(6))
     overlay(plot(6),tmpB(6))

     res0@cnLevels                    := cnwdlev
     tmpA(7)                           = gsn_csm_contour(wks,up3,res0)
     vecres@vcRefMagnitudeF            = vcRefMag            ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(7)                           = gsn_csm_vector(wks,vp3,wp3,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfp3(:,:)))
     res1@tiMainString                 = main_str(7)
     res1@cnLevels                    := cnsflev
     plot(7)                           = gsn_csm_contour(wks,sfp3,res1)
     overlay(plot(7),tmpA(7))
     overlay(plot(7),tmpB(7))

     res0@cnLevels                    := cnwdlevd
     tmpA(8)                           = gsn_csm_contour(wks,udp3,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd1             ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(8)                           = gsn_csm_vector(wks,vdp3,wdp3,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdp3(:,:)))
     res1@tiMainString                 = main_str(8)
     res1@cnLevels                    := cnsflevd
     plot(8)                           = gsn_csm_contour(wks,sfdp3,res1)
     overlay(plot(8),tmpA(8))
     overlay(plot(8),tmpB(8))

     res0@cnLevels                    := cnwdlevd
     tmpA(9)                           = gsn_csm_contour(wks,udd1,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(9)                           = gsn_csm_vector(wks,vdd1,wdd1,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdd1(:,:)))
     res1@tiMainString                 = main_str(9)
     res1@cnLevels                    := cnsflevd
     plot(9)                           = gsn_csm_contour(wks,sfdd1,res1)
     overlay(plot(9),tmpA(9))
     overlay(plot(9),tmpB(9))

     res0@cnLevels                    := cnwdlevd
     tmpA(10)                          = gsn_csm_contour(wks,upp1,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2             ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(10)                          = gsn_csm_vector(wks,vpp1,wpp1,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfpp1(:,:)))
     res1@tiMainString                 = main_str(10)
     res1@cnLevels                    := cnsflevd
     plot(10)                          = gsn_csm_contour(wks,sfpp1,res1)
     overlay(plot(10),tmpA(10))
     overlay(plot(10),tmpB(10))

     res0@cnLevels                    := cnwdlevd
     tmpA(11)                          = gsn_csm_contour(wks,uppdd1,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2       ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(11)                          = gsn_csm_vector(wks,vppdd1,wppdd1,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfppdd1(:,:)))
     res1@tiMainString                 = main_str(11)
     res1@cnLevels                    := cnsflevd
     plot(11)                          = gsn_csm_contour(wks,sfppdd1,res1)
     overlay(plot(11),tmpA(11))
     overlay(plot(11),tmpB(11))

     res0@cnLevels                    := cnwdlevd
     tmpA(12)                          = gsn_csm_contour(wks,udd2,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(12)                          = gsn_csm_vector(wks,vdd2,wdd2,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdd2(:,:)))
     res1@tiMainString                 = main_str(12)
     res1@cnLevels                    := cnsflevd
     plot(12)                          = gsn_csm_contour(wks,sfdd2,res1)
     overlay(plot(12),tmpA(12))
     overlay(plot(12),tmpB(12))

     res0@cnLevels                    := cnwdlevd
     tmpA(13)                          = gsn_csm_contour(wks,upp2,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2             ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(13)                          = gsn_csm_vector(wks,vpp2,wpp2,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfpp2(:,:)))
     res1@tiMainString                 = main_str(13)
     res1@cnLevels                    := cnsflevd
     plot(13)                          = gsn_csm_contour(wks,sfpp2,res1)
     overlay(plot(13),tmpA(13))
     overlay(plot(13),tmpB(13))

     res0@cnLevels                    := cnwdlevd
     tmpA(14)                          = gsn_csm_contour(wks,uppdd2,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2       ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(14)                          = gsn_csm_vector(wks,vppdd2,wppdd2,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfppdd2(:,:)))
     res1@tiMainString                 = main_str(14)
     res1@cnLevels                    := cnsflevd
     plot(14)                          = gsn_csm_contour(wks,sfppdd2,res1)
     overlay(plot(14),tmpA(14))
     overlay(plot(14),tmpB(14))

    ; panel plot only resources
     pres1                            = True                ; mods desired
     pres1@gsnFrame                   = True               ; save panel until both ready
    ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
    ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
    ;pres1@gsnPanelXF                 = (/0.07,0.57/)
     pres1@gsnPanelYWhiteSpacePercent = 8
     gsn_panel(wks,plot,(/5,3/),pres1)          ; create first panel
 
     delete([/ud1,ud2,ud3,up1,up2,up3,udp1,udp2,udp3,udd1,udd2,upp1,upp2,uppdd1,uppdd2/])
     delete([/vd1,vd2,vd3,vp1,vp2,vp3,vdp1,vdp2,vdp3,vdd1,vdd2,vpp1,vpp2,vppdd1,vppdd2/])
     delete([/wd1,wd2,wd3,wp1,wp2,wp3,wdp1,wdp2,wdp3,wdd1,wdd2,wpp1,wpp2,wppdd1,wppdd2/])
     delete([/sfd1,sfd2,sfd3,sfp1,sfp2,sfp3,sfdp1,sfdp2,sfdp3,sfdd1,sfpp1,sfdd2,sfpp2,sfppdd1,sfppdd2/])
  
 end do      


end

