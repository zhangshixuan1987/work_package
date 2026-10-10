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
  exp_name  = (/"ERA5", \
                "ERA_20C", \
                "NOAA_20C"/)
  nexps     = dimsizes(exp_name)

  time_tag  = (/"1981-2010"/)
  time_str  = (/"Year 1981-2010"/)

  main_str  = (/"ERA5", \
                "ERA_20C", \
                "NOAA_20C", \
                "ERA_20C - ERA5", \
                "NOAA_20C - ERA5", \
                "NOAA_20C - ERA_20C"/)
  
  d2r       = get_d2r("float")
  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
  facsfm    = 1e-9
  facwm     = 100.0 
  latS      = -90.0
  latN      =  90.0
 
 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    do k = 0,nexps-1,1

    fl   = data_dir+"/"+exp_name(k)+"_HC_analysis_"+seasons(i)+"_"+time_tag+".nc"
    f    = addfile(fl,"r")

    rad       = 4.0*atan(1.0)/180.0
    wgty      = cos(f->lat*rad)

    if (k.eq.0) then 
      ud1    = dim_avg_n_Wrap(f->U,0)
      vd1    = dim_avg_n_Wrap(f->V,0)
      wd1    = dim_avg_n_Wrap(f->W,0)
      sfd1   = dim_avg_n_Wrap(f->SF,0)
      sfd1   = sfd1 * facsfm
      wd1    = wd1 * facwm
    end if 

    if (k.eq.1) then
      ud2    = dim_avg_n_Wrap(f->U,0)
      vd2    = dim_avg_n_Wrap(f->V,0)
      wd2    = dim_avg_n_Wrap(f->W,0)
      sfd2   = dim_avg_n_Wrap(f->SF,0)
      sfd2   = sfd2 * facsfm
      wd2    = wd2 * facwm
    end if

    if (k.eq.2) then
      ud3    = dim_avg_n_Wrap(f->U,0)
      vd3    = dim_avg_n_Wrap(f->V,0)
      wd3    = dim_avg_n_Wrap(f->W,0)
      sfd3   = dim_avg_n_Wrap(f->SF,0)
      sfd3   = sfd3 * facsfm
      wd3    = wd3 * facwm
    end if

    end do  

    ;;difference between different products
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
    mstrdd3  = mstrx1 + "~C~" + mstrx2

    udd3     = ud2
    udd3     = ud3 - ud2
    vdd3     = vd2
    vdd3     = vd3 - vd2
    wdd3     = wd2
    wdd3     = wd3 - wd2
    sfdd3    = sfd2
    sfdd3    = sfd3 - sfd2
    glCordd3 = pattern_cor(sfd2,sfd3,1.0,0)
    glAvedd3 = wgt_arearmse(sfd2,sfd3, 1.0, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd3)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd3)
    mstrdd3  = mstrx1 + "~C~" + mstrx2

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype             = "png"
    wtype@wkWidth     = 2500
    wtype@wkHeight    = 2500
    fignam            = "Hadley_circulation_analysis_"+seasons(i)+"_"+time_tag
    wks               = gsn_open_wks(wtype,fignam)

    plot              = new((/6/),graphic)
    tmpA              = new((/6/),graphic)
    tmpB              = new((/6/),graphic)

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
     tmpA(1)                           = gsn_csm_contour(wks,ud2,res0)
     vecres@vcRefMagnitudeF            = vcRefMag              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(1)                           = gsn_csm_vector(wks,vd2,wd2,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfd2(:,:)))
     res1@tiMainString                 = main_str(1)
     res1@cnLevels                    := cnsflev
     plot(1)                           = gsn_csm_contour(wks,sfd2,res1)
     overlay(plot(1),tmpA(1))
     overlay(plot(1),tmpB(1))

     res0@cnLevels                    := cnwdlev
     tmpA(2)                           = gsn_csm_contour(wks,ud3,res0)
     vecres@vcRefMagnitudeF            = vcRefMag              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(2)                           = gsn_csm_vector(wks,vd3,wd3,vecres)
     res1@cnFillPalette                = Contour_color_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfd3(:,:)))
     res1@tiMainString                 = main_str(2)
     res1@cnLevels                    := cnsflev
     plot(2)                           = gsn_csm_contour(wks,sfd3,res1)
     overlay(plot(2),tmpA(2))
     overlay(plot(2),tmpB(2))


     res0@cnLevels                    := cnwdlevd
     tmpA(3)                           = gsn_csm_contour(wks,udd1,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(3)                           = gsn_csm_vector(wks,vdd1,wdd1,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdd1(:,:)))
     res1@tiMainString                 = main_str(3)
     res1@cnLevels                    := cnsflevd
     plot(3)                           = gsn_csm_contour(wks,sfdd1,res1)
     overlay(plot(3),tmpA(3))
     overlay(plot(3),tmpB(3))

     res0@cnLevels                    := cnwdlevd
     tmpA(4)                           = gsn_csm_contour(wks,udd2,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(4)                           = gsn_csm_vector(wks,vdd2,wdd2,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdd2(:,:)))
     res1@tiMainString                 = main_str(4)
     res1@cnLevels                    := cnsflevd
     plot(4)                           = gsn_csm_contour(wks,sfdd2,res1)
     overlay(plot(4),tmpA(4))
     overlay(plot(4),tmpB(4))

     res0@cnLevels                    := cnwdlevd
     tmpA(5)                           = gsn_csm_contour(wks,udd3,res0)
     vecres@vcRefMagnitudeF            = vcRefMagd2              ; define vector ref mag
     vecres@vcRefAnnoString1           = vecres@vcRefMagnitudeF
     tmpB(5)                           = gsn_csm_vector(wks,vdd3,wdd3,vecres)
     res1@cnFillPalette                = Diff_corlor_bar
     res1@gsnRightString               = "Mean: "+sprintf("%5.2f",avg(sfdd3(:,:)))
     res1@tiMainString                 = main_str(5)
     res1@cnLevels                    := cnsflevd
     plot(5)                           = gsn_csm_contour(wks,sfdd3,res1)
     overlay(plot(5),tmpA(5))
     overlay(plot(5),tmpB(5))

   ; panel plot only resources
     pres1                            = True                ; mods desired
     pres1@gsnFrame                   = True               ; save panel until both ready
    ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
    ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
    ;pres1@gsnPanelXF                 = (/0.07,0.57/)
     pres1@gsnPanelYWhiteSpacePercent = 8
     gsn_panel(wks,plot,(/2,3/),pres1)          ; create first panel
 
     delete([/ud1,ud2,udd1/])
     delete([/vd1,vd2,vdd1/])
     delete([/wd1,wd2,wdd1/])
     delete([/sfd1,sfd2,sfdd1/])
  
end do      


end

