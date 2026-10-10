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

  l_plot_test        = True
  statistical_siglev = 0.05

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/EOF_analysis/v2_4CO2/data"
  exp_name  = (/"v2.LR.piControl", \
                "v2.LR.piControl.CRYO", \
                "v2.LR.abrupt-4xCO2", \
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
 
  d2r       = get_d2r("float")
  left_str  = (/"v2.LR.piControl","v2.LR.abrupt-4xCO2", "4xCO2 - piControl (WC)", \ 
                "v2.LR.piControl.CRYO","v2.LR.abrupt-4xCO2.CRYO", "4xCO2 - piControl (CRYO)", \
                "CRYO - WC (piControl)", "CRYO - WC (4xCO2)", "CRYO - WC (4xCO2 - piControl)"/) 

  time_tag  = (/"0001-0050"/)
  time_str  = (/"Year 1-50"/)

 ;define regions for eof analysis
  mode       = ("Jet") 
  vars       = (/"LWCF"/)
  facmod     = (/1.0 /)
  varunt     = (/"W m~S~-2~N~"/)

do kk = 0, 2

  if(kk.eq.0)then
    time_tag  = (/"0001-0050"/)
    time_str  = (/"Year 1-50"/)
  end if

  if(kk.eq.1)then
    time_tag  = (/"0051-0100"/)
    time_str  = (/"Year 51-100"/)
  end if

  if(kk.eq.2)then
    time_tag  = (/"0101-0150"/)
    time_str  = (/"Year 101-150"/)
  end if

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))
    fl1     = "Jet_ind_U_"+seasons(i)+"_"+time_tag+".nc"
    f1      = addfile(fl1,"r")

    do k = 0,nexps-1,1

      fl2    = data_dir+"/"+exp_name(k)+"_SAM_"+vars+"_eof_"+seasons(i)+"_"+time_tag+".nc"
      f2     = addfile(fl2,"r")

      rad    = 4.0*atan(1.0)/180.0
      wgty   = cos(f2->lat*rad)
      varstr = vars+"_anom"
      var    = f2->$varstr$
      varstr = "jetd_pos_"+exp_name(k)
      vts    = f1->$varstr$
      rc1    = regCoef_n(vts, var, 0, 0)
      rc1!0  = "lat"
      rc1!1  = "lon"
      rc1&lat= var&lat
      rc1&lon= var&lon
       
      ;;;;;;;t-test;;;;;;;;
      drc1   = dimsizes(rc1)
      df1    = rc1@nptxy-2   ; degrees of freedom
      tval1  = rc1@tval      ; t-statistic
      xb1    = tval1
      xb1    = 0.50
      pb1    = betainc(df1/(df1+tval1^2),df1/2.0,xb1)
      prob1  = reshape(pb1,(/drc1(0),drc1(1)/))
      prob1  = where(prob1.le.statistical_siglev,rc1,rc1@_FillValue)
      rgprob = (/  tofloat(prob1) /)
      delete([/drc1,df1,tval1,xb1,pb1,prob1/])
      rgprob!0   = "lat"
      rgprob!1   = "lon"
      rgprob&lat = var&lat
      rgprob&lon = var&lon

      if (k.eq.0) then 
        eofd1  = rc1
        prod1  = rgprob 
      end if 

      if (k.eq.1) then
        eofd2  = rc1
        prod2  = rgprob
      end if

      if (k.eq.2) then
        eofp1  = rc1
        prop1  = rgprob
      end if

      if (k.eq.3) then
        eofp2  = rc1
        prop2  = rgprob
      end if

    delete ([/rgprob,var,vts,rc1/])

    end do 

    ;;difference due to 4CO2
    eofdp1   = eofd1
    eofdp1   = eofp1 - eofd1
    glCordp1 = pattern_cor(eofp1,eofd1,1.0,0)
    glAvedp1 = wgt_arearmse(eofp1,eofd1, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp1)
    mstrdp1  = mstrx1 + "~C~" + mstrx2

    eofdp2   = eofd2
    eofdp2   = eofp2 - eofd2
    glCordp2 = pattern_cor(eofp2,eofd2,1.0,0)
    glAvedp2 = wgt_arearmse(eofp2,eofd2, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp2)
    mstrdp2  = mstrx1 + "~C~" + mstrx2

    ;;difference due to configuration
    eofdd1   = eofd1
    eofdd1   = eofd2 - eofd1
    glCordd1 = pattern_cor(eofd1,eofd2,1.0,0)
    glAvedd1 = wgt_arearmse(eofd1,eofd2, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd1)
    mstrdd1  = mstrx1 + "~C~" + mstrx2

    eofpp1   = eofp1
    eofpp1   = eofp2 - eofp1
    glCorpp1 = pattern_cor(eofp1,eofp2,1.0,0)
    glAvepp1 = wgt_arearmse(eofp1,eofp2, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorpp1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvepp1)
    mstrpp1  = mstrx1 + "~C~" + mstrx2
   
    eofppdd  = eofpp1
    eofppdd  = eofpp1 - eofdd1
    glCorddpp= pattern_cor(eofpp1,eofdd1,1.0,0)
    glAveddpp= wgt_arearmse(eofpp1,eofdd1, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorddpp)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAveddpp)
    mstrppdd = mstrx1 + "~C~" + mstrx2

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fignam          = "Jet_pos_REG_"+vars+"_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)

    plot      = new((/9/),graphic)
    dum       = new((/4/),graphic)
    
    ;============================================================
    ; PLOTS
    ;============================================================
     res                       = True         
     res@gsnDraw               = False        ; don't draw yet
     res@gsnFrame              = False        ; don't advance frame yet
     res@gsnPolar              = "SH"
 
     res@mpFillOn              = False        ; turn off map fill
     res@mpMaxLatF             = -20
    ;res@mpCenterLonF          = 180

     FontHeightF = 0.018
     res@tiMainFontHeightF     = FontHeightF*1.5
     res@tmYLLabelFontHeightF  = FontHeightF
     res@tmXBLabelFontHeightF  = FontHeightF
     res@lbLabelFontHeightF    = FontHeightF
     res@gsnStringFontHeightF  = FontHeightF*1.2

     res@cnFillOn              = True         ; turn on color fill
     res@cnLinesOn             = False        ; True is default
     res@cnLineLabelsOn        = False        ; True is default

     res@lbLabelBarOn             = True
     res@lbOrientation            = "Vertical"
     res@lbTitleString            = "" 
     res@lbTitlePosition          = "Right"     ; title location
     res@lbTitleDirection         = "Across"    ; letter angle
     res@lbTitleAngleF            = 90.         ; title angle
     res@lbTitleFontHeightF       = FontHeightF ; font height
     res@lbBoxLineColor           = "black"  
    ;res@lbBoxLinesOn             = False
    ;res@pmLabelBarParallelPosF   = 0.61
     res@pmLabelBarOrthogonalPosF = 0.15

     res@tmXTOn                   = False
     res@tmYLLabelFontHeightF     = 0.015

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

     res@cnLevelSelectionMode = "ExplicitLevels"
     maxlev   = 11
     cnlevs1  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.0,0.1,0.5,1.0,2.0,4.0/)
     cnlevs2  = (/-4.0,-2.0,-1.0,-0.5,-0.1,0.1,0.5,1.0,2.0,4.0/)

    ;*******************************************
    ; first plot
    ;*******************************************
     res@gsnPolarLabelFontHeightF  = FontHeightF
     res@mpLabelFontHeightF        = FontHeightF*2.0
     res@mpGridLineThicknessF      = lthick
     res@mpLimbLineThicknessF      = lthick
     res@mpPerimLineThicknessF     = lthick
     res@mpGeophysicalLineThicknessF = lthick

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + "(" + varunt+")"
     res@gsnRightString    = "Jet position" 
     res@gsnCenterString   = "" 
     res@tiMainString      = left_str(0) 
     res@tiMainOffsetYF    = 0.01
     res@lbLabelBarOn      = True
     plot(0)   = gsn_csm_contour_map_polar(wks,eofd1,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + "(" + varunt+")"
     res@gsnRightString    = "Jet position" 
     res@gsnCenterString   = "" 
     res@tiMainOffsetYF    = 0.01
     res@tiMainString      = left_str(1)
     res@lbLabelBarOn      = True
     plot(1) = gsn_csm_contour_map_polar(wks,eofp1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = "~F33~D~N~~F~"+ vars + "(" + varunt+")"
     res@gsnRightString    = mstrdp1
     res@gsnCenterString   = "" 
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(2)
     res@lbLabelBarOn      = True
     plot(2) = gsn_csm_contour_map_polar(wks,eofdp1,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + "(" + varunt+")"     
     res@gsnRightString    = "Jet position" 
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.01
     res@tiMainString      = left_str(3)
     res@lbLabelBarOn      = True
     plot(3)   = gsn_csm_contour_map_polar(wks,eofd2,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + "(" + varunt+")" 
     res@gsnRightString    = "Jet position" 
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.01
     res@tiMainString      = left_str(4)
     res@lbLabelBarOn      = True
     plot(4) = gsn_csm_contour_map_polar(wks,eofp2,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = "~F33~D~N~~F~"+ vars + "(" + varunt+")"
     res@gsnRightString    = mstrdp2
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(5)
     res@lbLabelBarOn      = True
     plot(5) = gsn_csm_contour_map_polar(wks,eofdp2,res)
    
     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = "~F33~D~N~~F~"+ vars + "(" + varunt+")"
     res@gsnRightString    = mstrdd1
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(6)
     res@lbLabelBarOn      = True
     plot(6)   = gsn_csm_contour_map_polar(wks,eofdd1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = "~F33~D~N~~F~"+ vars + "(" + varunt+")"     
     res@gsnRightString    = mstrpp1
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(7)
     res@lbLabelBarOn      = True
     plot(7) = gsn_csm_contour_map_polar(wks,eofpp1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = "~F33~D~N~~F~"+ vars + "(" + varunt+")"
     res@gsnRightString    = mstrppdd
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(8)
     res@lbLabelBarOn      = True
     plot(8) = gsn_csm_contour_map_polar(wks,eofppdd,res)

     rescn2 = True
     rescn2@lbLabelBarOn    = False
     rescn2@gsnTickMarksOn  = False     ; no tickmarks
     rescn2@gsnDraw         = False     ; don't draw
     rescn2@gsnFrame        = False     ; don't advance frame
     rescn2@tiMainString    = ""
     rescn2@gsnLeftString   = ""
     rescn2@gsnRightString  = ""
     rescn2@cnMinLevelValF  = (1.0 - statistical_siglev)*100.0
     rescn2@cnMaxLevelValF  = 100.0
    ;rescn2@cnLevelSpacingF = 1.
     rescn2@cnMonoFillPattern = True
     rescn2@cnMonoFillColor = True
     rescn2@cnFillOn        = True      ; color fill
     rescn2@cnLinesOn       = False
     rescn2@cnLineLabelsOn  = False        ; True is default
     rescn2@cnInfoLabelOn   = False                ; turn off info label
     rescn2@cnFillMode      = "AreaFill" ; raster fill significance patterns
     rescn2@cnFillPattern   = 17
     rescn2@cnFillDotSizeF  = 0.0030
     rescn2@cnFillColor     = "black"
     rescn2@cnFillOpacityF  = 0.5
     dum(0) = gsn_csm_contour(wks,prod1,rescn2) ; contours are at 95,100
     dum(1) = gsn_csm_contour(wks,prop1,rescn2) ; contours are at 95,100
     dum(2) = gsn_csm_contour(wks,prod2,rescn2) ; contours are at 95,100
     dum(3) = gsn_csm_contour(wks,prop2,rescn2) ; contours are at 95,100

     overlay(plot(0),dum(0))
     overlay(plot(1),dum(1))
     overlay(plot(3),dum(2))
     overlay(plot(4),dum(3))

    ; panel plot only resources
     pres1                            = True                ; mods desired
     pres1@gsnFrame                   = True               ; save panel until both ready
    ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
    ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
    ;pres1@gsnPanelXF                 = (/0.07,0.57/)
     pres1@gsnPanelYWhiteSpacePercent = 8
     gsn_panel(wks,plot,(/3,3/),pres1)          ; create first panel
 
     delete([/eofd1,eofd2,eofdd1,eofp1,eofp2,eofpp1,eofppdd/])  
  
  end do      

end do 

end

