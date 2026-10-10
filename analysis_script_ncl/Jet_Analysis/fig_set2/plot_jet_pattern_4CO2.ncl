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

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
 
  d2r       = get_d2r("float")

  left_str  = (/"v2.LR.piControl            ", "v2.LR.abrupt-4xCO2           ", "4xCO2 - piControl (WC)             ", \
                "v2.LR.piControl.CRYOnoismf ", "v2.LR.CRYOnoismf.abrupt-4xCO2", "4xCO2 - piControl (CRYOnoismf)     ", \
                "v2.LR.piControl.CRYO       ", "v2.LR.abrupt-4xCO2.CRYO      ", "4xCO2 - piControl (CRYO)           ", \
                "CRYOnoismf - WC (piControl)", "CRYOnoismf - WC (4xCO2)      ", "CRYOnoismf - WC (4xCO2 - piControl)", \
                "CRYO - WC (piControl)      ", "CRYO - WC (4xCO2)            ", "CRYO - WC (4xCO2 - piControl)      "/)

 ;define regions for eof analysis
  mode       = ("JET") 
  vars       = (/"U"/)
  facmod     = (/1.0 /)
  varunt     = (/"m s~S~-1~N~"/)

  time_tag   = (/"0051-0150"/)
  time_str   = (/"Year 51-150"/)

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    do k = 0,nexps-1,1

      fl   = data_dir+"/"+exp_name(k)+"_"+mode+"_"+vars+"_eof_Monthly_"+time_tag+".nc"
      f    = addfile(fl,"r")

      rad       = 4.0*atan(1.0)/180.0
      wgty      = cos(f->lat*rad)
      utc_date  = cd_calendar(f->time,-1)
      yyyy      = yyyymm_to_yyyyfrac(utc_date, 0.0)

      if(seasons(i).eq."ANN")then
        xa  = month_to_annual (f->$vars$, 1)
      else
        xa  = month_to_season (f->$vars$, seasons(i))
      end if

      if (k.eq.0) then 
        eofd1  = dim_avg_n_Wrap(xa,0)
      end if 

      if (k.eq.1) then
        eofd2  = dim_avg_n_Wrap(xa,0)
      end if

      if (k.eq.2) then
        eofd3  = dim_avg_n_Wrap(xa,0)
      end if

      if (k.eq.3) then
        eofp1  = dim_avg_n_Wrap(xa,0)
      end if

      if (k.eq.4) then
        eofp2  = dim_avg_n_Wrap(xa,0)
      end if

      if (k.eq.5) then
        eofp3  = dim_avg_n_Wrap(xa,0)
      end if
      delete(xa)
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

    eofdp3   = eofd3
    eofdp3   = eofp3 - eofd3
    glCordp3 = pattern_cor(eofp3,eofd3,1.0,0)
    glAvedp3 = wgt_arearmse(eofp3,eofd3, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordp3)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedp3)
    mstrdp3  = mstrx1 + "~C~" + mstrx2

    ;;difference due to configuration
    eofdd1   = eofd1
    eofdd1   = eofd2 - eofd1
    glCordd1 = pattern_cor(eofd1,eofd2,1.0,0)
    glAvedd1 = wgt_arearmse(eofd1,eofd2, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd1)
    mstrdd1  = mstrx1 + "~C~" + mstrx2

    eofdd2   = eofd1
    eofdd2   = eofd3 - eofd1
    glCordd2 = pattern_cor(eofd1,eofd3,1.0,0)
    glAvedd2 = wgt_arearmse(eofd1,eofd3, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCordd2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvedd2)
    mstrdd2  = mstrx1 + "~C~" + mstrx2

    eofpp1   = eofp1
    eofpp1   = eofp2 - eofp1
    glCorpp1 = pattern_cor(eofp1,eofp2,1.0,0)
    glAvepp1 = wgt_arearmse(eofp1,eofp2, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorpp1)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvepp1)
    mstrpp1  = mstrx1 + "~C~" + mstrx2
   
    eofpp2   = eofp1
    eofpp2   = eofp3 - eofp1
    glCorpp2 = pattern_cor(eofp1,eofp3,1.0,0)
    glAvepp2 = wgt_arearmse(eofp1,eofp3, wgty, 1.0,0)
    mstrx1   = "corr. w/r " + sprintf("%.2f", glCorpp2)
    mstrx2   = "RMSD  w/r " + sprintf("%.2f", glAvepp2)
    mstrpp2  = mstrx1 + "~C~" + mstrx2

    eofppdd1  = eofpp1
    eofppdd1  = eofpp1 - eofdd1
    glCorddpp1= pattern_cor(eofpp1,eofdd1,1.0,0)
    glAveddpp1= wgt_arearmse(eofpp1,eofdd1, wgty, 1.0,0)
    mstrx1    = "CORR. w/r " + sprintf("%.2f", glCorddpp1)
    mstrx2    = "RMSD  w/r " + sprintf("%.2f", glAveddpp1)
    mstrppdd1 = mstrx1 + "~C~" + mstrx2

    eofppdd2  = eofpp2
    eofppdd2  = eofpp2 - eofdd2
    glCorddpp2= pattern_cor(eofpp2,eofdd2,1.0,0)
    glAveddpp2= wgt_arearmse(eofpp2,eofdd2, wgty, 1.0,0)
    mstrx1    = "CORR. w/r " + sprintf("%.2f", glCorddpp2)
    mstrx2    = "RMSD  w/r " + sprintf("%.2f", glAveddpp2)
    mstrppdd2 = mstrx1 + "~C~" + mstrx2

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fignam          = "Jet_pattern_4CO2_"+vars+"_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)

    plot      = new((/15/),graphic)
    
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

     FontHeightF = 0.020
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
     res@tmYLLabelFontHeightF     = 0.005

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
     cnlevs1  = (/-20,-15,-10,-5,-1,0,1,5,10,15,20/) 
     cnlevs2  = (/-8,-4,-2,-1,-0.5,-0.2,0.2,0.5,1,2,4,8/)

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
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = "" 
     res@gsnCenterString   = "" 
     res@tiMainString      = left_str(0) 
     res@tiMainOffsetYF    = 0.005
     res@lbLabelBarOn      = True
     plot(0)   = gsn_csm_contour_map_polar(wks,eofd1,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = "" 
     res@gsnCenterString   = "" 
     res@tiMainOffsetYF    = 0.005
     res@tiMainString      = left_str(1)
     res@lbLabelBarOn      = True
     plot(1) = gsn_csm_contour_map_polar(wks,eofp1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrdp1
     res@gsnCenterString   = "" 
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(2)
     res@lbLabelBarOn      = True
     plot(2) = gsn_csm_contour_map_polar(wks,eofdp1,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = ""
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.005
     res@tiMainString      = left_str(3)
     res@lbLabelBarOn      = True
     plot(3)   = gsn_csm_contour_map_polar(wks,eofd2,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = ""
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.005
     res@tiMainString      = left_str(4)
     res@lbLabelBarOn      = True
     plot(4) = gsn_csm_contour_map_polar(wks,eofp2,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrdp2
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(5)
     res@lbLabelBarOn      = True
     plot(5) = gsn_csm_contour_map_polar(wks,eofdp2,res)
 
     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = ""
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.005
     res@tiMainString      = left_str(6)
     res@lbLabelBarOn      = True
     plot(6)   = gsn_csm_contour_map_polar(wks,eofd3,res)

     res@cnFillPalette     = "BlueWhiteOrangeRed"
     res@cnLevels          := cnlevs1
     res@gsnLeftString     = vars + " ("+varunt+")"
     res@gsnRightString    = ""
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.005
     res@tiMainString      = left_str(7)
     res@lbLabelBarOn      = True
     plot(7) = gsn_csm_contour_map_polar(wks,eofp3,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrdp3
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(8)
     res@lbLabelBarOn      = True
     plot(8) = gsn_csm_contour_map_polar(wks,eofdp3,res)
   
     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrdd1
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(9)
     res@lbLabelBarOn      = True
     plot(9)   = gsn_csm_contour_map_polar(wks,eofdd1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrpp1
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(10)
     res@lbLabelBarOn      = True
     plot(10) = gsn_csm_contour_map_polar(wks,eofpp1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrppdd1
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(11)
     res@lbLabelBarOn      = True
     plot(11) = gsn_csm_contour_map_polar(wks,eofppdd1,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrdd2
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(12)
     res@lbLabelBarOn      = True
     plot(12)   = gsn_csm_contour_map_polar(wks,eofdd2,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrpp2
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(13)
     res@lbLabelBarOn      = True
     plot(13) = gsn_csm_contour_map_polar(wks,eofpp2,res)

     res@cnFillPalette     = "CBR_coldhot"
     res@cnLevels          := cnlevs2
     res@gsnLeftString     = ""
     res@gsnRightString    = mstrppdd2
     res@gsnCenterString   = ""
     res@tiMainOffsetYF    = 0.00
     res@tiMainString      = left_str(14)
     res@lbLabelBarOn      = True
     plot(14) = gsn_csm_contour_map_polar(wks,eofppdd2,res)

    ; panel plot only resources
     pres1                            = True                ; mods desired
     pres1@gsnFrame                   = True               ; save panel until both ready
    ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
    ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
    ;pres1@gsnPanelXF                 = (/0.07,0.57/)
     pres1@gsnPanelYWhiteSpacePercent = 8
     gsn_panel(wks,plot,(/5,3/),pres1)          ; create first panel
 
     delete([/eofd1,eofd2,eofd3,eofdd1,eofdd2,eofp1,eofp2,eofp3,eofpp1,eofpp2,eofppdd1,eofppdd2/])  
     delete([/utc_date,yyyy/])
  
end do 

end

