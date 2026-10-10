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
  obs_name  = (/"NOAA_20C","ERA_20C","ERA5"/)
  epi_name  = (/"v2.LR.piControl", \
                "v2.LR.piControl.CRYOnoismf", \
                "v2.LR.piControl.CRYO"/)
  epd_name  = (/"v2.LR.abrupt-4xCO2", \
                "v2.LR.CRYOnoismf.abrupt-4xCO2", \
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(epi_name)

  lfpi_strs = (/"WC (piControl)",\
                "CRYOnoismf (piControl)",\
                "CRYO (piControl)"/)
  lfpd_strs = (/"WC (4xCO2)",\
                "CRYOnoismf (4xCO2)",\
                "CRYO (4xCO2)"/)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
 
  d2r       = get_d2r("float")

 ;define regions for eof analysis
  mode       = ("JET") 
  vars       = (/"U"/)
  facmod     = (/1.0 /)
  varunt     = (/"m s~S~-1~N~"/)

  time_opd  = (/"1981-2010"/)
  time_opi  = (/"1901-1930"/)
  time_tag   = (/"0051-0150"/)
  time_str   = (/"Year 51-150"/)

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    ;read and process the reanalysis;;
    flopd     = data_dir+"/"+obs_name(0:2)+"_"+mode+"_"+vars+"_eof_Monthly_"+time_opd+".nc"
    fopd      = addfiles(flopd,"r")
    flopi     = data_dir+"/"+obs_name(0:1)+"_"+mode+"_"+vars+"_eof_Monthly_"+time_opi+".nc"
    fopi      = addfiles(flopi,"r")
    ListSetType (fopd, "join")       
    ListSetType (fopi, "join")       

    xpd    = dim_avg_n_Wrap(fopd[:]->U,0)
    xpi    = dim_avg_n_Wrap(fopi[:]->U,0)
    if(seasons(i).eq."ANN")then
      ypd  = month_to_annual (xpd, 1)
      ypi  = month_to_annual (xpi, 1)
    else
      ypd  = month_to_season (xpd, seasons(i))
      ypi  = month_to_season (xpi, seasons(i))
    end if

    eofopd = dim_avg_n_Wrap(ypd,0)
    eofopi = dim_avg_n_Wrap(ypi,0)

    delete([/xpd,xpi,ypd,ypi/])

    ;read and process the e3smv2;;
    do k = 0,nexps-1,1

      flpd   = data_dir+"/"+epd_name(k)+"_"+mode+"_"+vars+"_eof_Monthly_"+time_tag+".nc"
      fpd    = addfile(flpd,"r")
      flpi   = data_dir+"/"+epi_name(k)+"_"+mode+"_"+vars+"_eof_Monthly_"+time_tag+".nc"
      fpi    = addfile(flpi,"r")

      rad    = 4.0*atan(1.0)/180.0
      wgty   = cos(fpd->lat*rad)

      if(seasons(i).eq."ANN")then 
        xpd  = month_to_annual (fpd->U, 1)
        xpi  = month_to_annual (fpi->U, 1)
      else
        xpd  = month_to_season (fpd->U, seasons(i))
        xpi  = month_to_season (fpi->U, seasons(i))
      end if 
     
      eofpd    = dim_avg_n_Wrap(xpd,0)
      eofpi    = dim_avg_n_Wrap(xpi,0)
      delete([/xpd,xpi/])

      eofopdpi = eofopd
      eofopdpi = eofopd - eofopi
      eofpdpi  = eofpd
      eofpdpi  = eofpd - eofpi

      ;;difference due to configuration
      glCordd1 = pattern_cor(eofpd,eofopd,1.0,0)
      glAvedd1 = wgt_areaave(eofpd-eofopd,wgty, 1.0,0)
      glRmsdd1 = wgt_arearmse(eofpd,eofopd, wgty, 1.0,0)
      mstrx1   = "corr w/r obs.   " + sprintf("%5.2f", glCordd1)
      mstrx2   = "diff w/r obs.   " + sprintf("%5.2f", glAvedd1)
      mstrx3   = "rmsd w/r obs.   " + sprintf("%5.2f", glRmsdd1)
      mstrpd    = mstrx1 + "~C~" + mstrx2 + "~C~" + mstrx3
      delete([/glCordd1,glAvedd1,glRmsdd1,mstrx1,mstrx2,mstrx3/])

      glCordd1 = pattern_cor(eofpi,eofopi,1.0,0)
      glAvedd1 = wgt_areaave(eofpi-eofopi,wgty, 1.0,0)
      glRmsdd1 = wgt_arearmse(eofpi,eofopi, wgty, 1.0,0)
      mstrx1   = "corr w/r obs.   " + sprintf("%5.2f", glCordd1)
      mstrx2   = "diff w/r obs.   " + sprintf("%5.2f", glAvedd1)
      mstrx3   = "rmsd w/r obs.   " + sprintf("%5.2f", glRmsdd1)
      mstrpi    = mstrx1 + "~C~" + mstrx2 + "~C~" + mstrx3
      delete([/glCordd1,glAvedd1,glRmsdd1,mstrx1,mstrx2,mstrx3/])

      glCordd1 = pattern_cor(eofpdpi,eofopdpi,1.0,0)
      glAvedd1 = wgt_areaave(eofpdpi-eofopdpi,wgty, 1.0,0)
      glRmsdd1 = wgt_arearmse(eofpdpi,eofopdpi, wgty, 1.0,0)
      mstrx1   = "corr w/r obs.   " + sprintf("%5.2f", glCordd1)
      mstrx2   = "diff w/r obs.   " + sprintf("%5.2f", glAvedd1)
      mstrx3   = "rmsd w/r obs.   " + sprintf("%5.2f", glRmsdd1)
      mstrpdpi = mstrx1 + "~C~" + mstrx2 + "~C~" + mstrx3
      delete([/glCordd1,glAvedd1,glRmsdd1,mstrx1,mstrx2,mstrx3/])

     ;;;;;;;;;;;;;;;;;;;;;;;;;;
     ; Start to process plot 
     ;;;;;;;;;;;;;;;;;;;;;;;;;;
      if (k.eq.0) then
        wtype           = "png"
        wtype@wkWidth   = 2500
        wtype@wkHeight  = 2500
        fignam          = "Jet_pattern_obs_"+vars+"_"+seasons(i)+"_"+time_tag
        wks             = gsn_open_wks(wtype,fignam)
        plot            = new((/nexps+1,3/),graphic)
        tmp             = new((/nexps+1,3/),graphic)
        cnlevs1         = (/-20,-15,-10,-5,-1,0,1,5,10,15,20/)
        cnlevs2         = (/-8,-4,-2,-1,-0.1,0.1,1,2,4,8/)
      end if 

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

     res@cnLevelSelectionMode     = "ExplicitLevels"

    ;*******************************************
    ; first plot
    ;*******************************************
     res@gsnPolarLabelFontHeightF  = FontHeightF
     res@mpLabelFontHeightF        = FontHeightF*2.0
     res@mpGridLineThicknessF      = lthick
     res@mpLimbLineThicknessF      = lthick
     res@mpPerimLineThicknessF     = lthick
     res@mpGeophysicalLineThicknessF = lthick

     res@cnFillOn               = True         ; turn on color fill
     res@cnLinesOn              = False        ; True is default
     res@cnLineLabelsOn         = False        ; True is default
     res@gsnLeftString          = ""
     res@gsnRightString         = "" 
     res@gsnCenterString        = "" 

     res0 = res
     res0@cnFillPalette         = "BlueWhiteOrangeRed"
     res0@cnLevels             := cnlevs1
     res0@lbLabelBarOn          = True
     if (k.eq.0) then 
       res0@tiMainOffsetYF      = 0.015
       res0@gsnLeftString       = vars + " ("+varunt+")"
       res0@tiMainString        = "Obs_1901-1930"
       plot(k,0)   = gsn_csm_contour_map_polar(wks,eofopi,res0)
       res0@gsnLeftString       = vars + " ("+varunt+")"
       res0@tiMainString        = "Obs_1981-2010"
       plot(k,1)   = gsn_csm_contour_map_polar(wks,eofopd,res0)
     end if 
     res0@tiMainOffsetYF        = 0.0
     res0@gsnLeftString         = vars + " ("+varunt+")"
     res0@tiMainString          = lfpi_strs(k)
     res0@gsnRightString        = mstrpi
     plot(k+1,0)                = gsn_csm_contour_map_polar(wks,eofpi,res0)
     res0@tiMainString          = lfpd_strs(k)
     res0@gsnRightString        = mstrpd
     plot(k+1,1)                = gsn_csm_contour_map_polar(wks,eofpd,res0)

     res1 = res
     res1@cnFillPalette         = "CBR_coldhot"
     res1@cnLevels             := cnlevs2
     res1@lbLabelBarOn          = True
     if (k.eq.0) then
       res1@tiMainOffsetYF      = 0.015
       res1@gsnLeftString       = vars + " ("+varunt+")"
       res1@tiMainString        = "Obs_1981-2010 - Obs_1901-1930"
       plot(k,2)   = gsn_csm_contour_map_polar(wks,eofopdpi,res1)
     end if
     res1@tiMainOffsetYF        = 0.0
     res1@gsnLeftString         = vars + " ("+varunt+")"
     res1@tiMainString          = lfpd_strs(k)+ " - " + lfpi_strs(k)
     res1@gsnRightString        = mstrpdpi
     plot(k+1,2)                = gsn_csm_contour_map_polar(wks,eofpdpi,res1)

   end do 

   ; panel plot only resources
    pres1                            = True                ; mods desired
    pres1@gsnFrame                   = True               ; save panel until both ready
   ;pres1@gsnPanelBottom             = 0.45    ; space for label bar
   ;pres1@gsnPanelTop                = 0.85     ; only panel on lower half of page
   ;pres1@gsnPanelXF                 = (/0.07,0.57/)
    pres1@gsnPanelYWhiteSpacePercent = 8
   gsn_panel(wks,(/plot(0,0),plot(1,0),plot(2,0),plot(3,0),\
                   plot(0,1),plot(1,1),plot(2,1),plot(3,1),\
                   plot(0,2),plot(1,2),plot(2,2),plot(3,2)/),(/3,4/),pres1)
 
end do 

end

