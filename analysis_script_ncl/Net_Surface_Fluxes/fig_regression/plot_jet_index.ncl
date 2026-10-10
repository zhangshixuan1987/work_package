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

  data_dir  = "/lcrc/group/acme/ac.szhang/acme_scratch/data/"
  exp_name  = (/"v2.LR.piControl", \
                "v2.LR.abrupt-4xCO2", \
                "v2.LR.piControl.CRYO", \                 
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
 
  left_str  = (/"piControl (WC)", \
                "4xCO2 (WC)",\
                "piControl (CRYO)", \
                "4xCO2 (CRYO)"/) 

  cong_str  = (/"Water cycle Config.", "Cryosphere Config."/)
 ;define regions for eof analysis
  mode       = ("Jet")
  vars       = (/"U"/)
  facmod     = (/1.0 /)
  varunt     = (/"m s~S~-1~N~"/)

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
    fl        = "Jet_ind_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
    f         = addfile(fl,"r")

    rad       = 4.0*atan(1.0)/180.0
    utc_date  = cd_calendar(f->time,-1)
    yyyy      = yyyymm_to_yyyyfrac(utc_date, 0.0)

    plotdat1  = new ((/nexps,dimsizes(utc_date)/),double)
    plotdat2  = new ((/nexps,dimsizes(utc_date)/),double)

    do j = 0,nexps-1,1
      varstr  = "jeto_ins_"+exp_name(j)
      plotdat1(j,:)  = f->$varstr$
      varstr  = "jeto_pos_"+exp_name(j)
      plotdat2(j,:)  = f->$varstr$
    end do 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 1200
    wtype@wkHeight  = 1200
    fignam          = "Jet_intensity_"+vars+"_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)

    gsn_define_colormap(wks,"t2m_29lev")
    plot = new(4,graphic)

   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.018
    res@tiMainFontHeightF     = FontHeightF*1.5
    res@tmYLLabelFontHeightF  = FontHeightF
    res@tmXBLabelFontHeightF  = FontHeightF
    res@gsnStringFontHeightF  = FontHeightF*1.2

    if(wtype.eq."png") then
     lthick = 3.0
    else
     lthick = 1.0
    end if

    res@tmBorderThicknessF    = lthick
    res@tmXBMajorThicknessF   = lthick
    res@tmXBMinorThicknessF   = lthick
    res@tmYLMajorThicknessF   = lthick
    res@tmYLMinorThicknessF   = lthick
    res@tmYRMajorThicknessF   = lthick
    res@tmYRMinorThicknessF   = lthick

    res@vpHeightF             = 0.40        ; Changes the aspect ratio
    res@vpWidthF              = 0.85
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@tiYAxisString          = "Jet Index (standard)"          ; y-axis label
    res@tiXAxisString          = "Model time (year)"

    res@pmLegendDisplayMode    = "Always"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.20                   ; move units right
    res@pmLegendOrthogonalPosF = -0.45                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.10                  ; height of legend.
    res@lgLabelFontHeightF     = .016                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    ;*******************************************
    ; second plot
    ;*******************************************
    rts                       = res
    rts@gsnXYBarChart         = False ;True           ; create bar chart 
    rts@gsnLeftString         = ""
    rts@gsnRightString        = ""
    rts@tiMainString          = "" 
    rts@xyMonoLineColor       = False             ; want colored lines
    rts@xyLineColors          = (/1,27/)
    rts@xyLineThicknesses     = (/20.,20./)/10      ; line thicknesses
    rts@xyDashPatterns	      = (/0.,0./)      ; make all lines solid

    rts@tiYAxisString         = "Jet intensity ("+varunt+")"          ; y-axis label
    rts@trYMinF               = 10.0
    rts@trYMaxF               = 20.0
    rts@tmYLPrecision         = 3
    xyarr = wgt_runave_n_Wrap(plotdat1,(/1.,3,5,6,5,3,1/), 0, 1)
    vref = dim_avg_n(plotdat1,1)
    rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
    rts@gsnYRefLineDashPatterns := (/0.,0./)
    rts@gsnYRefLineThicknesses  := (/20.,20./)/4
    rts@gsnYRefLine             := (/vref(0),vref(1)/)
    rts@gsnLeftString            = cong_str(0)
    rts@gsnRightString           = "Jet intensity" 
    rts@xyExplicitLegendLabels   = "  "+left_str(0:1)
    plot(0) = gsn_csm_xy (wks,yyyy,xyarr(0:1,:),rts)

    rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
    rts@gsnYRefLineDashPatterns := (/0.,0./)
    rts@gsnYRefLineThicknesses  := (/20.,20./)/4
    rts@gsnYRefLine             := (/vref(2),vref(3)/)
    rts@gsnLeftString            = cong_str(1)
    rts@gsnRightString           = "Jet intensity"
    rts@xyExplicitLegendLabels   = "  "+left_str(2:3)
    plot(1) = gsn_csm_xy (wks,yyyy,xyarr(2:3,:),rts)

    rts@tiYAxisString         = "Jet position (degree)"          ; y-axis label
    rts@trYMinF               = -56.0
    rts@trYMaxF               = -38.0
    rts@tmYLPrecision         = 3
    xyarr = wgt_runave_n_Wrap(plotdat2,(/1.,3,5,6,5,3,1/), 0, 1)
    vref = dim_avg_n(plotdat2,1)
    rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
    rts@gsnYRefLineDashPatterns := (/0.,0./)
    rts@gsnYRefLineThicknesses  := (/20.,20./)/4
    rts@gsnYRefLine             := (/vref(0),vref(1)/)
    rts@gsnLeftString            = cong_str(0)
    rts@gsnRightString           = "Jet position"
    rts@xyExplicitLegendLabels   = "  "+left_str(0:1)
    plot(2) = gsn_csm_xy (wks,yyyy,xyarr(0:1,:),rts)

    rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
    rts@gsnYRefLineDashPatterns := (/0.,0./)
    rts@gsnYRefLineThicknesses  := (/20.,20./)/4
    rts@gsnYRefLine             := (/vref(2),vref(3)/)
    rts@gsnLeftString            = cong_str(1)
    rts@gsnRightString           = "Jet position"
    rts@xyExplicitLegendLabels   = "  "+left_str(2:3)
    plot(3) = gsn_csm_xy (wks,yyyy,xyarr(2:3,:),rts)

    ; panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/2,2/),pres1)          ; create first panel

    delete([/plotdat1,plotdat2,yyyy,xyarr,utc_date/])

  end do   

end do 

end

