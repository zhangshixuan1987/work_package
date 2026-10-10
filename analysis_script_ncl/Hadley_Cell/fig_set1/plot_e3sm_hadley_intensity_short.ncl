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
                "v2.LR.abrupt-4xCO2", \
                "v2.LR.piControl.CRYOnoismf", \
                "v2.LR.CRYOnoismf.abrupt-4xCO2", \
                "v2.LR.piControl.CRYO", \
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(exp_name)

  left_str  = (/"piControl (WC)", \
                "4xCO2 (WC)",\
                "piControl (WC+ICE_SHELF)", \
                "4xCO2 (WC+ICE_SHELF)",\
                "piControl (CRYO)", \
                "4xCO2 (CRYO)"/)

  cong_str  = (/"WC Config.", "WC+ICE_SHELF config.", "CRYO Config."/)
  ncong = dimsizes(cong_str)

  d2r       = get_d2r("float")
  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)

  minbin1   = (/  5,  10,   0/)
  maxbin1   = (/ 15,  25,  10/)
  minbin2   = (/-15, -10, -25/)
  maxbin2   = (/ -5,   0, -10/)

  varnam    = "SF"
  varfac    = 1e10
  varunt    = "10~S~10~N~ kg s~S~-1~N~" 

do k = 0, 2

  if(k.eq.0)then 
    time_tag  = (/"0001-0050"/)
    time_str  = (/"Year 1-50"/)
  end if 

  if(k.eq.1)then
    time_tag  = (/"0051-0100"/)
    time_str  = (/"Year 51-100"/)
  end if

  if(k.eq.2)then
    time_tag  = (/"0101-0150"/)
    time_str  = (/"Year 101-150"/)
  end if

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    do j = 0,nexps-1,1 

     fl        = data_dir+"/"+exp_name(j)+"_HC_analysis_ANN_"+time_tag+".nc"
     f         = addfile(fl,"r")

     lat       = f->lat
     rad       = 4.0*atan(1.0)/180.0
     wgty      = cos(lat*rad)
     utc_date  = cd_calendar(f->time(::12),-1)
     yyyy      = yyyymm_to_yyyyfrac(utc_date, 0.0)

     if (j.eq.0) then 
       plotdat1 = new ((/nexps,dimsizes(yyyy)/),double)
       plotdat2 = new ((/nexps,dimsizes(yyyy)/),double)
     end if 

    ;HC intensity is defined as the vertically averaged maximum/minimum value of SF 
    ;between 900 and 200 hPa for NH and SH, respectively. 
     tmp1 = dim_avg_n_Wrap(dim_max_n_Wrap(f->SF(:,{200:900},{0:90}),2),1)
     tmp2 = dim_avg_n_Wrap(dim_min_n_Wrap(f->SF(:,{200:900},{-90:0}),2),1)

     if(seasons(i).eq."ANN") then 
       plotdat1(j,:)  = month_to_annual(tmp1,1)
       plotdat2(j,:)  = month_to_annual(tmp2,1)
     else
       plotdat1(j,:)  = month_to_season(tmp1,seasons(i))
       plotdat2(j,:)  = month_to_season(tmp2,seasons(i))
     end if 
     delete([/tmp1,tmp2/])
    end do

    plotdat1 = plotdat1 / varfac
    plotdat2 = plotdat2 / varfac

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 1200
    wtype@wkHeight  = 1200
    fignam          = "HC_intensity_analysis_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)
    gsn_define_colormap(wks,"t2m_29lev")
    plot = new(ncong*2,graphic)

    LineColors      = (/1,27,9,18/)
    LineThicknesses = (/20.,20.,20.,10./)/4      ; line thicknesses
    DashPatterns    = (/0.,0.,0.,0./)      ; make all lines solid
    
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

    res@tiXAxisString          = "Model time (year)"

    res@pmLegendDisplayMode    = "Always"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of
    res@pmLegendParallelPosF   = 0.25                   ; move units right
    res@pmLegendOrthogonalPosF = -0.55                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.12                  ; height of legend.
    res@lgLabelFontHeightF     = .018                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    ;*******************************************
    ; first plot
    ;*******************************************
    rts                       = res
    rts@gsnXYBarChart         = False ;True           ; create bar chart
   ;rts@gsnAboveYRefLineColor = "red"           ; above ref line fill red
   ;rts@gsnBelowYRefLineColor = "blue"          ; below ref line fill blue

    rts@gsnLeftString         = ""
    rts@gsnRightString        = ""
    rts@tiMainString          = ""
    rts@xyMonoLineColor       = False             ; want colored lines
    rts@xyLineColors          = LineColors
    rts@xyLineThicknesses     = LineThicknesses
    rts@xyDashPatterns        = DashPatterns

    do kk = 0,ncong-1,1

      ;;NH branch
      xyarr  = plotdat1(kk*2:kk*2+1,:)
      vref   = dim_avg_n(xyarr,1)
      rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
      rts@gsnYRefLineDashPatterns := (/5.,5./)
      rts@gsnYRefLineThicknesses  := (/2.,2./)
      rts@gsnYRefLine             := vref
      rts@gsnLeftString            = cong_str(kk)
      rts@gsnRightString           = "HC intensity (NH branch)"
      rts@xyExplicitLegendLabels   = "  "+left_str(kk*2:kk*2+1)
      rts@trYMinF                  = minbin1(i)
      rts@trYMaxF                  = maxbin1(i)
      rts@tiYAxisString            = "~F33~Y~N~~S~*~N~~F~~B~N~N~ (x"+varunt+")"          ; y-axis label
      plot(kk)                     = gsn_csm_xy (wks,yyyy,xyarr,rts)

      ;;SH branch
      xyarr  = plotdat2(kk*2:kk*2+1,:)
      vref   = dim_avg_n(xyarr,1)
      rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
      rts@gsnYRefLineDashPatterns := (/5.,5./)
      rts@gsnYRefLineThicknesses  := (/2.,2./)
      rts@gsnYRefLine             := vref
      rts@gsnLeftString            = cong_str(kk)
      rts@gsnRightString           = "HC intensity (SH branch)"
      rts@xyExplicitLegendLabels   = "  "+left_str(kk*2:kk*2+1)
      rts@trYMinF                  = minbin2(i)
      rts@trYMaxF                  = maxbin2(i)
      rts@tiYAxisString            = "~F33~Y~N~~S~*~N~~F~~B~S~N~ (x"+varunt+")"          ; y-axis label
      plot(kk+ncong)               = gsn_csm_xy (wks,yyyy,xyarr,rts)
      delete([/xyarr,vref/])
    end do

    ; panel plot only resources
    pres1                 = True               ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/2,ncong/),pres1)      ; create first panel
    delete([/plotdat1,plotdat2,yyyy,utc_date/])
  
  end do      

end do 

end

