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

  data_dir  = "../fig_data/"
  exp_name  = (/"ERA_20C", \
                "NOAA_20C"/)
  nexps     = dimsizes(exp_name)

  time_tag  = (/"1901-1930"/)
  time_str  = (/"Year 1901-1930"/)
  time_unt  = "hours since 1901-01-01 00:00:00.0"

  yyyy      = ispan(1901,1930,1)
  yyyy!0    = "time"
  yyyy&time = yyyy

  left_str  = (/"HC intensity","HC extents", "HC width"/) 
  cong_str  = (/"NH", "SH"/)
  ncong     = dimsizes(cong_str)

  d2r       = get_d2r("float")
  season    = (/"ANN","DJF","JJA"/)
  nsea      = dimsizes(season)

  varnam    = (/"HC Intensity", "HC extent", "HC width"/)
  varunt    = (/"10~S~9~N~ kg s~S~-1~N~", \
                "degree" , \
                "degree" /)

  Xstring   = "Time (year)"

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nsea-1,1 

    print("working on season "+ season(i))

    if (season(i).eq."DJF") then
      minbin   = (/ 120., -60.,  26., -40,  36., -30/)
      maxbin   = (/ 240., -20.,  36., -30,  56., -20/)
    else if (season(i).eq."JJA") then
      minbin   = (/  10.,  -220.,  24., -36,  8., -60/)
      maxbin   = (/  30.,  -160.,  52., -26, 32., -36/)
    else 
      minbin   = (/  70., -130.,  28., -40, 26., -40/)
      maxbin   = (/ 130.,  -90.,  40., -28, 36., -30/)
    end if 
    end if 

    do j = 0,nexps-1,1 

     fl   = data_dir+"/"+exp_name(j)+"_HC_index_Monthly_"+time_tag+".nc"
     f    = addfile(fl,"r")

     time  = f->time 
     ntime = dimsizes(time)/12

     if(j.eq.0) then 
       hcin = new((/nexps,ntime/),double)
       hcis = new((/nexps,ntime/),double)
       hcei = new((/nexps,ntime/),double)
       hcen = new((/nexps,ntime/),double)
       hces = new((/nexps,ntime/),double)
       hcwn = new((/nexps,ntime/),double)
       hcws = new((/nexps,ntime/),double)
       hcwa = new((/nexps,ntime/),double)
     end if 

     if (season(i).eq."ANN") then
       hcin(j,:) = month_to_annual(f->hcin,1)
       hcis(j,:) = month_to_annual(f->hcis,1)
       hcei(j,:) = month_to_annual(f->hcei,1)
       hcen(j,:) = month_to_annual(f->hcen,1)
       hces(j,:) = month_to_annual(f->hces,1)
       hcwn(j,:) = month_to_annual(f->hcwn,1)
       hcws(j,:) = month_to_annual(f->hcws,1)
       hcwa(j,:) = month_to_annual(f->hcwa,1)
     else
       hcin(j,:) = month_to_season(f->hcin,season(i))
       hcis(j,:) = month_to_season(f->hcis,season(i))
       hcei(j,:) = month_to_season(f->hcei,season(i))
       hcen(j,:) = month_to_season(f->hcen,season(i))
       hces(j,:) = month_to_season(f->hces,season(i))
       hcwn(j,:) = month_to_season(f->hcwn,season(i))
       hcws(j,:) = month_to_season(f->hcws,season(i))
       hcwa(j,:) = month_to_season(f->hcwa,season(i))
     end if

    end do 

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 1200
    wtype@wkHeight  = 1200
    fignam          = "HC_intensity_analysis_"+season(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)

    plot = new(ncong*3,graphic)

    gsn_define_colormap(wks,"t2m_29lev")
    LineColors      = (/9,18/)
    LineThicknesses = (/20.,20./)/4      ; line thicknesses
    DashPatterns    = (/0.,0./)      ; make all lines solid
 
   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling

    FontHeightF = 0.028
    res@tiMainFontThicknessF     = 3.0
    res@tiMainFontHeightF        = FontHeightF*1.2
    res@tmYLLabelFontHeightF     = FontHeightF
    res@tmXBLabelFontHeightF     = FontHeightF
    res@tiXAxisFontHeightF       = FontHeightF
    res@tiYAxisFontHeightF       = FontHeightF
    res@gsnStringFontHeightF     = FontHeightF*1.1

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

    res@vpHeightF             = 0.50        ; Changes the aspect ratio
    res@vpWidthF              = 0.90
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@tiXAxisString         = Xstring
    res@tiYAxisString         = ""

    res@tmXTOn                = False       ; have tick marks

    res@pmLegendDisplayMode    = "Always"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of
    res@pmLegendParallelPosF   = 0.20                  ; move units right
    res@pmLegendOrthogonalPosF = -0.54                  ; more neg = down
    res@pmLegendWidthF         = 0.12                  ; Change width and
    res@pmLegendHeightF        = 0.15                  ; height of legend.
    res@lgLabelFontHeightF     = FontHeightF           ; change font height
    res@lgPerimOn              = False                 ; no box around

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

    rts@gsnYRefLineColors       := LineColors
    rts@gsnYRefLineDashPatterns := (/5.,5./)
    rts@gsnYRefLineThicknesses  := (/2.,2./)

    ;;NH and SH branch
    do ix = 0, 1

      if(ix.eq.0)then 
        hci = hcin
        hce = hcen
        hcw = hcwn
      else
        hci = hcis
        hce = hces
        hcw = hcws
      end if 

      miny = minbin(ix::2)
      maxy = maxbin(ix::2)

      if(ix.eq.0)then 
        rts@trYReverse             = False
      else
        rts@trYReverse             = True
      end if 

      rts@gsnRightString           = cong_str(ix)
      rts@xyExplicitLegendLabels   = "  "+exp_name

      vref = dim_avg_n(hci,1)
      rts0 = rts
      rts0@gsnLeftString            = left_str(0)+" (" + season(i) +")"
      rts0@gsnYRefLine             := (/vref/)
      rts0@tiYAxisString            = "("+varunt(0)+")"
      rts0@trYMinF                  = miny(0)
      rts0@trYMaxF                  = maxy(0)
      print(yyyy)
      print(hci)
      tmp1                          = gsn_csm_xy (wks,yyyy,hci,rts0)

      vref = dim_avg_n(hce,1)
      rts1 = rts
      rts1@gsnLeftString            = left_str(1)+" (" + season(i) +")"
      rts1@gsnYRefLine             := (/vref/)
      rts1@tiYAxisString            = "("+varunt(1)+")"
      rts1@trYMinF                  = miny(1)
      rts1@trYMaxF                  = maxy(1)
      tmp2                          = gsn_csm_xy (wks,yyyy,hce,rts1)

      vref = dim_avg_n(hcw,1)
      rts2 = rts
      rts2@gsnLeftString            = left_str(2)+" (" + season(i) +")"
      rts2@gsnYRefLine             := (/vref/)
      rts2@tiYAxisString            = "("+varunt(2)+")"
      rts2@trYMinF                  = miny(2)
      rts2@trYMaxF                  = maxy(2)
      tmp3                           = gsn_csm_xy (wks,yyyy,hcw,rts2)

     plot(ix::2) = (/tmp1,tmp2,tmp3/)
     delete([/hci,hce,hcw/])
    end do 
    delete([/hcin,hcis,hcei,hcen,hces,hcwn,hcws,hcwa/])

   ; panel plot only resources
   pres1                 = True               ; mods desired
   pres1@gsnFrame        = True               ; save panel until both ready
   gsn_panel(wks,plot,(/3,2/),pres1)      ; create first panel
  
 end do

end

