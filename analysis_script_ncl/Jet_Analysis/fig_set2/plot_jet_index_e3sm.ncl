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
  obs_name  = (/"NOAA_20C"/) ;"ERA_20C"/) ;,"ERA5"/)
  exp_name  = (/"v2.LR.piControl", \
                "v2.LR.abrupt-4xCO2", \
                "v2.LR.piControl.CRYOnoismf", \
                "v2.LR.CRYOnoismf.abrupt-4xCO2", \
                "v2.LR.piControl.CRYO", \
                "v2.LR.abrupt-4xCO2.CRYO"/)
  nexps     = dimsizes(exp_name)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
 
  left_str  = (/"WC (piControl)",\
                "WC (4xCO2)",\
                "CRYOnoismf (piControl)",\
                "CRYOnoismf (4xCO2)",\
                "CRYO (piControl)",\
                "CRYO (4xCO2)"/)

  cong_str  = (/"WC Config.", "CRYOnoismf config.", "CRYO Config."/)
  ncongs    = dimsizes(cong_str)

 ;define regions for eof analysis
  mode      = ("Jet")
  vars      = (/"U"/)
  facmod    = (/1.0 /)
  varunt    = (/"m s~S~-1~N~"/)

  time_obs  = (/"1901-1930","1981-2010"/)
  time_tag  = (/"0051-0150"/)
  time_str  = (/"Year 51-150"/)

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    ;read and process the reanalysis;;
    flopi     = "../fig_data/Jet_ind_"+vars+"_Monthly_"+time_obs(0)+".nc"
    fopi      = addfile(flopi,"r")
    flopd     = "../fig_data/Jet_ind_"+vars+"_Monthly_"+time_obs(1)+".nc"
    fopd      = addfile(flopd,"r")

    utc_date1 = cd_calendar(fopd->time,0)
    monthx1   = tointeger(utc_date1(:,1))
    if(seasons(i).eq."ANN")then
      indx1 = ind ((monthx1.ge.1).and.(monthx1.le.12))
    else if (seasons(i).eq."DJF")then
      indx1 = ind ((monthx1.eq.12).or.(monthx1.eq.1).or.(monthx1.eq.2))
    else if (seasons(i).eq."JJA") then
      indx1 = ind ((monthx1.eq.6).or.(monthx1.eq.7).or.(monthx1.eq.8))
    else if (seasons(i).eq."MAM") then
      indx1 = ind ((monthx1.eq.3).or.(monthx1.eq.4).or.(monthx1.eq.5))
    else if (seasons(i).eq."SON") then
      indx1 = ind ((monthx1.eq.9).or.(monthx1.eq.10).or.(monthx1.eq.11))
    else
      print("User did not specify the season, stop" )
      exit
    end if
    end if
    end if
    end if
    end if

    if(seasons(i).eq."ANN")then
      yobs   = ispan(0,dimsizes(indx1)-1,1)/12.0
    else
      yobs   = ispan(0,dimsizes(indx1)-1,1)*3.0/12.0
    end if

    plotobs1  = new ((/2,dimsizes(indx1)/),double)
    plotobs2  = new ((/2,dimsizes(indx1)/),double)

    varstr    = "jeto_ins_"+obs_name(0)
    plotobs1(0,:) = fopi->$varstr$(indx1)
    plotobs1(1,:) = fopd->$varstr$(indx1)
    ;varstr    = "jeto_ins_"+obs_name(1)
    ;plotobs1(2,:) = fopi->$varstr$(indx1)
    ;plotobs1(3,:) = fopd->$varstr$(indx1)

    varstr    = "jeto_pos_"+obs_name(0)
    plotobs2(0,:) = fopi->$varstr$(indx1)
    plotobs2(1,:) = fopd->$varstr$(indx1)
    ;varstr    = "jeto_pos_"+obs_name(1)
    ;plotobs2(2,:) = fopi->$varstr$(indx1)
    ;plotobs2(3,:) = fopd->$varstr$(indx1)

    vrefo1 = dim_avg_n(plotobs1,1)
    vrefo2 = dim_avg_n(plotobs2,1)

    mstro1 = sprintf("%4.1f", vrefo1)
    mstro2 = sprintf("%4.1f", vrefo2)

    delete([/indx1,utc_date1,monthx1/])
    
    ;read and process the e3smv2;;
    fl        = "../fig_data/Jet_ind_"+vars+"_Monthly_"+time_tag+".nc"
    f         = addfile(fl,"r")

    utc_date1 = cd_calendar(f->time,0)
    monthx1   = tointeger(utc_date1(:,1))
    if(seasons(i).eq."ANN")then
      indx1 = ind ((monthx1.ge.1).and.(monthx1.le.12))
    else if (seasons(i).eq."DJF")then
      indx1 = ind ((monthx1.eq.12).or.(monthx1.eq.1).or.(monthx1.eq.2))
    else if (seasons(i).eq."JJA") then
      indx1 = ind ((monthx1.eq.6).or.(monthx1.eq.7).or.(monthx1.eq.8))
    else if (seasons(i).eq."MAM") then
      indx1 = ind ((monthx1.eq.3).or.(monthx1.eq.4).or.(monthx1.eq.5))
    else if (seasons(i).eq."SON") then
      indx1 = ind ((monthx1.eq.9).or.(monthx1.eq.10).or.(monthx1.eq.11))
    else
      print("User did not specify the season, stop" )
      exit
    end if
    end if
    end if
    end if
    end if

    if(seasons(i).eq."ANN")then
      yyyy   = ispan(0,dimsizes(indx1)-1,1)/12.0
    else
      yyyy   = ispan(0,dimsizes(indx1)-1,1)*3.0/12.0
    end if

    plotdat1  = new ((/nexps,dimsizes(indx1)/),double)
    plotdat2  = new ((/nexps,dimsizes(indx1)/),double)

    do j = 0,nexps-1,1
      varstr  = "jeto_ins_"+exp_name(j)
      plotdat1(j,:)  = f->$varstr$(indx1)
      varstr  = "jeto_pos_"+exp_name(j)
      plotdat2(j,:)  = f->$varstr$(indx1)
    end do 
    vref1 = dim_avg_n(plotdat1,1)
    vref2 = dim_avg_n(plotdat2,1)

    mstr1 = sprintf("%4.1f", vref1)
    mstr2 = sprintf("%4.1f", vref2)

    delete([/indx1,utc_date1,monthx1/])

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fignam          = "Jet_intensity_"+vars+"_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)

    gsn_define_colormap(wks,"t2m_29lev")
    plot = new((ncongs+1)*2,graphic)

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
    res@pmLegendParallelPosF   = 0.25                   ; move units right
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

    do kk = 0, ncongs, 1

      rts@tiYAxisString         = "Jet intensity ("+varunt+")"          ; y-axis label
      rts@trYMinF               = 8.0
      rts@trYMaxF               = 20.0
      rts@tmYLPrecision         = 3
      rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
      rts@gsnYRefLineDashPatterns := (/0.,0./)
      rts@gsnYRefLineThicknesses  := (/20.,20./)/4

      if (kk.eq.0) then 
        rts@gsnLeftString            = "Obs.("+obs_name+")" 
        rts@gsnRightString           = "Jet intensity"
        rts@gsnYRefLine             := (/vrefo1(0),vrefo1(1)/)
        rts@xyExplicitLegendLabels   = "  "+time_obs + " (" + mstro1(0:1)+")"
        plot(kk) = gsn_csm_xy (wks,yobs,plotobs1,rts)
      else 
        kk0 = kk -1
        rts@gsnLeftString            = cong_str(kk0)
        rts@gsnRightString           = "Jet intensity" 
        rts@gsnYRefLine             := (/vref1(kk0*2),vref1(kk0*2+1)/)
        rts@xyExplicitLegendLabels   = "  "+left_str(kk0*2:kk0*2+1)+" (" + mstr1(kk0*2:kk0*2+1)+")"
        xyarr = wgt_runave_n_Wrap(plotdat1,(/1.,3,5,6,5,3,1/), 0, 1)
        plot(kk) = gsn_csm_xy (wks,yyyy,xyarr(kk0*2:kk0*2+1,:),rts)
      end if 

      rts@tiYAxisString         = "Jet position (degree)"          ; y-axis label
      rts@trYMinF               = -60.0
      rts@trYMaxF               = -36.0
      rts@tmYLPrecision         = 3
      rts@gsnYRefLineColors       := rts@xyLineColors(0:1)
      rts@gsnYRefLineDashPatterns := (/0.,0./)
      rts@gsnYRefLineThicknesses  := (/20.,20./)/4

      if (kk.eq.0) then
        rts@gsnLeftString            = "Obs.("+obs_name+")"
        rts@gsnRightString           = "Jet position"
        rts@gsnYRefLine             := (/vrefo2(0),vrefo2(1)/)
        rts@xyExplicitLegendLabels   = "  "+time_obs+" (" + mstro2(0:1)+")"
        plot(kk+ncongs+1) = gsn_csm_xy (wks,yobs,plotobs2,rts)
      else
        kk0 = kk -1
        rts@gsnLeftString            = cong_str(kk0)
        rts@gsnRightString           = "Jet position"
        rts@gsnYRefLine             := (/vref2(kk0*2),vref2(kk0*2+1)/)
        rts@xyExplicitLegendLabels   = "  "+left_str(kk0*2:kk0*2+1)+" (" + mstr2(kk0*2:kk0*2+1)+")"
        xyarr = wgt_runave_n_Wrap(plotdat2,(/1.,3,5,6,5,3,1/), 0, 1)
        plot(kk+ncongs+1) = gsn_csm_xy (wks,yyyy,xyarr(kk0*2:kk0*2+1,:),rts)
      end if 

    end do 

    delete([/plotobs1,plotobs2,plotdat1,plotdat2,yyyy,yobs,xyarr/])

    ; panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/2,4/),pres1)          ; create first panel

 end do   

end

