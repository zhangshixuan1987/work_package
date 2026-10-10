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

  left_strs = (/"WC (piControl)",\
                "WC (4xCO2)",\
                "CRYOnoismf (piControl)",\
                "CRYOnoismf (4xCO2)",\
                "CRYO (piControl)",\
                "CRYO (4xCO2)"/)

  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)

  time_tag  = (/"0051-0150"/)
  time_str  = (/"Year 51-150"/)
  time_prd  = (/"0051-0100", "0101-0150"/)
  nprd      = dimsizes(time_prd)

  yyyy      = ispan(51,150,1)
  yyyy!0    = "time"
  yyyy&time = yyyy

  cong_str  = (/"HC mean intensity (900-200hPa)"/)
  ncong = dimsizes(cong_str)

  latS              = -90.0
  latN              =  90.0
  YLValues          = ispan(toint(latS),toint(latN),15)
  YLMinorValues     = ispan(toint(latS),toint(latN), 5)
  YLLabels          = (/"90S","75S","60S","45S","30S","15S", "0",\
                        "15N","30N","45N","60N","75N","90N"/)
  Xstring           = "Time (year)"
  Ystring           = "Latitude"

  d2r       = get_d2r("float")
  seasons   = (/"ANN","DJF","JJA"/)
  nseas     = dimsizes(seasons)
  minbin    = (/-100, -180, -180/)
  maxbin    = (/ 100,  180,  180/)
  xbins     = (/  10,   20,  20/)

  varnam    = "SF"
  varfac    = 1e9
  varunt    = "10~S~9~N~ kg s~S~-1~N~" 

 ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1 

    print("working on season "+ seasons(i))

    do j = 0,nexps-1,1 

     fl        = data_dir+"/"+exp_name(j)+"_HC_analysis_Monthly_"+time_tag+".nc"
     f         = addfile(fl,"r")

     lat       = f->lat
     rad       = 4.0*atan(1.0)/180.0
     wgty      = cos(lat*rad)

    ;vertical average over 200-900hPa
     hc = f->SF
    ;printVarSummary(hc)

     if (seasons(i).eq."ANN") then
       hcm  = month_to_annual(hc,1)
     else
       hcm  = month_to_season(hc,seasons(i))
     end if
     hcm    = hcm / varfac
    ;printMinMax(hcm,False)

     hcm!0     = "time" 
     hcm&time  = yyyy
    ;printVarSummary(hcm)

     htmp   = hcm(lat|:,time|:,lev|:)
     hczm   = dim_avg_n_Wrap(htmp(:,:,{200:900}),2)
     hczm1  = dim_avg_n_Wrap(htmp(:,:,{400:600}),2)
    ;printVarSummary(hczm)

     tmpx1  = htmp({-40:0},:,{200:900})
     tmpx2  = htmp({0:40},:,{200:900})
     latx1  = lat({-40:0})
     latx2  = lat({0:40})
     nmin   = dim_minind(tmpx1,0)
     nmax   = dim_maxind(tmpx2,0)
    ;printVarSummary(nmax)
     latins = new ((/2,dimsizes(nmin(:,0))/),double)
     hcins  = new ((/2,dimsizes(nmin(:,0))/),double)
     do ii  = 0, dimsizes(nmax(:,0)) -1
       latins(0,ii) = avg(latx2(nmax(ii,:)))
       hcins(0,ii)  = avg(tmpx2(nmax(ii,:),ii,:))
     end do

     do ii  = 0, dimsizes(nmin(:,0)) -1 
       latins(1,ii) = avg(latx1(nmin(ii,:)))
       hcins(1,ii)  = avg(tmpx1(nmin(ii,:),ii,:))
     end do

     delete([/hcm,htmp,tmpx1,tmpx2,latx1,latx2,nmax,nmin/])
    ;printVarSummary(hczm)
    ;printMinMax(hcm,False)

    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Start to process plot 
    ;;;;;;;;;;;;;;;;;;;;;;;;;;
    if (j.eq.0) then 
      wtype           = "png"
      wtype@wkWidth   = 2500
      wtype@wkHeight  = 2500
      fignam          = "HC_mean_intensity_2d_"+seasons(i)+"_"+time_tag
      wks             = gsn_open_wks(wtype,fignam)
      plot = new(nexps,graphic)
      tmp1 = new(nexps,graphic)
      tmp2 = new(nexps,graphic)
    end if 

   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling

    res@vpHeightF             = 0.55        ; Changes the aspect ratio
    res@vpWidthF              = 0.90

    FontHeightF = 0.024
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
    res@tmXBValues               = ispan(min(yyyy)-1,max(yyyy),10)
    res@tmXBLabels               = tostring(res@tmXBValues)
    res@tmXBMinorValues          = ispan(min(yyyy)-1,max(yyyy),5)

    res@tmYLMode                 = "Explicit"
    res@tmYLValues               = YLValues
    res@tmYLLabels               = YLLabels
    res@tmYLMinorValues          = YLMinorValues

    res@tiXAxisString                = Xstring
    res@tiYAxisString                = "" ;Ystring

    res@tiMainString                 = ""
    res@gsnLeftString                = left_strs(j)
    res@gsnRightString               = "Stream function ("+varunt+")"
    res@gsnLeftStringOrthogonalPosF  = 0.02
    res@gsnRightStringOrthogonalPosF = 0.02

    res0 = res
    res0@cnLineLabelFontHeightF       = FontHeightF
    res0@lbLabelFontHeightF           = FontHeightF*0.92
    res0@lbTitleFontHeightF           = FontHeightF*0.92
    res0@lbLabelBarOn                 = True
    res0@lbOrientation                = "Vertical"
    res0@lbTitleString                = ""
    res0@lbTitlePosition              = "Right"     ; title location
    res0@lbTitleDirection             = "Across"    ; letter angle
    res0@lbTitleAngleF                = 90.         ; title angle
    res0@lbBoxLineColor               = "black"
   ;res0@lbBoxLinesOn                 = False
   ;res0@pmLabelBarParallelPosF       = 0.61
    res0@pmLabelBarOrthogonalPosF     = 0.0
    res0@pmLabelBarWidthF             = 0.15
    res0@pmLabelBarHeightF            = 0.58

    res0@cnFillOn                     = True                ; color on
    res0@cnFillPalette                = "BlueYellowRed" ;"MPL_RdYlBu" 
    res0@cnLinesOn                    = True                 ; no contour lines
    res0@cnLineLabelsOn               = True
    res0@cnMissingValFillColor        = "White" ;"gray76"               ; black = 1
    res0@cnLineColor                  = "black" ;"orangered";"Magenta";"Blue3"
    res0@cnInfoLabelOn                = False
    res0@cnLineLabelInterval          = 2
    res0@cnLineLabelFontColor         = res0@cnLineColor
    res0@cnLineLabelPlacementMode     = "constant"    ; choose constant label method
   ;res0@cnLevelSelectionMode         = "ExplicitLevels"       ; manual levels
    res0@cnLevelSelectionMode         = "ManualLevels" ; manual levels

    res0@cnMinLevelValF               = minbin(i)
    res0@cnMaxLevelValF               = maxbin(i)
    res0@cnLevelSpacingF              = xbins(i)
    res0@gsnContourZeroLineThicknessF = 10
    res0@gsnContourNegLineDashPattern = 5
    res0@gsnContourPosLineDashPattern = 0
    res0@cnLineThicknessF             = 1

    plot(j)  = gsn_csm_lat_time(wks, hczm({-45:45},:), res0) 

    res1                             = res
    res1@cnFillOn                    = False                ; color on
    res1@tiMainString                = ""
    res1@gsnLeftString               = ""
    res1@gsnRightString              = ""
    res1@cnLineColor                 = "Magenta"
    res1@cnLineLabelsOn              = False
    res1@cnLevelSelectionMode        = "ExplicitLevels"       ; manual levels
    res1@cnLevels                    = (/0.0/)
    res1@cnInfoLabelOn               = False
    res1@gsnContourZeroLineThicknessF= 10
    res1@cnLineThicknessF            = 5
    tmp1(j)    = gsn_csm_lat_time(wks, hczm1({-45:45},:), res1)
    overlay(plot(j),tmp1(j))

    res2                             = res
    res2@tiMainString                = ""
    res2@gsnLeftString               = ""
    res2@gsnRightString              = ""
    res2@xyMonoDashPattern           = False
    res2@xyDashPatterns              = (/2,2/)                  ; Make curves all solid
    res2@xyMarkLineMode              = "Lines"                  ; Markers *and* lines
    res2@xyMarkers                   = (/6,11,16/)              ; 3 different markers
    res2@xyMarkerColors             := (/"blue","red","green"/) ; 3 different colors
    res2@xyLineThicknesses           = (/  5.0,   5.0/) *2.0         ; make second line thicker
    res2@xyLineColors                = (/"Cyan","Cyan"/)      ; change line color
    tmp2(j)  = gsn_csm_xy (wks,yyyy,latins,res2) ; create plot
    overlay(plot(j),tmp2(j))
    delete([/res,res0,res1,res2/])
    delete([/rad,lat,wgty,hc,hczm,hczm1,latins,hcins/])
   end do 

    ; panel plot only resources
    pres1                 = True               ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/nexps/2,2/),pres1)      ; create first panel
    delete([/pres1,plot,tmp1,tmp2/])

  end do      

end
