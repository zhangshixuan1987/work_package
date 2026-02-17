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

  data_dir = "./"
  exp_name = (/"ERA5","CLIM_30yr","SSP245_30yr","SSP585_30yr"/)
  labList  = (/"ERA5","AMIP","AMIP+SSP245","AMIP+SSP585"/)
  nexps    = dimsizes(exp_name)
  ltrend   = new(nexps,string)

  season   = (/"DJF"/)
  ystr     = 1979
  yend     = 2014
  time_tag = ystr+"-"+yend
  region   = (/"Polar_Cap"/)

 ;declear the season or annual mean data to be used;;;;;
  do k = 0,nexps-1
  
   fl   = data_dir+"/"+exp_name(k)+"_HC_index_MON_"+time_tag+".nc"
   f    = addfile(fl,"r")

   time = ispan(ystr,yend,1)*1.0
   tmp  = month_to_seasonN(f->hces,season)
   hc  = tmp(0,:)
  ;hc  = dim_standardize_n_Wrap(hc,0,0)
   delete([/tmp/])

  ;3-yr smoothing
   wgt = (/0.25, 0.50, 0.25/)
   hc_smt = wgt_runave_Wrap(hc,wgt,0)
   hc_reg = hc

   ind0 = dim_minind(abs(time-1999),0)
   ind1 = ind(time.ge.1978.and.time.le.2001)
   ind2 = ind(time.ge.2000.and.time.le.2015)

   vmn  = (/avg(hc(ind1)),avg(hc(ind2)),avg(hc)/)

   jrc1 = regline(time(ind1),hc_smt(ind1))
   jp1  = trend_manken(hc(ind1), False, 0)
  ;hc_reg(0:ind0) = jrc1*(time(0:ind0)-jrc1@xave) + jrc1@yave
  ;hc_reg(0:ind0) = jp1(1)*(time(0:ind0)-jrc1@xave) + jrc1@yave
   
   jrc2 = regline(time(ind2),hc_smt(ind2))
   jp2  = trend_manken(hc(ind2), False, 0)
  ;hc_reg(ind0+1:) = jrc2*(time(ind0+1:)-jrc2@xave) + jrc2@yave
  ;hc_reg(ind0+1:) = jp2(1)*(time(ind0+1:)-jrc2@xave) + jrc2@yave
  ;hc_reg(ind0+1:) = hc_reg(ind0)  
   delete([/ind1,ind2,f,fl,f/])

  ;linear regression
   hc_reg = hc
   hc_reg = 0.01
   n       = 4
   jp3     = lspoly(time,hc, 1, n)
   do ii = 0,n-1,1
     hc_reg = hc_reg + jp3(ii)*time^(ii) 
   end do 
   hc_reg = hc@_FillValue 

   pltdat = new((/3,dimsizes(time)/),typeof(hc))
   pltdat(0,:) = (/hc/)
   pltdat(1,:) = (/hc_smt/)
   pltdat(2,:) = (/hc_reg/)
   trYReverse  = False
   minY        = -38
   maxY        = -32
   rstr        = season
   lstr        = "Hadley cell edge (SH)"
   ylab        = "Hadley cell edge (degree)"

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = "HC_loc_trend_"+season+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)

     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/16,29,26,18/)
     lthick = (/1.0,1.0,5.0,3.0/) *1.2
     ldash  = (/0,0,0,1/)
     opcity = (/0.5,0.5,1.0,1.0/)
     dum1   = new((/nexps/),graphic)
     dum2   = new((/nexps/),graphic)
     txid   = new((/nexps/),graphic)

   end if

   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.035
    res@tiMainFontHeightF     = FontHeightF
    res@tmYLLabelFontHeightF  = FontHeightF
    res@tmXBLabelFontHeightF  = FontHeightF
    res@gsnStringFontHeightF  = FontHeightF

    bthick = 3.0
    res@tmBorderThicknessF    = bthick
    res@tmXBMajorThicknessF   = bthick
    res@tmXBMinorThicknessF   = bthick
    res@tmYLMajorThicknessF   = bthick
    res@tmYLMinorThicknessF   = bthick
    res@tmYRMajorThicknessF   = bthick
    res@tmYRMinorThicknessF   = bthick

    res@vpHeightF             = 0.80        ; Changes the aspect ratio
    res@vpWidthF              = 0.85
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@pmLegendDisplayMode    = "Never"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.90                   ; move units right
    res@pmLegendOrthogonalPosF = -0.45                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.12                  ; height of legend.
    res@lgLabelFontHeightF     = .018                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    res@gsnLeftString           = ""
    res@gsnRightString          = ""
    res@tiMainString            = ""
    res@tiYAxisString           = "" 
    res@tiXAxisString           = ""

    res@trXMinF                 = 1978 ;floor(min(tsm03))
    res@trXMaxF                 = 2015 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************
    res1 = res
    res1@gsnYRefLine             = (/vmn(2)/)
    res1@gsnYRefLineColors       = colors(k)
    res1@gsnYRefLineThicknesses  = (/2.0/)
    res1@gsnYRefLineDashPatterns = (/2/)

    res1@gsnXRefLine             = (/2000/);(/2000,2000,2000,2000/)
    res1@gsnXRefLineColors       = "Black";colors(0:nexps-1) 
    res1@gsnXRefLineThicknesses  = (/2.0/)
    res1@gsnXRefLineDashPatterns = (/2/)

    res1@trYReverse              = trYReverse
    res1@trYMinF                 = minY
    res1@trYMaxF                 = maxY
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = (/colors(k),colors(k),colors(k)/)
    res1@xyLineThicknesses       = (/lthick(1),lthick(2),lthick(3)/)  ; line thicknesses
    res1@xyDashPatterns          = (/ldash(1),ldash(2),ldash(3)/)      ; make all lines solid
    if(k.eq.0) then
      res1@tiYAxisString         = ylab
      res1@tiXAxisString         = "Year"
      res1@gsnLeftString         = lstr 
      res1@gsnRightString        = rstr 
      res1@tiMainString          = ""
    else
      res1@gsnLeftString         = ""
      res1@gsnRightString        = ""
      res1@tiMainString          = ""
      res1@tiYAxisString         = ""
      res1@tiXAxisString         = ""
    end if

    dum1(k) = gsn_csm_xy (wks, time, pltdat(:,:), res1)

    if(k.eq.0) then
      plot = dum1(k)
    else
      overlay(plot,dum1(k))
    end if

    gres = True
    gres@YPosPercent          = 94.
    gres@XPosPercent          = 60.
    gres@Position             = "Top"
    gres@ItemSpacePercent     = 6
    gres@LineLabelWhiteSpacePercent = 2

    lineres = True
    lineres@lgLineColors       = colors(0:nexps-1) 
    lineres@lgLineOpacity      = opcity(1) 
    lineres@lgLineThicknesses  = lthick(1)*1.5
    lineres@lgDashIndexes      = ldash(1)  
    lineres@LineLengthPercent  = 6

    textres                    = True
    textres@lgPerimOn          = True
    textres@lgPerimFillColor   = "white"
    textres@lgLabels           = labList 
    textres@lgLabelFontHeights = FontHeightF*0.8
    textres@lgLabelColors      = colors(0:nexps-1) 
    textres@lgLabelOpacity     = opcity(1)         

    plot = simple_legend(wks,plot,gres,lineres,textres)
    delete([/res,res1,gres,lineres,textres/])
    delete([/time,pltdat,hc,hc_reg,hc_smt/])

   end do


   ;panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/1,1/),pres1)          ; create first panel

end

