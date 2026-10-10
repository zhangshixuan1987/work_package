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

  regions    = (/"Antarctic"/)
  latS       = -90.0
  latN       = -65.0
  nregs      = dimsizes(regions)

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/ozone_vint"
  exp_name   = (/"v1.LR.historical", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  ltrend     = new(nexps,string)

  labels    = (/"E3SMv1", "E3SMv2"/)
  nlabs     = dimsizes(labels)

  time_tag  = (/"1979-2014"/)
  seasons   = (/"SON"/)
  nseas     = dimsizes(seasons)
  maxlev    = 11

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1 

  do k = 0,nexps-1
  
   fl      = data_dir+"/"+regions+"_"+exp_name(k)+"_TCO_MON_"+time_tag+".nc"
   f       = addfile(fl,"r")

   tmp  = month_to_season(f->time,seasons(i))
   tmp@units = f->time@units
   time = yyyymm_to_yyyyfrac(cd_calendar(tmp,-1),0.0)
   delete(tmp)

   ;;ensemble mean
   tmp  = dim_avg_n_Wrap(f->TCO,0)
   tco  = month_to_season (tmp,seasons(i))
   delete(tmp)

   ;3-yr smoothing
   tco_smt = runave_Wrap(tco,3,0)
   tco_reg = tco

   ind0 = dim_minind(abs(time-1999),0)
   ind1 = ind(time.ge.1978.and.time.le.2001)
   ind2 = ind(time.ge.2000.and.time.le.2015)

   jrc1 = regline(time(ind1),tco_smt(ind1))
   jp1  = trend_manken(tco(ind1), False, 0)
  ;tco_reg(0:ind0) = jrc1*(time(0:ind0)-jrc1@xave) + jrc1@yave
  ;tco_reg(0:ind0) = jp1(1)*(time(0:ind0)-jrc1@xave) + jrc1@yave
   
   jrc2 = regline(time(ind2),tco_smt(ind2))
   jp2  = trend_manken(tco(ind2), False, 0)
  ;tco_reg(ind0+1:) = jrc2*(time(ind0+1:)-jrc2@xave) + jrc2@yave
  ;tco_reg(ind0+1:) = jp2(1)*(time(ind0+1:)-jrc2@xave) + jrc2@yave
  ;tco_reg(ind0+1:) = tco_reg(ind0)  
   delete([/ind1,ind2,f,fl,f/])

   if(k.eq.1) then
     n       = 4
     jp3     = lspoly(time,tco, 1, n)
     tco_reg = tco 
     tco_reg = 0.01
     do ii = 0,n-1,1
       tco_reg = tco_reg + jp3(ii)*time^(ii) 
     end do 
   else 
     tco_reg = tco 
     tco_reg = tco@_FillValue 
   end if 

   ;if(k.eq.0) then 
   ;  tco_reg = runave_Wrap(tco_reg,11,11)
   ;else
   ;  tco_reg = runave_Wrap(tco_reg,5,5)
   ;end if 

   pltdat = new((/3,dimsizes(time)/),typeof(tco))
   pltdat(0,:) = (/tco/)
   pltdat(1,:) = (/tco_smt/)
   pltdat(2,:) = (/tco_reg/)
   trYReverse  = False
   minY        = 200
   maxY        = 360
   rstr        = "SON" ; (1979-2014)"
   lstr        = "Total Column Ozone"
   ylab        = "TCO (DU)"

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = "TCO_ind_trend_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/18,29,10/)
     lthick = (/2.0,5.0,3.0/)*1.2
     ldash  = (/0,0,1/)
     opcity = (/1,1,1/)
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
 
    FontHeightF = 0.040
    res@tiMainFontHeightF     = FontHeightF
    res@tmYLLabelFontHeightF  = FontHeightF
    res@tmXBLabelFontHeightF  = FontHeightF
    res@gsnStringFontHeightF  = FontHeightF

    bthick = 3.5
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

    res@trYMinF                 = -6.0
    res@trYMaxF                 =  4.0
    res@trXMinF                 = 1978 ;floor(min(tsm03))
    res@trXMaxF                 = 2015 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************
    res1 = res
    res1@gsnYRefLine             = -999
    res1@gsnYRefLineColors       = "Black" ;colors(k)
    res1@gsnYRefLineThicknesses  = (/2.0/)
    res1@gsnYRefLineDashPatterns = (/2/)

    res1@gsnXRefLine             = 2000
    res1@gsnXRefLineColors       = "Black" 
    res1@gsnXRefLineThicknesses  = (/2.0/)
    res1@gsnXRefLineDashPatterns = (/2/)

    res1@trYReverse              = trYReverse
    res1@trYMinF                 = minY
    res1@trYMaxF                 = maxY
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = (/colors(k),colors(k),colors(k)/)
    res1@xyLineThicknesses       = (/lthick(0),lthick(1),lthick(2)/)      ; line thicknesses
    res1@xyDashPatterns          = (/ldash(0),ldash(1),ldash(2)/)      ; make all lines solid
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
    gres@YPosPercent          = 90.
    gres@XPosPercent          = 60.
    gres@Position             = "Top"
    gres@ItemSpacePercent     = 8.
    gres@LineLabelWhiteSpacePercent = 2

    lineres = True
    lineres@lgLineColors       = (/colors(0),colors(1)/)
    lineres@lgLineOpacity      = (/opcity(1),opcity(1)/)
    lineres@lgLineThicknesses  = (/lthick(1),lthick(1)/)
    lineres@lgDashIndexes      = (/ldash(1),ldash(1)/)
    lineres@LineLengthPercent  = 8

    textres                    = True
    textres@lgLabels           = labels
    textres@lgLabelFontHeights = FontHeightF
    textres@lgLabelColors      = (/colors(0),colors(1)/) 
    textres@lgLabelOpacity     = (/opcity(1),opcity(1)/)

    plot = simple_legend(wks,plot,gres,lineres,textres)

   end do

    delete([/time,pltdat,tco,tco_reg,tco_smt/])

   ;panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/1,1/),pres1)          ; create first panel

 end do   

end

