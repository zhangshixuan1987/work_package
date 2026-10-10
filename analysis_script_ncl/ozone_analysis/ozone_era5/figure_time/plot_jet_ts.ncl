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

  data_dir  = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/jet_ind"
  exp_name  = (/"v1.LR.historical", "v2.LR.historical"/)
  nexps     = dimsizes(exp_name)
  ltrend    = new(nexps+2,string)

  labels    = (/"E3SMv1", "E3SMv2"/)
  nlabs     = dimsizes(labels)

  time_tag  = (/"1979-2014"/)
  seasons   = (/"DJF"/)
  nseas     = dimsizes(seasons)
  maxlev    = 11

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1 

  do k = 0,nexps-1
  
   data_dir = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/jet_ind"
   fl      = systemfunc ("ls -1 " + data_dir+"/"+exp_name(k)+"_Jet_ind_U850_MON_*_"+time_tag+".nc")
   f       = addfiles(fl,"r")
   ListSetType (f, "join")       ; join

   tmp  = month_to_season(f[0]->time,seasons(i))
   tmp@units = f[0]->time@units
   time = yyyymm_to_yyyyfrac(cd_calendar(tmp,-1),0.0)
   delete(tmp)

   ;;ensemble mean
   tmp  = dim_avg_n_Wrap(f[:]->jeto_pos,0)
   jet  = month_to_season (tmp,seasons(i)) 
   delete(tmp)

   ;3-yr smoothing
   jet_smt = runave_Wrap(jet,3,0)
   jet_reg = jet

   ind0 = ind(time.eq.2000)
   ind1 = ind(time.ge.1979.and.time.le.2000)
   ind2 = ind(time.ge.2000.and.time.le.2015)

   jrc1 = regline(time(ind1),jet(ind1))
   jp1  = trend_manken(jet(ind1), False, 0)
   jet_reg(ind1) = jrc1*(time(ind1)-jrc1@xave) + jrc1@yave
  ;jet_reg(ind1) = jp1(1)*(time(ind1)-jrc1@xave) + jrc1@yave

   jrc2 = regline(time(ind2),jet(ind2))
   jp2  = trend_manken(jet(ind2), False, 0)
   jet_reg(ind2) = jrc2*(time(ind2)-jrc2@xave) + jrc2@yave
  ;jet_reg(ind2) = jp2(1)*(time(ind2)-jrc2@xave) + jrc2@yave
   jet_reg(ind2) = jet_reg(ind0)  
   delete([/ind1,ind2,f,fl,f/])

   pltdat = new((/3,dimsizes(time)/),typeof(jet))
   pltdat(0,:) = (/jet/)
   pltdat(1,:) = (/jet_smt/)
   pltdat(2,:) = (/jet_reg/)
   trYReverse  = False
   minY        = -54.0
   maxY        = -44.0
   rstr        = "Jet position"
   ylab        = "Jet position (~S~o~N~)"

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = "Jet_pos_trend_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,27,10,19,31/)
    ;colors = (/1,18,2,10,27/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/1,1,1,1,1/)
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
 
    FontHeightF = 0.032
    res@tiMainFontHeightF     = FontHeightF
    res@tmYLLabelFontHeightF  = FontHeightF
    res@tmXBLabelFontHeightF  = FontHeightF
    res@gsnStringFontHeightF  = FontHeightF

    if(wtype.eq."png") then
     bthick = 3.0
    else
     bthick = 1.0
    end if

    res@tmBorderThicknessF    = bthick
    res@tmXBMajorThicknessF   = bthick
    res@tmXBMinorThicknessF   = bthick
    res@tmYLMajorThicknessF   = bthick
    res@tmYLMinorThicknessF   = bthick
    res@tmYRMajorThicknessF   = bthick
    res@tmYRMinorThicknessF   = bthick

    res@vpHeightF             = 0.80        ; Changes the aspect ratio
    res@vpWidthF              = 0.80
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@pmLegendDisplayMode    = "Never"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.40                   ; move units right
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
    res@trXMinF                 = 1979 ;floor(min(tsm03))
    res@trXMaxF                 = 2010 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************
    res1 = res
    res1@gsnYRefLine             = -999
    res1@gsnYRefLineColors       = colors(k)
    res1@gsnYRefLineThicknesses  = (/1.0/)
    res1@gsnYRefLineDashPatterns = (/2/)

    res1@gsnXRefLine             = 2000
    res1@gsnXRefLineColors       = "Grey" 
    res1@gsnXRefLineThicknesses  = (/2.0/)
    res1@gsnXRefLineDashPatterns = (/1/)

    res1@trYReverse              = trYReverse
    res1@trYMinF                 = minY
    res1@trYMaxF                 = maxY
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = (/colors(k),colors(k),colors(k)/)
    res1@xyLineThicknesses       = (/1.0,lthick(k),lthick(k)/)      ; line thicknesses
    res1@xyDashPatterns          = (/ldash(k),ldash(k),1/)      ; make all lines solid
    if(k.eq.0) then
      res1@tiYAxisString         = ylab
      res1@tiXAxisString         = "Year"
      res1@gsnLeftString         = "DJF (1979-2014)"
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
    gres@XPosPercent          = 15.
    gres@Position             = "Top"
    gres@ItemSpacePercent     = 8.
    gres@LineLabelWhiteSpacePercent = 2

    lineres = True
    lineres@lgLineColors       = colors(0:nexps-1)
    lineres@lgLineOpacity      = opcity(0:nexps-1)
    lineres@lgLineThicknesses  = lthick(0:nexps-1)
    lineres@lgDashIndexes      = ldash(0:nexps-1)
    lineres@LineLengthPercent  = 8

    textres                    = True
    textres@lgLabels           = labels
    textres@lgLabelFontHeights = FontHeightF 
    textres@lgLabelColors      = colors(0:nexps-1)
    textres@lgLabelOpacity     = opcity(0:nexps-1)

    plot = simple_legend(wks,plot,gres,lineres,textres)

   end do

    delete([/time,pltdat,jet,jet_reg,jet_smt/])

   ;panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/1,1/),pres1)          ; create first panel

 end do   

end

