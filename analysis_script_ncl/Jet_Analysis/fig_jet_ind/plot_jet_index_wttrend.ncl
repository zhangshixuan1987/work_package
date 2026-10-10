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

  data_dir   = "../data"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  labels     = (/"ERA5","v2 AMIP", "v2 Coupled"/)

  l_ins_trend = new(nexps,string)
  l_pos_trend = new(nexps,string)


 ;define regions for eof analysis
  mode       = ("Jet")
  vars       = (/"U850"/)
  facmod     = (/1.0/)
  varunt     = (/"m s~S~-1~N~"/)

  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

   do k = 0,nexps-1

    print("working on season "+ seasons(i))
    fl = data_dir+"/"+exp_name(k)+"_Jet_ind_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
    f  = addfile(fl,"r")

    time  = yyyymm_to_yyyyfrac(cd_calendar(f->time,-1),0.0)

    if(seasons(i).eq."DJF")
     time = time + 1
    end if 

    vins  = dim_avg_n_Wrap(f->jeto_ins,0)
    vpos  = dim_avg_n_Wrap(f->jeto_pos,0)

    rc1   = regline(time(:),vins(:)) 
    rc2   = regline(time(:),vpos(:))
    vireg = rc1*time + rc1@yintercept
    vpreg = rc2*time + rc2@yintercept
    rc1   = rc1 * 10.0
    rc2   = rc2 * 10.0

    rcx1 = regCoef(time,vins)             
    tval = rcx1@tval
    df   = rcx1@nptxy -2 
    b    = tval
    b    = 0.5        
    prb1 = betainc(df/(df+tval^2),df/2.0,b)

    rcx2 = regCoef(time,vpos)
    tval = rcx2@tval
    df   = rcx2@nptxy -2
    b    = tval
    b    = 0.5  
    prb2 = betainc(df/(df+tval^2),df/2.0,b)
    print(prb1 + " " +prb2)
    l_ins_trend(k) = sprintf("%.2f",avg(rc1)) 
    l_pos_trend(k) = sprintf("%.2f",avg(rc2)) 
    
    if(prb1.le.0.05) then 
      l_ins_trend(k) = l_ins_trend(k) +"~S~**~N~"
    else if (prb1.le.0.1) then 
      l_ins_trend(k) = l_ins_trend(k) +"~S~*~N~"
    else
      l_ins_trend(k) = l_ins_trend(k) +""
    end if 
    end if 

    if(prb2.le.0.05) then
      l_pos_trend(k) = l_pos_trend(k) +"~S~**~N~"
    else if (prb2.le.0.1) then
      l_pos_trend(k) = l_pos_trend(k) +"~S~*~N~"
    else
      l_pos_trend(k) = l_pos_trend(k) +""
    end if
    end if 

    delete([/rc1,rc2,rcx1,tval,df,b,rcx2,prb1,prb2/])
          
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = mode+"_trend_"+vars+"_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,10,27,19,31/)
    ;colors = (/1,18,2,10,27/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/1,1,1,1,1/)
     dum1   = new((/2,nexps/),graphic)
     dum2   = new((/2,nexps/),graphic)
     txid   = new((/2,nexps/),graphic)                  
     plot   = new(2,graphic)
   end if 

   ;============================================================
   ; PLOTS
   ;============================================================
    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.024
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

    res@vpHeightF             = 0.40        ; Changes the aspect ratio
    res@vpWidthF              = 0.90
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

    res@trXMinF                 = 1978 ;floor(min(tsm03))
    res@trXMaxF                 = 2015 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************
   do ii = 0,1

    if (ii.eq.0) then
      pltdat = new((/2,dimsizes(time)/),typeof(vins))
      pltdat(0,:) = (/vins/)
      pltdat(1,:) = (/vireg/)
      minY = 10.0
      maxY = 18.0
      ystr = vars + " (m s~S~-1~N~)"      ; y-axis label
      ylab = "Jet intensity (" + vars+ ")"
    else
      pltdat = new((/2,dimsizes(time)/),typeof(vins))
      pltdat(0,:) = (/vpos/)
      pltdat(1,:) = (/vpreg/)
      minY = -54.0
      maxY = -44.0
      ystr = "Latitude (degree)"      ; y-axis label
      ylab = "Jet position (" + vars+ ")"
    end if 

    res1 = res
    res1@gsnYRefLine             = -9999;avg(pltdat(0,:))              ; reference line
    res1@gsnYRefLineColors       = colors(k)
    res1@gsnYRefLineThicknesses  = (/2.0/)
    res1@gsnYRefLineDashPatterns = (/5/)
    res1@xyMonoLineColor         = False             ; want colored lines

    res1@trYMinF                 = minY
    res1@trYMaxF                 = maxY
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = (/colors(k),colors(k)/)
    res1@xyLineThicknesses       = (/lthick(k),lthick(k)/)      ; line thicknesses
    res1@xyDashPatterns          = (/ldash(k),1/)      ; make all lines solid
    if(k.eq.0) then
      res1@tiYAxisString          = ystr
      res1@tiXAxisString          = "Year"
      res1@gsnLeftString          = ylab
      res1@gsnRightString         = seasons(i)
      res1@tiMainString           = ""
    else
      res1@gsnLeftString          = ""
      res1@gsnRightString         = ""
      res1@tiMainString           = ""
      res1@tiYAxisString          = ""
      res1@tiXAxisString          = ""
    end if

    dum1(ii,k)     = gsn_csm_xy (wks, time, pltdat(:,:), res1)

    if(k.eq.0) then
      plot(ii) = dum1(ii,k)
    else
      overlay(plot(ii),dum1(ii,k))
    end if

    gres = True
    gres@YPosPercent          = 95.     
    gres@XPosPercent          = 15.    
    gres@Position             = "Top"        
    gres@ItemSpacePercent     = 8.          
    gres@LineLabelWhiteSpacePercent = 2 

    lineres = True
    lineres@lgLineColors       = colors(0:nexps-1)
    lineres@lgLineOpacity      = opcity(0:nexps-1)
    lineres@lgLineThicknesses  = lthick(0:nexps-1)
    lineres@lgDashIndexes      = ldash(0:nexps-1)
    lineres@LineLengthPercent  = 5      

    textres                    = True
    textres@lgLabels           = labels
    textres@lgLabelFontHeights = FontHeightF*0.8
    textres@lgLabelColors      = colors(0:nexps-1)
    textres@lgLabelOpacity     = opcity(0:nexps-1)

   ;plot(ii) = simple_legend(wks,plot(ii),gres,lineres,textres)

   ;;add text;;
   txres                  = True
  ;txres@txFont           = "helvetica-bold"
   txres@txFontHeightF    = FontHeightF*0.8
   txres@txFontThicknessF = 4.0
   txres@txJust           = "BottomLeft"
   if(ii.eq.0) then 
     ltrend = l_ins_trend(k)
     xpos   = 1982
     ypos   = 17.0
   else
     ltrend = l_pos_trend(k)
     xpos   = 1982
     ypos   = -45.0
   end if 
  
   txres@txFontColor   = colors(k)
   txtstr = labels(k) + " ("+ltrend +")"
   txid(ii,k) = gsn_add_text (wks,plot(ii),txtstr,xpos+k*8,ypos ,txres)

   end do

   delete([/time,pltdat,vireg,vins,vpreg,vpos/])
   end do 

  ;panel plot only resources
   pres1                 = True                ; mods desired
   pres1@gsnFrame        = True               ; save panel until both ready
   gsn_panel(wks,plot,(/2,1/),pres1)          ; create first panel

 end do   

end

