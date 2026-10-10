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

  data_dir   = "../data/hc_ind/"
  exp_name   = (/"ERA5", "v2.LR.historical","v1.LR.historical"/)
  labels     = (/"ERA5", "v2 Coupled",      "v1 Coupled"/)
  nexps      = dimsizes(exp_name)
  ltrend     = new(nexps,string)

  regions    = (/"Tropical"/)
  nregs      = dimsizes(regions)

 ;define regions for eof analysis
  varnam     = (/"HC"/)
  varnamtr   = (/"HC edge (SH)"/)
  facmod     = (/1.0/)
  varunt     = (/"Degree"/)
  thresh     = 220.0 ;unit: Du 

  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  Min1       = (/120,240,250,220,250,150/)
  Max1       = (/360,320,300,320,300,360/)

  Min2       = (/0,   0,  0,  0,  0, 0/)
  Max2       = (/20, 20, 20, 20, 20, 20/)

  Min3       = (/260,280,280,280,270,240/)
  Max3       = (/400,350,360,370,320,400/)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1 

  print("working on season "+ seasons(i))

  do k = 0,nexps-1
 
   fl  = systemfunc("ls " + data_dir+exp_name(k)+"_HC_index_MON"+"_en*_"+time_tag+".nc")
   f   = addfiles(fl,"r")
   ListSetType (f, "join")

   vdat0 = f[:]->hces
   vdat0!0 = "ens"
   vdat  = vdat0(time|:,ens|:)

   if(seasons(i).eq."MON") then
     var  = vdat
   else if(seasons(i).eq."ANN") then 
     var  = month_to_annual(vdat,0)
   else
     var  = month_to_annual(vdat,0)
     do ii = 0,dimsizes(vdat(0,:))-1
       tmp = month_to_seasonN(vdat(:,ii),seasons(i))
       var (:,ii) = (/tmp(0,:)/)
       delete(tmp)
     end do 
   end if 
   end if 

   nens  = dimsizes(var(0,:))
   if(k.eq.0) then 
     date   = cd_calendar(f[0]->time,0)
     ystr   = toint(min(date(:,0)))
     yend   = toint(max(date(:,0)))
     if(seasons(i).eq."MON") then
       time = yyyymm_to_yyyyfrac(cd_calendar(f[0]->time,-1),0.0)
     else if(seasons(i).eq."ANN") then
       time = ispan(ystr,yend,1)*1.0
     else
       time = ispan(ystr,yend,1)*1.0
     end if
     end if
     pltdat = new((/nexps,nens+1,dimsizes(var(:,0))/),double)
     delete([/date,ystr,yend/])
   end if 

   do ie = 0,nens,1
     if(ie.eq.0) then 
       pltdat(k,ie,:)  = dim_avg_n(var,1)
     else
       pltdat(k,ie,:)  = var(:,ie-1)
     end if 
   end do 

   delete([/var,vdat0,vdat,fl,f/])

   end do ;end of ensembles 

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
    wtype           = "png"
    wtype@wkWidth   = 1200
    wtype@wkHeight  = 1200
    fignam          = regions+"_"+varnam+"_"+seasons(i)+"_"+time_tag
    wks             = gsn_open_wks(wtype,fignam)
    gsn_define_colormap(wks,"t2m_29lev")
    colors = (/1,27,10,19,31/)
    lthick = (/3.0,3.0,3.0,5.0,5.0/)
    ldash  = (/0,0,0,0,0/)
    opcity = (/0.2,0.2,0.2,0.2,0.2/)
    dum1   = new(nexps,graphic)

    res                       = True         
    res@gsnDraw               = False        ; don't draw yet
    res@gsnFrame              = False        ; don't advance frame yet
    res@gsnScale              = True        ; force text scaling
 
    FontHeightF = 0.030
    res@tiMainFontHeightF     = FontHeightF*1.5
    res@tmYLLabelFontHeightF  = FontHeightF*1.5
    res@tmXBLabelFontHeightF  = FontHeightF*1.5
    res@gsnStringFontHeightF  = FontHeightF*1.5

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
    res@vpWidthF              = 0.90
    res@vpXF                  = 0.10        ; change start locations
    res@vpYF                  = 0.75        ; the plot

    res@pmLegendDisplayMode    = "Never"              ; turn on legend
    res@pmLegendSide           = "Top"                 ; Change location of 
    res@pmLegendParallelPosF   = 0.40                   ; move units right
    res@pmLegendOrthogonalPosF = -0.45                  ; more neg = down
    res@pmLegendWidthF         = 0.10                  ; Change width and
    res@pmLegendHeightF        = 0.12                  ; height of legend.
    res@lgLabelFontHeightF     = .024                   ; change font height
    res@lgPerimOn              = False                 ; no box around

    res@tmXTOn                 = False       ; have tick marks

    res@gsnLeftString           = ""
    res@gsnRightString          = ""
    res@tiMainString            = ""
    res@tiYAxisString           = "" 
    res@tiXAxisString           = ""

    res@gsnYRefLine             = (/-9999./)              ; reference line
    res@gsnYRefLineColors       = (/"Grey25"/)
    res@gsnYRefLineThicknesses  = (/0.5/)
    res@gsnYRefLineDashPatterns = (/0/)
    res@xyMonoLineColor         = False             ; want colored lines

    res@trXMinF                 = 1978 ;floor(min(tsm03))
    res@trXMaxF                 = 2015 ;ceil(max(tsm03))

    do k = 0,nexps-1,1

    ;*******************************************
    ; First plot
    ;*******************************************

    res1  = res
    res1@xyMonoLineColor         = False             ; want colored lines
    res1@xyLineColors            = new(nens+1,integer)
    res1@xyLineThicknesses       = new(nens+1,float)
    res1@xyDashPatterns          = new(nens+1,integer)
    res1@xyLineOpacities         = new(nens+1,float)
    res1@xyLineColors(:)         = colors(k)
    if(seasons(i).eq."MON") then 
      res1@xyLineThicknesses(:)    = lthick(k)*0.1      ; line thicknesses
    else
      res1@xyLineThicknesses(:)    = lthick(k)
    end if 
    res1@xyDashPatterns(:)       = ldash(k)      ; make all lines solid
    res1@xyLineOpacities(:)      = 1.0
    res1@xyLineThicknesses(1:)   = lthick(k)/2.0
    res1@xyLineOpacities(1:)     = 0.5

   ;*******************************************
   ; First plot
   ;*******************************************
    if(k.eq.0) then
      res1@trYMinF                =  min( pltdat) - abs(max( pltdat))/100.0
      res1@trYMaxF                =  max( pltdat) + abs(max( pltdat))/100.0
    end if
    if(k.eq.0) then
      res1@tiYAxisString          = "HC edge ("+ varunt+")"      ; y-axis label
      res1@tiXAxisString          = "Year"
      res1@gsnLeftString          = "Hadley cell edge (SH)"
      res1@gsnRightString         = seasons(i)
      res1@tiMainString           = ""
    else
      res1@gsnLeftString          = ""
      res1@gsnRightString         = ""
      res1@tiMainString           = ""
      res1@tiYAxisString          = ""
      res1@tiXAxisString          = ""
    end if

    dum1(k)                      = gsn_csm_xy (wks, time, pltdat(k,:,:), res1)

    ;res2 = res
    ;res2@gsnLeftString          = ""
    ;res2@gsnRightString         = ""
    ;res2@tiMainString           = ""
    ;res2@tiYAxisString          = ""
    ;res2@tiXAxisString          = ""
    ;res2@gsnXYFillColors        = colors(k)
    ;res2@gsnXYFillOpacities     = opcity(k)
    ;res2@xyLineColor            = -1          
    ;dum2(k) = gsn_csm_xy (wks, time, pltdat(k,1:,:), res2)
    ;overlay(dum1(k),dum2(k))

    if(k.eq.0) then
      plot = dum1(k)
    else
      overlay(plot,dum1(k))
    end if 
 
   delete([/res1/])

   end do 

   gres = True
   gres@YPosPercent          = 95.     
   gres@XPosPercent          = 65.    
   gres@Position             = "Top"        
   gres@ItemSpacePercent     = 8.          
   gres@LineLabelWhiteSpacePercent = 2 

   lineres = True
   lineres@lgLineColors       = colors(0:nexps-1)
  ;lineres@lgLineOpacity      = opcity(0:nexps-1)
   lineres@lgLineThicknesses  = lthick(0:nexps-1)
   lineres@lgDashIndexes      = ldash(0:nexps-1)
   lineres@LineLengthPercent  = 5      

   textres                    = True
   textres@lgLabels           = labels(0:nexps-1)
   textres@lgLabelFontHeights = FontHeightF*1.2
   textres@lgLabelColors      = colors(0:nexps-1)
  ;textres@lgLabelOpacity     = opcity(0:nexps-1)

   plot = simple_legend(wks,plot,gres,lineres,textres)

   ;;add text;;
   txid = new(3*(nexps+1),graphic)                  ; Nine text strings
   txres                  = True             
  ;txres@txFont           = "helvetica-bold" 
   txres@txFontHeightF    = FontHeightF*2.0
   txres@txFontThicknessF = 4.0
   txres@txJust           = "BottomLeft"

   do m = 0,nexps-1,1
     txres@txFontColor   = colors(m)
     txtstr = labels(m) + " ("+ltrend(m) +")"
     txid(m) = gsn_add_text (wks,plot,txtstr, 1982.0, -15-4*m ,txres)
   end do 
  
   ;panel plot only resources
   pres1                 = True                ; mods desired
   pres1@gsnFrame        = True               ; save panel until both ready
   gsn_panel(wks,plot,(/3,1/),pres1)          ; create first panel

   delete([/time,pltdat,res,dum1,plot,txid,txres,pres1,lineres,gres/])

 end do   

end

