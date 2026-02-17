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

  obs_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/sam_ind"
  obs_name   = (/"M03","GW99"/)
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/sam_ind"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  ltrend     = new(nexps+2,string)

  labels     = (/"M03","GW99","ERA5", "v2 AMIP", "v2 Coupled"/)
  nlabs      = dimsizes(labels)

 ;define regions for eof analysis
  mode       = ("SAM")
  ieof       = 1
  vars       = (/"PSL"/)
  facmod     = (/1.0/)
  varunt     = (/"Pa / std dev"/)

  time_tag   = (/"1979-2014"/)
  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1 

  print("working on season "+ seasons(i))

  fl0   = data_dir+"/Marshall_SAM_AOI_"+seasons(i)+"_1957-2020.nc"
  f0    = addfile(fl0,"r")
  tsm03 = yyyymm_to_yyyyfrac(cd_calendar(f0->time,-1),0.0)
  dsm03 = f0->SAM_AOI 
  rc    = regline(tsm03,dsm03)
  rsm03 = rc*tsm03 + rc@yintercept
  indx  = ind(tsm03.ge.1981.and.tsm03.le.2010)
  x     = tsm03(indx)
  y     = dsm03(indx)
  ltrend(0) = sprintf("%.2f",regline(x,y)*10.0)
  delete([/x,y,rc,indx,f0,fl0/])

  fl1   = data_dir+"/ERA5_SAM_AOI_"+seasons(i)+"_1979-2014.nc"
  f1    = addfile(fl1,"r")
  time  = yyyymm_to_yyyyfrac(cd_calendar(f1->time,-1),0.0)
  indx  = ind(time.ge.1981.and.time.le.2010)
  tgw99 = time(indx)
  sgw99 = f1->SAM_AOI(:,indx) 
  nens  = dimsizes(sgw99(:,0))
  rc    = new(nens,typeof(sgw99))
  do ie = 0,nens-1,1
    rc(ie) = regline(tgw99(:),sgw99(ie,:)) * 10.0
  end do 
  ltrend(1) = sprintf("%.2f",avg(rc)) + " ~F34~1~N~~F~ " + sprintf("%.2f",stddev(rc)) 
  delete([/time,rc,indx,f1,fl1/])

  do k = 0,nexps-1
  
   fl   = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars+"_"+seasons(i)+"_"+time_tag+".nc"
   f    = addfile(fl,"r")
   time = yyyymm_to_yyyyfrac(cd_calendar(f->time,-1),0.0)
   indx = ind(time.ge.1981.and.time.le.2010)
   teof = time(indx)
   seof = f->eof_ts(:,ieof-1,indx) 
   nens  = dimsizes(seof(:,0))
   rc    = new(nens,typeof(seof))
   do ie = 0,nens-1,1
     rc(ie) = regline(teof(:),seof(ie,:)) * 10.0
   end do

   ltrend(k+2) = sprintf("%.2f",avg(rc)) + " ~F34~1~N~~F~ " + sprintf("%.2f",stddev(rc))
   delete([/time,rc,indx,f,fl/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
     wtype           = "pdf"
    ;wtype@wkWidth   = 2500
    ;wtype@wkHeight  = 2500
     fignam          = mode+"_EOF_ts_"+vars+"_"+seasons(i)+"_"+time_tag
     wks             = gsn_open_wks(wtype,fignam)
     gsn_define_colormap(wks,"t2m_29lev")
     colors = (/1,27,10,19,31/)
     colors = (/1,18,2,10,27/)
     lthick = (/3.0,3.0,3.0,5.0,5.0/)
     ldash  = (/0,0,0,0,0/)
     opcity = (/1,1,1,1,1/)
     dum1   = new(nexps+2,graphic)
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

    res@gsnYRefLine             = (/0./)              ; reference line
    res@gsnYRefLineColors       = (/"Grey25"/)
    res@gsnYRefLineThicknesses  = (/0.5/)
    res@gsnYRefLineDashPatterns = (/0/)
    res@xyMonoLineColor         = False             ; want colored lines

    res@trYMinF                 = -6.0
    res@trYMaxF                 =  4.0
    res@trXMinF                 = 1957 ;floor(min(tsm03))
    res@trXMaxF                 = 2020 ;ceil(max(tsm03))

   ;*******************************************
   ; First plot
   ;*******************************************
   if (k.eq.0) then 
      res1  = res
      res1@xyLineColors             = (/colors(0),colors(0)/)
      res1@xyLineThicknesses        = (/lthick(0),2.0/)      ; line thicknesses
      res1@xyDashPatterns           = (/ldash(0),1.0/)      ; make all lines solid
      pltdat = new((/2,dimsizes(rsm03)/),typeof(rsm03))
      pltdat(0,:) = (/dsm03/)
      pltdat(1,:) = (/rsm03/)
      plot = gsn_csm_xy (wks, tsm03, pltdat, res1)
      delete(pltdat)
      
      res2 = res
      res2@tiYAxisString          = "Standardized index"      ; y-axis label
      res2@tiXAxisString          = "Year"
      res2@gsnLeftString          = "SAM index (" + vars+ ")"
      res2@gsnRightString         = seasons(i) 
      res2@xyLineColors           = colors(1)
      res2@xyDashPatterns         = ldash(1)      ; make all lines solid
      res2@xyLineThicknesses      = lthick(1)      ; line thicknesses
      pltdat = dim_avg_n_Wrap(sgw99,0)
      dum1(0) =  gsn_csm_xy (wks, tgw99, pltdat, res2)
      overlay(plot,dum1(0))
      delete(pltdat)

    end if 

    ;*******************************************
    ; second plot
    ;*******************************************
    rts                        = res
    rts@gsnLeftString          = ""
    rts@gsnRightString         = ""
    rts@tiMainString           = "" 
    rts@tiYAxisString          = ""
    rts@tiXAxisString          = ""
    rts@xyLineColors           = colors(k+2)
    rts@xyDashPatterns         = ldash(k+2)      ; make all lines solid
    rts@xyLineThicknesses      = lthick(k+2)      ; line thicknesses

    pltdat    = dim_avg_n_Wrap(seof,0)
    dum1(k+2) = gsn_csm_xy (wks,teof,pltdat,rts)           ; create plot
    overlay(plot,dum1(k+2))
    delete(pltdat)
    delete([/teof,seof/])
   end do 

   gres = True
   gres@YPosPercent          = 95.     
   gres@XPosPercent          = 15.    
   gres@Position             = "Top"        
   gres@ItemSpacePercent     = 8.          
   gres@LineLabelWhiteSpacePercent = 2 

   lineres = True
   lineres@lgLineColors       = colors
   lineres@lgLineOpacity      = opcity
   lineres@lgLineThicknesses  = lthick
   lineres@lgDashIndexes      = ldash
   lineres@LineLengthPercent  = 5      

   textres                    = True
   textres@lgLabels           = labels
   textres@lgLabelFontHeights = FontHeightF*0.8
   textres@lgLabelColors      = colors
   textres@lgLabelOpacity     = opcity

  ;plot = simple_legend(wks,plot,gres,lineres,textres)

   delete([/tsm03,tgw99,dsm03,rsm03,sgw99/])

   ;;add text;;
   txid = new(nlabs+1,graphic)                  ; Nine text strings
   txres                  = True             
  ;txres@txFont           = "helvetica-bold" 
   txres@txFontHeightF    = FontHeightF*0.8
   txres@txFontThicknessF = 4.0
   txres@txJust        = "BottomLeft"
   do m = 0,nlabs,1
     if(m.eq.0) then
       txres@txFontColor   = colors(m)
       txtstr = "slope (1979-2014):"
     else 
       txres@txFontColor   = colors(m-1)
       txtstr = labels(m-1) + " ("+ltrend(m-1) +")"
     end if 
     if(m.le.1) then 
       txid(m) = gsn_add_text (wks,plot,txtstr, 1958.0, -4.5 - m ,txres)
     else if(m.le.3) then
       txid(m) = gsn_add_text (wks,plot,txtstr, 1978.0, -4.5 - (m-2) ,txres)
     else 
       txid(m) = gsn_add_text (wks,plot,txtstr, 1998.0, -4.5 - (m-4),txres)
     end if 
     end if 
   end do 

   ;panel plot only resources
    pres1                 = True                ; mods desired
    pres1@gsnFrame        = True               ; save panel until both ready
    gsn_panel(wks,plot,(/1,1/),pres1)          ; create first panel

 end do   

end

