;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Figure 5 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/sam_eof"
  time_tag   = (/"1979-2014"/)
  mode       = ("SAM")
  ieof       = 1

  vars       = (/"PSL","Z500"/)
  varunt     = (/"hPa", "m"/)
  nvars      = dimsizes(vars)
  miny       = (/-6, -60/)
  maxy       = (/ 3,  30/)

  fig_set    = (/"no_mbias", "no_mbias", "no_mbias", "wt_mbias", "wt_mbias"/)
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical", "v2.LR.amip", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  exp_str    = (/"ERA5", "v2 AMIP", "v2 Coupled", "v2 AMIP (bias)", "v2 Coupled (bias)"/)
  seasons    = (/"MON"/)
  nseas      = dimsizes(seasons)

  xtype      = "frequency" ; or "period"
  plot       = new(nexps,graphic)

  do j = 0, nseas-1,1

   do i = 0,nvars-1,1

    plotFileLabel = "fig_spectrum_anl_"+vars(i) +"_"+ seasons(j) 
    load  "./fig_setup.ncl"

   ;;prepare the data for plot;;
   do k = 0,nexps-1,1
    fl     = systemfunc("ls "+data_dir+"/"+fig_set(k)+"/"+exp_name(k)+"_EOF_PCS_"+vars(i)+"_"+seasons(j)+"_en"+"*_"+time_tag+".nc")
    f        = addfiles(fl,"r")
    ListSetType (f, "join")

    eof_ts = f[:]->eof_ts(:,ieof-1,:)

    delete([/fl,f/])
   ;************************************************
   ;calculate mean spectrum spectrum and lag1 auto cor
   ;************************************************
   nens = dimsizes(eof_ts(:,0))
   ntim = dimsizes(eof_ts(0,:))

   ;loop over each segment of length ntim
   spcavg = new ( ntim/2, typeof(eof_ts))
   spcavg = 0.0
   r1zsum = 0.0
   
   d   = 0    ; detrending opt: 0=>remove mean 1=>remove mean + detrend
   sm  = 5    ; <3, Do no smoothing. spcx contains raw spectra estimates (periodogram).
              ; >3, average utilizing modified Daniell smoothing, should be odd number
   pct = 0.10 ; percent taper: (0.0 <= pct <= 1.0) 0.10 common. 

   do n=0,nens-1
      sdof   = specx_anal(eof_ts(n,:),d,sm,pct) ; current segment spc
      spcavg = spcavg + sdof@spcx               ; sum spc of each segment
      r1     = sdof@xlag1                       ; extract segment lag-1
      r1zsum = r1zsum  + 0.5*log((1+r1)/(1-r1)) ; sum the Fischer Z
   end do

   r1z    = r1zsum/nens                         ; average r1z
   r1     = (exp(2*r1z)-1)/(exp(2*r1z)+1)       ; transform back

   spcavg = spcavg/nens                         ; average spectrum
  ;smooth the spectral power by a three-point running mean filter 
  ;in the frequency domain to reduce the uncertainty of energy peaks
   spcavg = runave_Wrap (spcavg, 3, 1)
   
   delete([/eof_ts/])

  ;************************************************
  ; Assign mean spectrum to data object
  ;************************************************
   df       = 2.0*nens                           ; deg of freedom
   df@spcx  = spcavg                             ; assign the mean spc
   df@frq   = sdof@frq
   if(xtype.eq."frequency") then 
     df@sca   = sdof@frq
   else 
     df@sca   = 1.0/sdof@frq / 12.0
   end if 

   df@xlag1 = r1

  ;************************************************
  ; calculates a red noise confidence interval.
  ; This function returns four curves:
  ; the calculated spectrum
  ; the "red noise" curve
  ; the curves indicating the upper and lower confidence bounds
  ;************************************************
   splt0 = specx_ci(df,  0.1, 0.90)                 ; confidence interval
   splt1 = specx_ci(df, 0.05, 0.95)
   nx    = dimsizes(splt1(:,0))
   ny    = dimsizes(splt1(0,:))
   splt  = new((/nx+1,ny/),typeof(splt0))
   splt(0:2,:) = splt1(0:2,:)
   splt(3,:)   = splt0(3,:)
   splt(4,:)   = splt1(3,:)
   delete([/splt0,splt1,nx,ny/])

  labels       = (/"spectrum", "red noise", "5% bound", "90% bound", "95% bound"/)
  leftstring   = exp_str(k)
  right_string = "PC1 of "+vars(i)
  ystring      = "Variance"
  if(xtype.eq."frequency") then
    xstring      = "Frequency (cycles/month)" ;"Period (years)" ;"Frequency (cycles/month)"
  else 
    xstring      = "Period (years)"           ;"Frequency (cycles/month)"
  end if 

  ;;;;;start to plot;;;;;;
  res               = True           ; plot mods desired
  res@vpXF          = 0.15
  res@vpYF          = 0.8
  res@vpWidthF      = 0.8                      ; set width of plot
  res@vpHeightF     = 0.6                      ; set height of plot

  res@gsnDraw       = False
  res@gsnFrame      = False
  res@tmYROn        = False            ; Turn off right tickmarks.
  res@tmXTOn        = False            ; Turn off top tickmarks

  FontHeight = 0.028
  res@gsnStringFontHeightF  = FontHeight ;*0.93
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight
  res@tiYAxisFontHeightF    = FontHeight
  res@tiMainFontThicknessF  = 8.0
  res@tiYAxisFontThicknessF = 5.0
  res@tiXAxisFontThicknessF = 5.0

  res@tmBorderThicknessF       = bthick
  res@tmXBMajorThicknessF      = bthick
  res@tmXBMinorThicknessF      = bthick/2.0
  res@tmYLMajorThicknessF      = bthick
  res@tmYLMinorThicknessF      = bthick/2.0
  res@tmYRMajorThicknessF      = bthick
  res@tmYRMinorThicknessF      = bthick/2.0

  res@tmYMajorGrid                = False        ; Add horizontal grid lines
 ;res@tmYMinorGrid                = False
 ;res@tmYMajorGridLineColor       = "LightGrey"
 ;res@tmYMinorGridLineColor       = "White"
  res@tmYMajorGridLineDashPattern = 0
 ;res@tmYMinorGridLineDashPattern = 5
  res@tmYMajorGridThicknessF      = 0.3
  res@tmGridDrawOrder             = "PreDraw"    ; new resource added in NCL V6.5.0

  res@gsnYRefLine             = 0.0             ; create a reference line   
 ;res@gsnYRefLineDashPattern  = 2             
 ;res@gsnYRefLineColor        = "Black"
 ;res@gsnYRefLineThicknessF   = 1.0

  res@tiMainString            = ""
  res@gsnLeftString           = "" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor  = ""
  res@gsnRightString          = ""
  res@gsnRightStringFontColor = ""

  if(xtype.eq."frequency") then
    res@trXMinF               = 0.0 ;min(latrgd)
    res@trXMaxF               = 0.5 ;max(latrgd)
    res@tmXBMode              = "Explicit"
    res@tmXBValues            = fspan(res@trXMinF,res@trXMaxF,6)
    res@tmXBLabels            = "" + sprintf("%.1f",res@tmXBValues)
    res@tmXBMinorValues       = fspan(res@trXMinF,res@trXMaxF,16)
  else
    res@trXReverse            = True         
    res@trXMinF               = log(0.15)   
    res@trXMaxF               = log(20.0) 
    res@tmXBMode              = "Explicit"
    res@tmXBValues            = log((/30,20,10,5,2,1,0.5,0.2,0.1/))
    res@tmXBLabels            = "" + (/30,20,10,5,2,1,0.5,0.2,0.1/)
   ;res@tmXBMinorValues       = fspan(res@trXMinF,res@trXMaxF,16)
  end if 

  res@trYMinF               = 0.0
  res@trYMaxF               = 12.0
  res@tiMainString          = "" ;mainstr
  res@tiXAxisString         = xstring
  res@tiYAxisString         = ystring
  res@xyLineDashSegLenF     = 0.4
 
  res0                         = res
  res0@gsnLeftString           = leftstring ;"Solid: 1800s       Dash: 300s"
  res0@gsnLeftStringFontColor  = "Black"
  res0@gsnRightString          = ""
  res0@gsnRightStringFontColor = "Black"

  res0@xyMarkLineModes      = "Lines"
  res0@xyLineThicknessF     = lnthick            ; make thicker
  res0@xyDashPatterns       = dashind
  res0@xyMarkerSizes        = 0.012
  res0@xyLineColors         = colors
  res0@xyMarkers            = markers
  res0@xyMarkerColors       = colors
  if(xtype.eq."frequency") then
    plot(k) = gsn_csm_xy(wks,df@sca, splt,res0)
  else
    plot(k) = gsn_csm_xy(wks,log(df@sca), splt,res0)
  end if 

  gres = True

  if(xtype.eq."frequency") then

    gres@YPosPercent                = 60.  
    gres@XPosPercent                = 60.   
  else 
    gres@YPosPercent                = 60.
    gres@XPosPercent                = 20.
  end if 
  gres@Position                   = "Bottom" 
  gres@ItemSpacePercent           = 8.   
  gres@LineLabelWhiteSpacePercent = 4 

  lineres = True
  lineres@lgLineColors            = colors
  lineres@lgLineOpacity           = 1
  lineres@lgLineThicknesses       = lnthick
  lineres@lgDashIndexes           = dashind
  lineres@LineLengthPercent       = 8    

  textres = True
  textres@lgLabels           = labels
  textres@lgLabelFontHeights = FontHeight
  textres@lgLabelColors      = colors
  textres@lgLabelOpacity     = 1
  plot(k) = simple_legend(wks,plot(k),gres,lineres,textres)

 end do 

  panres = True                                         ; panel resource list
  panres@gsnMaximize                = False
  panres@gsnPaperMargin             = 0.3
  panres@gsnPanelBottom             = 0.38   
  panres@gsnPanelYWhiteSpacePercent = 4.0
  gsn_panel(wks,(/plot/),(/1,nexps/),panres)
 ;frame(wks)

 end do

end do 

end



