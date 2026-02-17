;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Figure 5 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

  data_dir   = "../data"
  time_tag   = (/"1979-2014"/)
  mode       = ("SAM")
  ieof       = 1

  vars       = (/"PSL","Z500"/)
  varunt     = (/"hPa", "m"/)
  nvars      = dimsizes(vars)
  miny       = (/-6, -60/)
  maxy       = (/ 3,  30/)

  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)
  labels     = (/"ERA5          ", \
                 "v2 AMIP    ", \
                 "v2 Coupled"/)

  seasons    = (/"MON","ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)

  plot       = new(nvars,graphic)
  dum1       = new((/nvars,nexps/),graphic)
  dum2       = new((/nvars,nexps/),graphic)

  do j = 0, nseas-1,1

    plotFileLabel = "fig_zonal_mean_anomaly_" + seasons(j)
    load  "./fig_setup.ncl"

   do i = 0,nvars-1,1

   print("working on season "+ seasons(j))

   ;;prepare the data for plot;;
   do k = 0,nexps-1,1
   ;fl   = data_dir+"/"+exp_name(k)+"_anom_"+vars(i)+"_"+seasons(j)+"_"+time_tag+".nc"
    fl   = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars(i)+"_"+seasons(j)+"_"+time_tag+".nc"
    f    = addfile(fl,"r")

    vd1  = dim_avg_n_Wrap(f->eof_reg(:,0,:,:),2)

    if(vars(i).eq."PSL") then 
      vd1 = vd1/100.0 ; to hPa
    end if 

    mean = dim_avg_n_Wrap(vd1,0)
    mstd = dim_stddev_n_Wrap(vd1,0)
    printVarSummary(vd1)
    if(k.eq.0) then 
      lat       = f->lat
      nlat      = dimsizes(lat)
      plotmn    = new((/nexps,nlat/),float)
      plotmx    = new((/2,nexps,nlat/),float)
      plotxlab  = new((/nexps,nlat/),float)
    end if 

    plotmn(k,:)   = mean 
    plotmx(0,k,:) = mean - mstd
    plotmx(1,k,:) = mean + mstd
    delete([/fl,f,vd1,mean,mstd/])
   end do 

   leftstring   = "Zonal mean " +vars(i)+ " anomaly"
   right_string = ""
   xlabnam      = seasons
 
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

  bthick = 1.0
  res@tmBorderThicknessF    = bthick
  res@tmXBMajorThicknessF   = bthick
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMajorThicknessF   = bthick
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0

   res@tmYMajorGrid                = False        ; Add horizontal grid lines
  ;res@tmYMinorGrid                = False
  ;res@tmYMajorGridLineColor       = "LightGrey"
  ;res@tmYMinorGridLineColor       = "White"
   res@tmYMajorGridLineDashPattern = 0
  ;res@tmYMinorGridLineDashPattern = 5
   res@tmYMajorGridThicknessF      = 0.3
   res@tmGridDrawOrder             = "PreDraw"    ; new resource added in NCL V6.5.0

  res@gsnYRefLine             = 0.0             ; create a reference line   
 ;res@gsnYRefLineDashPattern = 2             
 ;res@gsnYRefLineColor       = "Black"
 ;res@gsnYRefLineThicknessF  = 1.0

  res@tiMainString            = ""
  res@gsnLeftString           = "" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor  = ""
  res@gsnRightString          = ""
  res@gsnRightStringFontColor = ""


  res@trXMinF               = -90 ;min(latrgd)
  res@trXMaxF               = -20 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = ispan(-90,-20,10)
  res@tmXBLabels            = "" + res@tmXBValues
  res@tmXBMinorValues       = ispan(-90,90,10)

  res@trYMinF               = miny(i)
  res@trYMaxF               = maxy(i)
  res@tiMainString          = "" ;mainstr
  res@tiXAxisString         = "Latitude (~S~o~N~S)"
  res@tiYAxisString         = vars(i)+ "  (" + varunt(i) + ")"
  res@xyLineDashSegLenF     = 0.4
 
  do k = 0,nexps-1,1
    res0                         = res
    res0@gsnLeftString           = leftstring ;"Solid: 1800s       Dash: 300s"
    res0@gsnLeftStringFontColor  = "Black"
    res0@gsnRightString          = ""
    res0@gsnRightStringFontColor = "Black"

    res0@xyMarkLineModes      = "Lines"
    res0@xyLineThicknessF     = lnthick(k)/5.0            ; make thicker
    res0@xyDashPatterns       = dashind(k)
    res0@xyMarkerSizes        = 0.012
    res0@xyLineColors         = colors(k)
    res0@xyMarkers            = markers(k)
    res0@xyMarkerColors       = colors(k)
    dum1(i,k) = gsn_csm_xy(wks,lat,plotmn(k,:),res0)
   
    gsres0                    = res                        ; poly res
    gsres0@tmYROn = False            ; Turn off right tickmarks.
    gsres0@tmYLOn = False            ; Turn off top tickmarks
    gsres0@tmXTOn = False            ; Turn off top tickmarks
    gsres0@tmXBOn = False            ; Turn off top tickmarks
    gsres0@gsnXYFillColors    = colors(k)                ; color chosen
    gsres0@gsnXYFillOpacities = 0.1
    gsres0@xyLineColors       = (/-1,-1/)
    dum2(i,k)   = gsn_csm_xy (wks,lat,plotmx(:,k,:),gsres0)
    
    if(k.eq.0) then 
      plot(i) = dum1(i,k)
      overlay(plot(i),dum2(i,k))
    else 
      overlay(plot(i),dum1(i,k))
      overlay(plot(i),dum2(i,k))
    end if 

  end do 

  gres = True
  gres@YPosPercent                = 15.  
  gres@XPosPercent                = 60.   
  gres@Position                   = "Bottom" 
  gres@ItemSpacePercent           = 8.   
  gres@LineLabelWhiteSpacePercent = 4 

  lineres = True
  lineres@lgLineColors            = colors(0:nexps-1)
  lineres@lgLineOpacity           = 1
  lineres@lgLineThicknesses       = lnthick(0:nexps-1)/5.0
  lineres@lgDashIndexes           = dashind(0:nexps-1)
  lineres@LineLengthPercent       = 8    

  textres = True
  textres@lgLabels           = labels
  textres@lgLabelFontHeights = FontHeight
  textres@lgLabelColors      = colors(0:nexps-1)
  textres@lgLabelOpacity     = 1
  plot(i) = simple_legend(wks,plot(i),gres,lineres,textres)

 end do 

  panres = True                                         ; panel resource list
  panres@gsnMaximize                = False
  panres@gsnPaperMargin             = 0.3
  panres@gsnPanelBottom             = 0.38   
  panres@gsnPanelYWhiteSpacePercent = 4.0
  gsn_panel(wks,(/plot(0),plot(1)/),(/1,2/),panres)
 ;frame(wks)

 end do

end



