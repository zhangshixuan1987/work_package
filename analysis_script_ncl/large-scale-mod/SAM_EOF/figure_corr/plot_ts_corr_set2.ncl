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
  mlag       = 6 ; lag correlation 
  psig       = 0.10445

  vars       = (/"PSL","Z500"/)
  nvars      = dimsizes(vars)

  exp_name   = (/"v2.LR.amip", \
                 "v2.LR.historical",\
                 "v2.LR.historical",\
                 "v2.LR.historical",\
                 "v2.LR.historical",\
                 "v2.LR.historical"/)
  nexps      = dimsizes(exp_name)

  labels     = (/"v2 AMIP            ", \
                 "v2 Coupled         ", \
                 "v2 Coupled (No DJF)", \
                 "v2 Coupled (No SON)", \
                 "v2 Coupled (No JJA)", \
                 "v2 Coupled (No MAM)"/)

  seasons    = (/"MON"/)
  nseas      = dimsizes(seasons)

  plot       = new(nvars,graphic)
  dum1       = new((/nvars,nexps/),graphic)
  dum2       = new((/nvars,nexps/),graphic)

  plotFileLabel = "fig_ts_corr_set2"
  load  "./fig_setup_set2.ncl"

  do i = 0,nvars-1,1

    ;;prepare the data for plot;;
    plotmn   = new((/nexps,mlag+1/),float)
    plotmx   = new((/2,nexps,mlag+1/),float)
    
    do k = 0,nexps-1,1

       fl    = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars(i)+"_"+seasons+"_"+time_tag+".nc"
       f     = addfile(fl,"r")
       time  = yyyymm_to_yyyyfrac(cd_calendar(f->time,-1),0.0)
       indx  = ind(time.ge.1981.and.time.lt.2011)
       eoft  = f->eof_ts(:,ieof-1,indx)
       tutc  = cd_calendar(f->time(indx),0)       
      ;print(time(indx))
      ;exit

       if(k.eq.2) then 
         indm = ind(tutc(:,1).ne.12.and.tutc(:,1).ne.1.and.tutc(:,1).ne.2) 
         rc    = esacr_n(eoft(:,indm),mlag,1)
         delete(indm)
       else if(k.eq.3) then
         indm = ind(tutc(:,1).ne.9.and.tutc(:,1).ne.10.and.tutc(:,1).ne.11)
         rc    = esacr_n(eoft(:,indm),mlag,1)
         delete(indm)
       else if(k.eq.4) then
         indm = ind(tutc(:,1).ne.6.and.tutc(:,1).ne.7.and.tutc(:,1).ne.8)
         rc    = esacr_n(eoft(:,indm),mlag,1)
         delete(indm)
       else if(k.eq.5) then
         indm = ind(tutc(:,1).ne.3.and.tutc(:,1).ne.4.and.tutc(:,1).ne.5)
         rc    = esacr_n(eoft(:,indm),mlag,1)
         delete(indm)
       else 
         rc    = esacr_n(eoft,mlag,1)
       end if 
       end if 
       end if
       end if 

       plotmn(k,:) = dim_avg_n(rc,0)
       plotmx(0,k,:) = dim_avg_n(rc,0) - dim_stddev_n(rc,0)
       plotmx(1,k,:) = dim_avg_n(rc,0) + dim_stddev_n(rc,0)
       delete([/time,indx,eoft,rc/])

   end do 

   if(i.eq.0) then 
     leftstring   = "Auto corr. of SAM index" 
     right_string = vars(i)+ " EOF"
   else
     leftstring   = "Auto corr. of SAM index"
     right_string = vars(i)+ " EOF"
   end if 
   
   xlabnam  = seasons
   xx       = ispan(0,mlag,1)*1.0

  ;;;;;start to plot;;;;;;
  res               = True           ; plot mods desired
  res@vpWidthF      = 0.8                     ; set width of plot
  res@vpHeightF     = 0.8                      ; set height of plot
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@tmYROn        = False            ; Turn off right tickmarks.
  res@tmXTOn        = False            ; Turn off top tickmarks

  FontHeight = 0.032
  res@gsnStringFontHeightF  = FontHeight
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight
  res@tiYAxisFontHeightF    = FontHeight

  ;res@tiMainFontThicknessF  = 8.0
  ;res@tiYAxisFontThicknessF = 5.0
  ;res@tiXAxisFontThicknessF = 5.0

  bthick = 1.0
  res@tmBorderThicknessF    = bthick
  res@tmXBMajorThicknessF   = bthick
  res@tmXBMinorThicknessF   = res@tmXBMajorThicknessF/2.0
  res@tmYLMajorThicknessF   = bthick
  res@tmYLMinorThicknessF   = res@tmYLMajorThicknessF/2.0

  res@gsnYRefLine            = psig            ; create a reference line   
  res@gsnYRefLineDashPattern = 0             
  res@gsnYRefLineColor       = "Grey25"
  res@gsnYRefLineThicknessF  = 1.0

  res@tiMainString            = ""
  res@gsnLeftString           = "" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor  = ""
  res@gsnRightString          = ""
  res@gsnRightStringFontColor = ""

  res@trXMinF               = 0    ;min(latrgd)
  res@trXMaxF               = mlag+1 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = ispan(res@trXMinF,res@trXMaxF,1)*1.0
  res@tmXBMinorValues       = ispan(res@trXMinF,res@trXMaxF,1)*1.0
  res@tmXBLabels            = ispan(res@trXMinF,res@trXMaxF,1)*1.0

  res@tmYMajorGrid                = False        ; Add horizontal grid lines
 ;res@tmYMinorGrid                = False
 ;res@tmYMajorGridLineColor       = "LightGrey"
 ;res@tmYMinorGridLineColor       = "White"
  res@tmYMajorGridLineDashPattern = 0
 ;res@tmYMinorGridLineDashPattern = 5
  res@tmYMajorGridThicknessF      = 0.3
  res@tmGridDrawOrder             = "PreDraw"    ; new resource added in NCL V6.5.0

 ;res@tmXTAutoPrecision    = False
 ;res@tmYLPrecision        =  2
  res@tmXBOn               = True            ; Turn off top tickmarks
  res@tmXBLabelsOn         = True            ; Turn off top tickmarks

   ;;;;plot global mean biases
   res0 = res
   res0@trYMinF            = -0.1 
   res0@trYMaxF            =  1.1 ;ceil(max(pltobs)) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
   res0@tmYLMode           = "Explicit"
   res0@tmYLValues        := fspan(res0@trYMinF,res0@trYMaxF,7)
   res0@tmYLMinorValues   := fspan(res0@trYMinF,res0@trYMaxF,13)
   res0@tmYLLabels        := sprintf("%.1f",res0@tmYLValues )

   res0@tiMainString                = ""
   res0@tiYAxisString               = "Correlation coefficient"
   res0@tiXAxisString               = "Laged months"
   res0@gsnLeftString               = leftstring
   res0@gsnLeftStringOrthogonalPosF = 0.0
   res0@gsnLeftStringFontColor      = "Black"
   res0@gsnRightString              = right_string
   res0@gsnRightStringFontColor     = "Black"

   res0@xyMarkLineMode      := "MarkLines"
   res0@xyLineThicknessF     = lnthick(0:nexps-1)
   res0@xyDashPatterns       = dashind(0:nexps-1)
   res0@xyMarkerSizes        = mkszf(0:nexps-1)
   res0@xyMarkerThicknessF   = lnthick(0:nexps-1)
   res0@xyLineColors         = colors(0:nexps-1)
   res0@xyMarkers            = markers(0:nexps-1)
   res0@xyMarkerColors       = colors(0:nexps-1)

   plot(i) = gsn_csm_xy(wks,xx,plotmn(:,:),res0)

 ; add uncertainty bars
   xx1    = xx
   nx1    = dimsizes(xx1)
   dumx1  = new((/nx1/),"graphic")
   plres  = True
   plres@gsLineThicknessF = 1.5
   plres@tfPolyDrawOrder = "Draw"
   do iexp = 1,nexps-1,1
    plres@gsLineColor = colors(iexp)
    do ix = 0,nx1-1
      barxx = (/xx1(ix),xx1(ix)/)
      baryy = (/plotmx(0,iexp,ix),plotmx(1,iexp,ix)/)
      str1  = unique_string("bar")
      dumx1@$str1$ =  gsn_add_polyline( wks, plot(i), barxx, baryy, plres )
      barxx = (/xx1(ix)-0.08,xx1(ix)+0.08/)
      baryy = (/plotmx(0,iexp,ix),plotmx(0,iexp,ix)/)
      str1  = unique_string("bar")
      dumx1@$str1$ =  gsn_add_polyline( wks, plot(i), barxx, baryy, plres )
      barxx = (/xx1(ix)-0.08,xx1(ix)+0.08/)
      baryy = (/plotmx(1,iexp,ix),plotmx(1,iexp,ix)/)
      str1  = unique_string("bar")
      dumx1@$str1$ =  gsn_add_polyline( wks, plot(i), barxx, baryy, plres )
   end do
   end do 

  gres = True
  gres@YPosPercent                = 50.
  gres@XPosPercent                = 30.
  gres@Position                   = "Bottom"
  gres@ItemSpacePercent           = 8.
  gres@LineLabelWhiteSpacePercent = 4

  lineres = True
  lineres@lgLineColors            = colors(0:nexps-1)
  lineres@lgLineOpacity           = 1
  lineres@lgLineThicknesses       = lnthick(0:nexps-1)
  lineres@lgDashIndexes           = dashind(0:nexps-1)
  lineres@LineLengthPercent       = 8

  textres = True
  textres@lgLabels           = labels(0:nexps-1)
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
 panres@gsnPanelXWhiteSpacePercent = 4.0
 gsn_panel(wks,(/plot(0),plot(1)/),(/1,2/),panres)
;frame(wks)

end



