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
  nvars      = dimsizes(vars)

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

  plotFileLabel = "fig_eof_variance"
  load  "./fig_setup.ncl"

  do i = 0,nvars-1,1

    ;;prepare the data for plot;;
    plotvars = new((/3,nseas,nexps/),float)
    do j = 0, nseas-1,1
    do k = 0,nexps-1,1
     fl      = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars(i)+"_"+seasons(j)+"_"+time_tag+".nc"
     f       = addfile(fl,"r")
     eof_pcv = f->eof_pcv(:,ieof-1)
     mean    = avg(eof_pcv)
     mstd    = stddev(eof_pcv)
     plotvars(0,j,k) = mean 
     plotvars(1,j,k) = mean - mstd
     plotvars(2,j,k) = mean + mstd
     delete(eof_pcv)
    end do 
   end do 
 
   if(i.eq.0) then 
     leftstring   = "EOF variance explained for " +vars(i)
     right_string = ""
   else
     leftstring   = "EOF variance explained for " +vars(i)
     right_string = ""
   end if 
   
   xlabnam  = seasons
 
  ;;prepare the data for plot;;
   plotmn    = new((/nexps,nseas/),float)
   plotmx    = new((/2,nexps,nseas/),float)
   plotxlab  = new((/nexps,nseas/),float)
   do j = 0, nseas-1,1
    do k = 0,nexps-1,1
     fl      = data_dir+"/"+exp_name(k)+"_EOF_PCS_"+vars(i)+"_"+seasons(j)+"_"+time_tag+".nc"
     f       = addfile(fl,"r")
     eof_pcv = f->eof_pcv(:,ieof-1)
     mean    = avg(eof_pcv)
     mstd    = stddev(eof_pcv)
     plotmn(k,j)   = mean
     plotmx(0,k,j) = mean - mstd
     plotmx(1,k,j) = mean + mstd
     plotxlab(k,j) = j+1+0.5*(k-nexps/2.0)
     delete([/fl,f,eof_pcv,mean,mstd/])
    end do
   end do

  ;;;;;start to plot;;;;;;
  res               = True           ; plot mods desired
  res@vpXF          = 0.15
  res@vpYF          = 0.8
  res@vpWidthF      = 0.9                      ; set width of plot
  res@vpHeightF     = 0.3                      ; set height of plot
  res@gsnDraw       = False
  res@gsnFrame      = False
  res@tmYROn        = False            ; Turn off right tickmarks.
  res@tmXTOn        = False            ; Turn off top tickmarks

  FontHeight = 0.020
  res@gsnStringFontHeightF  = FontHeight*1.15
  res@tmXBLabelFontHeightF  = FontHeight
  res@tmYLLabelFontHeightF  = FontHeight
  res@tmYRLabelFontHeightF  = FontHeight
  res@tiMainFontHeightF     = FontHeight
  res@tiXAxisFontHeightF    = FontHeight*0.3
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

  res@gsnYRefLine             = 0.0             ; create a reference line   
  ;res@gsnYRefLineDashPattern = 2             
  ;res@gsnYRefLineColor       = "Black"
  ;res@gsnYRefLineThicknessF  = 1.0

  res@tiMainString            = ""
  res@gsnLeftString           = "" ;"Solid: 1800s       Dash: 300s"
  res@gsnLeftStringFontColor  = ""
  res@gsnRightString          = ""
  res@gsnRightStringFontColor = ""

  
  res@trXMinF               = 0.0 ;min(latrgd)
  res@trXMaxF               = 2*nseas+1.25 ;max(latrgd)
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = fspan(1.075,2*nseas+.075,nseas) 
  res@tmXBMinorValues       = fspan(1.075,2*nseas+.075,nseas) 
  res@tmXBLabels            = xlabnam
 ;res@tmXBLabelDirection    = "Across"
 ;res@tmXBLabelAngleF       = 90.0
 ;res@tmXBLabelJust         = "CenterRight"
 ;res@tmYLMinorPerMajor     = 2

   res@tmYMajorGrid                = False        ; Add horizontal grid lines
  ;res@tmYMinorGrid                = False
  ;res@tmYMajorGridLineColor       = "LightGrey"
  ;res@tmYMinorGridLineColor       = "White"
   res@tmYMajorGridLineDashPattern = 0
  ;res@tmYMinorGridLineDashPattern = 5
   res@tmYMajorGridThicknessF      = 0.3
   res@tmGridDrawOrder             = "PreDraw"    ; new resource added in NCL V6.5.0

  res@gsnXYBarChart                   = True ; create bar chart
  res@gsnXYBarChartOutlineThicknessF  = 0.5;1.008
  res@gsnXYBarChartFillLineThicknessF = 1.0

  ;res@xyLineThicknessF     = 10.0            ; make thicker
  ;res@xyDashPatterns       = dashind(0:nseas)
  ;res@xyMarkerSizes        = 0.012
  ;res@xyLineColors         = colors(0:nseas)
  ;res@xyMarkers            = markers(0:nseas)
  ;res@xyMarkerColors       = colors(0:nseas)

  ;res@tmXTAutoPrecision    = False
  ;res@tmYLPrecision        =  2
   res@tmXBOn               = True            ; Turn off top tickmarks
   res@tmXBLabelsOn         = True            ; Turn off top tickmarks

  do iexp = 0,nexps-1,1

   res@xyMarkLineModes           := "Lines"
   res@xyLineColors              = (/-1,-1/);(/colors(0),colors(0)/)
   res@gsnXYBarChartBarWidth     = barWidth(iexp)          ; change bar widths
   res@gsnXYBarChartColors       = colors(iexp)
   res@gsnXYBarChartFillOpacityF = bopc(iexp)
   res@gsnXYBarChartPatterns     = 0
   res@gsnXYBarChartFillScaleF   = 0.5

   xx  = fspan(1.25,2*nseas+0.25,nseas)+0.61*(iexp-nexps/2.0)

   ;;;;plot global mean biases
   res0 = res
   res0@trYMinF              =  0.0 ;floor(min(pltobs)) ;min(ynflx)-0.1 ;You could also just use
   res0@trYMaxF              =  70.0   ;ceil(max(pltobs)) ;max(ynflx)+0.1 ;min/max of u0,u1,u2.
   res0@tmYLMode             = "Explicit"
   res0@tmYLValues          := fspan(res0@trYMinF,res0@trYMaxF,8)
   res0@tmYLMinorValues     := fspan(res0@trYMinF,res0@trYMaxF,15)
   res0@tmYLLabels          := sprintf("%.0f",res0@tmYLValues )
   if(iexp.eq.0)then
     res0@tiMainString                = ""
     res0@tiYAxisString               = "Percent (%)"
     res0@tiXAxisString               = ""
     res0@gsnLeftString               = leftstring
     res0@gsnLeftStringOrthogonalPosF = 0.04
     res0@gsnLeftStringFontColor      = "Black"
     res0@gsnRightString              = right_string
     res0@gsnRightStringFontColor     = "Black"
   else
     res0@tiMainString                = ""
     res0@tiYAxisString               = ""
     res0@tiXAxisString               = ""
     res0@gsnLeftString               = ""
     res0@gsnRightString              = ""
   end if

   dum1(i,iexp) = gsn_csm_xy(wks,xx,plotmn(iexp,:),res0)

   ;; for bars with fill pattern 
   if(i.eq.1) then 
    res0@gsnXYBarChartPatterns         = pattns(i)  
    res0@gsnXYBarChartColors           = 1
    res0@gsnXYBarChartFillScaleF       = BarChartFillScale *2.0
    res0@gsnXYBarChartFillOpacityF     = 0.5 
    res0@gsnXYBarChartFillLineThicknessF = 0.01
    dum2(i,iexp) = gsn_csm_xy(wks,xx,plotmn(iexp,:),res0)
    overlay(dum1(i,iexp),dum2(i,iexp))
   end if 

   if(iexp.eq.0)then
    plot(i) = dum1(i,iexp)
   else
    overlay (plot(i),dum1(i,iexp))
   end if

 ; add uncertainty bars
   xx1    = xx
   nx1    = dimsizes(xx1)
   dumx1  = new((/nx1/),"graphic")
   plres  = True
   plres@gsLineThicknessF = 1.5
   plres@gsLineColor = "Black"
   plres@tfPolyDrawOrder = "Draw"
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

  ;----------------------
  ;Attach a legend
  ;----------------------
  lgres                    = True
  lgres@vpWidthF           = 0.08          ; width of legend (NDC)
  lgres@vpHeightF          = 0.05           ; height of legend (NDC)
  lgres@lbBoxMajorExtentF  = 0.12          ; puts space between color boxes
  lgres@lbBoxMinorExtentF  = 0.18
  lgres@lbMonoFillPattern  = False          ; Solid fill pattern
  lgres@lgLabelFontHeightF = FontHeight*0.2     ; legend label font thickness
  lgres@lgMonoDashIndex    = False
  lgres@lgMonoLineThickness= False
  lgres@lgPerimOn          = True
  lgres@lgPerimThicknessF  = 1.0
  lgres@lbLabelJust        = "CenterLeft"  ; left justify labels

  do iexp = 0,nexps-1,1
    ;lgres@lbBoxLineDashPattern = 2 ;dashind(iexp) 
     lgres@lbBoxLineColor       = "black"
     lgres@lbBoxLineThicknessF  = 0.5
     lgres@lbFillColors         = colors(iexp)
     lgres@lbFillOpacityF       = bopc(1)*3.0
     lgres@lbFillPatterns       = 0;pattns(iexp)
    ;lgres@gsnXYBarChartFillScaleF = 0.4
    ;lgres@lbFillPatterns         = pattns(iexp)
     print("ndguv pattern : " + iexp + "   " + pattns(iexp)) 
     lgres@vpWidthF           = 0.12          ; width of legend (NDC)
     ypos = 0.90 + (nexps-iexp-1)*0.024
     gsn_labelbar_ndc(wks,1,labels(iexp),0.825,ypos,lgres)

     ypos = 0.60 + (nexps-iexp-1)*0.024
     gsn_labelbar_ndc(wks,1,labels(iexp),0.825,ypos,lgres)

     if(i.eq.1)  then 
       lgres@lbFillPatterns         = pattns(i)
       lgres@lbFillColors           = 1 
       lgres@lbFillScaleF           = BarChartFillScale 
       lgres@lbFillLineThicknessF   = 0.01
       ypos = 0.60 + (nexps-iexp-1)*0.024
       gsn_labelbar_ndc(wks,1,labels(iexp),0.825,ypos,lgres)
     end if 

   end do 

 end do

 panres = True                                         ; panel resource list
 panres@gsnMaximize                = True
 panres@gsnPaperMargin             = 0.3
 panres@gsnPanelBottom             = 0.38   
 panres@gsnPanelYWhiteSpacePercent = 4.0
 gsn_panel(wks,(/plot(0),plot(1)/),(/2,1/),panres)
;frame(wks)

end



