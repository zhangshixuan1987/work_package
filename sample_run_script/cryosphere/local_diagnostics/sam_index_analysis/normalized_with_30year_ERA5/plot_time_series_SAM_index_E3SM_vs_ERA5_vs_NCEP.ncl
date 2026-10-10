;*****************************************************
; xy_18.ncl
;
; Concepts illustrated:
;   - Filling the area between two curves in an XY plot
;   - Labeling the bottom X axis with years
;   - Drawing a main title on three separate lines
;   - Calculating a weighted average
;   - Changing the size/shape of an XY plot using viewport resources
;   - Manually creating a legend
;   - Overlaying XY plots on each other
;   - Maximizing plots after they've been created
;
;*****************************************************
;
; These files are loaded by default in NCL V6.2.0 and newer
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
; load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;*****************************************************
begin

 Datdir   = "./"

 ExpNames = (/"NCEP","ERA5","E3SM_NE30", "E3SM_RRM", "E3SM_NE120"/)
 nexps    = dimsizes(ExpNames)
 Labels   = (/"NCEP","ERA5","E3SM (NE30)", "E3SM (RRM)", "E3SM (NE120)"/)

 maxy               = 8.0
 miny               = -8.0
 yint               = 1.0
 amJust             = "BottomLeft"    ; Use bottom right corner of box
 amParallelPosF     = -0.48             ; Move legend to right
 amOrthogonalPosF   =  0.48             ; Move legend down.



 ;;;;;only have year 2010 to 2012 data;;;;
 nyear  = 3
 nmonth = 12
 ncol   = nmonth + 1
 ntime = nyear*nmonth

 yyyy   = new(nyear,float)
 mm     = new(nmonth,float)
 time   = new(ntime,float)
 yyyy   = ispan(2010,2012,1)
 mm     = ispan(1,12,1)
 time  = ispan(1,ntime,1)
 dateF = time

 do iy = 0,nyear-1,1
  do im = 0, nmonth-1,1
    n = iy*nmonth + im 
    dateF(n) = yyyy(iy) + (mm(im)-1)*1.0/nmonth  
  end do
 end do 

 plotdata = new ((/nexps,ntime/),float)

 ;*****************************************************
 ; open files and read in monthly data
 ;*****************************************************

 do iexp = 0,nexps-1,1

  File             = "SAM_Norm_AOI_Monthly_2010-2012_"+ExpNames(iexp)+".txt"
  data             = asciiread(File,(/nyear,ncol/),"float")     ; 2 x 70 x 70
  plotdata(iexp,:) = ndtooned(data(:,1:)) 
 
 end do 

 ;plotdata = runave (plotdata, 3, 0)

;;;;correlation;;;;
 r      = escorc_n(plotdata(0,:),plotdata(1:,:),0,1)   ; r is a scalar

 plotFileLabel ="monthly_variation_2010-2012"
;***************************************************** 
; create plot
;***************************************************** 
  OUTDir  = "./"
  wtype = "png"
 ;wtype@wkPaperSize     = "A4"
 ;wtype@wkOrientation   = "landscape"    ;;pictures displayed along the length
  wks = gsn_open_wks(wtype,OUTDir+"SAM_indices_"+plotFileLabel)

  ;;;;set up line color, thickness, style etc ;;;;;;;;;;;;;;;;;;;;;;;
  gsn_define_colormap(wks,"cb_9step")
  colors  = (/2,6,28,3,1/)
  markers = (/16,11,8,2,7/)       ;;line markers
  mkszf   = (/0.8,0.8,0.8,0.8/)/2.0   ;;line marker size
  lnthick = (/3,3,3,3,3,3/)*2.0   ;;line thichness
  dashind = (/0,0,0,0,0,0/)

  res                    = True              ; plot mods desired
  res@gsnDraw            = False             ; don't draw yet
  res@gsnFrame           = False            ; don't advance frame yet
  res@vpHeightF 	 = 0.4               ; change aspect ratio of plot
  res@vpWidthF 	         = 0.8

  res@tmXTOn                 = False
  res@tmYROn                 = False

  FontHeight = 0.022
  res@gsnStringFontHeightF      = FontHeight
  res@gsnRightStringFontHeightF = FontHeight
  res@gsnLeftStringFontHeightF  = FontHeight*1.2
  res@tmXBLabelFontHeightF      = FontHeight
  res@tmYLLabelFontHeightF      = FontHeight
  res@tmYRLabelFontHeightF      = FontHeight
  res@tiMainFontHeightF         = FontHeight
  res@tiXAxisFontHeightF        = FontHeight*1.2
  res@tiYAxisFontHeightF        = FontHeight*1.2

  ;res@gsnLeftStringParallelPosF    = -0.15
  res@gsnLeftStringOrthogonalPosF  =  0.04
  res@gsnRightStringOrthogonalPosF =  0.04

  res@xyMarkLineMode         = "Lines" ;"MarkLines"
  res@xyLineColors           = colors
  res@xyMarkers              = markers
  res@xyMarkerColors         = colors
  res@xyMarkerSizes          = 0.015*mkszf
  res@xyDashPatterns         = dashind
  res@xyLineThicknesses      = lnthick


  res@tmXBMode          = "Explicit"
  res@tmXBValues        = time(::3)
  res@tmXBLabels        = "" + res@tmXBValues
  res@tmXBMinorValues   = time

  res@tiXAxisString      = "Months science January 2010"  
  res@tiYAxisString      = "SAM Index"    ; y-axis label
  res@tiMainString       = "Variations of SAM Index (monthly)" ; title


  res@trYMaxF          =  maxy
  res@trYMinF          =  miny
  res@tmYLTickSpacingF =  yint
 ;res@trXMinF          = 1890              ; set x-axis minimum

  res@gsnYRefLine           = 0.0             ; create a reference line
  res@gsnYRefLineDashPattern= 0
  res@gsnYRefLineColor      = "Grey";"Black"
  res@gsnYRefLineThicknessF = 1.0

  plot = gsn_csm_xy (wks,time,plotdata,res)        ; create plot

  lgres                     = True
  lgres@lgItemType          = "Lines"
  lgres@lgLineLabelFontThicknessF= 3
  lgres@lgBoxMinorExtentF   = 0.9
  lgres@lgLineDashSegLenF   = 0.15
  lgres@vpWidthF            = 0.22          ; width of legend (NDC)
  lgres@vpHeightF           = 0.16           ; height of legend (NDC)
  lgres@lgLabelFontHeightF  = .24            ; legend label font thickness
  lgres@lgMonoDashIndex     = False
  lgres@lgMonoLineThickness = False
  lgres@lgPerimOn           = False
  lgres@lgPerimThicknessF   = 2.0

  lgres1 = lgres
  lgres1@lgLineColors       = colors ;(/"Black","Black","Black"/)
  lgres1@lgDashIndexes      = dashind
  lgres1@lgLineThicknesses  = lnthick
  Labels                    = "  "+Labels
  legend1 = gsn_create_legend (wks,dimsizes(Labels), Labels,lgres1)

  ;;;;;;;;;;;;
  amres                      = True
  amres@amJust               = amJust    ; Use bottom right corner of box
  amres@amParallelPosF       = amParallelPosF             ; Move legend to right
  amres@amOrthogonalPosF     = amOrthogonalPosF           ; Move legend down.
  annoid0 = gsn_add_annotation(plot,legend1,amres)  ; add legend to the top figure
  delete([/lgres, lgres1, amres/])

  draw(plot)
  frame(wks)
 
end
