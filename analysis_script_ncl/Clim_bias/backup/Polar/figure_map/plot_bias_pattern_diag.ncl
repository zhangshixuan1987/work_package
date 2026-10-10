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

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/reg_analysis_era5"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"v2 AMIP","v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  seasons    = (/"ANN","DJF","JJA","MAM","SON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions = (/"Global"/)
  nregs   = dimsizes(Regions)
  latS    = -90.
  latN    = -20.
  lonW    =  0.
  lonE    = 360.

  load "./share_var_info.ncl"
  nvar = dimsizes(varList)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do j = 0,nvar-1,1

  var   = varList(j)
  varunt = untList(j)

  do k = 0,nexps-1

  print(data_dir)
  print(exp_name(k))
  print(seasons(i))
   
  fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + var+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
  f1   = addfiles(fl1,"r")
  ListSetType (f1, "join")
  
 ;print(fl)

  vnm  = var+"_c"
  vd1  = f1[:]->$vnm$(:,{latS:latN},{lonW:lonE})
  time = f1[0]->time

  vdm  = dim_avg_n_Wrap(vd1,0)
  rad  = 4.0*atan(1.0)/180.0
  wgty = cos(f1[0]->lat({latS:latN})*rad)
  nens = dimsizes(vd1(:,0,0))
  tmp1 = dim_max_n(dim_max_n(vd1,2),1) 
  tmp2 = dim_min_n(dim_min_n(vd1,2),1) 
  tmp3 = wgt_areaave(vd1,wgty,1.0,1)

  glm1  = sprintf("%3.1f", avg(tmp1)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp1))
  glm2  = sprintf("%3.1f", avg(tmp2)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp2))
  glm3  = sprintf("%3.1f", avg(tmp3)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp3))
  mstr1 = "Min:  " + glm1
  mstr2 = "Max:  " + glm2
  mstr3 = "Mean: " + glm3
  right_str = mstr1 + "~C~" + mstr2 + "~C~" + mstr3

  delete([/tmp1,tmp2,tmp3,glm1,glm2,glm3,mstr1,mstr2,mstr3/])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Start to process plot 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (k.eq.0) then 
    ;;;start to process the plot;;;;
    load "./share_color_bar_for_polar_map.ncl"
    wtype           = "png"
    wtype@wkWidth   = 2500
    wtype@wkHeight  = 2500
    fe              = var+"_bias_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)
    cmap1           = read_colormap_file("BlGrYeOrReVi200")
    cmap2           = read_colormap_file("BlueDarkRed18")
    plot            = new(nexps,graphic)
  end if 

  if (seasons(i).eq."MON") then
    main_str = "Monthly mean"
  else
    main_str = seasons(i)+" mean"
  end if

  exp_str  = exp_strs(k)
  left_str = ""

  ;============================================================
  ; PLOTS
  ;============================================================
   res                       = True         
   res@gsnDraw               = False        ; don't draw yet
   res@gsnFrame              = False        ; don't advance frame yet
   res@gsnPolar              = "SH"
 
   res@mpFillOn              = False        ; turn off map fill
   res@mpMaxLatF             = -20
  ;res@mpCenterLonF          = 180

   FontHeightF = 0.018
   res@tiMainFontHeightF     = FontHeightF*1.5
   res@tmYLLabelFontHeightF  = FontHeightF
   res@tmXBLabelFontHeightF  = FontHeightF
   res@lbLabelFontHeightF    = FontHeightF
   res@gsnStringFontHeightF  = FontHeightF*1.2

   res@cnFillOn              = True         ; turn on color fill
   res@cnLinesOn             = False        ; True is default
   res@cnLineLabelsOn        = False        ; True is default

   res@lbLabelBarOn             = False
   res@lbOrientation            = "Vertical"
   res@lbTitleString            = "" 
   res@lbTitlePosition          = "Right"     ; title location
   res@lbTitleDirection         = "Across"    ; letter angle
   res@lbTitleAngleF            = 90.         ; title angle
   res@lbTitleFontHeightF       = FontHeightF ; font height
   res@lbBoxLineColor           = "black"  
  ;res@lbBoxLinesOn             = False
  ;res@pmLabelBarParallelPosF   = 0.61
   res@pmLabelBarOrthogonalPosF = 0.15

   res@tmXTOn                   = False

   if(wtype.eq."png") then
    lthick = 3.0
   else
    lthick = 1.0
   end if
  
   res@tmBorderThicknessF       = lthick
   res@tmXBMajorThicknessF      = lthick
   res@tmXBMinorThicknessF      = lthick
   res@tmYLMajorThicknessF      = lthick
   res@tmYLMinorThicknessF      = lthick
   res@tmYRMajorThicknessF      = lthick
   res@tmYRMinorThicknessF      = lthick

  ;*******************************************
  ; first plot
  ;*******************************************
   res@gsnPolarLabelFontHeightF    = FontHeightF
   res@mpLabelFontHeightF          = FontHeightF*2.0
   res@mpGridLineThicknessF        = lthick* 0.5
   res@mpLimbLineThicknessF        = lthick* 0.5
   res@mpPerimLineThicknessF       = lthick* 0.5
   res@mpGeophysicalLineThicknessF = lthick

   res@gsnLeftString     = ""
   res@gsnRightString    = "" 
   res@gsnCenterString   = "" 
   res@tiMainString      = ""
   res@tiMainOffsetYF    = 0.005

   res0                      = res
   res0@cnFillPalette        = "temp_19lev" ;"BlueDarkRed18" ;"BlueWhiteOrangeRed" ; "CBR_coldhot"
   res0@cnLevelSelectionMode = "ExplicitLevels"
   res0@tiMainString         = exp_str
   res0@gsnLeftString        = left_str
   res0@gsnRightString       = right_str
   res0@cnLevels            := s2
   res0@cnLinesOn            = True        ; True is default
   res0@cnLineColor          = "Black"                   ; color of second contours
   res0@cnLineThicknessF     = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1

   plot(k) = gsn_csm_contour_map_polar(wks,vdm,res0)

  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012

   pres1@gsnPanelLabelBar           = True                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = main_str + " of " + var + " (" + varunt +")"
   pres1@lbTitleFontHeightF         = FontHeightF*0.65
   pres1@lbLabelFontHeightF         = FontHeightF*0.6 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
  ;pres1@pmLabelBarParallelPosF     = 0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01
   pres1@pmLabelBarWidthF           = 0.95
   pres1@pmLabelBarHeightF          = 0.08

  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 8
   pres1@gsnPanelXWhiteSpacePercent = 4
   gsn_panel(wks,plot,(/1,nexps/),pres1)          ; create first panel
  
 end do 

end do 

end

