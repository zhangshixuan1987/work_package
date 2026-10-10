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

  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim/"
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "v2 AMIP"   , "v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"


  vars       = (/"TREFHT"/)
  varunt     = (/"K"/)
  varfac     = (/1/)

  seasons    = (/"ANN","DJF","JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions = (/"Global"/)
  nregs   = dimsizes(Regions)
  latS    = -90.
  latN    =  90.
  lonW    =  0.
  lonE    = 360.
 
  l_plot_latlon_box  = True

  rad    = 4.0*atan(1.0)/180.0
  re     = 6371220.0
  rr     = re*rad

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  print("working on " + vars)

  do k = 0,nexps-1

  fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
  f1   = addfiles(fl1,"r")
  ListSetType (f1, "join")
  
  time = f1[0]->time
  vnm  = vars+"_o"
  var0 = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})
  var0!0 = "ens"  
  if(seasons(i).eq."ANN") then 
    var1 = month_to_annual(var0(time|:,ens|:,lat|:,lon|:),1) 
    printVarSummary(var1) 
  else
    var1 = month_to_season(var0(time|:,ens|:,lat|:,lon|:),seasons(i))
    printVarSummary(var1)
  end if 

  gradLatLon  = grad_latlon_cfd (var1, var1&lat, var1&lon, True, False)
  grad_lat    = gradLatLon[0] 
  grad_lon    = gradLatLon[1] 

  delete([/var0,gradLatLon/])

  grad_lat = grad_lat * rr
  grad_lon = grad_lon * rr
  printMinMax(grad_lat,False)
  printMinMax(grad_lon,False)

 ;replace the generic attributes with more specific attributes
  grad_lon@long_name = "longitudinal gradient (derivative)"
  grad_lat@long_name = "latitudinal gradient (derivative)"
  grad_lat@units     = "K/km"
  grad_lon@units     = "K/km"

  ;;long-term mean 
  vd1 = dim_avg_n_Wrap(grad_lat,0)

  ;;ensemble mean  
  vdm  = dim_avg_n_Wrap(dim_avg_n_Wrap(grad_lat,0),0)

  ;;global annual mean
  rad  = 4.0*atan(1.0)/180.0
  wgty = cos(f1[0]->lat({latS:latN})*rad)
  nens = dimsizes(vd1(:,0,0))
  tmp1 = dim_max_n(dim_max_n(vd1,2),1) 
  tmp2 = dim_min_n(dim_min_n(vd1,2),1) 
  tmp3 = wgt_areaave(vd1,wgty,1.0,0)

  glm1  = sprintf("%3.1f", avg(tmp1)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp1))
  glm2  = sprintf("%3.1f", avg(tmp2)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp2))
  glm3  = sprintf("%3.1f", avg(tmp3)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp3))
  right_str = "Min / Max / Mean:  " + glm1+" / "+glm2+" / "+glm3

  delete([/time,var1,grad_lat,grad_lon,fl1,f1,vd1,tmp1,tmp2,tmp3,glm1,glm2,glm3/])

  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  ; Start to process plot 
  ;;;;;;;;;;;;;;;;;;;;;;;;;;
  if (k.eq.0) then 
    ;;;start to process the plot;;;;
    load "./share_fig_info.ncl"
    wtype           = "png"
    if(wtype.eq."png") then
      lthick          = 3.0
      BorderThick     = 3.0
      wtype@wkWidth   = 2500
      wtype@wkHeight  = 2500
    else
      lthick          = 1.0
      BorderThick     = 1.0
    end if

    fe              = vars+"_bias_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)

    lndcolor        = "Grey25"
    plot            = new(nexps,graphic)
    color2          = read_colormap_file("cmocean_balance")
    ncols           = dimsizes(color2(:,0))
    color2(ncols/2-5:ncols/2+5,:) = 0

   end if 

   main_str = "Mean bias"
   exp_str  = exp_strs(k)
   left_str = ""

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   BorderThick = 1.0
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.024
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF

   mpres     = res
   resx      = res
   resx1     = res

   if(wtype.eq."png") then
    lthick = 3.0
   else
    lthick = 1.0
   end if

;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 30
   mpres@gsnMajorLonSpacing     = 45
   mpres@gsnMinorLatSpacing     = 15
   mpres@gsnMinorLonSpacing     = 15
  ;mpres@mpProjection           = "Robinson"       ; choose projection
   mpres@mpLabelsOn             = False
   mpres@mpPerimOn              = False
   mpres@mpFillOn               = False
   mpres@mpGeophysicalLineColor      = lndcolor  ; color of cont. outlines
   mpres@mpGeophysicalLineThicknessF = lthick          ; thickness of outlines
  ;mpres@mpOutlineOn            = True
  ;mpres@mpOutlineDrawOrder     = "PostDraw"
  ;mpres@mpFillDrawOrder        = "Predraw"
  ;mpres@mpGridAndLimbOn        = True
  ;mpres@mpGridLatSpacingF      = 90.0
  ;mpres@mpGridLonSpacingF      = 180.0
  ;mpres@mpGridLineThicknessF   = 3.0
  ;mpres@mpGridAndLimbDrawOrder ="PreDraw"
   mpres@mpCenterLonF           = 180.         ; default is 0 [GM]
   mpres@mpMinLatF              = latS
   mpres@mpMaxLatF              = latN
   mpres@mpMinLonF              = lonW
   mpres@mpMaxLonF              = lonE
  ;mpres@mpOceanFillColor       = "lightskyblue1"
  ;mpres@mpLandFillColor        = "gray33"

   mpres@tmYLLabelsOn     = True
   mpres@tmYLOn           = True
   mpres@tmYROn           = False

   plot(k) = gsn_csm_map(wks,mpres)

  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnAddCyclic          = True
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005
   res0@gsnStringFontHeightF  = FontHeightF*1.1
   res0@lbLabelFontHeightF    = FontHeightF*0.6
   res0@lbTitleFontHeightF    = FontHeightF*0.8

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = main_str + " of " + vars + " (" + varunt +")"
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.2
  ;res0@pmLabelBarWidthF         = 0.95
  ;res0@pmLabelBarHeightF        = 0.08

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color2
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = right_str
   res0@cnLevels             := cnlevs2
   res0@cnLinesOn             = False    ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1

   pltrgd  = gsn_csm_contour(wks,vdm,res0)
   overlay(plot(k) ,pltrgd)

   ;;;;;add a box for selected region;;;;

   if(l_plot_latlon_box) then
      regx = new((/4/),graphic)
      resb                  = True ; polyline mods desired
      resb@gsLineColor      = "Magenta" ; color of lines
      resb@gsLineThicknessF = 2.5 ; thickness of lines
      do ii = 0,3
       regx(ii)=gsn_add_polyline(wks,plot(k),xptsa(ii:ii+1),yptsa(ii:ii+1),resb)
      end do
   end if


  end do ; end of exps

 ; panel plot only resources
   pres1                            = True      ; mods desired
   pres1@gsnFrame                   = True      ; save panel until both ready
   pres1@gsnPanelMainString         = "" ;main_str
   pres1@gsnPanelMainFontHeightF    = 0.012
   pres1@gsnPanelLabelBar           = False                ; add common colorbar
  ;pres1@lbOrientation              = "Vertical"
   pres1@lbTitleString              = main_str + " of " + vars + " (" + varunt +")"
   pres1@lbTitleFontHeightF         = FontHeightF*0.6
   pres1@lbLabelFontHeightF         = FontHeightF*0.55 ;  font height
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

end

