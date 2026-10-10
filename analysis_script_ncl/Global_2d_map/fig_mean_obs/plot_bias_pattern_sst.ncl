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

  l_plot_box = True
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim/"
  obs_name   = "HadSST" ;ERA5
  exp_name   = (/"ERA5", "v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"ERA5", "E3SM AMIP"   , "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  varList    = (/"SST"/)
  varUnit    = (/"~S~o~N~C"/)
  nvar       = dimsizes(varList)
  
  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  20.
  lonW       =   0.
  lonE       = 360.

 ;sub regions for analysis 
 ;western equatorial Pacific (Niño4 region, 16E-150W, 5S-5N)
 ;western and eastern equatorial Pacific (90W-150W,5-5N) 
 ;Southern Ocean 40-60S, 0-360
  subreg    = (/"Nino34",    "SO"/)
  slats     = (/    -5,      -60/)
  slatn     = (/     5,      -40/)
  slonw     = (/   270,      350/)
  slone     = (/   160,        0/)
  regcols   = (/"Black", "Cyan"/)
  reglthk   = (/3, 3/)
  nsreg     = dimsizes(subreg)

  do j = 0,nvar-1,1

  vars   = varList(j)
  vunt   = varUnit(j)
  print("working on " + vars)

  ;declear the season or annual mean data to be used;;;;;
  do i = 0, nseas-1,1

  do k = 0,nexps-1

   fl0  = systemfunc("ls -1 "+data_dir+"/"+obs_name+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
   f0   = addfiles(fl0,"r")
   ListSetType (f0, "join")

   fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
   f1   = addfiles(fl1,"r")
   ListSetType (f1, "join")
  
   vnm  = vars+"_o"
   vd0  = f0[:]->$vnm$(0,:,{latS:latN},{lonW:lonE})
   vd1  = f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE})

   vd0@_FillValue  = -99999
   vd1@_FillValue  = -99999

   vd1   = vd1 - conform(vd1,vd0,(/1,2,3/))
   vd1!0 = "ens"
   vd1!1 = "time"
   vd1!2 = "lat"
   vd1!3 = "lon"

   if(vars.eq."PSL") then
     vd1 = vd1/100.0
     vd1@units = "hPa"
     vunt    = "hPa"
   end if 

   if(vars.eq."TMQ") then
     vd1 = vd1/1000.0
     vd1@units = "kg m~S~-2~N~"
     vunt    = "kg m~S~-2~N~"
   end if

   if(vars.eq."TAUX".or.vars.eq."TAUY") then
     vd1 = vd1*100.0
     vd1@units = "10~S~-2~N~ N m~S~-2~N~"
     vunt    = "10~S~-2~N~ N m~S~-2~N~"
   end if

   print("working on season "+ seasons(i))
   if (seasons(i).eq."ANN") then 
     tmp = vd1(time|:,ens|:,lat|:,lon|:)
     va1 = month_to_annual(tmp,1)
   else 
     tmp = vd1(time|:,ens|:,lat|:,lon|:)
     va1 = month_to_season(tmp,seasons(i))
   end if 
   delete([/tmp/])
 
   vde  = dim_avg_n_Wrap(va1,0) 
   vdm  = dim_avg_n_Wrap(vde,0)
   rad  = 4.0*atan(1.0)/180.0
   wgty = cos(f1[0]->lat({latS:latN})*rad)
   nens = dimsizes(vde(:,0,0))
   tmp1 = dim_max_n(dim_max_n(vde,2),1) 
   tmp2 = dim_min_n(dim_min_n(vde,2),1) 
   tmp3 = wgt_areaave(vde,wgty,1.0,0)
   
   glm1  = sprintf("%3.1f", avg(tmp1)) + " ~F34~1~N~~F~ " + \
           sprintf("%3.2f", stddev(tmp1))
   glm2  = sprintf("%3.1f", avg(tmp2)) + " ~F34~1~N~~F~ " + \
           sprintf("%3.2f", stddev(tmp2))
   glm3  = sprintf("%3.1f", avg(tmp3)) + " ~F34~1~N~~F~ " + \
           sprintf("%3.2f", stddev(tmp3))
   right_str = "Max / Min / Mean:  "+ glm1+" / "+glm2+" / "+glm3
  
   delete([/fl0,f0,fl1,f1,tmp1,tmp2,tmp3,glm1,glm2,glm3/])
   delete([/vd0,vd1,vde,va1/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   if (k.eq.0) then 
    ;;;start to process the plot;;;;
    load "./share_fig_info.ncl"
    wtype           = "png"
    if(wtype.eq."png") then
      lthick          = 2.0
      BorderThick     = 2.0
      wtype@wkWidth   = 1200
      wtype@wkHeight  = 1200
    else
      lthick          = 1.0
      BorderThick     = 1.0
    end if

    fe              = vars+"_bias_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)

    lndcolor        = "Grey25"
    plot            = new(nexps,graphic)
    dum             = new((/nexps,nsreg/),graphic)
    color2          = read_colormap_file("cmocean_balance")
    ncols           = dimsizes(color2(:,0))
    color2(ncols/2-5:ncols/2+5,:) = 0

   end if 

   main_str = "Mean bias"
   exp_str  = exp_strs(k)
   left_str = exp_strs(k)
   lab_str  =  main_str + " of " + vars + " (" + vunt +")"

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False

   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeight = 0.018
   res@gsnStringFontHeightF = FontHeight*2.5
   res@tmXBLabelFontHeightF = FontHeight
   res@tmYLLabelFontHeightF = FontHeight
   res@tiMainFontHeightF    = FontHeight
   res@tiXAxisFontHeightF   = FontHeight
   res@tiYAxisFontHeightF   = FontHeight

   mpres     = res
   resx      = res
   resx1     = res
   
;---Control appearance of map.
   mpres@gsnMajorLatSpacing     = 20
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
 ;  mpres@mpLimitMode           = "LatLon"
   mpres@mpCenterLonF           = 180.0         ; default is 0 [GM]
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

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = main_str + " of " + vars + " (" + vunt +")"
  ;res0@lbTitlePosition          = "Right"     ; title location
  ;res0@lbTitleDirection         = "Across"    ; letter angle
  ;res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
   res0@pmLabelBarOrthogonalPosF = 0.25
   res0@pmLabelBarWidthF         = 0.62
   res0@pmLabelBarHeightF        = 0.15
   res0@lbLabelFontHeightF       = FontHeight*0.7
   res0@lbTitleFontHeightF       = FontHeight*0.75

   res0@cnFillOn              = True     ; color plot desired
  ;res0@cnFillMode            = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor = "White"
   res0@cnFillPalette         = color2
   res0@cnLevelSelectionMode  = "ExplicitLevels"
   res0@tiMainString          = "" ;exp_str
   res0@gsnLeftString         = left_str
   res0@gsnRightString        = right_str
   res0@cnLevels             := cnlevs2
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03

   pltrgd  = gsn_csm_contour(wks,vdm,res0)
   overlay(plot(k) ,pltrgd)

   printVarSummary(vdm)
   printMinMax(vdm&lat,False)
   printMinMax(vdm&lon,False)

   delete(vdm)

   if(l_plot_box) then 

    ;---Draw a box using gsn_add_polygon. See how the edges are curved.
    gnres                  = True
    gnres@gsFillColor      = "transparent"
    gnres@gsEdgesOn        = True
 
    txres                  = True
    txres@txFontHeightF    = 0.015
    txres@txJust           = "CenterLeft"

    do irg = 0,nsreg-2,1
      gnres@gsEdgeColor      = regcols(irg)
      gnres@gsLineColor      = regcols(irg)
      gnres@gsLineThicknessF = reglthk(irg)
      gnres@gsEdgeThicknessF = reglthk(irg)
      latbox = (/slats(irg),slats(irg),slatn(irg),slatn(irg),slats(irg)/)
      lonbox = (/slonw(irg),slone(irg),slone(irg),slonw(irg),slonw(irg)/)
      dum(k,irg)= gsn_add_polygon(wks, plot(k), lonbox, latbox, gnres)
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
   pres1@lbTitleString              = main_str + " of " + vars + " (" + vunt +")"
   pres1@lbTitleFontHeightF         = FontHeight*0.6
   pres1@lbLabelFontHeightF         = FontHeight*0.55 ;  font height
   pres1@lbBoxLineColor             = "black"
   pres1@lbBoxLinesOn               = True
  ;pres1@pmLabelBarParallelPosF     = 0.61
   pres1@pmLabelBarOrthogonalPosF   = 0.01
   pres1@pmLabelBarWidthF           = 0.95
   pres1@pmLabelBarHeightF          = 0.08
  ;pres1@gsnPanelBottom             = 0.45       ; space for label bar
  ;pres1@gsnPanelTop                = 0.85       ; only panel on lower half of page
  ;pres1@gsnPanelXF                 = (/0.07,0.57/)
   pres1@gsnPanelYWhiteSpacePercent = 2
   pres1@gsnPanelXWhiteSpacePercent = 2
   gsn_panel(wks,plot,(/nexps,1/),pres1)          ; create first panel
  
 end do 

end do 

end

