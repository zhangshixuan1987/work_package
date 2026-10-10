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
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"v2 AMIP"   , "v2 Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "1979-2014"

  seasons    = (/"MON"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions = (/"Global"/)
  nregs   = dimsizes(Regions)
  latS    = -90.
  latN    =  90.
  lonW    =  0.
  lonE    = 360.

  load "./share_var_info.ncl"
  nvar = dimsizes(varList)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do j = 0,nvar-1,1

  vars   = varList(j)
  varunt = untList(j)
  print("working on " + vars)

  do k = 0,nexps-1

  if(vars.eq."SST") then 
    fl0  = systemfunc("ls -1 "+data_dir+"/HadSST_anom_" + vars+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
    f0   = addfiles(fl0,"r")
    ListSetType (f0, "join")
  else
    fl0  = systemfunc("ls -1 "+data_dir+"/ERA5_anom_" + vars+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
    f0   = addfiles(fl0,"r")
    ListSetType (f0, "join")
  end if 

  fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_" +seasons(i)+"_en*_"+time_tag+".nc")
  f1   = addfiles(fl1,"r")
  ListSetType (f1, "join")
  
  time = f1[0]->time

  vnm  = vars+"_c"
  if(dimsizes(fl0).le.1) then 
    vd0  = dim_avg_n_Wrap(f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
  else 
    vd0  = dim_avg_n_Wrap(f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
  end if 

  if(dimsizes(fl1).le.1) then
    vd1  = dim_avg_n_Wrap(f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
  else
    vd1  = dim_avg_n_Wrap(f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
  end if 

  if(vars.eq."PSL") then
    vd0 = vd0/100.0
    vd1 = vd1/100.0
    vd0@units = "hPa"
    vd1@units = "hPa"
    varunt    = "hPa"
  end if 

  if(vars.eq."TMQ") then
    vd0 = vd0/1000.0
    vd1 = vd1/1000.0
    vd0@units = "kg m~S~-2~N~"
    vd1@units = "kg m~S~-2~N~"
    varunt    = "kg m~S~-2~N~"
  end if

  if(vars.eq."TAUX".or.vars.eq."TAUY") then
    vd0 = vd0*100.0
    vd1 = vd1*100.0
    vd0@units = "10~S~-2~N~ N m~S~-2~N~"
    vd1@units = "10~S~-2~N~ N m~S~-2~N~"
    varunt    = "10~S~-2~N~ N m~S~-2~N~"
  end if

  if (vars.eq."NETCF_SRF") then
    vd0 = vd0 * -1.0
  end if

  vdm  = vd1(0,:,:)
  vdm  = (/dim_avg_n_Wrap(vd1-vd0,0)/)
  rad  = 4.0*atan(1.0)/180.0
  wgty = cos(f1[0]->lat({latS:latN})*rad)
  nens = dimsizes(vd1(:,0,0))
  tmp1 = dim_max_n(dim_max_n(vd1-vd0,2),1) 
  tmp2 = dim_min_n(dim_min_n(vd1-vd0,2),1) 
  tmp3 = wgt_areaave(vd1-vd0,wgty,1.0,0)

  glm1  = sprintf("%3.1f", avg(tmp1)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp1))
  glm2  = sprintf("%3.1f", avg(tmp2)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp2))
  glm3  = sprintf("%3.1f", avg(tmp3)) + " ~F34~1~N~~F~ " + \
          sprintf("%3.2f", stddev(tmp3))
  right_str = "Max / Min / Mean:  "+ glm1+" / "+glm2+" / "+glm3

  delete([/fl0,f0,fl1,f1,vd0,vd1,tmp1,tmp2,tmp3,glm1,glm2,glm3/])

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
   res0@cnLinesOn             = False        ; True is default
   res0@cnInfoLabelOn         = False
   res0@cnLineLabelsOn        = False    ; turn off contour labels
   res0@cnLineColor           = "Black"                   ; color of second contours
   res0@cnLineThicknessF      = 2.                       ; line thickness
   res0@gsnContourZeroLineThicknessF = 5  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1

   pltrgd  = gsn_csm_contour(wks,vdm,res0)
   overlay(plot(k) ,pltrgd)
   delete(vdm)

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

end do 

end

