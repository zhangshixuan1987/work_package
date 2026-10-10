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

  flmsk      = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/lnd_msk/landmsk.nc"
  data_dir   = "/global/cfs/cdirs/e3sm/zhan391/local_diag/fig_data/clim_rad"
  obs_name   = "OAFlux"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP"   , "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  left_strs  = exp_name
  time_tag   = "2001-2014"

  seasons    = (/"ANN", "DJF", "JJA"/)
  nseas      = dimsizes(seasons)
  maxlev     = 11

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions = (/"Global"/)
  nregs   = dimsizes(Regions)
  latS    = -90.
  latN    =  10.
  lonW    =  0.
  lonE    = 360.

  varList = (/"LHFLX","SHFLX"/)
  untList = (/"W m~S~-2~N~", "W m~S~-2~N~"/)
  nvar    = dimsizes(varList)

 ;declear the season or annual mean data to be used;;;;;
 do i = 0, nseas-1,1

  print("working on season "+ seasons(i))

  do j = 0,nvar-1,1

  vars  = varList(j)
  vunt  = untList(j)
  print("working on " + vars)

  do k = 0,nexps-1

   fl0  = systemfunc("ls -1 "+data_dir+"/"+obs_name+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
   fl1  = systemfunc("ls -1 "+data_dir+"/"+exp_name(k)+"_anom_" + vars+ "_MON_en*_"+time_tag+".nc")
   fm0  = addfile(flmsk,"r")
   f0   = addfiles(fl0,"r")
   f1   = addfiles(fl1,"r")
   ListSetType (f0, "join")
   ListSetType (f1, "join")

   vnm  = vars+"_o"
   vd0  = dim_avg_n_Wrap(f0[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
   vd1  = dim_avg_n_Wrap(f1[:]->$vnm$(:,:,{latS:latN},{lonW:lonE}),0)
   lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})
 
   vdat = vd0
   vdat = vd1 - vd0
   vdat@_FillValue = -9999.0
   lms0 = conform(vdat,lmsk,(/1,2/))
   vdat = where(lms0.le.0,vdat,vdat@_FillValue)
   delete([/lmsk,lms0/])

   print("working on season "+ seasons(i))
   if (seasons(i).eq."ANN") then
     var1 = month_to_annual(vdat,1)
   else
     var1 = month_to_season(vdat,seasons(i))
   end if

   vdm  = dim_avg_n_Wrap(var1,0)
   rad  = 4.0*atan(1.0)/180.0
   wgty = cos(f1[0]->lat({latS:latN})*rad)
   glm1 = dim_max_n(dim_max_n(vdm,1),0)
   glm2 = dim_min_n(dim_min_n(vdm,1),0)
   glm3 = wgt_areaave(vdm,wgty,1.0,0)
   rstr = "Max / Min / Mean:  "+ sprintf("%.2f",glm1)+" / "+sprintf("%.2f",glm2)+" / "+sprintf("%.2f",glm3)

   print(rstr)
   delete([/fl0,f0,fl1,f1,glm1,glm2,glm3/])
   delete([/vdat,var1,vd0,vd1/])

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

    fe              = vars+"_bias_vs_"+obs_name+"_"+seasons(i)
    wks             = gsn_open_wks(wtype,fe)
    lndcolor        = "Grey25"
    plot            = new(nexps,graphic)
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
   res0                          = res
   res0@gsnAddCyclic             = True
   res0@gsnLeftString            = ""
   res0@gsnRightString           = ""
   res0@gsnCenterString          = ""
   res0@tiMainString             = ""
   res0@tiMainOffsetYF           = 0.005

   res0@lbLabelBarOn             = True
  ;res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = lab_str 
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

   res0@cnFillOn                 = True     ; color plot desired
  ;res0@cnFillMode               = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor    = "White"
   res0@cnFillPalette            = color2
   res0@cnLevelSelectionMode     = "ExplicitLevels"
   res0@tiMainString             = "" 
   res0@gsnLeftString            = left_str
   res0@gsnRightString           = rstr
   res0@cnLevels                := cnlevs2
   res0@cnLinesOn                = False        ; True is default
   res0@cnInfoLabelOn            = False
   res0@cnLineLabelsOn           = False    ; turn off contour labels
   res0@cnLineColor              = "Black"                   ; color of second contours
   res0@cnLineThicknessF         = 2.                       ; line thickness
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
   pres1@lbTitleString              = lab_str 
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
   delete([/plot,res,res0,mpres,pres1/])

 end do 

end do 

end

