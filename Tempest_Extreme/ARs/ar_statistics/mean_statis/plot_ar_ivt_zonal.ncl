;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;This script used to calculate and plot the zonal mean variables from CAM;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;The terms zonal and meridional are used to describe directions on a globe;;;;;;;
;;;;;;;Zonal means "along a latitudinal circle" or "in the west–east direction";;;;;;;;
;;;;;;;while meridional means "along a longitudinal circle" (a.k.a. meridian) ;;;;;;;;;
;;;;;;;or "in the north–south direction" [2];;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

 ;;; flags to control the plot style
  l_plot_latlon_box = False
  datdir      = "/global/cscratch1/sd/zhan391/DARPA_project/darpa_project/TempestExtremes/ARs/global_mean/"
  outdir      = "./"
  wgtFileDir  = "./"
  wgtFileName = "map_ne30pg2_to_cmip6_180x360_aave.20200201.nc"

  ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
  Casdir  = (/"ERA5", "CLIM", "NDGUV3_tau06", "NDGUVTQ3_tau06"/)
  Groups  = (/"ERA5", "CLIM", "NDGUV3_tau06", "NDGUVTQ3_tau06"/)
  Labstr  = (/"ERA5", "E3SM (CLIM)", "E3SM (NDGUV)", "E3SM (NDGUVTQ)"/)
  ngrps   = dimsizes(Groups)
  year    = 2017
  timplt  = 2017

  do i = 0,ngrps-1,1

   fna = systemfunc("ls -1 " + datdir+Casdir(i)+"_zonal_mean_"+"*.nc")
   fla = addfiles(fna, "r")
   print(fna)

   lat    = fla[0]->lat
   var0   = fla[:]->VIWVN_PW_AR_zon(:,:)  
   var1   = fla[:]->VIWVN_PW_NONAR_zon(:,:)  
   
   data = new((/2,dimsizes(lat)/),typeof(var0))
   data(0,:) = dim_avg_n(var0,0)
   data(1,:) = dim_avg_n(var1,0)
   data!0    = "ngrp"
   data!1    = "lat"
   data&lat  = lat

   labels    = "  " +(/"AR points", "Non-AR points"/)
   printMinMax(data,False)
   delete([/fna,fla,var0,var1/])

   ;----------------------------------------------------------------------
   ; Plotting section
   ;----------------------------------------------------------------------

   if(i.eq.0)then
    ;;;;;;;open a wks for plot;;;;;;;;;;;;;;;;;;;
    wtype = "pdf"
    wtype@wkPaperSize  = "A4"
   ;wtype@wkOrientation   = "landscape" ;;pictures displayed along the length
 
    plotname = "ar_ivt_zonal_mean"
    wks1 = gsn_open_wks(wtype,plotname)
    plot = new (ngrps, graphic)
    dum1 = new (ngrps, graphic)
    dum2 = new (ngrps, graphic)
    print(" ")
    print("........................... ")
    print("plot : " + plotname+".pdf")
    print("........................... ")
   end if 

   ;;;set up the plot information;;;;
   load "./var_share_colorbar.ncl"

   Figstr = "" 

;---Resources to share between both plots
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False
   res@vpWidthF       = 0.8
   res@vpHeightF      = 0.6

   BorderThick = 0.5
   res@tmBorderThicknessF    = BorderThick
   res@tmXBMajorThicknessF   = BorderThick
   res@tmXBMinorThicknessF   = BorderThick*0.5
   res@tmYLMajorThicknessF   = BorderThick
   res@tmYLMinorThicknessF   = BorderThick*0.5
   res@tmXTMajorThicknessF   = BorderThick
   res@tmXTMinorThicknessF   = BorderThick*0.5
   res@tmYRMajorThicknessF   = BorderThick
   res@tmYRMinorThicknessF   = BorderThick*0.5

   FontHeightF = 0.028
   res@tiMainFontHeightF    = FontHeightF
   res@tmYLLabelFontHeightF = FontHeightF
   res@tmXBLabelFontHeightF = FontHeightF
   res@gsnStringFontHeightF  = FontHeightF*1.1

   res@tmXTOn                 = False
   res@tmYROn                 = False

   res@tiMainString       = ""
   res@gsnLeftString                 = "";leftstring ;var_name0    ; long_name is too long!
   res@gsnRightString                = "";rightstring
   res@gsnCenterString               = ""
   res@gsnRightStringOrthogonalPosF  = 0.04
   res@gsnLeftStringOrthogonalPosF   = 0.04
   res@gsnCenterStringOrthogonalPosF = 0.04

  ;add a legend
   if(i.eq.0) then 
     res@pmLegendDisplayMode    = "Always"              ; turn on legend
   else
     res@pmLegendDisplayMode    = "Never"              ; turn on legend
   end if 
   res@pmLegendSide             = "Top"                 ; Change location of
   res@pmLegendParallelPosF     = .30                   ; move units right
   res@pmLegendOrthogonalPosF   = -1.2                  ; more neg = down
   res@pmLegendWidthF           = 0.14                  ; Change width and
   res@pmLegendHeightF          = 0.16                  ; height of legend.
   res@lgLabelFontHeightF       = FontHeightF*0.9       ; change font height
   res@lgPerimOn                = False                 ; no box around
   res@xyExplicitLegendLabels   = "  "+labels


   res@xyMarkLineMode         = "Lines"               ; Markers *and* lines
   res@xyMarkerColors         = (/colors(i),colors(i)/)
   res@xyMarkers              = (/marker(i),marker(i)/)
   res@xyMarkerSizes          = (/mksizf(i),mksizf(i)/)
   res@xyMarkerThicknesses    = (/mkthik(i),mkthik(i)/)
   res@xyMarkerOpacities      = (/mkopcy(i),mkopcy(i)/)
   res@xyLineColors           = (/colors(i),colors(i)/)
  ;res@xyDashPatterns         = (/lnpatt(i),lnpatt(i)/)
   res@xyLineThicknesses      = (/lnthik(i),lnthik(i)/)
   res@xyDashPatterns         = (/0,2/)
   res@xyLineOpacities        = (/lnopcy(i),lnopcy(i)/)

   res1 = res 
   res1@trYMaxF                =  10.                    ; axis max
   res1@trYMinF                = -40.                    ; axis min
   res1@tmYLMode               = "Explicit"                ; explicit labels
   res1@tmYLValues             = fspan(res1@trYMinF,res1@trYMaxF,6)
   res1@tmYLMinorValues        = fspan(res1@trYMinF,res1@trYMaxF,16)
   res1@tmYLLabels             = ""+sprintf("%5.0f",res1@tmYLValues)
   res1@trXMaxF               = -25                  ; axis max
   res1@trXMinF               = -90                  ; axis min
   res1@tmXBMode              = "Explicit"                ; explicit labels
   res1@tmXBValues            = (/-90, -75, -60, -45, -30/)
   res1@tmXBMinorValues       = fspan(res1@trXMinF,res1@trXMaxF,14)
   res1@tmXBLabels            = ""+abs(res1@tmXBValues)+"S"
   res1@tiXAxisString         = "Latitude"
   res1@gsnLeftString         = "Southern Hemisphere" 
   res1@tiYAxisString         = "Zonal IVT (kg m~S~-2~N~ s~S~-1~N~)"
   dum1(i) = gsn_csm_xy(wks1,lat({-90:-25}),data(:,{-90:-25}), res1)

   res2 = res
   res2@trYMaxF                =  30.                    ; axis max
   res2@trYMinF                = -10.                    ; axis min
   res2@tmYLMode               = "Explicit"                ; explicit labels
   res2@tmYLValues             = fspan(res2@trYMinF,res2@trYMaxF,5)
   res2@tmYLMinorValues        = fspan(res2@trYMinF,res2@trYMaxF,13)
   res2@tmYLLabels             = ""+sprintf("%5.0f",res2@tmYLValues)
   res2@trXMaxF               = 90                  ; axis max
   res2@trXMinF               = 25                  ; axis min
   res2@tmXBMode              = "Explicit"                ; explicit labels
   res2@tmXBValues            = (/30, 45, 60, 75, 90/)
   res2@tmXBMinorValues       = fspan(res1@trXMinF,res1@trXMaxF,14)
   res2@tmXBLabels            = ""+abs(res2@tmXBValues)+"N"
   res2@tiXAxisString         = "Latitude"
   res2@tiYAxisString         = "Zonal IVT (kg m~S~-2~N~ s~S~-1~N~)"
   res2@gsnLeftString         = "Northern Hemisphere" 
   dum2(i) = gsn_csm_xy(wks1,lat({25:90}),data(:,{25:90}), res2)


  if(i.eq.0) then 
    plot(0) = dum1(i)
    plot(1) = dum2(i)
  else
    overlay(plot(0),dum1(i))
    overlay(plot(1),dum2(i))
  end if 

  delete([/res,res1/])
  delete([/lat,data/])
end do ; end of experiment group loop

  gres = True
  gres@YPosPercent = 95.               ; expressed as %, 0->100, sets position of bottom border of legend
  gres@XPosPercent = 50.               ; expressed as %, 0->100, sets position of left border of legend
  gres@Position    = "Top"            ; YPosPercent setting reflects the "Top" or "Bottom" of legend
  gres@ItemSpacePercent = 8.          ; expressed as %, 0->100, space between legend items (Default = 5.)
  gres@LineLabelWhiteSpacePercent = 2 ; expressed as %, 0->100, space between line and label (Default = 3.)

  lineres = True
  lineres@lgLineColors      = colors
  lineres@lgLineOpacity     = lnopcy
  lineres@lgLineThicknesses = lnthik
  lineres@lgDashIndexes     = lnpatt
  lineres@LineLengthPercent = 8                        ; expressed as %, 0->100, length of line

  textres = True
  textres@lgLabels           = Labstr
  textres@lgLabelFontHeights = FontHeightF*0.9 
  textres@lgLabelColors      = colors
  textres@lgLabelOpacity     = lnopcy

 ;plot(0) = simple_legend(wks1,plot(0),gres,lineres,textres)
  plot(1) = simple_legend(wks1,plot(1),gres,lineres,textres)

;---Draw both plots in a panel
;***************************************
; panel first set of plots
;***************************************
  pres1                          = True
  pres1@gsnFrame                 = True      ; don't advance frame yet
  pres1@gsnMaximize              = False
  pres1@gsnPanelLabelBar         = False
  pres1@pmLabelBarParallelPosF   = -0.010
  pres1@pmLabelBarOrthogonalPosF = 0.0
  pres1@lbLabelFontHeightF       = FontHeightF * 0.45
  pres1@lbOrientation            = "vertical"  ; vertical label bar
  pres1@pmLabelBarHeightF        = 0.21
  pres1@pmLabelBarWidthF         = 0.05
  pres1@gsnPanelMainString       = " "
  pres1@gsnPanelMainFontColor    = "Black"
  pres1@gsnPanelMainFontHeightF  = 0.028
  pres1@gsnPanelYWhiteSpacePercent = 0.
  pres1@gsnPanelXWhiteSpacePercent = 3.
  gsn_panel(wks1,(/plot(0:1)/),(/1,2/),pres1)

end
