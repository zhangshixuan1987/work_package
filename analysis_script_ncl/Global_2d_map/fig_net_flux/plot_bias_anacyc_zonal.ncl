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
  data_dir   = "./"
  exp_name   = (/"v2.LR.amip", "v2.LR.historical"/)
  exp_strs   = (/"E3SM AMIP", "E3SM Coupled"/)
  nexps      = dimsizes(exp_name)
  time_tag   = "2001-2014"

  ylabels    = (/"Jan","Feb","Mar","Apr","May","Jun",\
                 "Jul","Aug","Sep","Oct","Nov","Dec"/)
  ystring    = "Zonally Averaged Annual Cycle"
  maxlev     = 11

  varList    = (/"dTs","dFQn","dSAF","dCRE","dnSAF","dFLW","dFLWC","dLH","dSH"/)
  vstList    = (/"~F33~D~N~~F~TS","~F33~D~N~~F~Qn","~F33~D~N~~F~SAF(SW)","~F33~D~N~~F~CRE", \
                 "~F33~D~N~~F~no_SAF(SW)","~F33~D~N~~F~FLW(down)","~F33~D~N~~F~FLWC(down,clear-sky)", \
                 "-~F33~D~N~~F~(LH)","-~F33~D~N~~F~(SH)"/)
  varUnit    = (/"~S~o~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C","W m~S~-2~N~C"/)
  nvar       = dimsizes(varList)

  ;;;;;;;;;;select the regions for plot;;;;;;;;;;;;;
  Regions    = (/"SH"/)
  nregs      = dimsizes(Regions)
  latS       = -90.
  latN       =  10.
  lonW       =   0.
  lonE       = 360.

  do j = 0,nvar-1,1

   vars   = varList(j)
   vunt   = varUnit(j)
   print("working on " + vars)

   do k = 0,nexps-1

    fm0  = addfile(flmsk,"r")
    lmsk = fm0->LANDFRAC(0,{latS:latN},{lonW:lonE})

    fl0  = data_dir+exp_name(k)+"_energy_MON_"+time_tag+".nc"
    f0   = addfile(fl0,"r") 

    if(vars.eq."dSUM") then 
      vdat = f0->dSAF1(:,:,{latS:latN},{lonW:lonE})
      vdat = vdat + f0->dSAF2(:,:,{latS:latN},{lonW:lonE})
      do i = 3,nvar-2,1
       vnm  = varList(i)
       vdat = vdat + f0->$vnm$(:,:,{latS:latN},{lonW:lonE})
      end do
    else if(vars.eq."dSAF") then
      vdat = f0->dSAF1(:,:,{latS:latN},{lonW:lonE})
      vdat = vdat + f0->dSAF2(:,:,{latS:latN},{lonW:lonE})
    else 
      vnm  = varList(j)
      vdat = f0->$vnm$(:,:,{latS:latN},{lonW:lonE})
    end if 
    end if 

    lms0 = conform(vdat,lmsk,(/2,3/))
    vdat = where(lms0.le.0,vdat,vdat@_FillValue)
    delete([/lmsk,lms0/])

    vmn = dim_avg_n_Wrap(vdat,0) 
    vsd = dim_stddev_n_Wrap(vdat,0)
    vmn1 = clmMonTLL(vmn)
    vsd1 = clmMonTLL(vsd)
    vdm  = dim_avg_n_Wrap(vmn1,2)
    vds  = dim_avg_n_Wrap(vsd1,2)
    delete([/fl0,f0,vdat,vmn,vsd,vmn1,vsd1/])

   ;;;;;;;;;;;;;;;;;;;;;;;;;;
   ; Start to process plot 
   ;;;;;;;;;;;;;;;;;;;;;;;;;;
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

   fe              = vars+"_bias_zonal_"+exp_name(k)+"_anacyc"
   wks             = gsn_open_wks(wtype,fe)
   lndcolor        = "Grey25"
   plot            = new(nexps,graphic)
   dumx            = new(nexps,graphic)
   dum0            = new(nexps,graphic)
   dum1            = new(nexps,graphic)

   main_str = "mean bias"
   exp_str  = exp_strs(k)
   left_str = exp_strs(k)
   lab_str  = "~F33~D~N~~F~T (~S~o~N~C)" 
   rstr     = vstList(j)

   color          = read_colormap_file("cmp_b2r") 
   color(31:32,:) = 0
   cnlevs1        = (/0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0/)
   cnlevs2        = (/-6.0,-5.0,-4.0,-3.0,-2.0,-1.0,-0.2,0.2,1.0,2.0,3.0,4.0,5.0,6.0/)
   minlev         =  0
   maxlev         =  2
   intlev         =  0.2

  ;============================================================
  ; PLOTS
  ;============================================================
   res                = True     ; Plot modes desired.
   res@gsnDraw        = False
   res@gsnFrame       = False
   res@gsnMaximize    = False     ; Maximize plot
   res@tmXTOn         = False
   res@tmYROn         = False
   res@vpWidthF       = 0.6
   res@vpHeightF      = 0.6

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
   res@gsnStringFontHeightF  = FontHeight*1.1
   res@tmXBLabelFontHeightF  = FontHeight
   res@tmYLLabelFontHeightF  = FontHeight
   res@tiMainFontHeightF     = FontHeight
   res@tiXAxisFontHeightF    = FontHeight*1.1
   res@tiYAxisFontHeightF    = FontHeight*1.1
   res@lbLabelFontHeightF    = FontHeight*0.95
   res@lbTitleFontHeightF    = FontHeight*1.05
   
   res@tiXAxisString    = "" ;"Latitude"
   res@tiXAxisOffsetYF  = 0.11
   res@tiYAxisString    = "Month"

   res@tmXBMode         = "Explicit"
   res@tmXBValues       = ispan(-90,0,15)
   res@tmXBMinorValues  = ispan(-90,0,5)
   res@tmXBLabels       = (/"90S","75S","60S","45S","30S","15S","0"/)

   res@tmYLMode         = "Explicit"      
   res@tmYLValues       = ispan(0,11,1)
   res@tmYLLabels       = ylabels
   res@tmYLLabelJust    = "CenterLeft"
  ;*******************************************
  ; first plot
  ;*******************************************
   res0                       = res
   res0@gsnLeftString         = ""
   res0@gsnRightString        = ""
   res0@gsnCenterString       = ""
   res0@tiMainString          = ""
   res0@tiMainOffsetYF        = 0.005

   res0@lbLabelBarOn             = True
   res0@lbOrientation            = "Vertical"
   res0@lbTitleString            = lab_str 
   res0@lbTitlePosition          = "Right"     ; title location
   res0@lbTitleDirection         = "Across"    ; letter angle
   res0@lbTitleAngleF            = 90.         ; title angle
   res0@lbBoxLineColor           = "black"
  ;res0@lbBoxLinesOn             = False
  ;res0@pmLabelBarParallelPosF   = 0.61
  ;res0@pmLabelBarOrthogonalPosF = 0.10
   res0@pmLabelBarWidthF         = 0.08
   res0@pmLabelBarHeightF        = 0.62
   res0@lbTitlePosition          = "Right"

   res0@cnFillOn                 = True     ; color plot desired
  ;res0@cnFillMode               = "RasterFill"   ; turn raster on
   res0@cnMissingValFillColor    = "White"
   res0@cnFillPalette            = color
   res0@cnLevelSelectionMode     = "ExplicitLevels"
   res0@tiMainString             = "" ;exp_str
   res0@gsnLeftString            = left_str
   res0@gsnRightString           = rstr
   res0@cnLevels                := cnlevs2
   res0@cnLinesOn                = False        ; True is default
   res0@cnInfoLabelOn            = False
   res0@cnLineLabelsOn           = False    ; turn off contour labels
   res0@cnLineColor              = "Black"               ; color of second contours
   res0@cnLineThicknessF         = lthick           ; line thickness
   res0@gsnContourZeroLineThicknessF = lthick*2.0  ; eliminates zero contour
   res0@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res0@gsnRightStringOrthogonalPosF = 0.03
   res0@gsnLeftStringOrthogonalPosF  = 0.03
   plot(k)  = gsn_csm_contour(wks,vdm,res0)

   res1 = res0
   res1@cnFillOn                     = False     ; color plot desired
   res1@tiMainString                 = "" ;exp_str
   res1@gsnLeftString                = ""
   res1@gsnRightString               = ""
  ;res1@cnLevelSelectionMode         = "ManualLevels"
  ;res1@cnMinLevelValF               = minlev
  ;res1@cnMaxLevelValF               = maxlev
  ;res1@cnLevelSpacingF              = intlev
   res1@cnLevelSelectionMode         = "ExplicitLevels"
   res1@cnLevels                    := cnlevs1
   res1@cnLinesOn                    = True
   res1@cnInfoLabelOn                = False
   res1@cnLineLabelsOn               = True    ; turn off contour labels
   res1@cnLineLabelFontHeightF       = FontHeight*0.8
   res1@cnConstFLabelPerimOn         = False
  ;res1@cnLineLabelBackgroundColor   = "Transparent"
   res1@cnLineLabelFormat            = "0@;*.2f"                      ; Only one value after decimal point
  ;res1@cnLineLabelDensityF          = 2.0
   res1@cnLineColor                  = "Black"                   ; color of second contours
   res1@cnLineThicknessF             = lthick*2.0                ; line thickness
   res1@gsnContourZeroLineThicknessF = lthick*2.5
   res1@gsnContourNegLineDashPattern = 1  ; sets negative contours to dash pattern 1
   res1@gsnRightStringOrthogonalPosF = 0.03
   res1@gsnLeftStringOrthogonalPosF  = 0.03
   dumx(k)  = gsn_csm_contour(wks,vds,res1)
   overlay(plot(k),dumx(k))

  ;==============================================
  ; create xy plot
  ;=============================================
   xyres                  = True         ; xy plot mods desired
   xyres@vpWidthF         = .20          ; set width of second plot
   xyres@vpWidthF         = 0.6
   xyres@vpHeightF        = 0.2

   xyres@gsnDraw          = False        ; don't draw yet
   xyres@gsnFrame         = False        ; don't advance frame yet

   xyres@tiXAxisString    = "" ;"Latitude"
   xyres@tiXAxisOffsetYF  = 0.11
   xyres@tiYAxisString    = lab_str

   xyres@gsnCenterString  = ""  ; add title
   xyres@txFontHeightF    = .015         ; change font height

   xyres@gsnYRefLine             = 0.0          ; create a reference line   
   xyres@gsnYRefLineThicknessF   = 3.0
   xyres@gsnYRefLineDashPatterns = 5
   xyres@gsnYRefLineColors       = "Black"

   xyres@trXMinF          = min(vdm&lat)
   xyres@trXMaxF          = max(vdm&lat)
   xyres@tmXBMode         = "Explicit"
   xyres@tmXBValues       = ispan(-90,0,15)
   xyres@tmXBMinorValues  = ispan(-90,0,5)
   xyres@tmXBLabels       = (/"90S","75S","60S","45S","30S","15S","0"/)

   xyres@trYMinF          = -12
   xyres@trYMaxF          =  12
   xyres@tmYLMode         = "Explicit"
   xyres@tmYLValues       = ispan(-20,20,4)
   xyres@tmYLMinorValues  = ispan(-20,20,1)
   xyres@tmYLLabels       = (/"","","","-8","-4","0","4","8","","",""/)

   xyres@xyLineThicknessF = 2.0
   xyres@xyLineColor      = "Black"
   xyres@xyMarkLineMode   = "MarkLines"               ; Markers *and* lines
   xyres@xyMarkers        = 16
   xyres@xyMarkerColors  := xyres@xyLineColor
   xyres@xyDashPattern    = 0                         ; Make curves all solid

   y       = dim_avg_n_Wrap(vdm,0)           ; average chi
   x       = vdm&lat
   dum0(k) = gsn_csm_xy(wks, x,y,xyres)

   xyres1  = xyres
   xyres1@gsnAttachPlotsXAxis = True 
   dum1(k) = gsn_attach_plots(plot(k),(/dum0(k)/),res0,xyres1)
  ;overlay(plot(k),dum(k))

   draw(plot(k))
   frame(wks)

   delete([/vdm,vds,dumx,dum0,dum1/])

  end do ; end of exps

end do 

end

