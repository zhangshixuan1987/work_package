;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;This script is used to plot the wrf simulations vs NEXRAD data ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  LOAD FUNCTIONS AND PROCEDURES
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

 datdir = "/pscratch/sd/z/zhan391/DARPA_project/post_process/SCREAM/"
 ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
 figkey    = "FutureALL_ARs"
 Explist   = (/"ERA5", "E3SM_ARNEP_721x1440_NDG","E3SM_ARNEP_721x1440_SSP245", \
               "E3SM_ARNEP_721x1440_SSP585","E3SM_ARNEP_721x1440_SSP245ML", "E3SM_ARNEP_721x1440_SSP585ML"/)
 Labels    = (/"ERA5","RRM_CTL","RRM_SSP245","RRM_SSP585","RRM_SSP245ML","RRM_SSP585ML"/)
 ngrps     = dimsizes(Explist)
 plot      = new (ngrps, graphic)
 vector    = new (ngrps, graphic)
 countour  = new (ngrps, graphic)
 countourx = new (ngrps, graphic)

 year      = 2017
 date      = yyyymmddhh_time(2017, 2017, 1, "integer")
 timtag    = "2017020300-2017020712"
 timplt    = date({2017020300:2017020712})
 ntimes    = dimsizes(timplt)

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

 ;********************************
 ; select sub-regions
 ;********************************
 minlat =  12.0 ;0.0    ; deg N
 maxlat =  60.0 ;57.0    ; deg N
 minlon = 198.0 ;-100.0    ; - deg E = deg W
 maxlon = 265.0 ;-20.0    ; - deg E = deg W
 cenlon = 200.0 ;

 do i = 0,ntimes-1,1

  do j = 0,ngrps-1,1

   expnam = Explist(j)
   flnm0  = datdir + expnam +"_IVT_"+timtag+".nc"
   f0     = addfile (flnm0, "r")
   utc_date = toint(timplt(i))
   year     = utc_date/1000000
   month    = (utc_date - year*1000000)/10000
   day      = (utc_date - year*1000000-month*10000)/100
   hour     = utc_date - year*1000000-month*10000-day*100
   minute   = 0
   second   = 0
   timstr   = sprinti("%0.2iZ ", hour) + sprinti("%0.2i ", day) + \
             month_abbr(month) + " "  + sprinti("%0.4i", year)
   timcur   = timplt(i)
   delete([/year,month,day,hour,minute,second/])
   print(timcur + " " +timstr)

   tdate    = cd_calendar(f0->time,-3)
   indtim   = ind(tdate.eq.timcur)

   if(expnam.eq."ERA5") then 
     unam = "VAR_10U"
     vnam = "VAR_10V"
   else 
     unam = "UBOT"
     vnam = "VBOT"
   end if 

   if(isfilevar(f0,"latitude")) then
     lat    = f0->latitude
     lon    = f0->longitude
   else
     lat    = f0->lat
     lon    = f0->lon
   end if

   u850   = f0->TUQ(indtim(0),:,:)
   v850   = f0->TVQ(indtim(0),:,:) 
   qtavg  = f0->PRECT(indtim(0),:,:)
   if(expnam.eq."ERA5") then
    qtavg = qtavg*1000.0 
    qtavg@units = "mm hr~S~-1~N~"
   else
    qtavg = qtavg * 1000.0 * 3600.0
    qtavg@units = "mm hr~S~-1~N~"
   end if 
   psrf   = u850
   psrf   = sqrt(u850*u850+v850*v850)
   delete([/tdate,indtim/])

  if(j.eq.0) then 
    ;;;generate the zoom in information;;;;
    fig_type  = "png" ;;figure types: pdf or png
    fig_dir   = "./prec_sub_fig/"
    fignam    = "fig_prec_"+timcur
    fig_type@wkWidth  = 2500
    fig_type@wkHeight = 2500 
    gridthick         = 2.0
    lthick            = 2.0
    wks               = gsn_open_wks(fig_type, fig_dir + fignam)
    plot      = new (ngrps, graphic)
    vector    = new (ngrps, graphic)
    countour  = new (ngrps, graphic)
    countourx = new (ngrps, graphic)

   ;cmap = read_colormap_file("WhViBlGrYeOrRe")
    cmap = read_colormap_file("BlAqGrYeOrReVi200")
    cmap(0,:) = 0.0 
    cnlevs = (/0.1,1,2.5,5,7.5,10,12,14,16,18,20,22/) ; "mm/hr" 
  end if 

;---Set common resources for all plots 
  res                  = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@gsnMaximize      = False 
  res@tmXTOn           = False
  res@tmYROn           = False
  res@gsnLeftString    = ""
  res@gsnRightString   = ""
  res@tiMainString         = ""

  BorderThick = 2.0
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
  res@gsnStringFontHeightF = FontHeightF*1.2
  res@tiXAxisFontHeightF   = FontHeightF
  res@tiYAxisFontHeightF   = FontHeightF

;---Control appearance of map.
  mresP = res
  mresP@gsnMajorLatSpacing          = 10
  mresP@gsnMajorLonSpacing          = 10
  mresP@gsnMinorLatSpacing          = 5
  mresP@gsnMinorLonSpacing          = 5
 ;mresP@mpProjection                = "Robinson"       ; choose projection
  mresP@mpLabelsOn                  = False
  mresP@mpPerimOn                   = False
  mresP@mpDataBaseVersion           = "MediumRes"     ; Alias 'MediumRes'
  mresP@mpDataSetName               = "Earth..1"
  mresP@mpFillOn                    = True    ; False to turn off gray continents
  mresP@mpFillDrawOrder             = "Draw"
  mresP@mpOutlineOn                 = True    ; turn on continental outline
  mresP@mpOutlineBoundarySets       = "AllBoundaries" ;"Geophysical" ;"AllBoundaries"
  mresP@mpLandFillColor             = "Grey90" ;"Yellow2enRod1"
  mresP@mpInlandWaterFillColor      = "PaleTurquoise3"
 ;mresP@mpOceanFillColor            = "PaleTurquoise3"
  mresP@mpGeophysicalLineColor      = "Grey27"
  mresP@mpGeophysicalLineThicknessF = 0.5
  mresP@mpUSStateLineColor          = "Grey27"
  mresP@mpUSStateLineThicknessF     = 0.5
  mresP@mpNationalLineColor         = "Grey27"
  mresP@mpNationalLineThicknessF    = 0.5
  mresP@mpGridAndLimbOn             = "False"
  mresP@mpGridAndLimbDrawOrder      = "PostDraw"
  mresP@mpGridMaskMode              = "MaskLand"
; mresP@mpGridSpacingF              = 5.0
  mresP@mpGridLineColor             = "Grey27"
  mresP@pmTickMarkDisplayMode       = "Always"
; mresP@mpCenterLonF                = 180.0
; mresP@mpCenterLatF                = -90.     ; This is necessary to get the correct map

 ;---Zoom in on plot
  mresP@mpMinLatF              = minlat
  mresP@mpMaxLatF              = maxlat
  mresP@mpMinLonF              = minlon
  mresP@mpMaxLonF              = maxlon
  mresP@mpCenterLonF           = 0.5*(minlon+maxlon) 
  plot(j)                      = gsn_csm_map(wks,mresP)

  res1                                = res 
  res1@cnFillOn                       = False        ; Turn on contour fill
  res1@cnFillPalette                  = "WhBlGrYeRe" ;"BlAqGrYeOrReVi200"
  res1@cnLevelSelectionMode           = "ManualLevels"
  res1@cnMinLevelValF                 =  0
  res1@cnMaxLevelValF                 =  1500.0
  res1@cnLevelSpacingF                =  100.0
  res1@gsnLeftString                  = ""
  res1@gsnRightString                 = "IVT(kg m~S~-1~N~ s~S~-1~N~, contour)" ;timstr
  res1@cnLinesOn                      = True   ; turn off contour lines
  res1@cnLineLabelsOn                 = True   ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@cnLineLabelDensityF            = 2.0
  res1@cnLineLabelFontHeightF         = FontHeightF*0.8
  res1@cnLineLabelBackgroundColor     = -1
  res1@gsnContourLineThicknessesScale = 2.0
  res1@cnLineThicknessF               = 1.0
  res1@cnLineColor                    = "Black"
  res1@gsnContourZeroLineThicknessF   = 2.
  res1@gsnContourNegLineDashPattern   = 1

  countourx(j) = gsn_csm_contour(wks,psrf,res1)
  overlay(plot(j),countourx(j))

;;set vector;;
  res_vc                            = res
  res_vc@vcGlyphStyle               = "LineArrow" ;"CurlyVector"
  res_vc@vcLineArrowThicknessF      = 3
  res_vc@vcLineArrowHeadMaxSizeF    = 0.008
  res_vc@vcMinDistanceF             = 0.025
  res_vc@vcRefLengthF               = 0.04
  res_vc@vcRefAnnoOn                = True
  res_vc@vcRefAnnoSide              = "Bottom"
  res_vc@vcRefAnnoString2On         = False
  res_vc@vcRefAnnoPerimOn           = True
  res_vc@vcRefAnnoOrthogonalPosF    = -0.99
  res_vc@vcRefAnnoParallelPosF      = 0.10
  res_vc@vcRefAnnoBackgroundColor   = "White";"Transparent"; "Purple"
  res_vc@vcVectorDrawOrder          = "PostDraw"
  res_vc@gsnRightString             = "" 
  res_vc@vcRefAnnoFontHeightF       = FontHeightF*0.65
  res_vc@vcRefMagnitudeF            = 500
  res_vc@vcRefAnnoString1           = res_vc@vcRefMagnitudeF

  vector(j) = gsn_csm_vector(wks,u850,v850,res_vc)

;---Temperature filled contour plot
  ref_res                       = res
  ref_res@cnFillOn              = True  
  ref_res@cnLevelSelectionMode  = "ExplicitLevels"
 ;ref_res@cnFillMode            = "RasterFill"
  ref_res@cnInfoLabelOn         = False
  ref_res@cnLineLabelsOn        = False
  ref_res@cnLinesOn             = False
  ref_res@cnLineLabelBackgroundColor = -1    ; transparent
  ref_res@cnFillDrawOrder       = "Draw"

  ref_res@lbLabelBarOn             = False   ; Will turn on in panel later
  ref_res@lbLabelFontHeightF       = FontHeightF*0.6
  ref_res@lbTitleFontHeightF       = FontHeightF*0.8
  ref_res@pmLabelBarParallelPosF   = 0.50 ;61
  ref_res@pmLabelBarOrthogonalPosF = 0.20
 ;ref_res@lbOrientation            = "Horizonal" ; "vertical"  ; vertical label bar
  ref_res@pmLabelBarHeightF        = 0.12 ;0.60
  ref_res@pmLabelBarWidthF         = 0.65 ;0.06
 ;ref_res@lbTitleString            = "" ;Varunt;Varstr+ " ("+Varunt+")";
 ;ref_res@lbTitlePosition          = "Right"         ; title location
 ;ref_res@lbTitleDirection         = "Across"        ; letter angle
 ;ref_res@lbTitleAngleF            = 90.          ; title angle

  ref_res@cnFillPalette       := cmap 
  ref_res@cnLevels            := cnlevs 
  ref_res@cnLevelSelectionMode:= "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF    := lthick
  ref_res@gsnRightString      := "" 
  ref_res@gsnLeftString       := Labels(j)

  countour(j) = gsn_csm_contour(wks,qtavg,ref_res)

  overlay(plot(j),countour(j))
  overlay(plot(j),vector(j))

  delete([/lat,lon,u850,v850,psrf,qtavg,res,ref_res,res,res_vc/])
 end do 

  resP                     = True
  resP@gsnMaximize              = True
  resP@gsnPaperOrientation      = "portrait"
  resP@gsnPanelLabelBar         = True
  resP@pmLabelBarParallelPosF   =  0.01
  resP@pmLabelBarOrthogonalPosF = -0.02
  resP@lbLabelFontHeightF       = FontHeightF*0.4
  resP@lbTitleFontHeightF       = FontHeightF*0.45
  resP@lbTitleString            = "PRECT (mm hr~S~-1~N~)"
  resP@lbOrientation            = "horizontal"  ; vertical label bar
  resP@pmLabelBarHeightF        = 0.05
  resP@pmLabelBarWidthF         = 0.90
  resP@gsnPanelMainString       = "Current time: " + timstr
  resP@gsnPanelMainFontColor    = "Black"
  resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6
  resP@gsnPanelMainPosXF        = 0.85
 ;resP@txPosYF                  = 0.7
  resP@gsnPanelBottom           = 0.3
  resP@gsnPanelYWhiteSpacePercent = 0.
  resP@gsnPanelXWhiteSpacePercent = 3.
  gsn_panel(wks,plot,(/ngrps,2/),resP)

end do 

end
