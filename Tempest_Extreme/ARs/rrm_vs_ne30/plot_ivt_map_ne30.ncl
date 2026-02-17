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
 figkey    = "NDGALL_ARs"
 Explist   = (/"ERA5","E3SM_NE30_721x1440_CTL", "E3SM_NE30_721x1440_NDG"/)
 Labels    = (/"ERA5","NE30_CTL","NE30_NDG"/)
 ngrps     = dimsizes(Explist)

 plot      = new (ngrps, graphic)
 vector    = new (ngrps, graphic)
 countour  = new (ngrps, graphic)
 countourx = new (ngrps, graphic)

 year      = 2017
 date      = yyyymmddhh_time(2017, 2017, 1, "integer")
 timtag    = "2017020300-2017021112"
 timplt    = date({2017020600:2017021112})
 ntimes    = dimsizes(timplt)
 tptag     = "raw_0206-0211"

 month_abbr = (/"","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep", \
                   "Oct","Nov","Dec"/)

 ;********************************
 ; select sub-regions
 ;********************************
 minlat     =  10.0 ;0.0    ; deg N
 maxlat     =  60.0 ;57.0    ; deg N
 minlon     = 198.0 ;-100.0    ; - deg E = deg W
 maxlon     = 265.0 ;-20.0    ; - deg E = deg W
 cenlon     = 200.0 ;

 do j = 0,ngrps-1,1

  expnam   = Explist(j)
  flnm0    = datdir + expnam +"_IVT_"+timtag+".nc"
  f0       = addfile (flnm0, "r")
 
  tdate    = cd_calendar(f0->time,-3)
  indtim   = ind (tdate.ge.timplt(0).and.tdate.le.timplt(ntimes-1))

  if(isfilevar(f0,"latitude")) then 
    lat    = f0->latitude
    lon    = f0->longitude
  else 
    lat    = f0->lat
    lon    = f0->lon
  end if 

  nvars    = 3
  nlat     = dimsizes(lat)
  nlon     = dimsizes(lon)

  if(dimsizes(dimsizes(f0->TUQ)).ge.3)
   var     = new((/nvars,nlat,nlon/),typeof(lat))
   var!1   = "lat"
   var!2   = "lon"
   var&lat = lat
   var&lon = lon
   var(0,:,:) = (/dim_avg_n(f0->TUQ(indtim,:,:),0)/)     ; K  (time,lev,lat,lon)
   var(1,:,:) = (/dim_avg_n(f0->TVQ(indtim,:,:),0)/)     ; K  (time,lev,lat,lon)
   var(2,:,:) = (/sqrt(var(0,:,:)^2+var(0,:,:)^2)/)
  else
   var        = new((/nvars,nlat/),typeof(lat))
   var@lat    = lat 
   var@lon    = lon
   var(0,:)   = dim_avg_n(f0->TUQ(indtim,:),0)     ; K  (time,lev,lat,lon)
   var(1,:)   = dim_avg_n(f0->TVQ(indtim,:),0)     ; K  (time,lev,lat,lon)
   var(2,:)   = (/sqrt(var(0,:)^2+var(0,:)^2)/)
  end if 
  var@FillValue = -9999.0
  delete([/tdate,indtim/])

  if(j.eq.0) then 
    vob = var 
  end if 

  ;calculate the global/regional mean values
  rad   = 4.*atan(1.)/180.    ;
  wgty  = lat
  wgty  = cos(lat*rad)      ; cosine weights
  if(dimsizes(dimsizes(var)).ge.3)
    gomn  = wgt_areaave_Wrap(vob(2,{minlat:maxlat},{minlon:maxlon}),\
                             wgty({minlat:maxlat}), 1.0, 1)
    gosd  = stddev(vob(2,{minlat:maxlat},{minlon:maxlon}))
    gbmn  = wgt_areaave_Wrap(var(2,{minlat:maxlat},{minlon:maxlon}),\
                             wgty({minlat:maxlat}), 1.0, 1)
    grms  = wgt_arearmse(var(2,{minlat:maxlat},{minlon:maxlon}),\
                         vob(2,{minlat:maxlat},{minlon:maxlon}), \
                         wgty({minlat:maxlat}), 1.0, 1)
  else
    lon   = where(lon.le.0,lon+360.0,lon)
    indc  = ind(lat.ge.minlat.and.lat.le.maxlat.and.lon.ge.minlon.and.lon.le.maxlon)
    gomn  = dim_avg_wgt_n(vob(2,indc), wgty(indc), 0.0, 0)
    gosd  = stddev(vob(2,indc))
    gbmn  = dim_avg_wgt_n(var(2,indc), wgty(indc), 0.0, 0)
    grms  = dim_avg_wgt_n((var(2,indc)-vob(2,indc))^2, wgty(indc), 0.0, 0)
  end if 

  if(j.eq.0) then 
    rstr  = "Mean/RMS (IVT): "+ sprintf("%.2f",gomn) + "/" + sprintf("%.2f",gosd)
  else
    rstr  = "Bias/RMSE (IVT): "+ sprintf("%.2f",gbmn-gomn) + "/" + sprintf("%.2f",grms)
  end if 
  delete([/rad,wgty,gbmn,grms/])

 ;printVarSummary(var)
 ;wrf_smooth_2d(var, 3 )             ; smooth slp

 if(j.eq.0) then 
   ;;;generate the zoom in information;;;;
   fig_type  = "png" ;;figure types: pdf or png
   fig_dir   = "./"
   fignam    = "fig_ivt_ne30_mean_"+tptag
   fig_type@wkWidth  = 2500
   fig_type@wkHeight = 2500 
   gridthick         = 2.0
   lthick            = 2.0
   wks               = gsn_open_wks(fig_type, fig_dir + fignam)
   plot              = new (ngrps,graphic)

   cmap = read_colormap_file("precip4_11lev") 
   cmap(0,:) = 0.0 
  ;cnlevs = (/250,300,400,500,600,700,800,900,1000,1200,1400/)*1.0
   cnlevs = (/150,200,250,300,350,400,450,500,550,600/)*1.0

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
  res1@cnMinLevelValF                 =  960 ;989.0
  res1@cnMaxLevelValF                 =  1040.0
  res1@cnLevelSpacingF                =  3.0
  res1@gsnLeftString                  = ""
  res1@gsnRightString                 = "SLP(hPa, contour)" ;timstr
  res1@cnLinesOn                      = True   ; turn off contour lines
  res1@cnLineLabelsOn                 = True   ; turn off contour labels
  res1@cnInfoLabelOn                  = False
  res1@cnLineLabelDensityF            = 2.0
  res1@cnLineLabelFontHeightF         = FontHeightF*0.8
  res1@cnLineLabelBackgroundColor     = -1
  res1@gsnContourLineThicknessesScale = 2.0
  res1@cnLineThicknessF               = 1.5
  res1@cnLineColor                    = "Black"
  res1@gsnContourZeroLineThicknessF   = 2.
  res1@gsnContourNegLineDashPattern   = 1

;;set vector;;
  res_vc                            = res
  res_vc@vcGlyphStyle               = "LineArrow" ;"CurlyVector"
  res_vc@vcLineArrowThicknessF      = 3
  res_vc@vcLineArrowHeadMaxSizeF    = 0.008
  res_vc@vcMinDistanceF             = 0.02
  res_vc@vcRefLengthF               = 0.03
  res_vc@vcRefAnnoOn                = True
  res_vc@vcRefAnnoSide              = "Bottom"
  res_vc@vcRefAnnoString2On         = True
  res_vc@vcRefAnnoPerimOn           = True
  res_vc@vcRefAnnoOrthogonalPosF    = -1.05
  res_vc@vcRefAnnoParallelPosF      = 0.15
  res_vc@vcRefAnnoBackgroundColor   = "White";"Transparent"; "Purple"
  res_vc@vcVectorDrawOrder          = "PostDraw"
  res_vc@gsnRightString             = "" 
  res_vc@vcRefAnnoFontHeightF       = FontHeightF*0.65
  res_vc@vcRefMagnitudeF            = 250
  res_vc@vcRefAnnoString1           = res_vc@vcRefMagnitudeF
  res_vc@vcRefAnnoString2           = "(kg/m/s)"

  res_vc0 = res_vc
  tuq = var(0,:,:)
  tvq = var(1,:,:)
  tuq = where(tuq.lt.min(cnlevs),0,tuq)
  tvq = where(tvq.lt.min(cnlevs),0,tvq)
  vector(j) = gsn_csm_vector(wks,tuq,tvq,res_vc0)
  delete([/tuq,tvq/])

;---Temperature filled contour plot
  ref_res                       = res
  ref_res@cnFillOn              = True  
  ref_res@cnLevelSelectionMode  = "ExplicitLevels"
 ;ref_res@cnFillMode            = "RasterFill"
  ref_res@cnInfoLabelOn         = False
  ref_res@cnLineLabelsOn        = False
  ref_res@cnLinesOn             = True
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

  ref_res@cnFillPalette           := cmap 
  ref_res@cnLevels                := cnlevs 
  ref_res@cnLevelSelectionMode    := "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF        := lthick
  ref_res@gsnRightString          := rstr 
  ref_res@gsnLeftString           := Labels(j)

  ref_res0 = ref_res
  countour(j) = gsn_csm_contour(wks,var(2,:,:),ref_res0)

  overlay(plot(j),countour(j))
  overlay(plot(j),vector(j))

  delete([/lat,lon,var,res,ref_res,res,res_vc/])

 end do 

 resP                     = True
 resP@gsnMaximize              = True
 resP@gsnPaperOrientation      = "portrait"
 resP@gsnPanelLabelBar         = True
 resP@pmLabelBarParallelPosF   =  0.01
 resP@pmLabelBarOrthogonalPosF = -0.02
 resP@lbLabelFontHeightF       = FontHeightF*0.4
 resP@lbTitleFontHeightF       = FontHeightF*0.45
 resP@lbTitleString            = "IVT (kg m~S~-1~N~ s~S~-1~N~)"
 resP@lbOrientation            = "horizontal"  ; vertical label bar
 resP@pmLabelBarHeightF        = 0.05
 resP@pmLabelBarWidthF         = 0.90
 resP@gsnPanelMainString       = "" ;"Current time: " + timstr
 resP@gsnPanelMainFontColor    = "Black"
 resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6
 resP@gsnPanelMainPosXF        = 0.55
;resP@gsnPanelMainPosYF        = 0.4
;resP@txPosYF                  = 0.7
 resP@gsnPanelBottom           = 0.1
 resP@gsnPanelYWhiteSpacePercent = 0.	
 resP@gsnPanelXWhiteSpacePercent = 3.          
 gsn_panel(wks,plot,(/1,3/),resP)

end
