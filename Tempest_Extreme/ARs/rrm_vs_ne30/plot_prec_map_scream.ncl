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
 Explist   = (/"ERA5","E3SM_ARNEP_721x1440_CTL", "E3SM_ARNEP_721x1440_NDG"/)
 Labels    = (/"ERA5","RRM_CTL","RRM_NDG"/)
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
 minlat =  10.0 ;0.0    ; deg N
 maxlat =  60.0 ;57.0    ; deg N
 minlon = 198.0 ;-100.0    ; - deg E = deg W
 maxlon = 265.0 ;-20.0    ; - deg E = deg W
 cenlon = 200.0 ;

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

  nvars    = 1
  nlat     = dimsizes(lat)
  nlon     = dimsizes(lon)

  if(dimsizes(dimsizes(f0->PRECT)).ge.3)
   var      = dim_avg_n_Wrap(f0->PRECT(indtim,:,:),0)
  else
   var      = dim_avg_n_Wrap(f0->PRECT(indtim,:),0) 
   var@lat  = lat 
   var@lon  = lon
  end if 
  var@FillValue = -9999.0
  delete([/tdate,indtim/])

  if(expnam.eq."ERA5") then
   var = var*1000.0 * 24.0
   var@units = "mm day~S~-1~N~"
  else
   var = var * 1000.0 * 3600.0 * 24.0
   var@units = "mm day~S~-1~N~"
  end if

  if(j.eq.0) then 
    vob = var 
  end if 

  ;calculate the global/regional mean values
  rad   = 4.*atan(1.)/180.    ;
  wgty  = lat
  wgty  = cos(lat*rad)      ; cosine weights
  if(dimsizes(dimsizes(var)).ge.2)
    gomn  = wgt_areaave_Wrap(vob({minlat:maxlat},{minlon:maxlon}),\
                             wgty({minlat:maxlat}), 1.0, 0)
    gosd  = stddev(vob({minlat:maxlat},{minlon:maxlon}))
    gbmn  = wgt_areaave_Wrap(var({minlat:maxlat},{minlon:maxlon}),\
                             wgty({minlat:maxlat}), 1.0, 0)
    grms  = wgt_arearmse(var({minlat:maxlat},{minlon:maxlon}),\
                         vob({minlat:maxlat},{minlon:maxlon}), \
                         wgty({minlat:maxlat}), 1.0, 0)
  else
    lon   = where(lon.le.0,lon+360.0,lon)
    indc  = ind(lat.ge.minlat.and.lat.le.maxlat.and.lon.ge.minlon.and.lon.le.maxlon)
    gomn  = dim_avg_wgt_n(vob(indc), wgty(indc), 0.0, 0)
    gosd  = stddev(vob(indc))
    gbmn  = dim_avg_wgt_n(var(indc), wgty(indc), 0.0, 0)
    grms  = dim_avg_wgt_n((var(indc)-vob(indc))^2, wgty(indc), 0.0, 0)
    delete(indc)
  end if 

  if(j.eq.0) then 
    rstr  = "Mean/RMS (PRECT): "+ sprintf("%.2f",gomn) + "/" + sprintf("%.2f",gosd)
  else
    rstr  = "Bias/RMSE (PRECT): "+ sprintf("%.2f",gbmn-gomn) + "/" + sprintf("%.2f",grms)
  end if 
  delete([/rad,wgty,gbmn,grms/])

 ;printVarSummary(var)
 ;wrf_smooth_2d(var, 3 )             ; smooth slp

 if(j.eq.0) then 
   ;;;generate the zoom in information;;;;
   fig_type  = "png" ;;figure types: pdf or png
   fig_dir   = "./"
   fignam    = "fig_prec_rrm_mean_"+tptag
   fig_type@wkWidth  = 2500
   fig_type@wkHeight = 2500 
   gridthick         = 2.0
   lthick            = 2.0/10.0
   wks               = gsn_open_wks(fig_type, fig_dir + fignam)
   plot              = new (ngrps,graphic)

  ;cmap = read_colormap_file("WhViBlGrYeOrRe")
   cmap = read_colormap_file("BlAqGrYeOrReVi200")
   cmap(0,:) = 0.0
   cnlevs = (/0.1,1,2.5,5,7.5,10,12,14,16,18,20,25,30/) ; "mm/hr"

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
  countour(j) = gsn_csm_contour(wks,var,ref_res0)
  overlay(plot(j),countour(j))

  delete([/lat,lon,var,res,ref_res,res/])

 end do 

 resP                     = True
 resP@gsnMaximize              = True
 resP@gsnPaperOrientation      = "portrait"
 resP@gsnPanelLabelBar         = True
 resP@pmLabelBarParallelPosF   =  0.01
 resP@pmLabelBarOrthogonalPosF = -0.02
 resP@lbLabelFontHeightF       = FontHeightF*0.4
 resP@lbTitleFontHeightF       = FontHeightF*0.45
 resP@lbTitleString            = "PRECT (mm day~S~-1~N~)"
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
