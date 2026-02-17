;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;This script is used to plot the wrf simulations vs NEXRAD data ;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;  LOAD FUNCTIONS AND PROCEDURES
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRFUserARW.ncl"

begin

 datdir = "/pscratch/sd/z/zhan391/DARPA_project/post_process/AR_201702/"
 ;;;;;;;;;;select the regions to calculate the rmse;;;;;;;;;;;;;
 figkey    = "NDGALL_ARs"
 Explist   = (/"ERA5", "SCREAM"/)
 Labels    = (/"ERA5", "SCREAM"/)
 ngrps     = dimsizes(Explist)

 plot      = new (ngrps, graphic)
 vector    = new (ngrps, graphic)
 countour  = new (ngrps, graphic)
 countourx = new (ngrps, graphic)

 year      = 2017
 timtag    = "2017020100-2017022818"

 date      = yyyymmddhh_time(2007, 2017, 6, "integer")
 timplt    = date({2017020100:2017021600})
 ntimes    = dimsizes(timplt)

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

 do i = 0,ntimes-1,1

  do j = 0,ngrps-1,1

   expnam = Explist(j)
   flnm0  = datdir + expnam +"."+timtag+".nc"
   f0     = addfile (flnm0, "r")
  ;ListSetType (f1, "join")        ; concatenate (default)
  
   f1     = addfile(datdir+"E3SMv2_CLIM."+timtag+".nc","r")
   PHIS   = f1->PHIS(0,:,:)
 
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

   if(expnam.eq."ERA5") then 
     unam = "U_ref"
     vnam = "V_ref"
     tnam = "T_ref"
     qnam = "Q_ref"
     psnam= "PS_ref"
   else if (expnam.eq."E3SMv2_CLIM") then
     unam = "U"
     vnam = "V"
     tnam = "T"
     qnam = "Q"
     psnam= "PS"
   else 
     unam = "U_bf_ndg"
     vnam = "V_bf_ndg"
     tnam = "T_bf_ndg"
     qnam = "Q_bf_ndg"
     psnam= "PS_bf_ndg"
   end if 
   end if 

   lat    = f0->lat
   lon    = f0->lon

   g      = 9.81                                ; m/s
   hyai   = f0->hyai
   hybi   = f0->hybi
   p0     = f0->P0

   dum_ps = f0->$psnam$(i,:,:)
   dp     = dpres_hybrid_ccm (dum_ps,p0,hyai,hybi)  ; Pa [kg/(m s2)]

   tdate  = cd_calendar(f0->time,-3)
   indtim = ind(tdate.eq.timcur)
   U      = f0->$unam$(indtim,:,:,:)     ; K  (time,lev,lat,lon)
   V      = f0->$vnam$(indtim,:,:,:)
   T      = f0->$tnam$(indtim,:,:,:)
   Q      = f0->$qnam$(indtim,:,:,:)     ; K  (time,lev,lat,lon)
   delete([/tdate,indtim/])

   nlat   = dimsizes(lat)
   nlon   = dimsizes(lon)

   nvars = 3
   var   = new((/nvars,nlat,nlon/),typeof(lat))
   var!1 = "lat"
   var!2 = "lon"
   var@lat = lat
   var@lon = lon

  ;print(min(PS) +  " "+max(PS))
   vopt = 1     ; vertically weighted sum
   uqavg = wgt_vertical_n(U*Q, dp, vopt, 0)  ; (lat,lon)
   vqavg = wgt_vertical_n(V*Q, dp, vopt, 0)  ; (lat,lon)
   uqavg = uqavg/g
   vqavg = vqavg/g
   uqavg@long_name = "UQ*DP/g: vertically weighted sum"
   uqavg@units     = "kg/m2 m/s"
   vqavg@long_name = "VQ*DP/g: vertically weighted sum"
   vqavg@units     = "kg/m2 m/s"

   var(0,:,:) = (/uqavg/)
   var(1,:,:) = (/vqavg/)
   var(2,:,:) = (/sqrt(uqavg^2+vqavg^2)/)
   var!1      = "lat"
   var!2      = "lon"
   var&lat    = lat
   var&lon    = lon 
   var@FillValue = -9999.0

 ;printMinMax(var,False)
 ;printVarSummary(var)
 ;wrf_smooth_2d(var, 3 )             ; smooth slp

  ;psrf = dum_ps
   hyam = f0->hyam   ; read from a file the mid-layer coef
   hybm = f0->hybm   ; read from a file
   pres = pres_hybrid_ccm(dum_ps,p0,hyam,hybm)  ; (ntim,klvl,nlat,lmon), [Pa]
   nlev = dimsizes(hyam)
   psrf = pslec(T(nlev-1,:,:),PHIS,dum_ps,pres(nlev-1,:,:))
   psrf!0 = "lat"
   psrf!1 = "lon"
   psrf&lat = lat 
   psrf&lon = lon 

   psrf           = psrf/100.0
   psrf@units     = "hPa"
   psrf@FillValue = -9999.0

  ;printVarSummary(psrf)
  ;printMinMax(psrf,False)
  ;exit

   delete([/U,V,dum_ps,Q,dp,g,uqavg,vqavg, hyam, hybm, pres/])

 ;printMinMax(var,False)
 ;printVarSummary(var)
 ;wrf_smooth_2d(var, 3 )             ; smooth slp

  if(j.eq.0) then 
    ;;;generate the zoom in information;;;;
    fig_type  = "png" ;;figure types: pdf or png
    fig_dir   = "./ivt_sub_fig/"
    fignam    = "fig_ivt_"+timcur
    fig_type@wkWidth  = 2500
    fig_type@wkHeight = 2500 
    gridthick         = 2.0
    lthick            = 2.0
    wks               = gsn_open_wks(fig_type, fig_dir + fignam)
    plot              = new (ngrps,graphic)

    cmap = read_colormap_file("precip4_11lev") 
    cmap(0,:) = 0.0 
    cnlevs = (/250,300,400,500,600,700,800,900,1000,1200,1400/)*1.0

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

  countourx(j) = gsn_csm_contour(wks,psrf,res1)
  overlay(plot(j),countourx(j))

;;set vector;;
  res_vc                            = res
  res_vc@vcGlyphStyle               = "LineArrow" ;"CurlyVector"
  res_vc@vcLineArrowThicknessF      = 3
  res_vc@vcLineArrowHeadMaxSizeF    = 0.008
  res_vc@vcMinDistanceF             = 0.02
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

  tuq = var(0,:,:)
  tvq = var(1,:,:)
  tuq = where(tuq.lt.min(cnlevs),0,tuq)
  tvq = where(tvq.lt.min(cnlevs),0,tvq)
  vector(j) = gsn_csm_vector(wks,tuq,tvq,res_vc)
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

  ref_res@cnFillPalette       := cmap 
  ref_res@cnLevels            := cnlevs 
  ref_res@cnLevelSelectionMode:= "ExplicitLevels" ;"ManualLevels"
  ref_res@cnLineThicknessF    := lthick
  ref_res@gsnRightString      := "" ;"IVT (kg m~S~-1~N~ s~S~-1~N~)"
  ref_res@gsnLeftString       := Labels(j)

  countour(j) = gsn_csm_contour(wks,var(2,:,:),ref_res)

  overlay(plot(j),countour(j))
  overlay(plot(j),vector(j))

  delete([/lat,lon,var,psrf,res,ref_res,res,res_vc/])

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
  resP@gsnPanelMainString       = "Current time: " + timstr
  resP@gsnPanelMainFontColor    = "Black"
  resP@gsnPanelMainFontHeightF  = FontHeightF * 0.6
  resP@gsnPanelMainPosXF        = 0.55
 ;resP@gsnPanelMainPosYF        = 0.4
 ;resP@txPosYF                  = 0.7
  resP@gsnPanelBottom           = 0.1
  resP@gsnPanelYWhiteSpacePercent = 0.	
  resP@gsnPanelXWhiteSpacePercent = 3.          
  gsn_panel(wks,plot,(/1,ngrps/),resP)

end do 

end
